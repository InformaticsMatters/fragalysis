#!/usr/bin/env python
# coding=utf-8

"""standardise_molport_compounds.py

Processes MolPort vendor compound files, expected to contain pricing
information and generates a 'standard' tab-separated output.
We create a 'molport-standardised-compounds.tab' file that contains a 1st-line
'header' formed from the _OUTPUT_COLUMNS list.

Alan Christie
January 2019
"""

import argparse
import glob
import gzip
import logging
import os
import sys

from rdkit import RDLogger

from frag.utils import standardise_utils
from frag.std_utils import parser

from db.fairmolecules import MoleculeLoader

# The columns in our output file.
_OUTPUT_COLUMNS = parser.STANDARD_COLUMNS + \
                  ['PRICERANGE_1MG',
                   'PRICERANGE_5MG',
                   'PRICERANGE_50MG',
                   'BEST_LEAD_TIME']

# Configure basic logging
logger = logging.getLogger('molport')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The minimum number of columns in the input data and
# a map of expected column names indexed by column number.
#
# The molecule data is spread over a number of `txt.gz` files
# (i.e. files like `iis_smiles-000-000-000--000-499-999.txt.gz`)
# in a common directory where the files have the following header
# names and (0-based) positions:
#
# SMILES                0
# SMILES_CANONICAL      1
# MOLPORTID             2
# STANDARD_INCHI        3
# INCHIKEY              4
# PRICERANGE_1MG        5
# PRICERANGE_5MG        6
# PRICERANGE_50MG       7
# BEST_LEAD_TIME        8

expected_min_num_cols = 9
smiles_col = 0
compound_col = 2
cost_col = {1: 5, 5: 6, 50: 7}
blt_col = 8
expected_input_cols = {smiles_col: 'SMILES',
                       compound_col: 'MOLPORTID',
                       cost_col[1]: 'PRICERANGE_1MG',
                       cost_col[5]: 'PRICERANGE_5MG',
                       cost_col[50]: 'PRICERANGE_50MG',
                       blt_col: 'BEST_LEAD_TIME'}

# The output file.
# Which will be gzipped.
output_filename = 'standardised-compounds.tab'

# The compound identifier prefix
# the vendor uses in the the compound files...
supplier_prefix = 'MolPort-'
# The prefix we use in our fragment file
# and the prefix we use for our copy of the
molport_prefix = 'MOLPORT:'

# All the vendor compound IDs
vendor_compounds = set()

# Various diagnostic counts
num_vendor_mols = 0
num_vendor_molecule_failures = 0
num_inchi_failures = 0

# The line rate at which the process writes updates to stdout.
report_rate = 50000


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: %s', msg)
    sys.exit(1)


def standardise_vendor_compounds(file_name, source_id, limit):

    global vendor_compounds
    global num_vendor_mols
    global num_vendor_molecule_failures
    global num_inchi_failures

    logger.info('Standardising %s...', file_name)

    line_num = 0
    num_processed = 0

    loader = MoleculeLoader()
    session = loader.create_session()

    base_name = os.path.basename(file_name)

    with gzip.open(file_name, 'rt') as gzip_file:

        # Check first line (a space-delimited header).
        # This is a basic sanity-check to make sure the important column
        # names are what we expect.

        hdr = gzip_file.readline()
        field_names = hdr.split('\t')
        # Expected minimum number of columns...
        if len(field_names) < expected_min_num_cols:
            error('expected at least {} columns found {}'.
                  format(expected_input_cols, len(field_names)))
        # Check salient columns...
        for col_num in expected_input_cols:
            actual_name = field_names[col_num].strip().upper()
            if actual_name != expected_input_cols[col_num]:
                error('expected "{}" in column {} found "{}"'.
                      format(expected_input_cols[col_num],
                             col_num,
                             actual_name))

        # Columns look right...

        loader.create_input(session, base_name, source_id)

        for line in gzip_file:

            line_num += 1
            fields = line.split('\t')
            if len(fields) <= 1:
                continue

            if line_num % report_rate == 0:
                logger.info(' ...at compound {:,}'.format(line_num))

            osmiles = fields[smiles_col]
            compound_id = molport_prefix + fields[compound_col].split(supplier_prefix)[1]

            prices = {}
            cost_1mg = fields[cost_col[1]]
            cost_5mg = fields[cost_col[5]]
            cost_50mg = fields[cost_col[50]]

            if cost_1mg:
                p_1mg = parse_cost(cost_1mg)
                if p_1mg:
                    prices[1] = p_1mg
            if cost_5mg:
                p_5mg = parse_cost(cost_5mg)
                if p_5mg:
                    prices[5] = p_5mg
            if cost_50mg:
                p_50mg = parse_cost(cost_50mg)
                if p_50mg:
                    prices[50] = p_50mg

            # Add the compound (expected to be unique)
            # to our set of 'all compounds'.
            if compound_id in vendor_compounds:
                error('Duplicate compound ID ({})'.format(compound_id))
            vendor_compounds.add(compound_id)

            # Standardise and update global maps...
            # And try and handle and report any catastrophic errors
            # from dependent modules/functions.

            std_info = standardise_utils.standardise(osmiles)
            if not std_info.std:
                num_vendor_molecule_failures += 1
                continue
            num_vendor_mols += 1

            if std_info.inchik is None:
                num_inchi_failures += 1
            else:
                loader.insert_mol(osmiles, compound_id, prices=prices, std_info=std_info, session=session)

            # Enough?
            num_processed += 1
            if limit and num_processed >= limit:
                break

    loader.complete_input(session, num_processed, num_vendor_molecule_failures, num_inchi_failures)
    session.close()
    return num_processed

def parse_cost(str):
    if str == 'null':
        return None
    parts = str.split()
    if len(parts) == 1:
        return int(parts[0])
    if len(parts) == 2:
        if parts[0] == '<':
            return (None, int(parts[1]))
        if parts[0] == '>':
            return (int(parts[1]), None)
    elif len(parts) == 3:
        return (int(parts[0]), int(parts[2]))

if __name__ == '__main__':

    print(str(parse_cost('< 50')))
    print(str(parse_cost('> 50')))
    print(str(parse_cost('5 - 50')))
    print(str(parse_cost('null')))
    print(str(parse_cost('5')))
    sys.exit(0)

    parser = argparse.ArgumentParser('Vendor Compound Standardiser (MolPort)')
    parser.add_argument('--input', '-i',
                        help='The ChemSpace file to be processed.')
    parser.add_argument('--source-id', '-s',
                        help='The source ID in the database')
    parser.add_argument('--limit', '-l',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise')

    args = parser.parse_args()


    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Report any limiting...?
    if args.limit:
        logger.warning('Limiting processing to first {:,} molecules'.format(args.limit))

    num_processed = 0

    num_processed += standardise_vendor_compounds(args.input, args.source_id, args.limit)

    logger.info("Processed %s records", num_processed)

    # Summary
    logger.info('{:,} vendor molecules'.format(num_vendor_mols))
    logger.info('{:,} vendor molecule failures'.format(num_vendor_molecule_failures))
    logger.info('{:,} inchi generation failures'.format(num_inchi_failures))

