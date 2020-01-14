#!/usr/bin/env python

"""standardise_xchem_compounds.py

Processes xchem vendor compound files, and generates a 'standard'
tab-separated output.

We create a 'standardised-compounds.tab.gz' file that contains a 1st-line
'header' formed from the _OUTPUT_COLUMNS list.

Alan Christie
December 2019
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

# Configure basic logging
logger = logging.getLogger('xchem')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The columns in our output file.
# In this file we don't add any of our own.
_OUTPUT_COLUMNS = parser.STANDARD_COLUMNS

# The minimum number of columns in the input files and
# and a map of expected column names indexed by (0-based) column number.
#
# The 'standardised' files contain at least 2 columns...
#
# smiles    0
# ID        1

expected_min_num_cols = 2
smiles_col = 0
compound_col = 1
# User upper-case (regardless)
# this is how column header comparisons are made here
expected_input_cols = {compound_col: 'ID',
                       smiles_col: 'SMILES'}

# The output file.
# Which will be gzipped.
output_filename = 'standardised-compounds.tab'

# The prefix we use in our fragment file
prefix = 'XCHEM:'

# All the vendor compound IDs
vendor_compounds = set()
# A map of duplicate compounds and the number of duplicates.
# The index uses the vendor's original ID value, not our prefixed value.
duplicate_suffix = '-duplicate-'
vendor_duplicates = {}

# Various diagnostic counts
num_vendor_mols = 0
num_vendor_molecule_failures = 0
num_inchi_failures = 0

# The line rate at which the process writes updates to stdout.
report_rate = 100


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


def standardise_vendor_compounds(file_name, source_id, limit):
    """Process the given file and standardise the vendor
    information, writing it as tab-separated fields to the output.

    As we load the Vendor compounds we 'standardise' the SMILES and
    determine whether they represent an isomer or not.

    :param output_file: The tab-separated standardised output file
    :param file_name: The (compressed) file to process
    :param limit: Limit processing to this number of values (or all if 0)
    :returns: The number of items processed
    """
    global vendor_compounds
    global vendor_duplicates
    global duplicate_suffix
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
        # Check salient columns (ignoring case)...
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

            osmiles = fields[smiles_col].strip()
            vendor_id = fields[compound_col].strip()
            compound_id = prefix + vendor_id

            # Add the compound (expected to be unique)
            # to our set of 'all compounds'.
            if compound_id in vendor_compounds:
                # Get the number of duplicates (default of 1)
                # using the vendor's original ID as a key
                duplicate_count = vendor_duplicates.setdefault(vendor_id, 1)
                compound_id += '{}{}'.format(duplicate_suffix, duplicate_count)
                # Increment for any further duplicates
                vendor_duplicates[vendor_id] = duplicate_count + 1
            else:
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
                loader.insert_mol(osmiles, compound_id, std_info=std_info, session=session)

            # Enough?
            num_processed += 1
            if limit and num_processed >= limit:
                break

    loader.complete_input(session, num_processed, num_vendor_molecule_failures, num_inchi_failures)
    session.close()
    return num_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Standardiser (Enamine/CS-BB)')
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
