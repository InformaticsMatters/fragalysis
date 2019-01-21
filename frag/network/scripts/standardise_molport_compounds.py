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

from rdkit import Chem
from rdkit import RDLogger
from frag.utils.rdkit_utils import standardize

# The columns in our output file.
_OUTPUT_COLUMNS = ['O_SMILES',
                   'ISO_SMILES',
                   'NONISO_SMILES',
                   'CMPD_ID',
                   'PRICERANGE_1MG',
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
output_filename = 'molport-standardised-compounds.tab'

# The compound identifier prefix
# the vendor uses in the the compound files...
supplier_prefix = 'MolPort-'
# The prefix we use in our fragment file
# and the prefix we use for our copy of the
molport_prefix = 'MOLPORT:'

# All the vendor compound IDs
vendor_compounds = set()

# Various diagnostic counts
num_vendor_iso_mols = 0
num_vendor_mols = 0
num_vendor_molecule_failures = 0

# The line rate at which the process writes updates to stdout.
# Every 20 million?
augment_report_rate = 20000000


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: %s', msg)
    sys.exit(1)


def standardise_vendor_compounds(output_file, file_name):
    """Process the given file and standardise the vendor (and pricing)
    information, writing it as tab-separated fields to the output.

    As we load the Vendor compounds we 'standardise' the SMILES and
    determine whether they represent an isomer or not.

    :param output_file: The tab-separated standardised output file
    :param file_name: The (compressed) file to process
    """

    global vendor_compounds
    global num_vendor_iso_mols
    global num_vendor_mols
    global num_vendor_molecule_failures

    logger.info('Standardising %s...', file_name)

    num_lines = 0
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
            if field_names[col_num].strip() != expected_input_cols[col_num]:
                error('expected "{}" in column {} found "{}"'.
                      format(expected_input_cols[col_num],
                             col_num,
                             field_names[col_num]))

        # Columns look right...

        for line in gzip_file:

            num_lines += 1
            fields = line.split('\t')
            if len(fields) <= 1:
                continue

            o_smiles = fields[smiles_col]
            compound_id = molport_prefix + fields[compound_col].split(supplier_prefix)[1]

            # Add the compound (expected to be unique)
            # to our set of 'all compounds'.
            if compound_id in vendor_compounds:
                error('Duplicate compound ID ({})'.format(compound_id))
            vendor_compounds.add(compound_id)

            # Standardise and update global maps...
            # And try and handle and report any catastrophic errors
            # from dependent modules/functions.

            mol = None
            std = None
            iso = None
            noniso = None

            try:
                mol = Chem.MolFromSmiles(o_smiles)
            except Exception as e:
                logger.warning('MolFromSmiles(%s) exception: "%s"', o_smiles, e.message)
            if not mol:
                num_vendor_molecule_failures += 1
                logger.error('Got nothing from MolFromSmiles(%s).'
                             ' Skipping this Vendor compound'
                             ' (id=%s line=%s failures=%s)',
                             o_smiles, compound_id, num_lines,
                             num_vendor_molecule_failures)
                continue

            # Got a molecule.
            #
            # Try to (safely) standardise,
            # and create isomeric an non-isomeric representations.

            try:
                std = standardize(mol)
            except Exception as e:
                logger.warning('standardize(%s) exception: "%s"', o_smiles, e.message)
            if not std:
                num_vendor_molecule_failures += 1
                logger.error('Got nothing from standardize(%s).'
                             ' Skipping this Vendor compound'
                             ' (line=%s failures=%s)',
                             o_smiles, num_lines, num_vendor_molecule_failures)
                continue

            try:
                iso = Chem.MolToSmiles(std, isomericSmiles=True, canonical=True)
            except Exception as e:
                logger.warning('MolToSmiles(%s, iso) exception: "%s"', o_smiles, e.message)
            if not iso:
                num_vendor_molecule_failures += 1
                logger.error('Got nothing from MolToSmiles(%s, iso).'
                             ' Skipping this Vendor compound'
                             ' (line=%s failures=%s)',
                             o_smiles, num_lines, num_vendor_molecule_failures)
                continue

            try:
                noniso = Chem.MolToSmiles(std, isomericSmiles=False, canonical=True)
            except Exception as e:
                logger.warning('MolToSmiles(%s, noniso) exception: "%s"', o_smiles, e.message)
            if not noniso:
                num_vendor_molecule_failures += 1
                logger.error('Got nothing from MolToSmiles(%s, noniso).'
                             ' Skipping this Vendor compound'
                             ' (line=%s failures=%s)',
                             o_smiles, num_lines, num_vendor_molecule_failures)
                continue

            # Write the standardised data

            output = [o_smiles,
                      iso,
                      noniso,
                      compound_id,
                      fields[cost_col[1]].strip(),
                      fields[cost_col[5]].strip(),
                      fields[cost_col[50]].strip(),
                      fields[blt_col].strip()]

            output_file.write('\t'.join(output) + '\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Standardiser (MolPort)')
    parser.add_argument('vendor_dir',
                        help='The MolPort vendor directory,'
                             ' containing the ".gz" files to be processed.')
    parser.add_argument('vendor_prefix',
                        help='The MolPort vendor file prefix,'
                             ' i.e. "iis_smiles". Only files with this prefix'
                             ' in the vendor directory will be processed')
    parser.add_argument('output',
                        help='The output directory')

    args = parser.parse_args()

    # Create the output directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if not os.path.isdir(args.output):
        error('output ({}) is not a directory'.format(args.output))

    # -------
    # Stage 1 - Process Vendor Files
    # -------

    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    output_filename = os.path.join(args.output, '{}.gz'.format(output_filename))
    logger.info('Writing %s...', output_filename)
    with gzip.open(output_filename, 'wt') as output_gzip_file:

        # Write the header...
        output_gzip_file.write('\t'.join(_OUTPUT_COLUMNS) + '\n')

        # Process all the Vendor files...
        molport_files = glob.glob('{}/{}*.gz'.format(args.vendor_dir,
                                                     args.vendor_prefix))
        for molport_file in molport_files:
            standardise_vendor_compounds(output_gzip_file, molport_file)

    # Summary
    logger.info('{:,}/{:,} vendor molecules/iso'.format(num_vendor_mols, num_vendor_iso_mols))
    logger.info('{:,} vendor molecule failures'.format(num_vendor_molecule_failures))