#!/usr/bin/env python

"""standardise_chembl_compounds.py

Processes ChEMBL compound files, and generates a 'standard'
tab-separated output.
Download ChEMBL from here: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/

Unlike most other sources of molecules where the structures are specified as SMILES with ChEMBL
the structures are in SDF format which means that the osmiles data must be generated from the CTAB
section by RDKit. It is therefore not 100% guaranteed to exactly represent the structure.

We create a 'chembl-standardised-compounds.tab' file that contains a 1st-line
'header' consisting of OSMILES\tISO_SMILES\tNONISO_SMILES\tHAC\tCMPD_ID.

Tim Dudgeon
December 2019
"""

import argparse
import gzip
import logging
import sys

from rdkit import Chem, RDLogger

from frag.utils import standardise_utils
from frag.std_utils import parser

# Configure basic logging
logger = logging.getLogger('chembl')
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

# The output file.
# Which will be gzipped.
output_filename = 'standardised-compounds.tab'

# The prefix we use in our fragment file
vendor_prefix = 'CHEMBL:'

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
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


def standardise_vendor_compounds(root_dir, file_name, limit):
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
    global num_vendor_mols
    global num_vendor_molecule_failures
    global num_inchi_failures

    logger.info('Standardising %s...', file_name)

    count = 0
    num_processed = 0

    vendor = 'chembl'

    vendor_paths_file = open('/'.join([root_dir, 'datasets', vendor + '.txt']) , "wt")

    with gzip.open(file_name, 'rb') as input:

        suppl = Chem.ForwardSDMolSupplier(input)

        new_inchis = 0
        new_noniso = 0
        new_iso = 0

        for mol in suppl:

            count += 1

            if mol is None:
                logger.warning("Failed to handle record %s", count)
                continue

            chembl_id = mol.GetProp('chembl_id')

            if count % report_rate == 0:
                logger.info(' ...at compound {:,}'.format(count))

            osmiles = Chem.MolToSmiles(mol)
            # generate the ID. The data in the chembl_id field looks like CHEMBL153534
            # so we need to string off the CHEMBL bit as our generated ID looks like CHEMBL:153534
            compound_id = vendor_prefix + chembl_id[6:]

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
                x,y,z = standardise_utils.write_data(root_dir + '/tree', vendor, osmiles, std_info, compound_id, vendor_paths_file)
                new_inchis += x
                new_noniso += y
                new_iso += z

            # Enough?
            num_processed += 1
            if limit and num_processed >= limit:
                break

        logger.info("Created %s new InChi %s new noniso and %s new iso", new_inchis, new_noniso, new_iso)

    vendor_paths_file.close()
    return num_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Standardiser (ChEMBL)')
    parser.add_argument('--input', '-i',
                        help='The ChEMBL chembl_*.sdf.gz file to be processed.')
    parser.add_argument('--output', '-o',
                        help='The output directory')
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

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    logger.info('Writing to %s...', args.output)
    num_processed = 0

    num_processed += standardise_vendor_compounds(args.output,
                                                  args.input,
                                                  args.limit)

    logger.info("Processed %s records", num_processed)

    # Summary
    logger.info('{:,} vendor molecules'.format(num_vendor_mols))
    logger.info('{:,} vendor molecule failures'.format(num_vendor_molecule_failures))
    logger.info('{:,} inchi generation failures'.format(num_inchi_failures))
