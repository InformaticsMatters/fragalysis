#!/usr/bin/env python

"""standardise_chembl_compounds.py

Processes ChEMBL compound files, and generates a 'standard'
tab-separated output.
Download ChEMBL from here: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/

Unlike most other sources of molecules where the structures are specified as SMILES with ChEMBL
the structures are in SDF format which means that the osmiles data must be generated from the CTAB
section by RDKit. It is therefore not 100% guaranteed to exactly represent the structure.

Run like this:
python -m frag.network.scripts.standardise_chembl_compounds_3 -i ~/data/chembl/chembl_25.sdf.gz -s 1

Tim Dudgeon
December 2019
"""

import argparse
import gzip
import logging
import sys, os

from rdkit import Chem, RDLogger

from frag.utils import standardise_utils
from frag.std_utils import parser

from db import fairmolecules

from db.fairmolecules import MoleculeLoader

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
report_rate = 1000


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

    loader = MoleculeLoader()
    session = loader.create_session()

    base_name = os.path.basename(file_name)

    with gzip.open(file_name, 'rb') as input:

        suppl = Chem.ForwardSDMolSupplier(input)

        loader.create_input(session, base_name, source_id)

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
                loader.insert_mol(osmiles, compound_id, std_info=std_info, session=session)

            # Enough?
            num_processed += 1
            if limit and num_processed >= limit:
                break

    loader.complete_input(session, num_processed, num_vendor_molecule_failures, num_inchi_failures)
    session.close()
    return num_processed

if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Standardiser (ChEMBL)')
    parser.add_argument('--input', '-i',
                        help='The ChEMBL chembl_*.sdf.gz file to be processed.')
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
