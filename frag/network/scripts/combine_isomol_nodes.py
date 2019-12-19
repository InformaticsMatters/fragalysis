#!/usr/bin/env python3
# coding=utf-8

# A utility to combine (de-duplicate) the compressed isomol node CSV files
# from a prior combination of molecule builds.
#
# Alan Christie
# December 2019

import argparse
import glob
import gzip
import logging
import os
import sys

# Configure basic logging
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')

out_hdlr = logging.StreamHandler()
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(out_hdlr)

_RENAME_SUFFIX = '.combined'
_OUTPUT_FILENAME = 'combined-isomol-nodes.csv'

# The map of combined molecules.
# The index is the molecule,
# the content is a list of vendor compounds
_COMBINED_MOLS = {}
# The isomol header, used to start the output file.
# Assumed to be the value of the first line of the first file.
_ISOMOL_HEADER = ''
# The list of processed files
# Which will be renamed on success.
_COMBINED_FILES = []
# Number of (input) nodes encountered
_NODES_ENCOUNTERED = 0


def combine_isomols(isomol_filename: str) -> None:
    """Combines (de-duplicates) the isomol entries with the current set.
    The input file is assumed to be a compressed (gzip) file.

    The file contains 3 columns: -

    - smiles:ID string (and a possible namespace)
    - cmpd_ids:string[] (a list of strings)
    - :LABEL (a list of labels)

    The namespace is expected to be shared between for each IsoMol file
    """
    global _COMBINED_MOLS
    global _ISOMOL_HEADER
    global _NODES_ENCOUNTERED

    logger.info('Processing "{}"'.format(isomol_filename))
    with gzip.open(isomol_filename, 'rt') as this_isomol_file:

        # Skip (but keep) the header
        hdr = this_isomol_file.readline().strip()
        if not _ISOMOL_HEADER:
            _ISOMOL_HEADER = hdr
        else:
            # The header MUST be the same for each file
            if hdr != _ISOMOL_HEADER:
                logger.error('Mismatched header, expected "{}"'
                             ' from the first file but got "{}"'
                             ' from "{}"'.format(_ISOMOL_HEADER,
                                                 hdr,
                                                 isomol_filename))
                exit(1)

        # Now process the rest of the file
        for line in this_isomol_file:

            line_items = line.strip().split(',')
            # Must have at least 1 column value
            smiles = line_items[0]
            compounds = set()
            labels = set()

            if len(line_items) > 1:
                compounds.update(line_items[1].split(';'))
            if len(line_items) > 2:
                labels.update(line_items[2].split(';'))
            # Do we have this SMILES?
            if smiles in _COMBINED_MOLS:
                # Augment the current set of compounds
                compounds_set = _COMBINED_MOLS[smiles][0]
                compounds_set.update(compounds)
                _COMBINED_MOLS[smiles][0] = compounds_set
                # Augment the current set of labels
                labels_set = _COMBINED_MOLS[smiles][1]
                labels_set.update(labels)
                _COMBINED_MOLS[smiles][1] = labels_set
            else:
                _COMBINED_MOLS[smiles] = [compounds, labels]

            _NODES_ENCOUNTERED += 1

    # Remember the file (for renaming)
    _COMBINED_FILES.append(isomol_filename)


def write_combination(output_dir: str) -> int:
    """Writes the _COMBINED_MOLS map to the combined file.
    The output consists of three columns; a SMILES string,
    list of compounds and list of labels.
    """

    global _COMBINED_MOLS
    assert _ISOMOL_HEADER

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    output_filename = os.path.join(output_dir,
                                   '{}.gz'.format(_OUTPUT_FILENAME))
    logger.info('Writing %s...', output_filename)
    with gzip.open(output_filename, 'wt') as output_gzip_file:

        # Write the header...
        num_written = 0
        output_gzip_file.write(_ISOMOL_HEADER + '\n')
        for smiles in _COMBINED_MOLS:
            compounds = list(_COMBINED_MOLS[smiles][0])
            labels = list(_COMBINED_MOLS[smiles][1])
            output_gzip_file.write('{},{},{}\n'.format(smiles,
                                                       ';'.join(compounds),
                                                       ';'.join(labels)))
            num_written += 1

    return num_written


def rename_combined_files() -> None:

    for filename in _COMBINED_FILES:
        logger.info('Renaming {}'.format(filename))
        os.rename(filename, filename + _RENAME_SUFFIX)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='IsoMol Node Combiner',
                                     epilog="If successful, i.e. if there are"
                                            " files to combine, the original"
                                            " files are renamed so they do not"
                                            " appear in any subsequent"
                                            " load-neo4j script that would"
                                            " be generated by the"
                                            " graph_create_loader script.")
    parser.add_argument('isomol_dir',
                        help='The directory of isomol files'
                             '(more than 1 expected).')

    args = parser.parse_args()

    # Input must exists as a directory
    if not os.path.exists(args.isomol_dir) and \
            not os.path.isdir(args.isomol_dir):
        logger.error('The IsoMol directory does not exist')
        sys.exit(1)

    # Before we open the output file
    # get a list of all the input files
    # If there's only one - we have nothing to do.
    isomol_files = glob.glob('{}/*-isomol-nodes.csv.gz'.format(args.isomol_dir))

    if len(isomol_files) == 0:
        logger.error("Couldn't find any IsoMol node files")
        sys.exit(1)

    if len(isomol_files) > 1:

        # Process all the IsoMol files (there will be more than one)...
        for isomol_file in isomol_files:
            combine_isomols(isomol_file)
        # Write the combination and rename the inputs...
        num_mols = write_combination(args.isomol_dir)
        rename_combined_files()

        # Summary
        logger.info('{:,} files combined'.format(len(isomol_files)))
        logger.info('{:,} nodes encountered'.format(_NODES_ENCOUNTERED))
        logger.info('{:,} nodes written'.format(num_mols))
