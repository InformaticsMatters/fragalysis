#!/usr/bin/env python
#
# Based on build_db.py, this module builds the graph network from
# the Informatics Matters 'standard' (uncompressed) file representation.
#
# Alan Christie
# February 2019

import argparse
import os
import sys
import collections
import time

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network, write_data_as_csv
from frag.std_utils.parser import verify_header, get_standard_items

node_count = 0
edge_count = 0

nodes_f = None
edges_f = None

def write_node(node):
    global node_count
    # print("writing node", node.SMILES)
    nodes_f.write(node.as_csv() + ',,F2\n')
    node_count += 1

def write_edge(edge):
    global edge_count
    # print("writing edge", edge.get_label())
    edges_f.write(edge.as_csv() + ',FRAG\n')
    edge_count += 1

def write_data(node_holder):
    for node in node_holder.node_list:
        write_node(node)
    for edge in node_holder.edge_list:
        write_edge(edge)

def fragment_and_write(chunk, chunk_number):

    t0 = int(round(time.time() * 1000))
    node_holder = fragment_mols(chunk)
    t1 = int(round(time.time() * 1000))
    size = node_holder.size()
    write_data(node_holder)
    print("Handled chunk {0} with {1} nodes and {2} edges".format(chunk_number, size[0], size[1]))
    node_holder = None


def fragment_mols(chunk, verbosity=0):

    attrs = []
    for smiles in chunk:
        attr = Attr(smiles, ["EM"])
        attrs.append(attr)

    # Build the network
    node_holder = NodeHolder(iso_flag=False)
    max_frags = 0
    node_holder = build_network(attrs, node_holder, max_frags, None, verbosity=verbosity, recurse=True)
    # Write the data out
    # print(str(node_holder.size()))
    # for node in node_holder.node_list:
    #     print(str(node))
    # for edge in node_holder.get_edges():
    #         print(str(edge))

    return node_holder

def main():
    """Read in a 'standard' file - then write out into a specified directory
    """
    parser = argparse.ArgumentParser(
        description="Convert un-compressed standard SMILES"
                    " for Astex Fragment network."
    )
    parser.add_argument("--input")
    parser.add_argument("--base_dir")
    parser.add_argument('-c', '--chunk-size',
                        type=int, default=0,
                        help='Chunk size')
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise')
    parser.add_argument('-s', '--skip',
                        type=int, default=0,
                        help='Number of molecules to skip molecules'
                             ' in the input file')
    parser.add_argument('--min-hac',
                        type=int, default=0,
                        help='Limit processing to molecules with at least this'
                             ' number of heavy atoms')
    parser.add_argument('--max-hac',
                        type=int, default=0,
                        help='Limit processing to molecules with no more than'
                             ' this number of heavy atoms')
    parser.add_argument('--max-frag',
                        type=int, default=0,
                        help='Limit processing to molecules with no more than'
                             ' this number of initial fragment (no limit if 0)')
    parser.add_argument("--isomeric", dest="iso_flag", action="store_true")
    parser.add_argument("--non_isomeric", dest="iso_flag", action="store_false")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", dest="verbosity", action="store_const", const=1)
    group.add_argument("-vv", dest="verbosity", action="store_const", const=2)

    parser.set_defaults(verbosity=0)
    parser.set_defaults(iso_flag=True)
    args = parser.parse_args()

    # Do we have an input and base directory?
    if not args.input:
        print('ERROR: Must specify an input')
        sys.exit(1)
    if not os.path.isfile(args.input):
        print('ERROR: input (%s) does not exist' % args.input)
        sys.exit(2)
    if not args.base_dir:
        print('ERROR: Must specify a base directory')
        sys.exit(3)

    global nodes_f
    global edges_f
    nodes_f = open(os.path.join(args.base_dir, "nodes.csv"), "w")
    edges_f = open(os.path.join(args.base_dir, "edges.csv"), "w")

    with open(args.input, 'r') as standard_file:

        # Read (and verify) the header...
        hdr = standard_file.readline()
        verify_header(hdr)

        # Process the rest of the file...
        num_skipped = 0
        num_processed = 0

        chunk = []
        chunk_count = 0

        for line in standard_file:

            std = get_standard_items(line)

            # Do we need to skip molecules before processing?
            if num_skipped < args.skip:
                num_skipped += 1
                continue

            # HAC within range?
            # If not, skip this line.
            if std.hac < args.min_hac or args.max_hac > 0 and std.hac > args.max_hac:
                continue

            smiles = std.noniso
            id = std.cmpd_id

            chunk.append(smiles)

            num_processed += 1
            if num_processed % args.chunk_size == 0:
                chunk_count += 1
                fragment_and_write(chunk, chunk_count)
                chunk = []
                # print("Processed chunk", chunk_count)

            # Enough?
            # print("Processed mol", num_processed)
            if args.limit and num_processed >= args.limit:
                break

        if chunk:
            fragment_and_write(chunk, chunk_count)
            chunk_count += 1

    nodes_f.close()
    edges_f.close()
    print("Wrote {0} nodes and {1} edges in {2} chunks of {3}".format(node_count, edge_count, chunk_count, args.chunk_size))


if __name__ == "__main__":
    main()
