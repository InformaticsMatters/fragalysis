#!/usr/bin/env python
#
# Based on build_db.py, this module builds the graph network from
# the Informatics Matters 'standard' (uncompressed) file representation.
#
# Alan Christie
# January 2020

import argparse
import os
import sys
import collections

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network
from frag.std_utils.parser import verify_header, get_standard_items

cache = set()
node_count = 0
edge_count = 0

nodes_f = None
edges_f = None


def write_node(node):
    global node_count
    # print("writing node", node.SMILES)
    nodes_f.write(node.as_csv() + ',,F2\n')
    node_count += 1
    cache.add(node.SMILES)


def write_edge(edge):
    global edge_count
    # print("writing edge", edge.get_label())
    edges_f.write(edge.as_csv() + ',FRAG\n')
    edge_count += 1


def write_data(node_holder):
    need_further_processing = []
    # if no edges then just set status to complete
    if node_holder.size()[1] == 0:

        if node_holder.node_list:
            node = node_holder.node_list.pop()
            smiles = node.SMILES
            if smiles not in cache:
                cache.add(smiles)
                write_node(node)

    else:
        # so we need to process the edges

        # group the edges by parent and prepare the nonisomols
        grouped_parent_edges = collections.OrderedDict()
        parent_nodes = {}
        r_smiles = None
        for edge in node_holder.get_edges():
            # handle the parent nonisomol
            p_node = edge.NODES[0]
            p_smiles = p_node.SMILES
            parent_nodes[p_smiles] = p_node
            if p_smiles not in cache:
                if not r_smiles:
                    r_smiles = p_smiles
                if p_smiles in grouped_parent_edges:
                    grouped_parent_edges[p_smiles].append(edge)
                else:
                    # it's new to this group of edges
                    grouped_parent_edges[p_smiles] = [edge]

        for p_smiles in grouped_parent_edges:
            p_node = parent_nodes[p_smiles]
            write_node(p_node)
            edges = grouped_parent_edges[p_smiles]
            grouped_child_edges = collections.OrderedDict()
            # child_nodes = {}
            for edge in edges:
                c_node = edge.NODES[1]
                c_smiles = c_node.SMILES
                # child_nodes[c_smiles] = c_node
                if c_smiles in grouped_child_edges:
                    grouped_child_edges[c_smiles].append(edge)
                else:
                    grouped_child_edges[c_smiles] = [edge]

            for c_smiles in grouped_child_edges:
                if c_smiles not in cache:
                    need_further_processing.append(c_smiles)
                edges = grouped_child_edges[c_smiles]
                for edge in edges:
                    write_edge(edge)

    return need_further_processing


def fragment_and_write(smiles, base_dir, verbosity):
    node_holder = fragment_mol(smiles, base_dir, verbosity)
    need_further_processing = write_data(node_holder)
    for smiles in need_further_processing:
        fragment_and_write(smiles, base_dir, verbosity)


def fragment_mol(smiles, base_dir, verbosity):
    attrs = []
    attr = Attr(smiles, ["EM"])
    attrs.append(attr)

    # Build the network
    node_holder = NodeHolder(iso_flag=False)
    max_frags = 0
    node_holder = build_network(attrs, node_holder, max_frags, base_dir,
                                verbosity=verbosity, recurse=False)
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
    parser.add_argument('-r', '--report-interval', type=int, default=0,
                        help='Reporting interval')

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", dest="verbosity", action="store_const", const=1)
    group.add_argument("-vv", dest="verbosity", action="store_const", const=2)

    parser.set_defaults(verbosity=0)
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

    if not os.path.isdir(args.base_dir):
        os.mkdir(args.base_dir)

    global nodes_f
    global edges_f
    nodes_f = open(os.path.join(args.base_dir, 'nodes.csv'), 'w')
    edges_f = open(os.path.join(args.base_dir, 'edges.csv'), 'w')

    with open(args.input, 'r') as standard_file:

        # Read (and verify) the header...
        hdr = standard_file.readline()
        verify_header(hdr)

        # Process the rest of the file...
        num_skipped = 0
        num_processed = 0

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

            fragment_and_write(std.noniso, args.base_dir, args.verbosity)

            # Enough?
            num_processed += 1
            if args.report_interval > 0 and num_processed % args.report_interval == 0:
                print("Processed mol", num_processed)
            if args.limit and num_processed >= args.limit:
                break

    nodes_f.close()
    edges_f.close()
    print("Wrote {0} nodes and {1} edges".format(node_count, edge_count))


if __name__ == "__main__":
    main()
