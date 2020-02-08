#!/usr/bin/env python
#
# Based on build_db.py, this module builds the graph network nodes and edges from
# the Informatics Matters 'standard' (uncompressed) file representation.
# The output is these files:
# nodes.csv containing the molecules. The fields are: SMILES, HAC, RAC, NUM_CHILDREN, NUM_EDGES, TIME_MS
# edges.csv containing the edges. The fields are: PARENT_SMILES, CHILD_SMILES, LABEL
# rejects.smi containing the SMILES that were rejected because of the fragment count filter.
#
# Tim Dudgeon
# February 2020

import argparse
import os
import sys
import collections
import time

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network

cache = set()
node_count = 0
edge_count = 0
rejects_count = 0

base_dir = None
nodes_f = None
edges_f = None
rejects_f = None

def write_node(node, time_ms, num_children, num_edges):
    global node_count
    # print("writing node", node.SMILES)
    nodes_f.write(','.join([node.SMILES, str(node.HAC), str(node.RAC), str(num_children), str(num_edges), str(time_ms)]) + '\n')
    node_count += 1
    cache.add(node.SMILES)

def write_edge(edge):
    global edge_count
    # print("writing edge", edge.get_label())
    edges_f.write(edge.as_csv() + '\n')
    edge_count += 1

def write_reject(smiles):
    global rejects_f

    if not rejects_f:
        rejects_f = open(os.path.join(base_dir, "rejects.smi"), "w")
    rejects_f.write(smiles + '\n')


def write_data(node_holder, time_ms, num_children, num_edges):
    need_further_processing = []
    # if no edges then just set status to complete
    if node_holder.size()[1] == 0:
        if node_holder.node_list:
            node = node_holder.node_list.pop()
            smiles = node.SMILES
            if smiles not in cache:
                write_node(node, time_ms, 0, 0)
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

        fcc1 = len(grouped_parent_edges)
        for p_smiles in grouped_parent_edges:
            p_node = parent_nodes[p_smiles]
            write_node(p_node, time_ms, num_children, num_edges)
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

def fragment_and_write(smiles, max_frags=0, verbosity=0):

    global rejects_count

    t0 = time.time()
    node_holder = fragment_mol(smiles, verbosity=verbosity)
    t1 = time.time()
    time_ms = int(round((t1 - t0) * 1000))
    size = node_holder.size()
    # the number of children is the number of nodes minus one (the parent)
    num_children = size[0] - 1
    if 0 < max_frags < num_children:
        #print("Skipping", smiles, "with edge count of", size[1])
        write_reject(smiles)
        rejects_count += 1
        return
    # print("Handling mol {0} with {1} nodes and {2} edges".format(smiles, size[0], size[1]))
    need_further_processing = write_data(node_holder, time_ms, num_children, size[1])
    reprocess_count = len(need_further_processing)
    for smiles in need_further_processing:
        # print("Recursing for ", child.smiles)
        fragment_and_write(smiles, max_frags=max_frags, verbosity=verbosity)
    node_holder = None

def fragment_mol(smiles, verbosity=0):

    attrs = []
    attr = Attr(smiles, ["EM"])
    attrs.append(attr)

    # Build the network
    node_holder = NodeHolder(iso_flag=False)
    node_holder = build_network(attrs, node_holder, base_dir=None, verbosity=verbosity, recurse=False)
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
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise')
    parser.add_argument('-s', '--skip',
                        type=int, default=0,
                        help='Number of molecules to skip molecules'
                             ' in the input file')
    parser.add_argument('--max-frag',
                        type=int, default=0,
                        help='Limit processing to molecules with no more than'
                             ' this number of initial fragment (no limit if 0)')
    parser.add_argument('-r', '--report-interval', type=int, default=1000, help='Reporting interval')


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

    global base_dir
    global nodes_f
    global edges_f
    global rejects_f
    base_dir = args.base_dir

    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)


    nodes_f = open(os.path.join(base_dir, "nodes.csv"), "w")
    edges_f = open(os.path.join(base_dir, "edges.csv"), "w")

    with open(args.input, 'r') as standard_file:

        # Process the rest of the file...
        num_skipped = 0
        num_processed = 0

        for line in standard_file:

            # Do we need to skip molecules before processing?
            if num_skipped < args.skip:
                num_skipped += 1
                continue

            line = line.strip()
            if line:
                fragment_and_write(line, max_frags=args.max_frag, verbosity=args.verbosity)

            # Enough?
            num_processed += 1
            if num_processed % args.report_interval == 0:
                print("Processed mol", num_processed)
            if args.limit and num_processed >= args.limit:
                break

    nodes_f.close()
    edges_f.close()
    if rejects_f:
        rejects_f.close()

    print("Processed {0} molecules, wrote {1} nodes and {2} edges, {3} rejects".format(num_processed, node_count, edge_count, rejects_count))


if __name__ == "__main__":
    main()
