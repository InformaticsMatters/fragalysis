#!/usr/bin/env python

import os
from datetime import datetime
import argparse

# Issue a reminder to stdout after every N lines...
# (just to re-assure the user we're doing something on very large files)
# On a 2.7GHz i7 we should see around 13Million node lines/min
# and about 9Million edge lines/min.
# Set to 0 for quiet operation.
LINE_REMINDER = 13000000


def do_for_dir(input_dir):
    os.chdir(input_dir)
    prop_dict = {
        "nodes.txt": [None, "smiles:ID(F2)", "hac:INT", "chac:INT", "osmiles"],
        "edges.txt": [None, ":START_ID(F2)", ":END_ID(F2)", "label"],
    }
    # build a map of SMILES to chemical IDs
    # from the attributes file...
    num_supplier_identities_found = 0
    f_name = "attributes.txt"
    if LINE_REMINDER:
        print('%s Processing %s/%s...' % (datetime.now(), input_dir, f_name))
    attrs = {}
    with open(f_name) as attr_f:
        for line in attr_f:

            line_parts = line.split()
            # Translate supplier names (in line_parts[3]).
            #  'Z' -> 'REAL:Z'
            #  'MolPort-' -> 'MolPort:'
            # Ignore others
            supplier = None
            if line_parts[3].startswith('Z'):
                supplier = 'REAL:' + line_parts[3]
            elif line_parts[3].startswith('MolPort-'):
                supplier = 'MOLPORT:' + line_parts[3].split('MolPort-', 1)[1]

            if supplier:
                num_supplier_identities_found += 1
                # Each attribute is a potential list of suppliers
                # for a given molecule (line_parts[1]).
                # We separate suppliers using ';'.
                if line_parts[1] in attrs:
                    attrs[line_parts[1]] += ';' + supplier
                else:
                    attrs[line_parts[1]] = supplier
    if LINE_REMINDER:
        if not num_supplier_identities_found:
            print('%s WARNING: Found no supplier identities (!)' %
                  (datetime.now(), num_supplier_identities_found))
        else:
            print('%s Found %d supplier identities' %
                  (datetime.now(), num_supplier_identities_found))

    for f_name in prop_dict:

        if LINE_REMINDER:
            print('%s Processing %s/%s...' %
                  (datetime.now(), input_dir, f_name))
        line_num = 0

        out_f = open(f_name.replace(".txt", ".csv"), "w")
        with open(f_name) as in_f:
            for line in in_f:

                if LINE_REMINDER:
                    line_num += 1
                    if line_num % LINE_REMINDER == 0:
                        print('%s ...%d' % (datetime.now(), line_num))

                line_spl = line.split()
                out_l = []
                for i, x in enumerate(prop_dict[f_name]):
                    if x:
                        out_l.append(line_spl[i])
                if f_name == "nodes.txt":
                    attr_cmp_id = attrs.get(line_spl[1], None)
                    if attr_cmp_id:
                        cmpd_id = attr_cmp_id
                        out_l.append(cmpd_id)
                        # This is where we can add tags - like CHEAP - EXPENSIVE
                        out_l.append("MOL;F2")
                    else:
                        out_l.append("")
                        out_l.append("F2")
                    out_f.write(",".join(out_l) + "\n")
                elif f_name == "edges.txt":
                    out_f.write(",".join(out_l) + "\n")

        if LINE_REMINDER:
            print('%s Processed (%d)' % (datetime.now(), line_num))
        out_f.flush()
        out_f.close()


def main():
    parser = argparse.ArgumentParser(
        description="Convert text files output from build_db into CSV files for neo4j."
    )
    parser.add_argument("--input_dir", required=True)
    args = parser.parse_args()

    node_list = ["nodes-header.csv"]
    edge_list = ["edges-header.csv"]
    with open("edges-header.csv", "w") as out_f:
        out_f.write(",".join([x for x in [":START_ID(F2)", ":END_ID(F2)", "label"]]))
    with open("nodes-header.csv", "w") as out_f:
        out_f.write(
            ",".join(
                [
                    x
                    for x in [
                        "smiles:ID(F2)",
                        "hac:INT",
                        "chac:INT",
                        "osmiles",
                        "cmpd_id",
                        ":LABEL",
                    ]
                    if x
                ]
            )
            + "\n"
        )
    do_for_dir(args.input_dir)
    node_list.append("nodes.csv")
    edge_list.append("edges.csv")
    # Now this is the command to run
    print(
        " ".join(
            [
                "/var/lib/neo4j/bin/neo4j-admin import",
                "--database",
                "new.db",
                "--nodes",
                '"' + ",".join(node_list) + '"',
                "--relationships:F2EDGE",
                '"' + ",".join(edge_list) + '"',
            ]
        )
    )


if __name__ == "__main__":
    main()
