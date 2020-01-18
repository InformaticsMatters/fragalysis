import argparse, time
from db.fairmolecules import MoleculeLoader

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network

from rdkit import RDLogger

# Suppress basic RDKit logging...
RDLogger.logger().setLevel(RDLogger.ERROR)

loader = MoleculeLoader()

def run_many(limit, source_id, hac_max, dump_cache):

    reporting_interval = 500

    session = loader.create_session()
    mols = read_mols(session, limit, hac_max, source_id)
    session.commit()
    inserted_nonisomol_total = 0
    inserted_edge_total = 0
    reprocess_total = 0
    print("Read {0} mols".format(len(mols)))
    t0 = time.time()
    t1 = t0
    count = 0
    for mol in mols:
        count += 1
        # print("Handling SMILES {0} {1}".format(count, mol.smiles))
        inserted_nonisomol_count, inserted_edge_count, reprocess_count = fragment_and_load_mol(session, mol, count)
        inserted_nonisomol_total += inserted_nonisomol_count
        inserted_edge_total += inserted_edge_count
        reprocess_total += reprocess_count

        if count % reporting_interval == 0:
            t2 = time.time()
            dur = int(t2 - t1)
            t1 = t2
            cache_size, cache_found, cache_miss_not_exists, cache_miss_exists = loader.get_cache_stats()
            print("{0} molecules processed in {1}s. Cache now has {2} entries, {3} hits, {4} not in db, {5} in db"\
                  .format(count, dur, cache_size, cache_found, cache_miss_not_exists, cache_miss_exists))

    if count % reporting_interval != 0:
        t2 = time.time()
        dur = int(t2 - t1)
        t1 = t2
        cache_size, cache_found, cache_miss_not_exists, cache_miss_exists = loader.get_cache_stats()
        print("{0} mols in {1}s. Cache has {2} entries, {3} hits, {4} not in db, {5} in db" \
              .format(count, dur, cache_size, cache_found, cache_miss_not_exists, cache_miss_exists))

    print("Total counts: smiles={0} edges={1} reprocess={2}".format(inserted_nonisomol_total, inserted_edge_total, reprocess_total))

    if dump_cache:
        filename = "/tmp/fragmenter-{}.txt".format(loader.frag_id)
        print("Writing cache to {}".format(filename))
        loader.dump_cache(filename)
        print("Finished writing cache")


def fragment_and_load_mol(session, mol, count):

    t0 = int(round(time.time() * 1000))
    node_holder = fragment_mol(mol)
    t1 = int(round(time.time() * 1000))
    size = node_holder.size()
    # print("Handling mol {0} {1} with {2} nodes and {3} edges".format(count, mol.smiles, size[0], size[1]))
    added_child_nonisos, inserted_nonisomol_count, inserted_edge_count = loader.insert_frags(session, node_holder, t1 - t0)
    reprocess_count = len(added_child_nonisos)
    for child in added_child_nonisos:
        # print("Recursing for ", child.smiles)
        n, e, r = fragment_and_load_mol(session, child, count)
        inserted_nonisomol_count += n
        inserted_nonisomol_count += e
        reprocess_count += r
    return inserted_nonisomol_count, inserted_edge_count, reprocess_count

def run_single(smiles):

    session = loader.create_session()
    mol = loader.read_smiles_for_fragmentation(session, smiles)
    session.commit()
    if mol:
        print("Read mol")
        fragment_and_load_mol(session, mol)
    else:
        print("SMILES not found")
    session.commit()

def load_mol(mol, session=None):
    if not session:
        session = loader.create_session()
    node_holder = fragment_mol(mol)
    inserted_nonisomol_count, inserted_edge_count = loader.insert_frags(session, mol, node_holder)
    session.commit()
    return inserted_nonisomol_count, inserted_edge_count

def read_mols(session, limit, hac_max, source_id):
    return loader.read_mols_for_fragmentation(session, limit=limit, hac_max=hac_max, frag_status=None, source_id=source_id)

def fragment_mol(isomol, verbosity=0):

    attrs = []

    attr = Attr(isomol.smiles, ["EM", isomol.id])
    attrs.append(attr)
    # Build the network
    node_holder = NodeHolder(iso_flag=False)
    max_frags = 0
    node_holder = build_network(attrs, node_holder,
                                max_frags, isomol.smiles, verbosity=verbosity, recurse=False)
    # Write the data out
    # print(str(node_holder.size()))
    # for node in node_holder.node_list:
    #     print(str(node))
    # for edge in node_holder.get_edges():
    #         print(str(edge))

    return node_holder


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Fragmenter')
    parser.add_argument('--limit', '-l', type=int, default=10, help='Number of molecules to be processed.')
    parser.add_argument('--hac-max', type=int, default=36, help='Max HAC to be processed.')
    parser.add_argument('--source-id', '-i', type=int, help='Source ID to be processed.')
    parser.add_argument('--smiles', '-s', help='SMILES to process (testing only)')
    parser.add_argument('--dump-cache', '-d', action='store_true', help='Dump cache details at end of run')

    args = parser.parse_args()


    if args.smiles:
        run_single(args.smiles)
    else:
        run_many(args.limit, args.source_id, args.hac_max, args.dump_cache)