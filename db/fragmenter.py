import argparse
from db.fairmolecules import MoleculeLoader

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network

loader = MoleculeLoader()

def run_many(limit):

    session = loader.create_session()
    mols = read_mols(session, limit)
    session.commit()
    inserted_nonisomol_total = 0
    inserted_edge_total = 0
    print("Read {0} mols".format(len(mols)))
    count = 0
    for mol in mols:
        count += 1
        print("Handling SMILES {0} {1}".format(count, mol.smiles))
        node_holder = fragment_mol(mol)
        inserted_nonisomol_count, inserted_edge_count = loader.insert_frags(session, mol, node_holder)
        inserted_nonisomol_total += inserted_nonisomol_count
        inserted_edge_total += inserted_edge_count
        session.commit()
    print("Total counts: smiles={0} edges={1}".format(inserted_nonisomol_total, inserted_edge_total))

def run_single(smiles):

    session = loader.create_session()
    mol = loader.read_smiles_for_fragmentation(session, smiles)
    session.commit()
    if mol:
        print("Read mol")
        node_holder = fragment_mol(mol)
        loader.insert_frags(session, mol, node_holder)
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

def read_mols(session, limit):
    return loader.read_mols_for_fragmentation(limit=limit, session=session)

def fragment_mol(isomol, verbosity=0):

    attrs = []

    attr = Attr(isomol.smiles, ["EM", isomol.id])
    attrs.append(attr)
    # Build the network
    node_holder = NodeHolder(iso_flag=False)
    max_frags = 0
    node_holder = build_network(attrs, node_holder,
                                max_frags, isomol.smiles, verbosity=verbosity, recurse=True)
    # Write the data out
    print(str(node_holder.size()))
    # for node in node_holder.node_list:
    #     print(str(node))
    # for edge in node_holder.get_edges():
    #         print(str(edge))

    return node_holder


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Fragmenter')
    parser.add_argument('--limit', '-l', type=int, default=10, help='Number of molecules to be processed.')
    parser.add_argument('--smiles', '-s', help='SMILES to process (testing only)')

    args = parser.parse_args()


    if args.smiles:
        run_single(args.smiles)
    else:
        run_many(args.limit)