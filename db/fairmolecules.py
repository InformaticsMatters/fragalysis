import os, sys, argparse
from sqlalchemy import Column, PrimaryKeyConstraint, ForeignKey, UniqueConstraint, Integer, String, Text, SmallInteger
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy import create_engine

from frag.utils import standardise_utils

from rdkit import Chem

Base = declarative_base()

class Inchi(Base):
    __tablename__ = 'inchi'
    id = Column(Integer, primary_key=True)

    inchik = Column(Text(), nullable=False, index=True)
    inchis = Column(Text(), nullable=False)

class NonIsomol(Base):
    __tablename__ = 'nonisomol'
    id = Column(Integer, primary_key=True)

    smiles = Column(Text(), nullable=False, unique=True)
    inchik = Column(Text())
    inchis = Column(Text())
    hac = Column(SmallInteger())
    rac = Column(SmallInteger())
    fc = Column(SmallInteger())
    fs = Column(SmallInteger(), index=True)
    inchi_id = Column(Integer, ForeignKey('inchi.id'), nullable=False)
    inchi = relationship(Inchi)

class Isomol(Base):
    __tablename__ = 'isomol'
    id = Column(Integer, primary_key=True)

    smiles = Column(Text(), nullable=False, unique=True)
    inchik = Column(Text())
    inchis = Column(Text())
    nonisomol_id = Column(Integer, ForeignKey('nonisomol.id'), nullable=False)
    nonisomol = relationship(NonIsomol)

class Source(Base):
    __tablename__ = 'source'
    id = Column(Integer, primary_key=True)

    name = Column(Text(), nullable=False)
    version = Column(Text(), nullable=False)
    currency = Column(Text())
    UniqueConstraint('name', 'version', name='uq_source')

class MolSource(Base):
    __tablename__ = 'mol_source'
    id = Column(Integer, primary_key=True)

    smiles = Column(Text(), nullable=False)
    code = Column(Text(), nullable=False)
    source_id = Column(Integer, ForeignKey('source.id', ondelete='CASCADE'), nullable=False)
    nonisomol_id = Column(Integer, ForeignKey('nonisomol.id'))
    isomol_id = Column(Integer, ForeignKey('isomol.id'))
    source = relationship(Source)
    nonisomol = relationship(NonIsomol)
    isomol = relationship(Isomol)

class Price(Base):
    __tablename__ = 'price'
    id = Column(Integer, primary_key=True)
    quantity_mg = Column(Integer)
    price = Column(Integer)
    price_min = Column(Integer)
    price_max = Column(Integer)
    molsource_id = Column(Integer, ForeignKey(MolSource.id))
    molsource = relationship(MolSource)

class Edge(Base):
    __tablename__ = 'edge'
    id = Column(Integer, primary_key=True)
    parent_id = Column(Integer, ForeignKey(NonIsomol.id))
    child_id = Column(Integer, ForeignKey(NonIsomol.id))
    label = Column(Text(), nullable=False)
    parent = relationship(NonIsomol, foreign_keys=[parent_id])
    child = relationship(NonIsomol, foreign_keys=[child_id])

class MoleculeLoader:
    """
    This is the loader class.
    The database connection URL need to be defined, either by passing in the url parameter or using the FM_DB_URL
    environment variable. The URL will be something like this: postgres://user:password@localhost:5432/database
    """

    def __init__(self, db_url=None):
        if db_url:
            self.url = db_url
        else:
            self.url = os.environ.get('FM_DB_URL')
        if self.url == None:
            raise('ERROR: Must define database URL using the db_url parameter or the FM_DB_URL environment variable')

        self.engine = create_engine(self.url)
        self.DBSession = sessionmaker(bind=self.engine)
        self.frag_cache = {}

    def gen_std_info(self, nonisosmiles):
        mol = Chem.MolFromSmiles(nonisosmiles)
        hac = mol.GetNumHeavyAtoms()
        rac = 0
        for atom in mol.GetAtoms():
            if atom.IsInRing():
                rac += 1
        inchis, inchik = standardise_utils.gen_inchi(mol, '')
        return inchis, inchik, hac, rac

    def create_session(self):
        return self.DBSession()

    def create_tables(self):
        # Create all tables in the engine.
        Base.metadata.create_all(self.engine)
        print("Tables created")

    def drop_tables(self):
        Base.metadata.drop_all(self.engine)
        print("Tables dropped")

    def insert_inchi(self, session, inchis, inchik):
        # print("Handling InChi")

        i = session.query(Inchi).filter(Inchi.inchik == inchik, Inchi.inchis == inchis).first()
        if i == None:
            i = Inchi(inchis=inchis, inchik=inchik)
            session.add(i)
            #print("Inserted Inchi")
            return i, True
        else:
            #print("Found existing InChi")
            return i, False


    def insert_noniso(self, session, inchi, smiles, inchis, inchik, hac, rac, fs=None):
        # print("Handling noniso SMILES")

        noniso = session.query(NonIsomol).filter(NonIsomol.smiles == smiles).first()
        if noniso == None:
            noniso = NonIsomol(smiles=smiles,
                               inchis=inchis, inchik=inchik, inchi=inchi,
                               hac=hac, rac=rac, fs=fs)
            session.add(noniso)
            return noniso, True
            #print("Inserted NonIso SMILES")
        else:
            #print("Found existing NonIso SMILES")
            return noniso, False



    def insert_iso(self, session, std_info, noniso):
        # print("Handling iso SMILES")

        iso = session.query(Isomol).filter(Isomol.smiles == std_info.iso).first()
        if iso == None:
            iso = Isomol(smiles=std_info.iso, inchis=std_info.iso_inchis, inchik=std_info.iso_inchik, nonisomol=noniso)
            session.add(iso)
            #print("Inserted Iso SMILES")
        #else:
            #print("Found existing Iso SMILES")

        return iso

    def insert_source_mol(self, session, osmiles, source_id, source_code, nonisomol, isomol):
        ms = MolSource(smiles=osmiles, code=source_code, source_id=source_id, nonisomol=nonisomol, isomol=isomol)
        session.add(ms)
        #print("Inserted MolSource")
        return ms

    def insert_price(self, session, quantity, price, mol_source):
        p = Price(quantity_mg=quantity, price=price, molsource=mol_source)
        session.add(p)
        return p

    def insert_mol(self, osmiles, source_id, source_code, prices=None, std_info=None, session=None):
        if not std_info:
            std_info = standardise_utils.standardise(osmiles)
        if not session:
            session = self.create_session()

        try:
            inchi, inchi_added = self.insert_inchi(session, std_info.inchis, std_info.inchik)
            #inchi, smiles, inchis, inchik, hac, rac
            noniso, noniso_added = self.insert_noniso(session, inchi, std_info.noniso,
                                                      std_info.noniso_inchis, std_info.noniso_inchik,
                                                      std_info.hac, std_info.rac)
            isomol = None
            nonisomol = None
            if std_info.noniso == std_info.iso:
                nonisomol = noniso
            else:
                isomol = self.insert_iso(session, std_info, noniso)

            mol_source = self.insert_source_mol(session, osmiles, source_id, source_code, nonisomol, isomol)
            if prices:
                # TODO - handle price ranges as well (min max values)
                for q in prices:
                    self.insert_price(session, q, prices[q], mol_source)

            session.commit()

        except:
            print("Failed to handle molecule {0} {1}. inchi: {2} noniso: {3} iso: {4}".format(source_code, osmiles, std_info.noniso, std_info.iso))
            session.rollback()


    def insert_source(self, name, version, currency):
        session = self.create_session()
        s = Source(name=name, version=version, currency=currency)
        session.add(s)
        session.commit()
        print("Inserted source " + str(s.id))
        return s


    def list_sources(self):
        session = self.create_session()
        sources = session.query(Source).all()
        session.commit()
        return sources


    def delete_source(self, id):
        session = self.create_session()
        session.query(Source).filter(Source.id == int(id)).delete()
        session.commit()


    def read_mols_for_fragmentation(self, session, frag_status=None, limit=100):
        mols = session.query(NonIsomol).filter(NonIsomol.fs == frag_status).limit(limit).all()
        for mol in mols:
            mol.fs = 1
        return mols

    def read_smiles_for_fragmentation(self, session, smiles):
        mol = session.query(NonIsomol).filter(NonIsomol.smiles == smiles).first()
        if mol:
            mol.fs = 1
        return mol

    def insert_frags(self, session, nonisomol, node_holder):
        cache = self.frag_cache

        noniso_added = set()
        inserted_nonisomol_count = 0
        inserted_edge_count = 0
        for edge in node_holder.get_edges():
            p = edge.NODES[0]
            c = edge.NODES[1]
            p_smiles = p.SMILES
            c_smiles = c.SMILES
            label = edge.get_label()
            # print("Looking at edge {0} {1} {2}".format(p_smiles, c_smiles, label))
            p_noniso, p_added = self.find_or_insert_nonisofrag(session, p_smiles)
            c_noniso, c_added = self.find_or_insert_nonisofrag(session, c_smiles)
            # print("Parent/Child added: {0} {1}".format(p_added, c_added))
            if p_added:
                noniso_added.add(p_smiles)
                inserted_nonisomol_count += 1
            if c_added:
                noniso_added.add(c_smiles)
                inserted_nonisomol_count += 1

            if  p_smiles in noniso_added or c_smiles in noniso_added:
                # print("Edge added for {0} {1}".format(p_noniso.id, c_noniso.id))
                e = Edge(label=label, parent=p_noniso, child=c_noniso)
                session.add(e)
                inserted_edge_count += 1
            # else:
            #     print("Edges for {0} {1} already present".format(p_noniso.id, c_noniso.id))

        nonisomol.fs = 2
        print("Inserted {0} smiles and {1} edges".format(inserted_nonisomol_count, inserted_edge_count))
        print("Cache now has {0} entries".format(len(cache)))
        return inserted_nonisomol_count, inserted_edge_count

    def find_or_insert_nonisofrag(self, session, smiles):
        cache = self.frag_cache

        if smiles in cache:
            # print("Molecule {0} already present as {1}".format(smiles, cache[smiles].id))
            return cache[smiles], False
        else:
            # print("Looking to add Molecule {0}".format(smiles))
            inchis, inchik, hac, rac = self.gen_std_info(smiles)
            inchi, inchi_added = self.insert_inchi(session, inchis, inchik)
            # if inchi_added:
            #     print("Added InChi {}".format(inchis))
            # else:
            #     print("Found existing InChi {}".format(inchis))
            noniso, added = self.insert_noniso(session, inchi, smiles, inchis, inchik, hac, rac, fs=2)
            # if added:
            #     print("Added noniso {}".format(smiles))
            # else:
            #     print("Found existing noniso {}".format(smiles))
            cache[smiles] = noniso
            return noniso, added


def main():

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Fair molecules DB')
    parser.add_argument('--create-tables', '-c', action='store_true', help='Create tables')
    parser.add_argument('--drop-tables', '-d', action='store_true', help='Drop tables')
    parser.add_argument('--smiles', '-s', help='Molecule SMILES')
    parser.add_argument('--source-id', '-i', help='ID of Source (Primary key in source table)')
    parser.add_argument('--source-code', help="Source's code (vendor's ID for molecule)")
    parser.add_argument('--list-sources', action='store_true', help="List the currently defined sources")
    parser.add_argument('--source-name', help="Name for the source to be added")
    parser.add_argument('--source-version', help="Version for the source to be added")
    parser.add_argument('--source-currency', help="Currency (optional) for the source to be added")
    parser.add_argument('--delete-source', help="Delete the source with this ID")
    parser.add_argument('--db-url', help="URL of database to use")

    args = parser.parse_args()
    print("fairmolecules.py: ", args)

    loader = MoleculeLoader(args.db_url)

    if args.drop_tables:
        loader.drop_tables()

    if args.create_tables:
        loader.create_tables()

    if args.delete_source:
        loader.delete_source(args.delete_source)

    if args.source_name and args.source_version:
        loader.insert_source(args.source_name, args.source_version, args.source_currency)

    if args.list_sources:
        sources = loader.list_sources()
        for source in sources:
            print("Source id={0} name={1} version={2}".format(source.id, source.name, source.version))

    if args.smiles and args.source_id and args.source_code:
        loader.insert_mol(args.smiles, args.source_id, args.source_code)


if __name__ == "__main__":
    main()