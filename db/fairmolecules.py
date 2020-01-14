import os, sys, argparse
from sqlalchemy import Column, PrimaryKeyConstraint, ForeignKey, Index, UniqueConstraint, Integer, String, Text, SmallInteger, DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy import create_engine
from sqlalchemy import func

from frag.utils import standardise_utils

from rdkit import Chem

Base = declarative_base()

class Inchi(Base):
    __tablename__ = 'inchi'
    id = Column(Integer, primary_key=True)

    # inchik = Column(Text(), nullable=False, index=True)
    inchik = Column(Text(), nullable=False)
    inchis = Column(Text(), nullable=False)
    Index('ix_inchik', inchik, postgresql_using='hash')

class NonIsomol(Base):
    __tablename__ = 'nonisomol'
    id = Column(Integer, primary_key=True)

    # smiles = Column(Text(), nullable=False, unique=True)
    smiles = Column(Text(), nullable=False)
    inchik = Column(Text())
    inchis = Column(Text())
    hac = Column(SmallInteger())
    rac = Column(SmallInteger())
    fc = Column(SmallInteger())
    fs = Column(SmallInteger(), index=True)
    inchi_id = Column(Integer, ForeignKey('inchi.id'), nullable=False)
    inchi = relationship(Inchi)
    Index('uq_noniso_smiles', smiles, postgresql_using='hash')

class Isomol(Base):
    __tablename__ = 'isomol'
    id = Column(Integer, primary_key=True)

    # smiles = Column(Text(), nullable=False, unique=True)
    smiles = Column(Text(), nullable=False)
    inchik = Column(Text())
    inchis = Column(Text())
    nonisomol_id = Column(Integer, ForeignKey('nonisomol.id'), nullable=False)
    nonisomol = relationship(NonIsomol)
    Index('uq_iso_smiles', smiles, postgresql_using='hash')

class Source(Base):
    __tablename__ = 'source'
    id = Column(Integer, primary_key=True)

    name = Column(Text(), nullable=False)
    version = Column(Text(), nullable=False)
    currency = Column(Text())

    __table_args__ = (
        UniqueConstraint(name, version, name='uq_source'),
    )

class MolInput(Base):
    __tablename__ = 'mol_input'
    id = Column(Integer, primary_key=True)
    name = Column(Text(), nullable=False)
    started_date = Column(DateTime())
    finished_date = Column(DateTime())
    size = Column(Integer())
    molecule_failures = Column(Integer())
    inchi_failures = Column(Integer())
    inchi_hits = Column(Integer())
    inchi_miss = Column(Integer())
    nonisomol_hits = Column(Integer())
    nonisomol_miss = Column(Integer())
    isomol_hits = Column(Integer())
    isomol_miss = Column(Integer())
    source_id = Column(Integer, ForeignKey(Source.id, ondelete='CASCADE'), nullable=False)
    source = relationship(Source)

    __table_args__ = (
        UniqueConstraint(source_id, name, name='uq_input'),
    )


class MolSource(Base):
    __tablename__ = 'mol_source'
    id = Column(Integer, primary_key=True)

    smiles = Column(Text(), nullable=False)
    code = Column(Text(), nullable=False)
    source_id = Column(Integer, ForeignKey(Source.id, ondelete='CASCADE'), nullable=False)
    input_id = Column(Integer, ForeignKey(MolInput.id, ondelete='CASCADE'), nullable=False)
    nonisomol_id = Column(Integer, ForeignKey('nonisomol.id'))
    isomol_id = Column(Integer, ForeignKey('isomol.id'))
    source = relationship(Source)
    input = relationship(MolInput)
    nonisomol = relationship(NonIsomol)
    isomol = relationship(Isomol)


class Price(Base):
    __tablename__ = 'price'
    id = Column(Integer, primary_key=True)
    quantity_mg = Column(Integer)
    price = Column(Integer)
    price_min = Column(Integer)
    price_max = Column(Integer)
    molsource_id = Column(Integer, ForeignKey(MolSource.id, ondelete='CASCADE'))
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
        if self.url is None:
            raise 'ERROR: Must define database URL using the db_url parameter or the FM_DB_URL environment variable'

        self.engine = create_engine(self.url, isolation_level="READ UNCOMMITTED")
        self.DBSession = sessionmaker(bind=self.engine)

        # loader stats
        self.count = 0
        self.inchi_hits = 0
        self.inchi_miss = 0
        self.noniso_hits = 0
        self.noniso_miss = 0
        self.iso_hits = 0
        self.iso_miss = 0

        # cache stats
        self.frag_cache = {}
        self.cache_found = 0
        self.cache_miss_exists = 0
        self.cache_miss_not_exists = 0
        self.already_fragmented = 0

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
        if i is None:
            self.inchi_miss += 1
            i = Inchi(inchis=inchis, inchik=inchik)
            session.add(i)
            #print("Inserted Inchi")
            return i, True
        else:
            self.inchi_hits += 1
            #print("Found existing InChi")
            return i, False


    def insert_noniso(self, session, inchi, smiles, inchis, inchik, hac, rac, fs=None):
        # print("Handling noniso SMILES")

        noniso = session.query(NonIsomol).filter(NonIsomol.smiles == smiles).first()
        if noniso is None:
            self.noniso_miss += 1
            noniso = NonIsomol(smiles=smiles,
                               inchis=inchis, inchik=inchik, inchi=inchi,
                               hac=hac, rac=rac, fs=fs)
            session.add(noniso)
            #print("Inserted NonIso SMILES")
            return noniso, True
        else:
            self.noniso_hits += 1
            #print("Found existing NonIso SMILES")
            return noniso, False



    def insert_iso(self, session, std_info, noniso):
        # print("Handling iso SMILES")

        iso = session.query(Isomol).filter(Isomol.smiles == std_info.iso).first()
        if iso is None:
            self.iso_miss += 1
            iso = Isomol(smiles=std_info.iso, inchis=std_info.iso_inchis, inchik=std_info.iso_inchik, nonisomol=noniso)
            session.add(iso)
            #print("Inserted Iso SMILES")
        else:
            self.iso_hits += 1
            #print("Found existing Iso SMILES")

        return iso

    def insert_source_mol(self, session, osmiles, source_code, nonisomol, isomol):
        ms = MolSource(input=self.mol_input, source_id=self.mol_input.source_id, smiles=osmiles, code=source_code, nonisomol=nonisomol, isomol=isomol)
        session.add(ms)
        #print("Inserted MolSource")
        return ms

    def insert_price(self, session, quantity, price, mol_source):
        """"
        price will either be an integer value or a tuple of size 2 with a min and max integer value
        """
        if type(price) is tuple:
            p = Price(quantity_mg=quantity, price_min=price[0], price_max=price[1], molsource=mol_source)
        else:
            p = Price(quantity_mg=quantity, price=price, molsource=mol_source)
        session.add(p)
        return p

    def retry(self, func, *argv, **kwargs):
        try:
            return func(*argv, **kwargs)
        except:
            ex = sys.exc_info()[0]
            print("Retrying function {0} with args {1} and {2}".format(str(func), str(argv), str(kwargs)))
            print("Error was ", str(ex))
            return func(*argv, **kwargs)

    def insert_mol(self, osmiles, source_code, prices=None, std_info=None, session=None):

        self.count += 1

        if not std_info:
            std_info = standardise_utils.standardise(osmiles)
        if not session:
            session = self.create_session()

        try:
            inchi, inchi_added = self.retry(self.insert_inchi, session, std_info.inchis, std_info.inchik)

            noniso, noniso_added = self.retry(self.insert_noniso, session, inchi, std_info.noniso,
                                                      std_info.noniso_inchis, std_info.noniso_inchik,
                                                      std_info.hac, std_info.rac)
            isomol = None
            nonisomol = None
            if std_info.noniso == std_info.iso:
                nonisomol = noniso
            else:
                isomol = self.retry(self.insert_iso, session, std_info, noniso)

            mol_source = self.insert_source_mol(session, osmiles, source_code, nonisomol, isomol)
            if prices:
                # TODO - handle price ranges as well (min max values)
                for q in prices:
                    self.insert_price(session, q, prices[q], mol_source)

            session.commit()

        except Exception as e:
            session.rollback()
            print("Failed to handle molecule {0} {1}. inchi: {2} noniso: {3} iso: {4}".format(
                source_code, osmiles, std_info.inchis, std_info.noniso, std_info.iso))
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)


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

    def create_input(self, session, name, source_id):
        inp = MolInput(name=name, source_id=source_id, started_date=func.now())
        session.add(inp)
        session.commit()
        self.mol_input = inp

    def complete_input(self, session, count, molecule_failures, inchi_failures):
        inp = self.mol_input
        inp.molecule_failures = molecule_failures
        inp.inchi_failures = inchi_failures
        inp.size = count
        inp.inchi_hits = self.inchi_hits
        inp.inchi_miss = self.inchi_miss
        inp.nonisomol_hits = self.noniso_hits
        inp.nonisomol_miss = self.noniso_miss
        inp.isomol_hits = self.iso_hits
        inp.isomol_miss = self.iso_miss
        inp.finished_date = func.now()
        session.commit()

    def read_mols_for_fragmentation(self, session, hac_max=36, frag_status=None, limit=100, source_id=None):
        if source_id:
            mols = session.query(NonIsomol).join(MolSource)\
                .join(Source)\
                .filter(NonIsomol.fs == frag_status)\
                .filter(NonIsomol.hac <= hac_max)\
                .filter(Source.id == source_id)\
                .limit(limit)\
                .all()
        else:
            mols = session.query(NonIsomol)\
                .filter(NonIsomol.fs == frag_status)\
                .filter(NonIsomol.hac <= hac_max)\
                .limit(limit)\
                .all()
        for mol in mols:
            mol.fs = 1
            self.frag_cache[mol.smiles] = mol
        return mols

    def read_smiles_for_fragmentation(self, session, smiles):
        mol = session.query(NonIsomol).filter(NonIsomol.smiles == smiles).first()
        if mol:
            mol.fs = 1
        return mol

    def insert_frags(self, session, nonisomol, node_holder):

        need_processing = set()
        inserted_nonisomol_count = 0
        inserted_edge_count = 0

        # if no edges then just set status to complete
        if node_holder.size()[1] == 0:
            if nonisomol.fs != 2:
                # print("Marking as complete", nonisomol.smiles)
                nonisomol.fs = 2
        else:
            # so we need to process the edges

            cache = self.frag_cache

            # print("Handling", nonisomol.smiles)

            noniso_added = set()
            parents_encountered = set()
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
                parents_encountered.add(p_noniso)
                if p_added:
                    # print("Parent added", p_smiles)
                    noniso_added.add(p_smiles)
                    inserted_nonisomol_count += 1
                # else:
                #     print("Parent exists, status is {0}".format(p_noniso.fs))
                if c_added:
                    # print("Child added", c_smiles)
                    noniso_added.add(c_smiles)
                    need_processing.add(c_noniso)
                    inserted_nonisomol_count += 1
                else:
                    if c_noniso.fs != 2:
                        need_processing.add(c_noniso)
                    # print("Child exists {0}, status is {1}".format(c_smiles, c_noniso.fs))

                if  p_noniso.fs != 2:
                    # print("Edge added for {0} {1}".format(p_noniso.id, c_noniso.id))
                    e = Edge(label=label, parent=p_noniso, child=c_noniso)
                    session.add(e)
                    inserted_edge_count += 1
                # else:
                #     print("Edges for {0} {1} already present".format(p_noniso.id, c_noniso.id))

            for p in parents_encountered:
                if p.fs != 2:
                    p.fs = 2

            nonisomol.fs = 2
            # print("Inserted {0} smiles and {1} edges".format(inserted_nonisomol_count, inserted_edge_count))
            # print("Cache now has {0} entries. {1} need processing".format(len(cache), len(need_processing)))

        return need_processing, inserted_nonisomol_count, inserted_edge_count

    def find_or_insert_nonisofrag(self, session, smiles):
        cache = self.frag_cache

        if smiles in cache:
            # print("Molecule {0} already present as {1}".format(smiles, cache[smiles].id))
            self.cache_found += 1
            return cache[smiles], False
        else:
            # print("Looking to add Molecule {0}".format(smiles))
            inchis, inchik, hac, rac = self.gen_std_info(smiles)
            inchi, inchi_added = self.insert_inchi(session, inchis, inchik)
            # if inchi_added:
            #     print("Added InChi {}".format(inchis))
            # else:
            #     print("Found existing InChi {}".format(inchis))
            noniso, added = self.insert_noniso(session, inchi, smiles, inchis, inchik, hac, rac, fs=1)
            if added:
                self.cache_miss_not_exists += 1
            #     print("Added noniso {}".format(smiles))
            else:
                self.cache_miss_exists += 1
            #     print("Found existing noniso {}".format(smiles))
            cache[smiles] = noniso
            return noniso, added

    def get_cache_stats(self):
        return len(self.frag_cache), self.cache_found, self.cache_miss_not_exists, self.cache_miss_exists


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