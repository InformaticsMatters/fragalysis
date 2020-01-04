import os, sys, argparse
from sqlalchemy import Column, ForeignKey, UniqueConstraint, Integer, String, Text, SmallInteger
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy import create_engine

from frag.utils import standardise_utils

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
    inchi_id = Column(Integer, ForeignKey('inchi.id'), nullable=False)
    inchi = relationship(Inchi)

class Isomol(Base):
    __tablename__ = 'isomol'
    id = Column(Integer, primary_key=True)

    smiles = Column(Text(), nullable=False, unique=True)
    inchik = Column(Text())
    inchis = Column(Text())
    noniso_id = Column(Integer, ForeignKey('nonisomol.id'), nullable=False)
    noniso = relationship(NonIsomol)

class Source(Base):
    __tablename__ = 'source'
    id = Column(Integer, primary_key=True)

    name = Column(Text(), nullable=False)
    version = Column(Text(), nullable=False)
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
        #else:
            #print("Found existing InChi")

        return i

    def insert_noniso(self, session, std_info, inchi):
        # print("Handling noniso SMILES")

        noniso = session.query(NonIsomol).filter(NonIsomol.smiles == std_info.noniso).first()
        if noniso == None:
            noniso = NonIsomol(smiles=std_info.noniso,
                               inchis=std_info.noniso_inchis, inchik=std_info.noniso_inchik, inchi=inchi,
                               hac=int(std_info.hac), rac=int(std_info.rac))
            session.add(noniso)
            #print("Inserted NonIso SMILES")
        #else:
            #print("Found existing NonIso SMILES")

        return noniso

    def insert_iso(self, session, std_info, noniso):
        # print("Handling iso SMILES")

        iso = session.query(Isomol).filter(Isomol.smiles == std_info.iso).first()
        if iso == None:
            iso = Isomol(smiles=std_info.iso, inchis=std_info.iso_inchis, inchik=std_info.iso_inchik, noniso=noniso)
            session.add(iso)
            #print("Inserted Iso SMILES")
        #else:
            #print("Found existing Iso SMILES")

        return iso

    def insert_source_mol(self, session, osmiles, source_id, source_code, nonisomol, isomol):
        ms = MolSource(smiles=osmiles, code=source_code, source_id=source_id, nonisomol=nonisomol, isomol=isomol)
        session.add(ms)
        #print("Inserted MolSource")

    def insert_mol(self, osmiles, source_id, source_code, std_info=None, session=None):
        if not std_info:
            std_info = standardise_utils.standardise(osmiles)
        if not session:
            session = self.create_session()
        inchi = self.insert_inchi(session, std_info.inchis, std_info.inchik)
        noniso = self.insert_noniso(session, std_info, inchi)
        isomol = None
        nonisomol = None
        if std_info.noniso == std_info.iso:
            nonisomol = noniso
        else:
            isomol = self.insert_iso(session, std_info, noniso)

        self.insert_source_mol(session, osmiles, source_id, source_code, nonisomol, isomol)

        session.commit()


    def insert_source(self, name, version):
        session = self.create_session()
        s = Source(name=name, version=version)
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
        loader.insert_source(args.source_name, args.source_version)

    if args.list_sources:
        sources = loader.list_sources()
        for source in sources:
            print("Source id={0} name={1} version={2}".format(source.id, source.name, source.version))

    if args.smiles and args.source_id and args.source_code:
        loader.insert_mol(args.smiles, args.source_id, args.source_code)


if __name__ == "__main__":
    main()