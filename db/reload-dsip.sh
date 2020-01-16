#!/bin/bash

python -m db.fairmolecules -d -c --source-name chembl --source-version 25
python -m db.fairmolecules --source-name chemspace-bb --source-version 2019-07
python -m db.fairmolecules --source-name chemspace-bb --source-version 2019-12 --source-currency USD
python -m db.fairmolecules --source-name molport --source-version 2019-12 --source-currency USD
python -m db.fairmolecules --source-name molport --source-version 2020-01 --source-currency USD
python -m db.fairmolecules --source-name dsip --source-version 1

time python -m frag.network.scripts.standardise_xchem_compounds_3 -s 6 -i ~/data/xchem/dsip.txt.gz