#!/bin/bash

for f in $@
do
  time python -m frag.network.scripts.standardise_chemspace_bb_compounds_3 -i $f -s 2 &
  echo "Loading file $f"
done
