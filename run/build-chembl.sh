#!/bin/bash

echo "Using DATA_DIR of $DATA_DIR"
echo "Using RESULTS_DIR of $RESULTS_DIR"

find $RESULTS_DIR/tree -name chembl.yml -exec rm {} \;
rm -rf $RESULTS_DIR/datasets/chembl*
python -m frag.network.scripts.standardise_chembl_compounds_2 --input $DATA_DIR/chembl_25.sdf.gz --output $RESULTS_DIR $@
