#!/bin/bash

echo "Using DATA_DIR of $DATA_DIR"
echo "Using RESULTS_DIR of $RESULTS_DIR"

find $RESULTS_DIR/tree -name chemspace_bb.yml -exec rm {} \;
rm -rf $RESULTS_DIR/datasets/chemspace*
python -m frag.network.scripts.standardise_chemspace_bb_compounds_2 --input $DATA_DIR/Jul2019_CS_BB_MOD_selected_set_SMILES.txt.gz --output $RESULTS_DIR $@
