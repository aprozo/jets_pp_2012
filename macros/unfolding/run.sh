#!/bin/bash
# source this file to set up unfolding environment
export ROOUNFOLD_HOME="/home/prozorov/install/RooUnfold"
export LD_LIBRARY_PATH="$ROOUNFOLD_HOME/lib:$ROOTSYS/lib"   # keep conda OFF here

radii=(0.2 0.3 0.4 0.5 0.6)


for r in "${radii[@]}"; do
    filename=/home/prozorov/dev/star/jets_pp_2012/output/jets_embedding_R"$r.root"
    echo "Processing file: $filename"
    root -l -b -q "unfold.cxx+(\"$filename\")"
done