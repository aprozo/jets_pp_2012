#!/bin/bash
# source this file to set up unfolding environment
export ROOUNFOLD_HOME="/home/prozorov/install/RooUnfold"
# 1) Headers for ACLiC
export ROOT_INCLUDE_PATH="${ROOUNFOLD_HOME}/include:${ROOT_INCLUDE_PATH}"

# 2) Runtime linker path for the .so
export LD_LIBRARY_PATH="${ROOUNFOLD_HOME}/lib:${ROOTSYS}/lib:${LD_LIBRARY_PATH}"
 # keep conda OFF here

TRIGGERS=(JP2 HT2)
RADII=(0.2 0.3 0.4 0.5 0.6)

# TRIGGERS=(JP2)
# RADII=(0.6)

for TRIGGER in "${TRIGGERS[@]}"; do
  for R in "${RADII[@]}"; do
    echo "Unfolding for trigger: $TRIGGER and R=$R"
    filename=/home/prozorov/dev/star/jets_pp_2012/output/merged_matching_"$TRIGGER"_R"$R.root"
    root -l -b -q \
  -e 'gSystem->Load("libRooUnfold"); gSystem->AddIncludePath(Form("-I%s/include", getenv("ROOUNFOLD_HOME")));' \
  "unfold.cxx+(\"$filename\")"
  done
done