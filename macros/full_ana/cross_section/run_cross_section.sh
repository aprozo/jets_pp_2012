#!/bin/bash

# Wrapper script to run cross-section analysis with conda + RooUnfold
# Usage: ./run_cross_section.sh

set -e

echo "=========================================="
echo "Cross-Section Analysis"
echo "=========================================="
echo ""

working_dir=/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/cross_section
cd $working_dir

# Check if physics conda environment exists
if ! conda env list | grep -q "physics"; then
    echo "Error: 'physics' conda environment not found"
    echo "Please create it first: conda create -n physics root -c conda-forge"
    exit 1
fi

# Check if RooUnfold exists
if [ ! -d "/home/prozorov/install/RooUnfold" ]; then
    echo "Error: RooUnfold not found at /home/prozorov/install/RooUnfold"
    exit 1
fi

echo "Setting up environment..."
echo "  - Activating conda physics environment"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate physics

echo "  - Setting up RooUnfold"
export ROOUNFOLD_HOME="/home/prozorov/install/RooUnfold"
export ROOT_INCLUDE_PATH="${ROOUNFOLD_HOME}/include:${ROOT_INCLUDE_PATH}"
export LD_LIBRARY_PATH="${ROOUNFOLD_HOME}/lib:${ROOTSYS}/lib:${LD_LIBRARY_PATH}"

echo ""
echo "Environment setup complete:"
echo "  ROOT version: $(root-config --version)"
echo "  ROOUNFOLD_HOME: $ROOUNFOLD_HOME"
echo ""

# run matching first
echo "Running matching step..."
cd $working_dir/../matching
root -l -b -q matching.cpp
echo "Matching step complete."
echo ""

#run unfolding 
echo "Running unfolding step..."
cd $working_dir/../unfolding
./run.sh
echo "Unfolding step complete."
echo ""

cd $working_dir

# Run the cross-section macro
root -l -b -q \
  -e 'gSystem->Load("libRooUnfold"); gSystem->AddIncludePath(Form("-I%s/include", getenv("ROOUNFOLD_HOME")));' \
cross_section.cpp+
echo ""
echo "=========================================="
echo "Analysis complete!"
echo "=========================================="
echo ""
