#!/bin/bash

set -e

# Always operate relative to this script's directory so it works from anywhere.
working_dir="$(cd "$(dirname "$0")" && pwd)"


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


# echo "matching..."
# cd $working_dir/matching
# root -l -b -q matching.cpp
# echo "Matching step complete."
# echo ""

# echo "unfolding..."
# cd $working_dir/unfolding

# root -l -b -q \
#   -e 'gSystem->Load("libRooUnfold"); gSystem->AddIncludePath(Form("-I%s/include", getenv("ROOUNFOLD_HOME")));' \
#   unfold.cxx+
# echo "Unfolding step complete."
# echo ""


# Run the cross-section macro

cd $working_dir/cross_section
root -l -b -q \
  -e 'gSystem->Load("libRooUnfold"); gSystem->AddIncludePath(Form("-I%s/include", getenv("ROOUNFOLD_HOME")));' \
cross_section.cpp+


echo ""
echo "=========================================="
echo "Analysis complete!"
echo "=========================================="
echo ""

