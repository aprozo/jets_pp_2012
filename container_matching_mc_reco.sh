#!/bin/bash
source /usr/local/root/bin/thisroot.sh
export FASTJETDIR=/usr/local/fastjet
#export PYTHIA8DIR=/gpfs01/star/pwg/elayavalli/pythia8303
export STARPICOPATH=/usr/local/eventStructuredAu
export JETREADER=/usr/local/jetreader_build
export LD_LIBRARY_PATH=/usr/local/jetreader_build/lib:/lib/:/usr/local/eventStructuredAu:/usr/local/RooUnfold:/usr/local/fastjet/lib:/usr/local/root/lib::/.singularity.d/libs

input_file1=${1}
output_file=${2}

# if input_file1 and output_file are not provided, just run macro without arguments else run with provided arguments
if [ -z "$input_file1" ] && [ -z "$output_file" ]; then
    root -l ./macros/matching_mc_reco.cxx++
else
    root -l "./macros/matching_mc_reco.cxx++(\"${input_file1}\",\"${output_file}\")"
fi  
