#!/bin/bash
source /usr/local/root/bin/thisroot.sh
export FASTJETDIR=/usr/local/fastjet
#export PYTHIA8DIR=/gpfs01/star/pwg/elayavalli/pythia8303
export STARPICOPATH=/usr/local/eventStructuredAu
export JETREADER=/usr/local/jetreader_build
export LD_LIBRARY_PATH=/usr/local/jetreader_build/lib:/lib/:/usr/local/eventStructuredAu:/usr/local/RooUnfold:/usr/local/fastjet/lib:/usr/local/root/lib::/.singularity.d/libs

input_file=${1}
# if input_file is not provided, just run macro without arguments
if [ -z "$input_file" ]; then
    root -l ./macros/matching_mc_reco.cxx+
else
    root -l "./macros/matching_mc_reco.cxx+(\"${input_file}\")"
fi  
