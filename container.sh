#!/bin/bash
source /usr/local/root/bin/thisroot.sh
export FASTJETDIR=/usr/local/fastjet
#export PYTHIA8DIR=/gpfs01/star/pwg/elayavalli/pythia8303
export STARPICOPATH=/usr/local/eventStructuredAu
export JETREADER=/usr/local/jetreader_build
export LD_LIBRARY_PATH=/usr/local/jetreader_build/lib:/lib/:/usr/local/eventStructuredAu:/usr/local/RooUnfold:/usr/local/fastjet/lib:/usr/local/root/lib::/.singularity.d/libs
input_file=${1}
output_file=$(basename $input_file)
data_type=${2}

make

treeType=JetTree
picoType=pico
trigger=pp
geantnum=0

if [[ $data_type = HT2 ]]; then
    trigger=ppHT
fi

if [[ $data_type = JP2 ]]; then
    trigger=ppJP
fi

if [[ $data_type = geant ]]; then
    geantnum=1
fi

if [[ $data_type = mc ]]; then
    treeType=JetTreeMc
    picoType=mcpico
    trigger=All
fi

./bin/RunppAna \
    -i $input_file \
    -intype $picoType \
    -c $treeType \
    -trig $trigger \
    -o tree_$output_file \
    -N -1 \
    -pj 0.001 50 \
    -pc 0.2 30 \
    -lja antikt \
    -ec 1 \
    -R 0.4 \
    -hadcorr 0.9999999 \
    -geantnum $geantnum
