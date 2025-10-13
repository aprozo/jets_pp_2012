#!/bin/bash

# HT2 JP2 MB
# geant_HT2 geant_JP2 
# mc_HT2 mc_JP2

mode=geant_JP2 

filename=""

if [ $mode == HT2 ]; then
    filename=$PWD/data/picoDst/ppHT2Run12/split/pp12Pico_pass4_hadded_part0.root

elif [ $mode == JP2 ]; then
    filename=$PWD/data/picoDst/ppJP2Run12/split/sum9_part1.root

elif [ $mode == geant_HT2 ] || [ $mode == mc_HT2 ] || [ $mode == geant_JP2 ] || [ $mode == mc_JP2 ]; then
    filename=/gpfs01/star/pwg/youqi/run12/embedding/P12id/picos/20235003/out/pt-hat55999_027.root
    #filename=/gpfs01/star/pwg/prozorov/TStarJetPicoMaker/submit/production/pt_hat55999_027.root
else
    echo "Unknown mode: $mode"
    exit 1
fi

# test 
singularity exec -e -B /gpfs01 star_star.simg bash ./container.sh $filename $mode

