#!/bin/bash

# HT2 JP2 MB
# geant_HT2 geant_JP2 
# mc_HT2 mc_JP2

is_embedding=true

trigger=JP2

filename=""

if [ $is_embedding == false ] && [ $trigger == HT2 ]; then
    filename=$PWD/data/picoDst/ppHT2Run12/split/pp12Pico_pass4_hadded_part0.root

elif [  $is_embedding == false ] && [ $trigger == JP2 ]; then
    filename=$PWD/data/picoDst/ppJP2Run12/split/sum9_part1.root

elif [ $is_embedding == true ]; then
    # filename=/gpfs01/star/pwg/youqi/run12/embedding/P12id/picos/20235003/out/pt-hat55999_027.root
    filename=/gpfs01/star/pwg/prozorov/TStarJetPicoMaker/submit/production/pt_hat911_250.root
else
    echo "Unknown trigger: $trigger"
    exit 1
fi



if [ $is_embedding == true ]; then
    echo "Running embedding"
    singularity exec -e -B /gpfs01 star_star.simg bash ./container.sh $filename mc_$trigger
    singularity exec -e -B /gpfs01 star_star.simg bash ./container.sh $filename geant_$trigger
# find mc_file from current directory
    mc_file=$(ls . | grep "mc_${trigger}.*\.root" | head -n 1)
    echo "Matching MC file: $mc_file"
    if [ -z "$mc_file" ]; then
        echo "MC file not found for matching!"
        exit 1
    fi
    singularity exec -e -B /gpfs01 star_star.simg bash ./container_matching_mc_reco.sh $mc_file true

else
    echo "Running data"
    singularity exec -e -B /gpfs01 star_star.simg bash ./container.sh $filename $trigger
fi
