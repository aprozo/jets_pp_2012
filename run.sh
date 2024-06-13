#!/bin/bash
singularity exec -e \
    -B /gpfs01 -B /star\
    star_star.simg \
    bash /gpfs01/star/pwg/prozorov/jets_pp_2012/container.sh ./ppRun12Datapicos/ppMBRun12/split/sum0_part4.root
