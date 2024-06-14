#!/bin/bash
# singularity exec -e \
#     -B /gpfs01 -B /star star_star.simg \
#     bash /gpfs01/star/pwg/prozorov/jets_pp_2012/container.sh ./ppRun12Datapicos/ppMBRun12/split/sum0_part4.root MB

singularity exec -e \
    -B /gpfs01 -B /star star_star.simg \
    bash /gpfs01/star/pwg/prozorov/jets_pp_2012/container.sh /gpfs01/star/pwg/robotmon/ppRun12_analysis_code/embedding/P12id/picos/pp2012embed_2023/out/pt-hat911_71.root mc
