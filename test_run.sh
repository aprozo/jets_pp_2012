#!/bin/bash
# singularity exec -e \
#     -B /gpfs01 -B /star star_star.simg \
#     bash /gpfs01/star/pwg/prozorov/jets_pp_2012/container.sh ./ppRun12Datapicos/ppMBRun12/split/sum0_part4.root MB

# ./bin/RunppAna -i  ppRun12Datapicos/test_embedding_pt-hat911_71.root -intype mcpico -c JetTreeMc -trig all -o test_embedding_pt-hat911_71.root -N 100 -pj 1 2000 -pc 0.2 1000 -lja antikt -ec 1 -R 0.4

# ./bin/RunppAna -i /gpfs01/star/pwg/prozorov/jets_pp_2012/ppRun12Datapicos/ppJP2Run12/sum0.root -intype pico -c JetTree -trig ppJP2 -o output/test.root -N 100000 -pj 1 2000 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -hadcorr 0.9999999

# ./bin/RunppAna -i /gpfs01/star/pwg/elayavalli/ppRun12Datapicos/ppRun12Embpicos/Cleanpp12Pico_pt11_15_g4.root -intype pico -c JetTree -trig ppJP2 -o Results/test_geant_pt11_15_g4.root -N 100 -pj 1 2000 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -hadcorr 0.9999999

# ./bin/RunppAna -i /gpfs01/star/pwg/elayavalli/ppRun12Datapicos/ppRun12Embpicos/Cleanpp12Pico_pt11_15_g4.root -intype mcpico -c JetTreeMc -trig all -o Results/test_pythia_pt11_15_g4.root -N 100 -pj 1 2000 -pc 0.2 1000 -lja antikt -ec 1 -R 0.4
# singularity exec -e \
#     -B $PWD star_star.simg \
#     bash ./container.sh output/test/pt-hat2535_35.root geant

# singularity exec -e \
#     -B $PWD star_star.simg \
#     bash ./container.sh output/test/pt-hat2535_35.root geant

# cp tree_pt-hat2535_35.root tree_pt-hat2535_35_geant.root

#test embedding
# singularity exec -e \
#     -B $PWD star_star.simg \
#     bash gpfs01/star/pwg/youqi/run12/embedding/P12id/picos/20235003/out/pt-hat34_047.root geant

# singularity exec -e \
#     -B $PWD -B /gpfs01 star_star.simg \
#     bash ./container.sh /gpfs01/star/pwg/youqi/run12/embedding/P12id/picos/20235003/out/pt-hat34_047.root geant

singularity exec -e \
    -B $PWD star_star.simg \
    bash ./container_matching_mc_reco.sh

# test HT2 data
# singularity exec -e \
#     -B $PWD star_star.simg \
#     bash ./container.sh ppRun12Datapicos/ppHT2Run12/split/pp12Pico_pass4_hadded_part0.root HT2

# test JP2 data
# singularity exec -e \
#     -B $PWD star_star.simg \
#     bash ./container.sh ppRun12Datapicos/ppJP2Run12/split/sum0_part1.root JP2

# test MB data
# singularity exec -e \
#     -B $PWD star_star.simg \
#     bash ./container.sh ppRun12Datapicos/ppMBRun12/split/sum4_part24.root MB

# test matching in embedding
# singularity exec -e \
#     -B $PWD star_star.simg \
#     bash ./container_matching_mc_reco.sh tree_pt-hat1115_066.root matching_test_tree.root
