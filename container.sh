#!/bin/bash
#source runimage.sh
source setup.sh
#make
./bin/RunppTestAna -i /gpfs01/star/pwg/prozorov/jets_pp_2012/ppRun12Datapicos/ppJP2Run12/sum0.root -intype pico -c JetTree -trig ppJP2 -o /gpfs01/star/pwg/prozorov/jets_pp_2012/output/out/out_$1 -N -1 -pj 10 40 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -hadcorr 0.9999999 -geantnum 0
