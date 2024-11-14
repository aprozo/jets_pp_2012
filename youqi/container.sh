#!/bin/bash

cd /gpfs01/star/pwg/youqi/run12/final_0628/
echo $1
source setup.sh
make
./bin/RunppTestAna -i /gpfs01/star/pwg/youqi/run12/final_0628/raw_data/$1 -intype pico -c JetTree -trig ppJP2 -o /gpfs01/star/pwg/youqi/run12/final_0628/results/data/out_$1 -N -1 -pj 15 2000 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -geantnum 0 -hadcorr 1 -towunc 0 -fakeeff 1
