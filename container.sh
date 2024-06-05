#!/bin/bash
export FASTJETDIR=/usr/local/fastjet
#export PYTHIA8DIR=/gpfs01/star/pwg/elayavalli/pythia8303
export STARPICOPATH=/usr/local/eventStructuredAu
export JETREADER=/usr/local/jetreader_build
export LD_LIBRARY_PATH=/usr/local/jetreader_build/lib:/lib/:/usr/local/eventStructuredAu:/usr/local/RooUnfold:/usr/local/fastjet/lib:/usr/local/root/lib::/.singularity.d/libs

make
./bin/RunppTestAna -i /gpfs01/star/pwg/prozorov/jets_pp_2012/ppRun12Datapicos/ppJP2Run12/sum0.root -intype pico -c JetTree -trig ppJP2 -o /gpfs01/star/pwg/prozorov/jets_pp_2012/output/out.root -N 10000 -pj 10 40 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -hadcorr 0.9999999 -geantnum 0
