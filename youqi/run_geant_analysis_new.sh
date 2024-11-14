#!/bin/bash

for file in `cat geant_new.list`
	do
	./bin/RunppTestAna -i '/gpfs01/star/pwg/youqi/run12/embedding/P12id/picos/pp2012embed_2023/out/'$file -intype pico -c JetTree -trig ppJP2 -o 'results/geant/Geant'$file -N -1 -pj 15 2000 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -geantnum 1 -hadcorr 1 -towunc 0 -fakeeff 1
done
