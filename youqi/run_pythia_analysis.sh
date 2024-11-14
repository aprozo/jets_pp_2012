#!/bin/bash

for file in `cat 20235003.list`
	do
	./bin/RunppTestAna -i '/gpfs01/star/pwg/youqi/run12/embedding/P12id/picos/20235003/out'$file -intype mcpico -c JetTreeMc -trig all -o 'results_20235003/pythia/Pythia'$file -N -1 -pj 5 2000 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -geantnum 1 -jetnef 1
	#./bin/RunppTestAna -i '/gpfs01/star/pwg/youqi/run12/embedding/P12id/picos/pp2012embed_2023/out/'$file -intype mcpico -c JetTreeMc -trig all -o 'results/pythia/Pythia'$file -N -1 -pj 5 2000 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -geantnum 1 -jetnef 1 
	#./bin/RunppTestAna -i '/gpfs01/star/pwg/elayavalli/ppRun12Embpicos/'$file -intype mcpico -c JetTreeMc -trig all -o 'Results/pythia/Pythia'$file -N -1 -pj 5 2000 -mj 0.0 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -geantnum 1 -jetnef 1

done
