#!/bin/bash

# input files sum0.root to sum9.root
# approximate nevents = 54500000
# split each file into 100 files
# each file will have 54500000/100 = 545000 events
# last file will have what is left after nsplit-1

nevents=54500000 # approximate number of events
nsplit=25        # number of files to split into
nevents_per_file=$((nevents / nsplit))
input_file=$1

# path_to_data="/gpfs01/star/pwg/prozorov/jets_pp_2012/ppRun12Datapicos/ppMBRun12"

# cd $path_to_data
# mkdir -p splitFiles
# output_path=$path_to_data/splitFiles
# for i in {0..9}; do
input_tree=$input_file:JetTree

filename=$(basename $input_file)
filename=${filename%.root}

echo "============================="

for part in $(seq 0 $((nsplit - 2))); do
    echo "splitting ${input_tree} into sum${i}_${part}.root"
    start_event=$((part * nevents_per_file))
    last_event=$((start_event + nevents_per_file - 1))
    rooteventselector -f $start_event -l $last_event $input_tree ${filename}_part${part}.root
done
last_event=$((last_event + 1)) # last event in the last file)
echo "splitting last part because don't know what is total number of events in file"
rooteventselector -f $last_event $input_tree ${filename}_part$((nsplit - 1)).root
# done
