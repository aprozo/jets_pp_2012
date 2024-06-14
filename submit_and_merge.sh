#!/bin/bash
echo "cleaning that directory"
rm *session.xml
rm -r run12data*.package
rm run12*.zip
rm schedTemplateExp.xml
rm -r sched*.package

# echo "cleaning previous output"
# rm output/HT/*.root

echo "additional out/ report/ csh/"
rm sumbit/scheduler/csh/*
rm sumbit/scheduler/list/*
rm sumbit/scheduler/report/*

echo "sending jobs"
# python submit/submit.py # one job scheduler per file (1 hour to send all 1000 jobs, done almost immediately after that)

data_types=(geant mc)

for data_type in ${data_types[@]}; do
    echo "submitting $data_type"
    mkdir -p output/$data_type
    # check if data_file_$data_type.list exists
    if [ ! -f lists/data_file_$data_type.list ]; then
        echo "data_file_$data_type.list does not exist"
        break
    fi
    star-submit-template -template mysubmit.xml -entities trigger=$data_type
done

# wait for jobs to finish using condor_q
# then merge the output files
./wait_for_job.csh

for data_type in ${data_types[@]}; do
    singularity exec -e -B /gpfs01 -B /star star_star.simg \
        bash merge.sh $data_type
done
