#!/bin/bash
echo "cleaning that directory"
rm *session.xml
rm -r run12data*.package
rm run12*.zip
rm schedTemplateExp.xml
rm -r sched*.package

echo "cleaning previous output"
rm output/out/*.root
rm output/log/*

echo "additional out/ report/ csh/"
rm sumbit/scheduler/csh/*
rm sumbit/scheduler/list/*
rm sumbit/scheduler/report/*

echo "sending jobs"
# python submit/submit.py # one job scheduler per file (1 hour to send all 1000 jobs, done almost immediately after that)
star-submit mysubmit.xml # one scheduler for all files

# wait for jobs to finish using condor_q
# then merge the output files
./wait_for_job.csh

singularity exec -e \
    -B /gpfs01 -B /star \
    star_star.simg \
    bash merge.sh
