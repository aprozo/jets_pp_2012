#!/bin/bash
echo "cleaning that directory"
rm *session.xml
rm -r run12data*.package
rm run12*.zip
rm schedTemplateExp.xml
rm -r sched*.package

echo "cleaning previous output"
rm output/HT/*.root

echo "additional out/ report/ csh/"
rm sumbit/scheduler/csh/*
rm sumbit/scheduler/list/*
rm sumbit/scheduler/report/*

echo "sending jobs"
# python submit/submit.py # one job scheduler per file (1 hour to send all 1000 jobs, done almost immediately after that)

star-submit-template -template mysubmit.xml -entities trigger=MB
star-submit-template -template mysubmit.xml -entities trigger=HT2
star-submit-template -template mysubmit.xml -entities trigger=embedding

# wait for jobs to finish using condor_q
# then merge the output files
./wait_for_job.csh

singularity exec -e -B /gpfs01 -B /star star_star.simg \
    bash merge.sh MB
singularity exec -e -B /gpfs01 -B /star star_star.simg \
    bash merge.sh HT2
singularity exec -e -B /gpfs01 -B /star star_star.simg \
    bash merge.sh embedding



