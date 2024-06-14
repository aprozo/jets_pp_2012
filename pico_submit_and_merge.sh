#!/bin/bash
echo ""
echo "========================================"
echo "=====cleaning previous iteration========"
echo "========================================"
echo ""

rm *session.xml
rm -r run12data*.package
rm run12*.zip
rm schedTemplateExp.xml
rm -r sched*.package
rm sumbit/scheduler/csh/*
rm sumbit/scheduler/list/*
rm sumbit/scheduler/report/*

echo ""
echo "========================================"
echo "=============sending jobs==============="
echo "========================================"
echo ""

data_types=(geant mc HT2 MB)

for data_type in ${data_types[@]}; do
    echo "submitting $data_type"
    mkdir -p output/$data_type
    # check if data_file_$data_type.list exists
    if [ ! -f lists/data_file_$data_type.list ]; then
        echo "data_file_$data_type.list does not exist"
        break
    fi
    star-submit-template -template submit/mysubmit.xml -entities trigger=$data_type
done

# wait for jobs to finish using condor_q
# then merge the output files
./wait_for_job.csh

echo ""
echo "========================================"
echo "=========trees are finished============="
echo "===========sending analysis============="
echo "========================================"
echo ""

for data_type in ${data_types[@]}; do
    ls -d $PWD/output/$data_type/*.root >lists/trees_$data_type.list
    mkdir -p output/$data_type/hists
    star-submit-template -template submit/ana_trees.xml -entities trigger=$data_type
done

echo ""
echo "========================================"
echo "========merging tree files=============="
echo "========================================"
echo ""

for data_type in ${data_types[@]}; do
    singularity exec -e -B /gpfs01 star_star.simg \
        bash scripts/merge_trees.sh $data_type
done

./wait_for_job.csh

echo ""
echo "========================================"
echo "========merging hists files============="
echo "========================================"
echo ""

for data_type in ${data_types[@]}; do
    singularity exec -e -B /gpfs01 star_star.simg \
        bash scripts/merge_hists.sh $data_type
done
