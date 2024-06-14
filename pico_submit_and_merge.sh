#!/bin/bash

# Function to clean up files
cleanup() {
    rm *session.xml
    rm -r run12data*.package
    rm run12*.zip
    rm schedTemplateExp.xml
    rm -r sched*.package
    rm sumbit/scheduler/csh/*
    rm sumbit/scheduler/list/*
    rm sumbit/scheduler/report/*
}
# Function to submit jobs
submit_tree_create() {
    local data_type=$1
    echo "submitting $data_type"
    mkdir -p output/$data_type
    # check if data_file_$data_type.list exists
    if [ ! -f lists/data_file_$data_type.list ]; then
        echo "data_file_$data_type.list does not exist"
        break
    fi
    star-submit-template -template submit/mysubmit.xml -entities trigger=$data_type
}

# Function to send analysis
submit_analysis() {
    local data_type=$1
    ls -d $PWD/output/$data_type/*.root >lists/trees_$data_type.list
    mkdir -p output/$data_type/hists
    star-submit-template -template submit/ana_trees.xml -entities trigger=$data_type
}
# Function to merge hists
merge_hists() {
    local data_type=$1
    singularity exec -e -B /gpfs01 star_star.simg \
        bash scripts/merge_hists.sh $data_type

}
# Function to merge trees
merge_trees() {
    local data_type=$1
    singularity exec -e -B /gpfs01 star_star.simg \
        bash scripts/merge_trees.sh $data_type
}

rerun_trees() {
    echo ""
    echo "========================================"
    echo "=============sending jobs==============="
    echo "========================================"
    echo ""

    for data_type in ${data_types[@]}; do
        submit_tree_create $data_type
    done

    ./wait_for_job.csh
    echo ""
    echo "========================================"
    echo "=========trees are finished============="
}

run_analysis() {
    echo ""
    echo "========================================"
    echo "===========sending analysis============="
    echo "========================================"
    echo ""
    for data_type in ${data_types[@]}; do
        submit_analysis $data_type
    done

    ./condor_control.sh
    for data_type in ${data_types[@]}; do
        merge_hists $data_type
    done
    wait
    echo ""
    echo "========================================"
    echo "=========analysis is finished==========="

}

####################################################################################################
####################################################################################################
####################################################################################################
# Main script execution

data_types=(MB HT2 mc geant)
cleanup
# rerun_trees
run_analysis
