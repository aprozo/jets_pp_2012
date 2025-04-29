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

# Function to submit analysis jobs for each data type
rerun_trees() {
    echo ""
    echo "========================================"
    echo "=============sending jobs==============="
    echo "========================================"
    echo ""

    for data_type in ${data_types[@]}; do
        echo "Submitting jobs for data type: $data_type"
        # Check if the list file exists before submitting
        mkdir -p output/$data_type
        # check if data_file_$data_type.list exists
        if [ ! -f lists/data/data_file_$data_type.list ]; then
            echo "data_file_$data_type.list does not exist"
            break
        fi
        star-submit-beta-template -template submit/mysubmit.xml -entities trigger=$data_type
    done

    ./scripts/condor_control.sh
    echo ""
    echo "========================================"
    echo "=========trees are finished============="
    echo "========================================"
    echo ""

    for data_type in ${data_types[@]}; do
        ls -d $PWD/output/$data_type/*.root >lists/trees/trees_$data_type.list
        # if JP2 or MB or HT2, merge trees
        if [ $data_type == "JP2" ] || [ $data_type == "MB" ] || [ $data_type == "HT2" ]; then
            echo "========================================"
            echo "Merging trees for data type: $data_type"
            echo "========================================"

            singularity exec -e -B /gpfs01 star_star.simg \
                bash scripts/merge_trees.sh $data_type
        fi
    done
}


matching_mc_geant() {
    echo ""
    echo "========================================"
    echo "=========matching mc and geant=========="
    echo "========================================"
    echo ""
    mkdir -p output/matching_mc_reco
    star-submit-beta submit/matching_mc_reco.xml
    ./scripts/condor_control.sh
    setup 64b
    hadd -f -k -j output/jets_embedding.root output/matching_mc_reco/*.root
    echo ""
    echo "========================================"
    echo "=========matching is finished==========="
}

####################################################################################################

# data_types=(JP2 HT2 MB mc geant)
# data_types=(mc geant)
cleanup
# rerun_trees
matching_mc_geant
