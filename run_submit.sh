#!/bin/bash
# Function to clean up files
cleanup() {
    echo "cleaning that directory"
    rm *session.xml
    rm -r *.package/

    rm *.zip
    rm schedTemplateExp.xml
    rm sched*.package
    rm *.dataset

    echo "additional out/ report/ csh/"
    rm -rf submit/scheduler/csh/
    rm -rf submit/scheduler/list/
    rm -rf submit/scheduler/report/
    rm -rf submit/scheduler/gen
    rm -rf submit/log/

    mkdir -p submit/scheduler/gen
    mkdir -p submit/scheduler/csh
    mkdir -p submit/scheduler/list
    mkdir -p submit/scheduler/report
    mkdir -p submit/log/

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
            # add here loop merge by radius
            for R in ${radii[@]}; do
                echo "Merging for R=$R"
                singularity exec -e -B /gpfs01 star_star.simg \
                    hadd -f -k -j output/jets_${data_type}_R${R}.root output/${data_type}/*_R${R}.root
            done
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
    echo "merging trees"

    for R in ${radii[@]}; do
        echo "Merging for R=$R"
        singularity exec -B /gpfs01 star_star.simg \
            hadd -f -k -j output/matching_mc_reco/jets_embedding_R${R}.root output/matching_mc_reco/*_R${R}.root
    done
    echo ""
    echo "========================================"
    echo "=========matching is finished==========="
}

####################################################################################################

data_types=(JP2 HT2 MB mc geant)
radii=(0.2 0.3 0.4 0.5 0.6)

cleanup
rerun_trees
# if data_types contains mc or geant then run matching_mc_geant
if [[ " ${data_types[@]} " =~ " mc " ]] || [[ " ${data_types[@]} " =~ " geant " ]]; then
    matching_mc_geant
fi
