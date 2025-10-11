#!/bin/bash
set -uo pipefail

################################################################################
# CONFIGURATION
################################################################################

# Select which triggers to process (HT2 JP2 MB)
TRIGGERS=(HT2 JP2)
# Select which types to process (data embedding)
TYPE=(data embedding)

RADII=(0.2 0.3 0.4 0.5 0.6)

cleanup() {
    echo "cleaning that directory"
    rm -f *session.xml schedTemplateExp.xml *.zip *.dataset
    rm -rf *.package/ submit/scheduler/{csh,list,report,gen} submit/log/
    mkdir -p submit/{scheduler/{gen,csh,list,report},log}
}

merge_trees() {
    local trigger=$1
    local type=$2
    local out_prefix=$3

    echo "Merging trees for trigger: $trigger and type: $type"
    for R in "${RADII[@]}"; do
        echo "  Merging for R=$R"
        singularity exec -e -B /gpfs01 star_star.simg \
            hadd -f -k -j output/$trigger/${out_prefix}_R${R}.root \
            output/$trigger/$type/*_R${R}.root
    done
}


submit_trees(){
    for trigger in "${TRIGGERS[@]}"; do
        for type in "${TYPE[@]}"; do

            echo "Processing trigger: $trigger"
            mkdir -p output/$trigger/$type/

            if [ "$type" == "embedding" ]; then
                filelist="lists/jet_pico_dst/embedding.list"
            else
                filelist="lists/jet_pico_dst/data_${trigger}.list"
            fi

            if [ ! -f "${filelist}" ]; then
                echo "WARNING: ${filelist} does not exist, skipping..."
                continue
            fi
        echo "Submitting jobs for trigger: $trigger and type: $type"
        star-submit-beta-template \
            -template submit/mysubmit.xml \
            -entities trigger="$trigger",type="$type",filelist="$filelist"
        done
    done
}

submit_matching() {
    mkdir -p lists/trees
    for trigger in "${TRIGGERS[@]}"; do
        type="embedding"
        ls -d $PWD/output/$trigger/$type/mc*.root > lists/trees/$type"_"$trigger.list
        mkdir -p output/$trigger/matching
        star-submit-beta-template \
            -template submit/matching_mc_reco.xml \
            -entities trigger=$trigger
    done
}


main() {

    # source ~/.cshrc_nfs4
    # source ~/.login_nfs4

    cleanup
    echo "Triggers: ${TRIGGERS[*]} | Types: ${TYPE[*]}"
    echo "=============sending data==============="
    echo ""

    submit_trees
    ./scripts/condor_control.sh
    echo ""
    echo "=========trees are finished============="
    echo ""
   

    # submit matching if type has embedding
    if [[ " ${TYPE[*]} " == *" embedding "* ]]; then

        echo "===========submitting matching=========="
        echo ""
        submit_matching
    fi

    # meanwhile after submitting matching jobs, merge data trees
    echo ""
    echo "===========merging data trees==========="
    for trigger in "${TRIGGERS[@]}"; do
         merge_trees "$trigger" "data" "merged_data_$trigger"
    done
    echo ""
    echo "===========merging is finished==========="

    ./scripts/condor_control.sh

    echo ""
    echo "==========merging embedding trees======="
    echo ""
    if [[ " ${TYPE[*]} " == *" embedding "* ]]; then
         for trigger in "${TRIGGERS[@]}"; do
            echo "Merging matching trees for trigger: $trigger"
            merge_trees "$trigger" "matching" "merged_matching_$trigger"
        done
    fi

    echo ""
    echo "=========all jobs are finished=========="

}

# Run main function
main