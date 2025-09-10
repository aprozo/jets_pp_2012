#!/bin/bash


source /usr/local/root/bin/thisroot.sh
export FASTJETDIR=/usr/local/fastjet
#export PYTHIA8DIR=/gpfs01/star/pwg/elayavalli/pythia8303
export STARPICOPATH=/usr/local/eventStructuredAu
export JETREADER=/usr/local/jetreader_build
export LD_LIBRARY_PATH=/usr/local/jetreader_build/lib:/lib/:/usr/local/eventStructuredAu:/usr/local/RooUnfold:/usr/local/fastjet/lib:/usr/local/root/lib::/.singularity.d/libs

# --- Input parameters ---
input_file=${1}
output_file=$(basename $input_file)
data_type=${2}

make clean
make

# Set default values for tree type, pico type, and trigger based on data type
treeType="JetTree"
picoType="pico"
jet_min_pt=5
trigger="All"
# Determine trigger type and tree/pico settings based on data_type

case "$data_type" in
    mc)
        treeType="JetTreeMc"
        picoType="mcpico"
        jet_min_pt=4
        ;;
    geant) ;;
    MB)  trigger="ppMB" ;;
    HT2) trigger="ppHT" ;;
    JP2) trigger="ppJP2" ;;
    *)   echo "Unknown data_type: $data_type"; usage ;;
esac

# Define arguments for the RunppAna command
args=(
    -i "$input_file"
    -intype "$picoType"
    -c "$treeType"
    -trig "$trigger"
    -o "tree_$output_file"
    -N -1
    -pj $jet_min_pt 2000
    -pc 0.2 30
    -lja "antikt"
    -ec 1
    -geantnum 1
)

# Append hadronic correction argument if running on Monte Carlo data
if [[ $data_type == "geant" ]]; then
    args+=(-hadcorr 0.9999999)
    args+=(-towunc 0)
    args+=(-fakeeff 1)
fi

if [[ $data_type == "mc" ]]; then
    args+=(-jetnef 1)
fi

# --- Sweep over jet radius ---
RADII=(0.2 0.3 0.4 0.5 0.6)

echo "Running with input:  $input_file"
echo "Data type:           $data_type"

# Execute the analysis command with specified arguments

for R in "${RADII[@]}"; do
   
    # add to ouput_file the radius info
    output_file_R=${output_file%.root}_R${R}.root
    args+=(-R "$R")
    args+=(-o "tree_$output_file_R")
    
    echo "Writing to output:   tree_${output_file_R}"
    printf ' %q' "${args[@]}"; echo
    ./bin/RunppAna "${args[@]}"
    # Remove the last 4 elements (-R and its value) to reset for the next iteration
    unset 'args[-1]'
    unset 'args[-1]'
    unset 'args[-1]'
    unset 'args[-1]'
done
