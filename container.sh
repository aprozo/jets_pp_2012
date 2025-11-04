#!/bin/bash
source /usr/local/root/bin/thisroot.sh
export FASTJETDIR=/usr/local/fastjet
export STARPICOPATH=/usr/local/eventStructuredAu
export JETREADER=/usr/local/jetreader_build
export LD_LIBRARY_PATH=/usr/local/jetreader_build/lib:/lib/:/usr/local/eventStructuredAu:/usr/local/RooUnfold:/usr/local/fastjet/lib:/usr/local/root/lib::/.singularity.d/libs

# --- Input parameters ---
input_file=${1}
data_type=${2}

output_file=${data_type}_$(basename $input_file)
# make clean
make


jet_min_pt=4.8
trigger=$data_type

# Set default values for tree type, pico type, and trigger based on data type
treeName="JetTree"
picoType="pico"
# if data_type contains "mc" like mc_HT2 or mc_JP2 or mc_MB
if [[ $data_type == mc_* ]]; then
    treeName="JetTreeMc"
    picoType="mcpico"
fi



# Define arguments for the RunppAna command
args=(
    -i "$input_file"
    -intype "$picoType"
    -c "$treeName"
    -trig "$trigger"
    -N -1
    -pj $jet_min_pt 200
    -pc 0.2 30
    -lja "antikt"
    -ec 1.0
    -geantnum 1
    -hadcorr 1
    -towunc 0
    -fakeeff 1
)

# remove 3 last rows if data_type is mc_*
if [[ $data_type == mc_* ]]; then
    unset 'args[-1]'
    unset 'args[-1]'
    unset 'args[-1]'
    unset 'args[-1]'
    unset 'args[-1]'
    unset 'args[-1]'
    args+=(-jetnef 1)
fi

# --- Sweep over jet radius ---
RADII=(0.2 0.3 0.4 0.5 0.6)
#  RADII=(0.6)
echo "Running with input:  $input_file"
echo "Data type:           $data_type"

# Execute the analysis command with specified arguments

for R in "${RADII[@]}"; do
   
    # add to ouput_file the radius info
    output_file_R=${output_file%.root}_R${R}.root
    args+=(-R "$R")
    args+=(-o "$output_file_R")
    
    echo "Writing to output:   ${output_file_R}"
    printf ' %q' "${args[@]}"; echo
    ./bin/RunppAna "${args[@]}"

    echo "command was ./bin/RunppAna ${args[@]}"
    # Remove the last 4 elements (-R and its value) to reset for the next iteration
    unset 'args[-1]'
    unset 'args[-1]'
    unset 'args[-1]'
    unset 'args[-1]'
done
