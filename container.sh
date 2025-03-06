#!/bin/bash
source /usr/local/root/bin/thisroot.sh
export FASTJETDIR=/usr/local/fastjet
#export PYTHIA8DIR=/gpfs01/star/pwg/elayavalli/pythia8303
export STARPICOPATH=/usr/local/eventStructuredAu
export JETREADER=/usr/local/jetreader_build
export LD_LIBRARY_PATH=/usr/local/jetreader_build/lib:/lib/:/usr/local/eventStructuredAu:/usr/local/RooUnfold:/usr/local/fastjet/lib:/usr/local/root/lib::/.singularity.d/libs
input_file=${1}
output_file=$(basename $input_file)
data_type=${2}

rm bin/RunppAna
make

# Set default values for tree type, pico type, and trigger based on data type
trigger="All"
treeType="JetTree"
picoType="pico"
jet_min_pt=5

# Determine trigger type and tree/pico settings based on data_type
case "$data_type" in
HT2)
    trigger="ppHT"
    ;;
JP2)
    trigger="ppJP"
    ;;
mc)
    treeType="JetTreeMc"
    picoType="mcpico"
    ;;
geant)
    jet_min_pt=5
    ;;
geant_JP2)
    trigger="ppJP"
    jet_min_pt=5
    ;;
geant_HT2)
    trigger="ppHT"
    jet_min_pt=5
    ;;

esac

# Display the input and output file information
echo "Running with input file: $input_file"
echo "Running with output file: $output_file"

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
    -R 0.4
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

# Execute the analysis command with specified arguments
./bin/RunppAna "${args[@]}"

echo "Command was: 
./bin/RunppAna ${args[@]}"
