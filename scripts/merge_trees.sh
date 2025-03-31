#!/bin/bash
data_type=$1
source ./scripts/setup.sh
hadd -f -k -j output/jets_${data_type}.root output/${data_type}/*.root
