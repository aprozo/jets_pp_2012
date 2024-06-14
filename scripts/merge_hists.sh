#!/bin/bash
data_type=$1
source scripts/setup.sh
hadd -f -k -j output/${data_type}_hists.root output/${data_type}/hists/*.root
