#!/bin/bash
trigger=$1
source setup.sh
hadd -f -k -j output/output_jets_${trigger}.root output/${trigger}/*.root
