#!/bin/bash
# Interactive shell inside the Stage-1/2 container (ROOT, fastjet, RooUnfold, the TStarJetPico reader),
# with the home directory and /gpfs01 bound. The image path comes from site.sh.
#   bash scripts/runimage.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/site.sh"
singularity shell -B "$HOME" -B /gpfs01 "$SIMG"
