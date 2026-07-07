#!/bin/bash
# Stage-2 driver (Bayes). Builds the per-trigger Miss/Fake response
# (unfolding/unfold.cxx), unfolds + normalizes + compares the data
# (cross_section.cpp), and overlays every trigger against Dmitry's Table III
# (plot_alltriggers.C).
#
# No environment flags: triggers, paths, iterations and thread count are all
# hardcoded in config.h. Edit config.h to change them.
#
# Usage:
#   bash run.sh                 # response + cross section + plot (all triggers)
#   bash run.sh response        # response only
#   bash run.sh cross_section   # unfold + normalize + plot only

set -eo pipefail

NEW_ANA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIMG=/gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg
step="${1:-all}"

# ---- re-exec inside star_star.simg -----------------------------------------
if [[ -z "${APPTAINER_NAME:-}${SINGULARITY_NAME:-}" ]]; then
    if [[ ! -f "$SIMG" ]]; then
        echo "Container image not found: $SIMG" >&2
        exit 1
    fi
    exec singularity exec -e -B /gpfs01 -B /gpfs/mnt/gpfs01 "$SIMG" bash "$0" "$@"
fi

# ---- inside the container ---------------------------------------------------
source /usr/local/root/bin/thisroot.sh
export ROOUNFOLD_HOME=/usr/local/RooUnfold
export ROOT_INCLUDE_PATH="${ROOUNFOLD_HOME}/src:${ROOT_INCLUDE_PATH:-}"
export LD_LIBRARY_PATH="${ROOUNFOLD_HOME}:${LD_LIBRARY_PATH:-}"

echo "ROOT: $(root-config --version)   ROOUNFOLD_HOME: $ROOUNFOLD_HOME"

run_response() {
    cd "$NEW_ANA_DIR/unfolding"
    root -l -b -q -e 'gSystem->Load("libRooUnfold");' 'unfold.cxx+'
}
run_cross_section() {
    cd "$NEW_ANA_DIR"
    root -l -b -q -e 'gSystem->Load("libRooUnfold");' 'cross_section.cpp+'
}
run_plot() {
    cd "$NEW_ANA_DIR"
    root -l -b -q "plot_alltriggers.C(\"${NEW_ANA_DIR}/\",\"comparison_alltriggers_R0.5.pdf\")"
}

case "$step" in
    response|unfold)  run_response ;;
    cross_section)    run_cross_section ; run_plot ;;
    all)              run_response ; run_cross_section ; run_plot ;;
    *) echo "Unknown step: $step (expected: all | response | cross_section)" >&2; exit 1 ;;
esac

echo ""
echo "Done. Outputs in $NEW_ANA_DIR :"
echo "  response_<T>_R0.5.root, xsec_<T>_R0.5.root, comparison_alltriggers_R0.5.pdf"
