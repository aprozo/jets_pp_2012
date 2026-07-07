#!/bin/bash
# Matrix-inversion driver — Dmitry's unregularized inversion, the second solver
# alongside the main Bayes pipeline (../run.sh). Builds the square etadet
# response (unfold.cxx) then unfolds + normalizes + compares the data
# (cross_section.cpp), on this repo's single-production trees.
#
# No environment flags: triggers, paths, iterations and thread count are all
# hardcoded in ../config.h. Same physical inputs as Bayes (etadet response +
# measured C(pt) trigger correction + per-trigger quote window); no fudge.
#
# Usage:
#   bash run_inverse.sh              # response + cross section (all triggers)
#   bash run_inverse.sh response     # square response only
#   bash run_inverse.sh cross_section

set -eo pipefail

CSI_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIMG=/gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg
step="${1:-all}"

if [[ -z "${APPTAINER_NAME:-}${SINGULARITY_NAME:-}" ]]; then
    if [[ ! -f "$SIMG" ]]; then
        echo "Container image not found: $SIMG" >&2
        exit 1
    fi
    exec singularity exec -e -B /gpfs01 -B /gpfs/mnt/gpfs01 "$SIMG" bash "$0" "$@"
fi

source /usr/local/root/bin/thisroot.sh
export ROOUNFOLD_HOME=/usr/local/RooUnfold
export ROOT_INCLUDE_PATH="${ROOUNFOLD_HOME}/src:${ROOT_INCLUDE_PATH:-}"
export LD_LIBRARY_PATH="${ROOUNFOLD_HOME}:${LD_LIBRARY_PATH:-}"

echo "ROOT: $(root-config --version)   ROOUNFOLD_HOME: $ROOUNFOLD_HOME"

run_response() {
    cd "$CSI_DIR"
    root -l -b -q -e 'gSystem->Load("libRooUnfold");' 'unfold.cxx+'
}
run_cross_section() {
    cd "$CSI_DIR"
    root -l -b -q -e 'gSystem->Load("libRooUnfold");' 'cross_section.cpp+'
}

case "$step" in
    response|unfold)  run_response ;;
    cross_section)    run_cross_section ;;
    all)              run_response ; run_cross_section ;;
    *) echo "Unknown step: $step (expected: all | response | cross_section)" >&2; exit 1 ;;
esac

echo ""
echo "Done. Outputs in ../ (config.h kWorkDir):"
echo "  response_<T>_R0.5_square.root, xsec_<T>_R0.5_invert.root, comparison_with_dmitriy_R0.5_invert.pdf"
