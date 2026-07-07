#!/bin/bash
# JPX promotion-combination driver — Dmitry's published nominal: the exclusive
# JP0+JP1+JP2 partition summed with prescale-recovery weights, unfolded ONCE by
# matrix inversion on the cell-filtered fine response.
#
# No environment flags: windows, weights, fine grid, floor and lambda are all
# hardcoded in promotion.h (+ shared physics in ../config.h).
#
# Usage:
#   bash run_jpx.sh [step] [systName]
#     step     = all | response | cross_section        (default all)
#     systName = a preset from config.h::Systematics() (default nominal)
#   bash run_jpx.sh cross_section       # filter + invert + compare only (fast)

set -eo pipefail

JPX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIMG=/gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg
step="${1:-all}"
syst="${2:-nominal}"

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

echo "ROOT: $(root-config --version)"

run_response() {
    cd "$JPX_DIR"
    root -l -b -q "response.cxx+(\"$syst\")"
}
run_cross_section() {
    cd "$JPX_DIR"
    root -l -b -q "cross_section.cpp+(\"$syst\")"
}

case "$step" in
    response)        run_response ;;
    cross_section)   run_cross_section ;;
    all)             run_response ; run_cross_section ;;
    *) echo "Unknown step: $step (expected: all | response | cross_section)" >&2; exit 1 ;;
esac

echo ""
echo "Done. Outputs in ../ (config.h kWorkDir):"
echo "  response_JPX_R0.5_fine.root, xsec_JPX_R0.5.root, comparison_with_dmitriy_R0.5_JPX.pdf"
