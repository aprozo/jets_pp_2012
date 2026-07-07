#!/bin/bash
# Measure the two in-situ trigger corrections that are frozen into
# config.h::TrigEffMeas as C(pt) = R(pt) x That(pt):
#   hw_ratio.C  -> R(pt)   = N(hardware-matched)/N(simulator-matched)
#   measure_T.C -> That(pt) = data/embedding turn-on shape
# These are DERIVATION tools: the production pipeline uses the frozen tables in
# config.h. Re-run these only to re-derive the tables (then hand-edit config.h).
#
# Usage:  bash corrections/run.sh
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATAPATH="$(cd "$HERE/../.." && pwd)/output/"   # <repo>/output/ (edit if it moves)
SIMG=/gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg

if [[ -z "${APPTAINER_NAME:-}${SINGULARITY_NAME:-}" ]]; then
    exec singularity exec -e -B /gpfs01 -B /gpfs/mnt/gpfs01 "$SIMG" bash "$0" "$@"
fi

source /usr/local/root/bin/thisroot.sh
cd "$HERE"

DATA="${DATAPATH%/}/merged_data_R0.5.root"
MATCH="${DATAPATH%/}/merged_matching_R0.5.root"

echo "=== R(pt): hardware vs simulator (hw_ratio.C) ==="
for T in JP0 JP1 JP2; do
    root -l -b -q "hw_ratio.C(\"$T\")"
done

echo "=== That(pt): data/embedding turn-on shape (measure_T.C) ==="
# adcRegisterPlus1 = DSM register + 1 : JP1 29, JP2 37
root -l -b -q "measure_T.C(\"$DATA\", \"$MATCH\", 29, \"JP1\")"
root -l -b -q "measure_T.C(\"$DATA\", \"$MATCH\", 37, \"JP2\")"

echo "Done. R -> hw_ratio_<T>.root/pdf ; That -> measure_T_<T>.pdf (stdout table)."
echo "Update config.h::TrigEffMeas with C = R x That per bin if re-deriving."
