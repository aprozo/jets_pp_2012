#!/bin/bash
# Build the systematic band. Config-as-code: the variation matrix lives in
# new_ana/config.h::Systematics() (typed presets), and this driver runs each and
# envelopes them. No environment flags, no variation table file.
#
#   1. build the nominal response, plus a rebuilt response per SHAPE variation
#      (jesShift/jerSmear -> response_<T>_R0.5_<name>.root);
#   2. run cross_section for every variation (reuses the nominal response unless
#      the variation is shape-level) -> xsec_<T>_R0.5_<name>.root, each stamped
#      with its full provenance;
#   3. combine.C takes the per-bin envelope -> xsec_<T>_R0.5_systband.root.
#
#   bash systematics/run_systematics.sh
#
# Triggers are whatever config.h lists. Add/remove sources by editing
# config.h::Systematics().
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEW_ANA="$(cd "$HERE/../new_ana" && pwd)"
SIMG=/gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg

if [[ -z "${APPTAINER_NAME:-}${SINGULARITY_NAME:-}" ]]; then
    exec singularity exec -e -B /gpfs01 -B /gpfs/mnt/gpfs01 "$SIMG" bash "$0" "$@"
fi

source /usr/local/root/bin/thisroot.sh
export ROOUNFOLD_HOME=/usr/local/RooUnfold
export ROOT_INCLUDE_PATH="${ROOUNFOLD_HOME}/src:${ROOT_INCLUDE_PATH:-}"
export LD_LIBRARY_PATH="${ROOUNFOLD_HOME}:${LD_LIBRARY_PATH:-}"

resp() { cd "$NEW_ANA/unfolding"; root -l -b -q -e 'gSystem->Load("libRooUnfold");' "unfold.cxx+(\"$1\")"; }
xsec() { cd "$NEW_ANA";           root -l -b -q -e 'gSystem->Load("libRooUnfold");' "cross_section.cpp+(\"$1\")"; }

# The committed variation matrix (name + needsResponse) from config.h.
mapfile -t VARS < <(cd "$HERE" && root -l -b -q list_systematics.C 2>/dev/null | sed -n 's/^SYST //p')

# 1 (nominal response) + 2 (nominal xsec) come from the "nominal" row below.
resp nominal
for row in "${VARS[@]}"; do
    name=${row%% *}
    needsResp=${row##* }
    if [[ "$name" != nominal && "$needsResp" == 1 ]]; then
        resp "$name"          # rebuild the shifted response for a shape variation
    fi
    xsec "$name"              # unfold + normalize + write xsec_<T>_R0.5<tag>.root
done

# 3: envelope band per trigger
cd "$HERE"
root -l -b -q combine.C

echo ""
echo "Done. Band -> $NEW_ANA/xsec_<T>_R0.5_systband.root"
