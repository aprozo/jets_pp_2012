#!/bin/bash
# Stage-1 worker (runs inside star_star.simg). ONE input pico -> ALL jet radii,
# ALL triggers in a single ResultTree. Triggers are NOT a production split: every
# per-jet trigger_match_JP0/JP1/JP2/HT2 bit and per-event fired_JP0/JP1/JP2 bit is
# always written, and the trigger is selected downstream at Stage-2. So there is
# exactly one data production and one embedding production.
#
#   container.sh <pico> data      [outdir]   -> data_<pico>_R<R>.root      (ResultTree)
#   container.sh <pico> embedding [outdir]   -> matched_<pico>_R<R>.root   (MatchedTree)
#                                               (+ intermediate mc_/geant_)
#
# The jet-find floor is 4.8 GeV for data, geant (reco) and mc (truth): the same
# low floor keeps low-pT reco jets for the response (matched <reco/mc> ~ Dmitry),
# and the per-trigger analysis floor (config.h::TrigPtFloor) is re-applied at
# Stage-2. Jets are found wide (|eta_phys|<1.0, hardcoded in ppAnalysis.cxx) and
# the |det_eta|<0.5 detector acceptance is applied downstream. The JP trigger
# match is hardware-first on the kOnline patches (ppAnalysis::match_jp) — no
# environment flags.
set -e
source /usr/local/root/bin/thisroot.sh
export FASTJETDIR=/usr/local/fastjet
export STARPICOPATH=/usr/local/eventStructuredAu
export JETREADER=/usr/local/jetreader_build
export LD_LIBRARY_PATH=/usr/local/jetreader_build/lib:/lib/:/usr/local/eventStructuredAu:/usr/local/RooUnfold:/usr/local/fastjet/lib:/usr/local/root/lib::/.singularity.d/libs

input_file=${1:?usage: container.sh <pico> <data|embedding> [outdir]}
mode=${2:?usage: container.sh <pico> <data|embedding> [outdir]}
output_dir=${3:-.}
mkdir -p "$output_dir"

# jet radii produced per input file
RADII=(${JETS_RADII:-0.5})

# Build the binary once if the sandbox did not ship it.
[ -x bin/RunppAna ] || make

# Collision-free output stem. A file list can mix productions with DUPLICATE
# basenames (e.g. 0.root in both the main and the supplementary pico dirs), so
# key the output name on a short hash of the FULL input path, not the basename
# alone — otherwise the duplicate jobs overwrite each other and ~7% of the data
# silently vanishes from the merge.
base=$(basename "$input_file")
tag=$(printf '%s' "$input_file" | cksum | cut -d' ' -f1)
stem="${base%.root}_${tag}"

# run_pass <type> : one RunppAna pass over $input_file for every radius.
#   type = data  -> detector level, real data       (pico / JetTree)
#   type = geant -> detector level, embedding reco   (pico / JetTree)
#   type = mc    -> particle level, embedding truth  (mcpico / JetTreeMc)
run_pass() {
    local type=$1
    local picoType="pico" treeName="JetTree"
    local args=(-i "$input_file" -N "${JETS_NEV:--1}" -pj 4.8 200 -pc 0.2 200 -lja "antikt" -ec 1.0 -jetnef 1)
    if [[ $type == mc ]]; then
        picoType="mcpico"; treeName="JetTreeMc"
        args+=(-geantnum 1)                        # truth: no hadronic/tower/track systematics
    else
        args+=(-geantnum 1 -hadcorr 1 -towunc 0 -fakeeff 1)
    fi
    args=(-intype "$picoType" -c "$treeName" -trig "JP2" "${args[@]}")
    for R in "${RADII[@]}"; do
        local out="${output_dir%/}/${type}_${stem}_R${R}.root"
        echo "[$type] R=$R -> $out"
        ./bin/RunppAna "${args[@]}" -R "$R" -o "$out"
    done
}

case "$mode" in
    data)
        run_pass data
        ;;
    embedding)
        run_pass mc
        run_pass geant
        for R in "${RADII[@]}"; do
            mc_name="mc_${stem}_R${R}.root"
            echo "[match] R=$R -> matched_${base%.root}_R${R}.root"
            root -l -b -q "macros/matching_mc_reco.cxx+(\"${mc_name}\", true, \"${output_dir%/}/\")"
        done
        ;;
    *)
        echo "unknown mode: $mode (expected data | embedding)" >&2
        exit 1 ;;
esac
