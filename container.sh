#!/bin/bash
# Stage-1 worker (runs inside star_star.simg). ONE input pico -> ALL jet radii,
# ALL trigger bits in a single tree; the trigger is selected downstream at
# Stage-2, so there is exactly one data and one embedding production.
#
#   container.sh <pico> data      [outdir] [radii]  -> data_<pico>_R<R>.root      (ResultTree)
#   container.sh <pico> embedding [outdir] [radii]  -> matched_<pico>_R<R>.root   (MatchedTree)
#                                                      (+ intermediate mc_/geant_)
#   (mbdata / mbembedding: the min-bias variants, separate measurement)
#
# Jet inputs are identical for data and embedding reco. Some of the selection
# lives upstream and is NOT visible here: the tower DB status mask and the
# >125 cm last-TPC-point cut are applied at Stage-0, the rest (track |eta|, DCA,
# nHits, pT window; tower Et and hadronic correction; anti-kT and the jet-find
# floor) by the RunppAna flags below. No event vetoes.
#
# Interactive-only conveniences (inert under condor: production.xml runs
# `singularity exec -e`, so the environment never reaches the worker):
#   JETS_RADII  default radius list when the 4th argument is absent (default 0.5)
#   JETS_NEV    event cap per pass for smoke tests (default -1 = all events)
set -e
source /usr/local/root/bin/thisroot.sh
export FASTJETDIR=/usr/local/fastjet
export STARPICOPATH=/usr/local/eventStructuredAu
export JETREADER=/usr/local/jetreader_build
export LD_LIBRARY_PATH=/usr/local/jetreader_build/lib:/lib/:/usr/local/eventStructuredAu:/usr/local/RooUnfold:/usr/local/fastjet/lib:/usr/local/root/lib::/.singularity.d/libs

input_file=${1:?usage: container.sh <pico> <data|embedding> [outdir] [radii]}
mode=${2:?usage: container.sh <pico> <data|embedding> [outdir] [radii]}
output_dir=${3:-.}
mkdir -p "$output_dir"

# Jet radii per input file. The 4th argument (underscore-separated, e.g.
# "0.2_0.3_0.4") wins: the environment does not survive `singularity exec -e`,
# so run_production.sh threads the radii through the XML as an argument.
RADII=(${JETS_RADII:-0.5})
[ -n "${4:-}" ] && RADII=(${4//_/ })
echo "[container] radii: ${RADII[*]}"

# Build the binary once if the sandbox did not ship it.
[ -x bin/RunppAna ] || make

# Collision-free output stem: a file list can mix productions with DUPLICATE
# basenames, so key the output name on a hash of the FULL input path. Using the
# basename alone lets duplicate jobs overwrite each other and silently drops
# data from the merge.
base=$(basename "$input_file")
tag=$(printf '%s' "$input_file" | cksum | cut -d' ' -f1)
stem="${base%.root}_${tag}"

# run_pass <type> : one RunppAna pass over $input_file for every radius.
#   type = data  -> detector level, real data       (pico / JetTree)
#   type = geant -> detector level, embedding reco   (pico / JetTree)
#   type = mc    -> particle level, embedding truth  (mcpico / JetTreeMc)
run_pass() {
    local type=$1
    local trig=${2:-JP2}                            # MB production passes "MB"; data/embedding use "JP2"
    local picoType="pico" treeName="JetTree"
    local args=(-i "$input_file" -N "${JETS_NEV:--1}" -pj ${PJMIN:-4.8} 200 -pc 0.2 200 -lja "antikt" -ec 2.5 -jetnef 1)
    if [[ $type == mc ]]; then
        picoType="mcpico"; treeName="JetTreeMc"
        args+=(-geantnum 1)                        # truth: no hadronic/tower/track systematics
    else
        args+=(-geantnum 1 -hadcorr 1 -towunc 0 -fakeeff 1)
    fi
    # Space-separated extra flag string, set by the mbdata / mbembedding modes.
    # -leadsdca 0.5 vetoes jets whose LEADING charged track is a fake (large
    # |sDCAxy|); a flat per-track sDCAxy cut instead removes ~a third of the
    # real bulk. Harmless on the mc pass (MC sDCAxy off).
    if [[ -n "${EXTRA_TRACK_CUTS:-}" ]]; then
        args+=(${EXTRA_TRACK_CUTS})
    fi
    args=(-intype "$picoType" -c "$treeName" -trig "$trig" "${args[@]}")
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
    mbdata)
        EXTRA_TRACK_CUTS="-leadsdca 0.5"
        run_pass mbdata "MB"
        ;;
    embedding)
        # Response input: truth jets from the pico MC tree (jet-find floor 1.5 GeV
        # so every reco jet above the analysis floor has a truth partner), reco jets
        # with the same settings as data; matching writes matched_* (MatchedTree).
        PJMIN=1.5 run_pass mc
        run_pass geant
        for R in "${RADII[@]}"; do
            mc_name="mc_${stem}_R${R}.root"
            echo "[match] R=$R -> matched_${stem}_R${R}.root"
            root -l -b -q "macros/matching_mc_reco.cxx+(\"${mc_name}\", true, \"${output_dir%/}/\")"
        done
        ;;
    mbembedding)
        # MB response: same dijet embedding, but the reco pass uses the tight MB
        # track cuts (as the MB data) and the matched trees carry vertex info.
        # Renamed to matchedMB_* so the JP embedding files are never overwritten.
        EXTRA_TRACK_CUTS="-leadsdca 0.5"
        PJMIN=1.5 run_pass mc          # truth: 1.5 GeV floor is REQUIRED — without it the truth
                                       # jets below 4.8 GeV are missing and their reco partners
                                       # are miscounted as fakes (MC forces the DCA cuts off)
        run_pass geant                 # reco: tight DCA applied
        for R in "${RADII[@]}"; do
            mc_name="mc_${stem}_R${R}.root"
            echo "[match] R=$R -> matchedMB_${stem}_R${R}.root"
            root -l -b -q "macros/matching_mc_reco.cxx+(\"${mc_name}\", true, \"${output_dir%/}/\")"
            mv -f "${output_dir%/}/matched_${stem}_R${R}.root" "${output_dir%/}/matchedMB_${stem}_R${R}.root"
        done
        ;;
    *)
        echo "unknown mode: $mode (expected data | mbdata | embedding | mbembedding)" >&2
        exit 1 ;;
esac
