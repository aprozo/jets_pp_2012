#!/bin/bash
# run_bits.sh — rebuild the per-event trigger-bit trees the min-bias chain efficiency is measured from.
#
# Runs vpd_trigger_eff_jp_pass1b.C over the jet-patch picos in parallel chunks: per JP-fired event it
# writes the hardware min-bias trigger ids and the offline VPD flags into bits4/, the input of
# vpd_trigger_eff_jp_pass2d.C. Only the event header is read, so this is I/O-light although it touches
# every pico. Finally copies the per-radius eps_chain3_R<R>.root deliverables (written by pass2d into
# bits4/) to new_ana/inputs/, where stage2/unfold.C reads them.
#
#   bash new_ana/eps_chain/run_bits.sh [nchunks] [njobs]     (default: 16 chunks, 6 at a time)
#   then, per radius, inside the container from new_ana/eps_chain/:
#   root -l -b -q 'vpd_trigger_eff_jp_pass2d.C("bits4","0.5")'
set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
source "$REPO/site.sh"
LIST="$REPO/lists/jet_pico_dst/data.list"
OUT="$HERE/bits4"
NCHUNK="${1:-16}"
NJOBS="${2:-6}"

# 1. split the pico list round robin, so every chunk sees a similar mix of runs
mkdir -p "$OUT/lists"
rm -f "$OUT/lists"/chunk_*.list
awk -v n="$NCHUNK" -v d="$OUT/lists" '{print > (d "/chunk_" (NR%n) ".list")}' "$LIST"
echo "[bits] $(wc -l < "$LIST") picos -> $NCHUNK chunks, $NJOBS at a time"

# One chunk: pass 1 over its picos inside the container, with the Stage-1 reader bound over the
# container's own copy so that the pico classes match the ones the trees were written with.
run_one() {
   local chunk=$1
   singularity exec -e -B /gpfs01 -B /gpfs/mnt/gpfs01 \
     -B "$REPO/lib/eventStructuredAu_final:/usr/local/eventStructuredAu" "$SIMG" bash -c \
     "source /usr/local/root/bin/thisroot.sh; cd $HERE; \
      root -l -b -q 'vpd_trigger_eff_jp_pass1b.C(\"$OUT/lists/chunk_${chunk}.list\",\"$OUT/bits4_chunk_${chunk}.root\")'" \
     > "$OUT/chunk_${chunk}.log" 2>&1
   echo "[bits] chunk $chunk done: $(grep -h '^\[pass1b\]' "$OUT/chunk_${chunk}.log" | tail -1)"
}

# 2. run the chunks, at most NJOBS at a time
for ((chunk = 0; chunk < NCHUNK; ++chunk)); do
   run_one "$chunk" &
   while (( $(jobs -rp | wc -l) >= NJOBS )); do
      wait -n
   done
done
wait
echo "[bits] all chunks finished"
ls -la "$OUT"/bits4_chunk_*.root | tail -3

# 3. publish whatever pass 2 has already produced next to the other Stage-2 inputs
cp -f "$OUT"/eps_chain3_R*.root "$HERE/../inputs/" 2>/dev/null && echo "eps_chain3_R*.root -> new_ana/inputs/"
