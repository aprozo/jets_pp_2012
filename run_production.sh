#!/bin/bash
# Stage-1 driver. ONE data production + ONE embedding production (triggers are a
# Stage-2 filter). Each condor job = one pico -> all radii, all trigger bits.
#
#   ./run_production.sh [build|submit|wait|merge|all]
#
# Outputs:
#   output/merged_data_R0.5.root      (ResultTree, all triggers)
#   output/merged_matching_R0.5.root  (MatchedTree, all triggers)
#
# star-submit-template runs on the HOST; the build + merge run inside star_star.simg.
# After any src/ change: rm -f run12prod.zip *.package (else the scheduler ships a
# stale sandbox binary).
set -u
RR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIMG=/gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg
DATALIST=lists/jet_pico_dst/data.list
EMBEDLIST=lists/jet_pico_dst/embedding.list
RADIUS=0.5
STAGE=${1:-all}
NDATA=$(wc -l < "$RR/$DATALIST")
NEMB=$(wc -l < "$RR/$EMBEDLIST")
say(){ echo "[$(date +%H:%M:%S)] $*"; }

stage_build(){
  say "building bin/RunppAna"
  singularity exec -e -B /gpfs01 "$SIMG" bash -c "cd $RR && source /usr/local/root/bin/thisroot.sh && \
     export FASTJETDIR=/usr/local/fastjet STARPICOPATH=/usr/local/eventStructuredAu && make"
  [ -x "$RR/bin/RunppAna" ] || { say "FATAL: build failed"; exit 1; }
}

stage_submit(){
  cd "$RR"
  rm -f run12prod.zip *.package
  mkdir -p output/data output/matching submit/log submit/scheduler/gen
  say "submitting data ($NDATA picos) + embedding ($NEMB picos)"
  star-submit-template -template submit/production.xml -entities type=data,filelist=$DATALIST
  star-submit-template -template submit/production.xml -entities type=embedding,filelist=$EMBEDLIST
}

stage_wait(){
  for i in $(seq 1 500); do
    local q nd nm
    q=$(condor_q -nobatch 2>/dev/null | grep -c "jets_" || true)
    nd=$(ls "$RR"/output/data/data_*_R${RADIUS}.root 2>/dev/null | wc -l)
    nm=$(ls "$RR"/output/matching/matched_*_R${RADIUS}.root 2>/dev/null | wc -l)
    say "poll $i: condor=$q data=$nd/$NDATA matched=$nm/$NEMB"
    [ "$nd" -ge "$NDATA" ] && [ "$nm" -ge "$NEMB" ] && break
    [ "$q" -eq 0 ] && [ "$i" -gt 5 ] && { say "queue drained (data=$nd matched=$nm)"; break; }
    sleep 180
  done
}

stage_merge(){
  local STG="$RR/scratch_merge"
  rm -rf "$STG"; mkdir -p "$STG"
  say "merging data -> merged_data_R${RADIUS}.root (chunked hadd)"
  singularity exec -e -B /gpfs01 "$SIMG" bash -c "
    source /usr/local/root/bin/thisroot.sh; cd $STG
    ls $RR/output/data/data_*_R${RADIUS}.root | split -l 600 - dchunk_
    k=0; for c in dchunk_*; do hadd -f dsub_\$k.root \$(cat \$c) >/dev/null 2>&1; k=\$((k+1)); done
    hadd -f $RR/output/merged_data_R${RADIUS}.root dsub_*.root 2>&1 | tail -1"
  say "merging matching -> merged_matching_R${RADIUS}.root"
  singularity exec -e -B /gpfs01 "$SIMG" bash -c "
    source /usr/local/root/bin/thisroot.sh
    hadd -f $RR/output/merged_matching_R${RADIUS}.root $RR/output/matching/matched_*_R${RADIUS}.root 2>&1 | tail -1"
  rm -rf "$STG"
  say "merged: $(du -h $RR/output/merged_data_R${RADIUS}.root 2>/dev/null | cut -f1) data, \
$(du -h $RR/output/merged_matching_R${RADIUS}.root 2>/dev/null | cut -f1) matching"
}

case "$STAGE" in
  build)  stage_build ;;
  submit) stage_submit ;;
  wait)   stage_wait ;;
  merge)  stage_merge ;;
  all)    stage_build; stage_submit; stage_wait; stage_merge ;;
  *) echo "usage: run_production.sh [build|submit|wait|merge|all]"; exit 1 ;;
esac
say "stage(s) [$STAGE] complete"
