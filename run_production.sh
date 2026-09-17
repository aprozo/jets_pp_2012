#!/bin/bash
# Stage-1 driver: build the binary, submit one condor job per pico, wait, merge.
# Each job produces the requested radii and all trigger bits (the trigger is a
# Stage-2 filter), so there is one data and one embedding production.
#
#   ./run_production.sh [build|submit|wait|merge|all|mb|mb-*|mbembed|mbembed-*] ["R1 R2 ..."]
#   build compiles the TStarJetPico reader from the maker repository into lib/ and then bin/RunppAna
#
# The optional 2nd argument is the space-separated radius list (default 0.5),
# e.g.  ./run_production.sh all "0.2 0.3 0.4"
#
# Outputs (per radius R):
#   output/merged_data_R<R>.root      (ResultTree, all triggers)
#   output/merged_matching_R<R>.root  (MatchedTree, all triggers)
#
# star-submit-template must run on the HOST; the build + merge run inside
# star_star.simg. After any src/ change: rm -f run12prod.zip *.package, else the
# scheduler ships a stale sandbox binary. Never rebuild the sandbox zip while
# previously submitted jobs are still starting up — they fail at unzip.
set -u
RR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$RR/site.sh"
DATALIST=lists/jet_pico_dst/data.list
EMBEDLIST=lists/jet_pico_dst/embedding2023.list  # 2023 request: response of the main deck (JP and MB)
EMBEDLIST2021=lists/jet_pico_dst/embedding.list  # 2021 request: response used to compare against the published analysis
DATALISTMB=lists/jet_pico_dst/dataMB.list        # VPDMB-nobsmd (min-bias) picos
STAGE=${1:-all}
RADII="${2:-0.5}"
RADII_ENT="${RADII// /_}"    # entity form for the XML (0.2_0.3_0.4)
FIRSTR="${RADII%% *}"        # progress polling uses the first radius
NDATA=$(wc -l < "$RR/$DATALIST")
NEMB=$(wc -l < "$RR/$EMBEDLIST")
NMB=$(wc -l < "$RR/$DATALISTMB" 2>/dev/null || echo 0)
say(){ echo "[$(date +%H:%M:%S)] $*"; }

# queued jobs matching a pattern, or -1 when condor_q failed (a timeout is NOT an empty queue)
qcount(){ local out; out=$(timeout 120 condor_q -nobatch 2>/dev/null) || { echo -1; return; }; echo "$out" | grep -c "$1" || true; }

# The TStarJetPico reader: the maker repository's StRoot/eventStructuredAu (the source of truth; it
# carries the hadronic correction over the analysis track list) compiled inside the analysis
# container into lib/eventStructuredAu (not tracked), which every job binds over the container's
# own copy. The ROOT dictionary is regenerated here, so the maker's own (built by the STAR ROOT 5)
# is left out of the copy.
stage_reader(){
  local src="$MAKER_REPO/StRoot/eventStructuredAu" dst="$RR/lib/eventStructuredAu"
  [ -d "$src" ] || { say "FATAL: reader sources not found at $src (MAKER_REPO in site.sh)"; exit 1; }
  mkdir -p "$dst"
  rsync -a --delete --exclude='TStarJetPicoRootDict*' --include='*.cxx' --include='*.h' --include='Makefile' \
        --include='Makefile.arch' --include='*.txt' --include='*.list' --exclude='*' "$src/" "$dst/"
  say "building the reader from $src"
  singularity exec -e -B /gpfs01 "$SIMG" bash -c "source /usr/local/root/bin/thisroot.sh && cd $dst && make -s" > "$dst/make.log" 2>&1
  [ -f "$dst/libTStarJetPico.so" ] || { say "FATAL: reader build failed (see $dst/make.log)"; exit 1; }
}

stage_build(){
  stage_reader
  say "building bin/RunppAna"
  singularity exec -e -B /gpfs01 -B "$RR/lib/eventStructuredAu:/usr/local/eventStructuredAu" "$SIMG" bash -c "cd $RR && source /usr/local/root/bin/thisroot.sh && \
     export FASTJETDIR=/usr/local/fastjet STARPICOPATH=/usr/local/eventStructuredAu && make"
  [ -x "$RR/bin/RunppAna" ] || { say "FATAL: build failed"; exit 1; }
}

stage_submit(){
  cd "$RR"
  rm -f run12prod.zip
  mkdir -p output/data output/matching submit/log submit/scheduler/gen
  say "submitting data ($NDATA picos) + embedding ($NEMB picos), radii: $RADII"
  star-submit-template -template submit/production.xml -entities type=data,filelist=$DATALIST,radii=$RADII_ENT
  star-submit-template -template submit/production.xml -entities type=embedding,filelist=$EMBEDLIST,radii=$RADII_ENT
}

stage_wait(){
  for i in $(seq 1 500); do
    local q nd nm
    q=$(qcount "jets_")
    nd=$(ls "$RR"/output/data/data_*_R${FIRSTR}.root 2>/dev/null | wc -l)
    nm=$(find "$RR"/output/matching -maxdepth 1 -name "matched_*_R${FIRSTR}.root" 2>/dev/null | wc -l)
    say "poll $i: condor=$q data=$nd/$NDATA matched=$nm/$NEMB (R=$FIRSTR)"
    [ "$nd" -ge "$NDATA" ] && [ "$nm" -ge "$NEMB" ] && break
    [ "$q" -eq 0 ] && [ "$i" -gt 5 ] && { say "queue drained (data=$nd matched=$nm)"; break; }
    sleep 180
  done
}

stage_merge(){
  local STG="$RR/scratch_merge"
  for R in $RADII; do
    rm -rf "$STG"; mkdir -p "$STG"
    say "merging data -> merged_data_R${R}.root (chunked hadd)"
    singularity exec -e -B /gpfs01 "$SIMG" bash -c "
      source /usr/local/root/bin/thisroot.sh; cd $STG
      ls $RR/output/data/data_*_R${R}.root | split -l 600 - dchunk_
      k=0; for c in dchunk_*; do hadd -f dsub_\$k.root \$(cat \$c) >/dev/null 2>&1; k=\$((k+1)); done
      hadd -f $RR/output/merged_data_R${R}.root dsub_*.root 2>&1 | tail -1"
    say "merging matching -> merged_matching_R${R}.root (chunked hadd; find, not a glob: the file count exceeds the argument limit)"
    singularity exec -e -B /gpfs01 "$SIMG" bash -c "
      source /usr/local/root/bin/thisroot.sh; cd $STG; rm -f mchunk_* msub_*.root
      find $RR/output/matching -maxdepth 1 -name 'matched_*_R${R}.root' | sort | split -l 600 - mchunk_
      k=0; for c in mchunk_*; do hadd -f msub_\$k.root \$(cat \$c) >/dev/null 2>&1; k=\$((k+1)); done
      hadd -f $RR/output/merged_matching_R${R}.root msub_*.root 2>&1 | tail -1"
    rm -rf "$STG"
    say "merged R=$R: $(du -h $RR/output/merged_data_R${R}.root 2>/dev/null | cut -f1) data, \
$(du -h $RR/output/merged_matching_R${R}.root 2>/dev/null | cut -f1) matching"
  done
}

# ---- MB (VPDMB-nobsmd) data-only production -----------------------------------
# The min-bias analog of the data flow: same Stage-1 jet pass on the VPDMB picos
# (container.sh mbdata -> RunppAna -trig MB), routed to DISTINCT paths so nothing
# overwrites the JP/HT production: output/dataMB/ and merged_dataMB_R<R>.root.
stage_submit_mb(){
  cd "$RR"
  if condor_q -nobatch 2>/dev/null | grep -q "jets_mb"; then
    say "WARNING: MB jobs still queued — rebuilding run12prod.zip WILL break them (Ctrl-C now, or backfill later)"
    sleep 10
  fi
  rm -f run12prod.zip; rm -rf run12prod.package     # ship the updated container.sh + dataMB.list
  mkdir -p output/dataMB submit/log submit/scheduler/gen
  say "submitting MB data ($NMB picos), radii: $RADII"
  star-submit-template -template submit/production.xml -entities type=mbdata,filelist=$DATALISTMB,radii=$RADII_ENT
}

stage_wait_mb(){
  for i in $(seq 1 500); do
    local q nd
    q=$(qcount "jets_mbdata")
    nd=$(ls "$RR"/output/dataMB/mbdata_*_R${FIRSTR}.root 2>/dev/null | wc -l)
    say "poll $i: condor=$q mbdata=$nd/$NMB (R=$FIRSTR)"
    [ "$nd" -ge "$NMB" ] && break
    [ "$q" -eq 0 ] && [ "$i" -gt 5 ] && { say "queue drained (mbdata=$nd)"; break; }
    sleep 180
  done
}

stage_merge_mb(){
  local STG="$RR/scratch_merge_mb"
  for R in $RADII; do
    rm -rf "$STG"; mkdir -p "$STG"
    # hadd does NOT de-duplicate: repeated production rounds leave several files
    # per pico index (mbdata_VPDMB_pp200_2012_<index>_<jobhash>_R<R>.root), so keep
    # only the newest per index (ls -t = newest first) or events are multi-counted.
    ls -t "$RR"/output/dataMB/mbdata_*_R${R}.root | awk -F'_' '!seen[$(NF-2)]++' > "$STG/dedup_list.txt"
    say "dedup: $(wc -l < "$STG/dedup_list.txt") unique picos (of $(ls "$RR"/output/dataMB/mbdata_*_R${R}.root | wc -l) files)"
    say "merging MB data -> merged_dataMB_R${R}.root (chunked hadd)"
    singularity exec -e -B /gpfs01 "$SIMG" bash -c "
      source /usr/local/root/bin/thisroot.sh; cd $STG
      split -l 600 dedup_list.txt dchunk_
      k=0; for c in dchunk_*; do hadd -f dsub_\$k.root \$(cat \$c) >/dev/null 2>&1; k=\$((k+1)); done
      hadd -f $RR/output/merged_dataMB_R${R}.root dsub_*.root 2>&1 | tail -1"
    rm -rf "$STG"
    say "merged MB R=$R: $(du -h $RR/output/merged_dataMB_R${R}.root 2>/dev/null | cut -f1) data"
  done
}

# ---- MB embedding (same dijet sample, tight DCA + vertex branches) -------------
# Produces the MB RESPONSE inputs. Distinct paths (output/matchingMB/,
# merged_matchingMB_R<R>.root) so the JP embedding files are NEVER overwritten.
stage_submit_mbembed(){
  cd "$RR"
  # This stage rebuilds run12prod.zip, so let the queue drain first (or plan a
  # backfill pass): jobs still starting up would fail at unzip.
  if condor_q -nobatch 2>/dev/null | grep -q "jets_mb"; then
    say "WARNING: MB jobs still queued — rebuilding run12prod.zip WILL break them (Ctrl-C now, or backfill later)"
    sleep 10
  fi
  # the merge has no de-duplication: a previous production must be moved away first
  if [ -n "$(ls output/matchingMB 2>/dev/null | head -1)" ]; then say "FATAL: output/matchingMB is not empty — move the previous production to TO_DELETE first"; exit 1; fi
  rm -f run12prod.zip; rm -rf run12prod.package
  mkdir -p output/matchingMB submit/log submit/scheduler/gen
  say "submitting MB embedding ($NEMB picos, same dijet sample), radii: $RADII"
  star-submit-template -template submit/production.xml -entities type=mbembedding,filelist=$EMBEDLIST,radii=$RADII_ENT
}

stage_wait_mbembed(){
  for i in $(seq 1 500); do
    local q nm
    q=$(qcount "jets_mbembedding")
    nm=$(find "$RR"/output/matchingMB -maxdepth 1 -name "matchedMB_*_R${FIRSTR}.root" 2>/dev/null | wc -l)
    say "poll $i: condor=$q matchedMB=$nm/$NEMB (R=$FIRSTR)"
    [ "$nm" -ge "$NEMB" ] && break
    [ "$q" -eq 0 ] && [ "$i" -gt 5 ] && { say "queue drained (matchedMB=$nm)"; break; }
    sleep 180
  done
}

stage_merge_mbembed(){
  for R in $RADII; do
    say "merging MB embedding -> merged_matchingMB_R${R}.root"
    singularity exec -e -B /gpfs01 "$SIMG" bash -c "
      source /usr/local/root/bin/thisroot.sh
      mkdir -p $RR/scratch_merge_mbe; cd $RR/scratch_merge_mbe; rm -f mchunk_* msub_*.root
      find $RR/output/matchingMB -maxdepth 1 -name 'matchedMB_*_R${R}.root' | sort | split -l 600 - mchunk_
      k=0; for c in mchunk_*; do hadd -f msub_\$k.root \$(cat \$c) >/dev/null 2>&1; k=\$((k+1)); done
      hadd -f $RR/output/merged_matchingMB_R${R}.root msub_*.root 2>&1 | tail -1; cd $RR; rm -rf $RR/scratch_merge_mbe"
    say "merged MB embedding R=$R: $(du -h $RR/output/merged_matchingMB_R${R}.root 2>/dev/null | cut -f1)"
  done
}

# ---- MB (min-bias mode) pass over the 2021 embedding ---------------------------
# Same container mode on the 2021 picos (per-run names e21_<sample>_<run>), distinct
# output (output/matchingMB2021, template submit/production_mb2021.xml) so the 2023
# min-bias trees are untouched; feeds the mb2021m step of new_ana/published/run.sh.
NEMB2021=$(wc -l < "$RR/$EMBEDLIST2021")
stage_submit_mbembed2021(){
  cd "$RR"
  if [ -n "$(ls output/matchingMB2021 2>/dev/null | head -1)" ]; then say "FATAL: output/matchingMB2021 is not empty"; exit 1; fi
  # The sandbox zip is REUSED (same container.sh/macros/bin) and deliberately not
  # rebuilt here: doing so while the 2023 jobs are still starting breaks them.
  [ -f run12prod.zip ] || { say "FATAL: run12prod.zip missing — run mbembed-submit first (it builds the sandbox)"; exit 1; }
  mkdir -p output/matchingMB2021 submit/log submit/scheduler/gen
  say "submitting MB embedding on the 2021 request ($NEMB2021 picos), radii: $RADII"
  star-submit-template -template submit/production_mb2021.xml -entities type=mbembedding,filelist=$EMBEDLIST2021,radii=$RADII_ENT
}
stage_wait_mbembed2021(){
  for i in $(seq 1 500); do
    local q nm
    q=$(qcount "jets_mbembedding2021")
    nm=$(find "$RR"/output/matchingMB2021 -maxdepth 1 -name "matchedMB_*_R${FIRSTR}.root" 2>/dev/null | wc -l)
    say "poll $i: condor=$q matchedMB2021=$nm/$NEMB2021 (R=$FIRSTR)"
    [ "$nm" -ge "$NEMB2021" ] && break
    [ "$q" -eq 0 ] && [ "$i" -gt 5 ] && { say "queue drained (matchedMB2021=$nm)"; break; }
    sleep 180
  done
}

case "$STAGE" in
  mbembed2021-submit) stage_submit_mbembed2021 ;;
  mbembed2021-wait)   stage_wait_mbembed2021 ;;
  build)  stage_build ;;
  submit) stage_submit ;;
  wait)   stage_wait ;;
  merge)  stage_merge ;;
  all)    stage_build; stage_submit; stage_wait; stage_merge ;;
  mb-submit) stage_submit_mb ;;
  mb-wait)   stage_wait_mb ;;
  mb-merge)  stage_merge_mb ;;
  mb)     stage_build; stage_submit_mb; stage_wait_mb; stage_merge_mb ;;
  mbembed-submit) stage_submit_mbembed ;;
  mbembed-wait)   stage_wait_mbembed ;;
  mbembed-merge)  stage_merge_mbembed ;;
  mbembed) stage_build; stage_submit_mbembed; stage_wait_mbembed; stage_merge_mbembed ;;
  *) echo "usage: run_production.sh [build|submit|wait|merge|all|mb|mb-*|mbembed|mbembed-*|mbembed2021-submit|mbembed2021-wait] [\"R1 R2 ...\"] (default radius 0.5)"; exit 1 ;;
esac
say "stage(s) [$STAGE] complete"
