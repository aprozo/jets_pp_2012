#!/bin/bash
# Shared shell functions of the Stage-2 drivers (published/run.sh, analysis/run.sh): ROOT inside the
# container from the stage2 directory, the data blocks, the response blocks of every systematic variant,
# one unfold per member and the band.
#
# Sourced with STAGE2 set to this directory; needs REPO, R (the jet radius) and site.sh's SIMG.
#   source "$STAGE2/lib.sh"
# Every step is skipped when its output is already there, so a driver can be re-run cheaply.

RES="$REPO/new_ana/results/R$R"
LOG="$RES/systematics_logs"
mkdir -p "$LOG"

# how many response builds may run at the same time: the per-run 2021 trees are small, the merged 2023
# tree is 23 GB and needs more memory per job
RESP_JOBS_2021=4
RESP_JOBS_2023=3

# a timestamped progress line
say() {
   echo "[$(date +%H:%M:%S)] $*"
}

# run ROOT inside the container, in the stage2 directory, with RooUnfold loaded
rt() {
   singularity exec -e -B /gpfs01 "$SIMG" bash -c "source /usr/local/root/bin/thisroot.sh; export ROOUNFOLD_HOME=/usr/local/RooUnfold; cd $STAGE2; root -l -b -q -e 'gSystem->Load(\"libRooUnfold\");' $*"
}

# block until fewer than $1 background jobs are running
wait_for_slot() {
   local limit=$1
   while [ "$(jobs -rp | wc -l)" -ge "$limit" ]; do
      wait -n
   done
}

# the systematic variants that need their own response build; the engine is the single source of the list
RESP_VARIANTS=$(rt "-e '.L systematics.C+' -e 'list_variants(\"response\")'" 2>/dev/null | sed -n 's/^VARIANTS //p')
[ -n "$RESP_VARIANTS" ] || { say "FATAL: variant list empty (systematics.C did not compile?)"; exit 1; }

# ---- data blocks: one pass writes the nominal and every data variant -----------------------------------
data_blocks() { # <jp|mb>
   local mode=$1
   local tag=""
   [ "$mode" = mb ] && tag="_mb"
   # the last name of the list is written last, so its presence means the whole pass finished
   local last
   last=$(rt "-e '.L systematics.C+' -e 'list_variants(\"databuilds\")'" 2>/dev/null | sed -n 's/^VARIANTS //p' | awk '{print $NF}')
   if [ -s "$RES/data_blocks_R${R}${tag}_${last}.root" ]; then
      say "data blocks ($mode) present"
      return
   fi
   say "data blocks ($mode): nominal + data variants"
   rt "-e '.L build_data.C+' -e 'build_data(\"$R\",\"$mode\")'" > "$LOG/build_data_$mode.log" 2>&1
}

# ---- response blocks from the per-run 2021 trees --------------------------------------------------------
# one variant of the 2021 response; "nominal" is the empty variant name
resp_2021_one() { # <variant|nominal> <file list>
   local name=$1
   local list=$2
   local log="$LOG/build_resp_2021_$name.log"
   local variant=$name
   [ "$variant" = nominal ] && variant=""
   rt "-e '.L build_resp.C+' -e 'build_resp(\"$R\",\"$list\",\"\",\"$variant\")'" > "$log" 2>&1
}

resp_2021() { # the nominal and every response variant from the per-run 2021 trees
   local list="$RES/matching_files.list"
   [ -s "$list" ] || find "$REPO/output/matching2021" -maxdepth 1 -name "matched_e21_*_R${R}.root" | sort > "$list"
   local todo=""
   local v
   for v in "" $RESP_VARIANTS; do
      [ -s "$RES/response_blocks_R${R}${v:+_$v}.root" ] || todo="$todo ${v:-nominal}"
   done
   [ -n "$todo" ] || { say "2021 responses present"; return; }
   say "2021 responses to build:$todo"
   rt "-e '.L build_resp.C+'" > "$LOG/compile_build_resp.log" 2>&1
   for v in $todo; do
      wait_for_slot $RESP_JOBS_2021
      resp_2021_one "$v" "$list" &
   done
   wait
}

# ---- response blocks from the merged 2023 tree -----------------------------------------------------------
# one variant of the 2023 response; the name carries the jet-definition prefix when there is one
resp_2023_one() { # <name|nominal> <e23|mb2023>
   local name=$1
   local which=$2
   local log="$LOG/build_resp_${which}_${name}.log"
   local variant=$name
   [ "$variant" = nominal ] && variant=""
   rt "-e '.L build_resp.C+' -e 'build_resp_merged(\"$R\",\"$which\",\"$variant\")'" > "$log" 2>&1
}

resp_2023() { # <e23|mb2023> [prefix]: the nominal and every response variant from the merged 2023 tree.
   # With a prefix (a jet definition of variants.h, e.g. noue) the builds are prefix and prefix_<variant>.
   local which=$1
   local prefix=${2:-}
   local todo=""
   local v
   for v in "" $RESP_VARIANTS; do
      local name="$v"
      [ -n "$prefix" ] && name="${prefix}${v:+_$v}"
      [ -s "$RES/response_blocks_R${R}_${which}${name:+_$name}.root" ] || todo="$todo ${name:-nominal}"
   done
   [ -n "$todo" ] || { say "$which responses present"; return; }
   say "$which responses to build:$todo"
   rt "-e '.L build_resp.C+'" > "$LOG/compile_build_resp.log" 2>&1
   for v in $todo; do
      wait_for_slot $RESP_JOBS_2023
      resp_2023_one "$v" "$which" &
   done
   wait
}

# ---- response blocks of the 2021 min-bias pass ------------------------------------------------------------
resp_2021_mb() { # the per-run trees of output/matchingMB2021, nominal only
   local list="$RES/matching_mb2021_files.list"
   [ -s "$list" ] || find "$REPO/output/matchingMB2021" -maxdepth 1 -name "matchedMB_*_R${R}.root" | sort > "$list"
   [ -s "$RES/response_blocks_R${R}_mb2021m.root" ] && { say "2021 min-bias response present"; return; }
   say "2021 min-bias response"
   rt "-e '.L build_resp.C+' -e 'build_resp(\"$R\",\"$list\",\"_mb2021m\")'" > "$LOG/build_resp_mb2021m.log" 2>&1
}

# ---- unfolding --------------------------------------------------------------------------------------------
unfold_one() { # <mode> <inv|bayes> <nIter>: skipped when the result is newer than every input block
   local md=$1 me=$2 n=$3
   local out="$RES/xsec_inversion_R${R}_${md}"; [ "$me" = bayes ] && out="${out}_bayes$n"; out="$out.root"
   [ "$md" = jp ] && [ "$me" = inv ] && out="$RES/xsec_inversion_R${R}.root"
   local newest; newest=$(ls -t "$RES"/response_blocks_R${R}*.root "$RES"/data_blocks_R${R}*.root 2>/dev/null | head -1)
   if [ -s "$out" ] && [ "$out" -nt "$newest" ]; then return; fi
   say "unfold $md $me $n"
   rt "-e '.L unfold.C+' -e 'unfold(\"$R\",\"$md\",\"$me\",$n)'" > "$LOG/unfold_${md}_${me}${n}.log" 2>&1 || say "  FAILED (see $LOG/unfold_${md}_${me}${n}.log)"
}

unfold_all() { # <mode> <inv|bayes> <nIter>: the nominal and every systematic member of the band
   local mode=$1
   local method=$2
   local nIter=$3
   local members
   members=$(rt "-e '.L systematics.C+' -e 'list_variants(\"members_$method\")'" 2>/dev/null | sed -n 's/^VARIANTS //p')
   local m
   for m in "" $members; do
      local md=$mode
      local n=$nIter
      # the regularisation members are the same mode at another iteration count, the rest are own modes
      case "$m" in
         iter2) n=2 ;;
         iter6) n=6 ;;
         "") ;;
         *) md="${mode}_$m" ;;
      esac
      unfold_one "$md" "$method" "$n"
   done
}

# ---- the systematic band -------------------------------------------------------------------------------------
combine() { # <mode> <inv|bayes> <nIter>: the band and its table
   local mode=$1
   local method=$2
   local nIter=$3
   local log="$LOG/combine_${mode}_${method}${nIter}.log"
   rt "-e '.L systematics.C+' -e 'combine(\"$R\",\"$mode\",\"$method\",$nIter)'" > "$log" 2>&1
   tail -n 14 "$log"
}
