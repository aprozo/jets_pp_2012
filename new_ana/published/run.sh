#!/bin/bash
# The published analysis, reproduced: exclusive jet-patch levels, the 2021 embedding, matrix inversion on the
# published bins, and the published systematic recipe. Compared with Table III (cross section) and Table I (band).
#
#   bash new_ana/published/run.sh [data|response|unfold|band|all|mb|mb2021m|ladder] [R]     (default: all 0.5)
#
#   data      detector spectra of the levels, nominal + data variants           -> results/R<R>/data_blocks_R<R>*.root
#   response  2021 eta-block responses, nominal + response variants (per run)   -> results/R<R>/response_blocks_R<R>*.root
#   unfold    the cross section, matrix inversion                               -> results/R<R>/xsec_inversion_R<R>.root
#   band      every member unfolded and combined                                -> results/R<R>/xsec_jp_R<R>_syst.{root,txt}
#   all       data, response, unfold, band
#   mb        the min-bias level with the 2023 min-bias response, inversion     -> xsec_inversion_R<R>_mb2023.root
#   mb2021m   the min-bias level with the 2021 min-bias response                -> xsec_inversion_R<R>_mb2021m*.root
#   ladder    every rung of the method / embedding / level ladder (needs analysis/run.sh bayes first) -> results/figures/published_ladder_*.pdf, copies in docs/talk/
#
# Figures: bash new_ana/plots.sh
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
STAGE2="$REPO/new_ana/stage2"
source "$REPO/site.sh"
STEP=${1:-all}
R=${2:-0.5}
source "$STAGE2/lib.sh"

case "$STEP" in
   # the detector-level spectra of the three exclusive jet-patch levels
   data)
      data_blocks jp
      ;;
   # the 2021 embedding response, per eta block, nominal and every systematic member
   response)
      resp_2021
      ;;
   # the nominal cross section only
   unfold)
      data_blocks jp
      resp_2021
      unfold_one jp inv 0
      ;;
   # the cross section plus the full systematic band
   band)
      data_blocks jp
      resp_2021
      unfold_all jp inv 0
      combine jp inv 0
      ;;
   # the default: data, response, unfold and band in one go
   all)
      data_blocks jp
      resp_2021
      unfold_all jp inv 0
      combine jp inv 0
      ;;
   # the min-bias level, unfolded with the 2023 min-bias response
   mb)
      data_blocks mb
      resp_2023 mb2023
      unfold_one mb2023 inv 0
      ;;
   # the min-bias level with the 2021 min-bias response, both unfolding methods
   mb2021m)
      data_blocks mb
      resp_2021_mb
      unfold_one mb2021m inv 0
      unfold_one mb2021m bayes 4
      ;;
   # every rung of the talk: method (inversion vs Bayes iterations), embedding (2021 vs 2023, calibrated),
   # and level (the single triggers), then the figures
   ladder)
      data_blocks jp
      data_blocks mb
      resp_2021
      resp_2023 e23
      resp_2023 mb2023
      resp_2021_mb
      unfold_one jp inv 0
      for n in 2 4 10 30 100; do unfold_one jp bayes $n; done   # 30 and 100 only for the error figure
      unfold_one jp_e23 inv 0
      unfold_one jp_e23cal inv 0
      for n in 2 4 10; do unfold_one jp_e23 bayes $n; done
      unfold_one jp_e23cal bayes 4
      for m in jp1_e23 jp2_e23 mb2021m mb2023; do
         unfold_one $m inv 0
         unfold_one $m bayes 4
      done
      rt "-e '.L ../published/plot_ladder.C' -e 'plot_ladder(\"$R\")'" > "$LOG/plot_ladder.log" 2>&1; tail -n 3 "$LOG/plot_ladder.log"
      ;;
   *)
      sed -n 2,16p "$0" >&2
      exit 1
      ;;
esac
say "done"
