#!/bin/bash
# The analysis: every trigger alone (MB, JP1, JP2, HT2), Bayesian unfolding (4 iterations) on the 2023
# embedding, with the systematic band, at R = 0.2-0.5; the ratios to JP1; UE on/off; the no-UE reference
# set on the 5-60 GeV bins; the radius ratios with a band.
#
#   bash new_ana/analysis/run.sh [bayes|michal|radii|all] [R]      (default: all 0.5)
#
#   bayes   per trigger: 2023 responses (nominal + variants), unfolds, bands, ratios to JP1
#             -> results/R<R>/xsec_<trigger>_R<R>_syst_bayes4.{root,txt}
#             -> results/R<R>/ratio_<trigger>_over_jp1_e23_R<R>_bayes4.{root,txt}
#   michal  the no-UE reference set: JP1 on the 5-60 GeV bins with its band, the combined and min-bias
#           columns, and the UE on/off table
#             -> results/Michal/xsec_R<R>_noUE_bins5-60.{root,txt}
#             -> results/R<R>/ue_onoff_R<R>_bayes4.txt
#   radii   sigma(R)/sigma(0.5) with a band, R = 0.2 0.3 0.4 (after bayes and michal at every radius)
#             -> results/R<R>/rratio_<mode>_R<R>_over_R0.5_bayes4.{root,txt}
#   all     bayes, michal (and radii when R = 0.5)
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

# The quoted result: every trigger unfolded with Bayes on the 2023 embedding, each systematic member
# unfolded and combined into a band, then the ratios of the other triggers to JP1.
bayes() {
   data_blocks jp
   data_blocks mb
   resp_2023 e23
   resp_2023 mb2023
   for T in jp1_e23 jp2_e23 ht2_e23 mb2023; do
      unfold_all "$T" bayes 4
      combine "$T" bayes 4
   done
   rt "-e '.L systematics.C+' -e 'ratios(\"$R\",\"jp1_e23\",\"mb2023 jp2_e23 ht2_e23\",\"bayes\",4)'" \
      > "$LOG/ratios.log" 2>&1
   tail -n 40 "$LOG/ratios.log"
}

# The no-UE reference set for the R_AA comparison: the raw-jet ("noue") jet definition on the coarse
# 5-60 GeV bins, with the JP1 band, plus the UE on/off table of the nominal bins.
michal() {
   data_blocks jp
   data_blocks mb
   resp_2023 e23
   resp_2023 e23 noue
   resp_2023 mb2023
   resp_2023 mb2023 noue
   unfold_all jp1_e23_noue_mbins bayes 4
   for md in jp_e23_noue_mbins mb2023_noue_mbins jp1_e23 jp2_e23 ht2_e23 mb2023 \
             jp1_e23_noue jp2_e23_noue ht2_e23_noue mb2023_noue; do
      unfold_one $md bayes 4
   done
   combine jp1_e23_noue_mbins bayes 4
   rt "-e '.L ../analysis/ue_onoff.C+' -e 'ue_onoff(\"$R\")'" > "$LOG/ue_onoff.log" 2>&1
   tail -n 2 "$LOG/ue_onoff.log"
   rt "-e '.L ../analysis/michal.C+' -e 'michal(\"$R\")'" > "$LOG/michal.log" 2>&1
   tail -n 3 "$LOG/michal.log"
}

# sigma(R) / sigma(0.5) with a systematic band, for the nominal and for the no-UE jet definition.
radii() {
   for md in jp1_e23 jp1_e23_noue_mbins; do
      rt "-e '.L systematics.C+' -e 'ratios_radii(\"$md\",\"0.2 0.3 0.4\",\"0.5\",\"bayes\",4)'" \
         > "$LOG/ratios_radii_$md.log" 2>&1
      tail -n 30 "$LOG/ratios_radii_$md.log"
   done
}

case "$STEP" in
   bayes)
      bayes
      ;;
   michal)
      michal
      ;;
   radii)
      radii
      ;;
   all)
      bayes
      michal
      [ "$R" = 0.5 ] && radii
      ;;
   *)
      sed -n 2,19p "$0" >&2
      exit 1
      ;;
esac
say "done"
