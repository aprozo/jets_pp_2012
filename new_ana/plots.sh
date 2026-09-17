#!/bin/bash
# Every figure of the deck, from the result files:  bash new_ana/plots.sh  -> new_ana/results/figures/
#   published/plot.C              Table III, the ladder, the systematic band against Table I
#   analysis/plot.C               triggers, ratios to JP1, radii, UE on/off, the no-UE set,
#                                 the generators at every radius, the hadronisation correction
#   published/embedding_diff_plot.C   what separates the 2023 embedding sample from the 2021 one
#   published/plot_ladder.C       every rung of the method / embedding / level ladder (also copied
#                                 to docs/talk/, which docs/talk/ladder.tex includes)
# Needs the results of published/run.sh (including its "ladder" step), analysis/run.sh (all radii)
# and models/run_models.sh compare; a macro whose inputs are missing says so and moves on.
# The models files are written by the LCG ROOT (6.34) and read here by the container ROOT; opening
# them prints harmless 'TList::Clear ... already deleted' messages from the streamer-info list,
# filtered below.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
source "$REPO/site.sh"

# run one macro of one deck inside the analysis container, with the noisy ROOT lines filtered out
run() {
   singularity exec -e -B /gpfs01 "$SIMG" bash -c \
      "source /usr/local/root/bin/thisroot.sh; cd $HERE/$1; root -l -b -q '$2' \
       2>&1 | grep -v 'TCanvas::Print\|has been created\|TList::Clear\|THashList::Delete'"
}

run published "plot.C+(\"0.5\")"
run analysis  "plot.C+"
run published "embedding_diff_plot.C+"
run published "plot_ladder.C+(\"0.5\")"

ls "$HERE/results/figures" | wc -l | xargs echo "figures:"
