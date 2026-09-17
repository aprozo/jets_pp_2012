#!/bin/bash
# Write the four Stage-1 input lists from the Stage-0 production folders (MAKER_REPO from site.sh).
#   bash scripts/make_pico_lists.sh <data-folder> <mb-folder> <emb2023-folder> <emb2021-folder> [<emb2021-folder> ...]
# folder = submit/<date>/job_<name> under MAKER_REPO; several 2021 folders are merged (one pico per run).
# Result: lists/jet_pico_dst/{data,dataMB,embedding2023,embedding}.list
set -eu

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$REPO/site.sh"
[ $# -ge 4 ] || { sed -n 2,5p "$0" >&2; exit 1; }

OUT="$REPO/lists/jet_pico_dst"
mkdir -p "$OUT"

# the three streams that come from a single production folder each
find "$MAKER_REPO/$1/production" -maxdepth 1 -name '*.root' | sort -V > "$OUT/data.list"
find "$MAKER_REPO/$2/production" -maxdepth 1 -name '*.root' | sort    > "$OUT/dataMB.list"
find "$MAKER_REPO/$3/production" -maxdepth 1 -name '*.root' | sort    > "$OUT/embedding2023.list"

# the 2021 embedding is one pico per run and may be spread over several production folders
shift 3
: > "$OUT/embedding.list"
for folder in "$@"; do
   find "$MAKER_REPO/$folder/production" -maxdepth 1 -name 'e21_*.root'
done | sort > "$OUT/embedding.list"

# every 2021 run must appear exactly once: a duplicate would double-count its generated events
nruns=$(sed 's#.*/##' "$OUT/embedding.list" | sort -u | wc -l)
[ "$nruns" = "$(wc -l < "$OUT/embedding.list")" ] || { echo "duplicate 2021 runs" >&2; exit 1; }

for name in data dataMB embedding2023 embedding; do
   echo "$name.list: $(wc -l < "$OUT/$name.list")"
done
