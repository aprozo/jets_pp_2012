# lists

| file | content | used by |
|---|---|---|
| `jet_pico_dst/data.list`, `dataMB.list` | the Stage-0 picos of the jet-patch/high-tower and the min-bias streams | `run_production.sh` |
| `jet_pico_dst/embedding2023.list`, `embedding.list` | the 2023 and the 2021 (per run) embedding picos | `run_production.sh` |
| `badrun_extras.list` | runs removed from data and embedding on top of `new_ana/config.h` | Stage-2 |
| `badtower_3407.list` | the one tower masked at Stage-1 (the DB status mask is applied at Stage-0) | `src/RunppAna.cxx` |
| `emb2021_events_per_run.txt` | generated events per (pt-hat sample, run) of the 2021 embedding: the sample normalisation | `stage2/build_resp.C` |
| `lum_perrun_VPDMB-nobsmd.txt` | per-run min-bias live time and prescale | `new_ana/eps_chain/` |

The pico lists are written by `scripts/make_pico_lists.sh` from the Stage-0 production folders.
Paths must start with `/gpfs01/`: the container binds that prefix only.
