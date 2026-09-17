# Stage-2: from the per-jet trees to the cross sections

Input: the Stage-1 trees in `../output/` (data: `merged_data_R<R>.root`, `merged_dataMB_R<R>.root`;
embedding: `merged_matching_R<R>.root`, `merged_matchingMB_R<R>.root` for the 2023 request,
`matching2021/` and `matchingMB2021/` per run for the 2021 request).
Output: `results/` (ignored by git), one folder per radius plus `Michal/`, `models/`, `figures/`.

## Folders, in the order they are used

| folder | what it does | driver |
|---|---|---|
| `inputs/` | luminosities, chain efficiency, the published table (tracked) | — |
| `eps_chain/` | measures the min-bias chain efficiency from the jet-patch data | `eps_chain/run_bits.sh` |
| `stage2/` | the engine: data spectra, responses, unfolding, systematic bands | called by the drivers below |
| `published/` | the published construction reproduced: exclusive levels, 2021 embedding, inversion, band vs Table I | `published/run.sh` |
| `analysis/` | the result: one trigger at a time, Bayes on the 2023 embedding, every radius, UE on/off, no-UE set, radius ratios | `analysis/run.sh` |
| `models/` | Pythia 6 / Pythia 8 / Herwig 7 particle-level spectra and the hadronisation correction | `models/run_models.sh` |
| `plots.sh` | every figure into `results/figures/` | `plots.sh` |

Shared headers at this level: `config.h` (paths, the published binning, bad runs),
`soft_reweight.h` and `vertex_reweight.h` (the two embedding event weights), `plot_style.h` (the figure style).
Every folder has its own README with its steps.

## The run, from scratch

```bash
bash new_ana/published/run.sh all 0.5                          # ~3 h: the quoted cross section and its band
for R in 0.5 0.2 0.3 0.4; do bash new_ana/analysis/run.sh all $R; done   # ~18 h: per trigger, bands, no-UE set, radius ratios
bash new_ana/models/run_models.sh all                          # 1-3 days: generators (host, LCG view)
bash new_ana/plots.sh                                          # ~2 min
```

The drivers start from the host and run ROOT inside the container (`../site.sh` holds the image path).
Every step skips work whose output already exists, so a driver can be rerun after a crash.
`eps_chain` is rerun only after a new Stage-1 data production; its result is tracked in `inputs/`.

## Names

Result files carry the *mode* of the unfolding: `<trigger>[_e23][_noue][_mbins][_<variant>][_bayes<n>]`.

- trigger: `jp` = the exclusive levels JP0 + JP1 + JP2 (the published construction), `jp1` / `jp2` = one
  jet-patch trigger alone, `ht2` = the high-tower trigger, `mb2023` / `mb2021m` = the min-bias level
  with the 2023 / 2021 min-bias response.
- `e23`: response from the 2023 embedding (default: 2021, per run). `noue`: raw jets at both levels
  (no underlying-event subtraction). `mbins`: result on the 5-60 GeV bins of the R_AA reference.
- variant: one member of the systematic band (`stage2/variants.h`). `bayes<n>`: Bayesian unfolding with
  n iterations (default: matrix inversion).

Files: `data_blocks_R<R>[_mb][_<variant>].root` (data), `response_blocks_R<R>[_e23|_mb2023|_mb2021m][_<variant>].root`
(embedding), `xsec_inversion_R<R>[_<mode>].root` (one unfolding), `xsec_<mode>_R<R>_syst[_bayes<n>].{root,txt}`
(the band), `ratio_<mode>_over_<mode>_R<R>...` and `rratio_<mode>_R<R>_over_R0.5...` (ratios with bands).
The `.txt` tables carry the quoted pT range in their header; bins outside it are in the `.root` only.
