# stage2: the engine

Five macros, one shared header and one shell library. They are not run by hand; `published/run.sh` and `analysis/run.sh`
call them through `lib.sh`. Every macro is a ROOT macro compiled with ACLiC inside the container.

## Steps

1. `build_data.C` — detector-level jet spectra of the data, per trigger level and eta block, from the
   merged Stage-1 tree. One pass writes the nominal spectra and every data variant of the band
   (`data_blocks_R<R>[_mb][_<variant>].root`, 550 bins of 0.1 GeV on 5-60 GeV).
   Levels: `jp0`, `jp1`, `jp2` (exclusive: fired by the hardware, decided by the trigger simulator,
   and not decided by the next higher patch), `jp1i`, `jp2i` (inclusive single triggers), `ht2`,
   and in the min-bias pass `mb`. Jets: matched to the firing patch, |eta| < 0.9 at detector and
   physics level, neutral fraction <= 0.95, raw-pT floor 6.0 (JP1) / 8.4 (JP2) GeV.
   Events: |vz| < 60 cm, run in the luminosity table and not in the bad-run lists.
2. `build_resp.C` — the response of the embedding, per level and eta block: particle spectrum
   (every generated event), detector spectrum, and the (detector, particle) matrix of the one-to-one
   pairs within dR < 0.2. Each pt-hat sample is weighted by sigma / N_generated (times the soft
   reweight and the vertex reweight), and stored separately so that `unfold.C` can apply the outlier
   filter per sample. `build_resp` reads the per-run 2021 trees (restricted to the runs kept in the
   data); `build_resp_merged` reads the merged 2023 tree. The response variants of the band
   (energy scale, thresholds, underlying event) are separate builds.
3. `unfold.C` — one cross section from one data file and one response file:
   the levels are summed with their luminosities, the outlier filter of the published analysis is
   applied per sample, the fine axes are rebinned to the analysis bins, and the result is obtained
   by matrix inversion (`inv`) or RooUnfoldBayes (`bayes`, n iterations). The HT2 level carries the
   trigger-efficiency factor data / embedding. Writes `xsec_inversion_R<R>[_<mode>].root`
   with `canonical` (|eta| < 0.5, pb/GeV) and its covariance.
4. `systematics.C` — the band: `combine` reads the nominal and every member of every component,
   forms the envelope per component, adds up and down in quadrature and writes the `_syst` file
   and table; `ratios` and `ratios_radii` do the same for trigger ratios and radius ratios, each
   member recomputed on both numerator and denominator so common shifts cancel.
5. `variants.h` — the single table of the band: member names, what each changes, and the components.

`common.h` — what the macros share: the pt-hat sample table and weights, the run selection (bad runs,
runs with a luminosity entry), the level table and detector windows, the fine binning, the matched-tree
reader and the one-to-one jet pairing (dR < 0.2).

`lib.sh` — the shell functions of the drivers: `data_blocks`, `resp_2021`, `resp_2023`,
`resp_2021_mb`, `unfold_one`, `unfold_all`, `combine`. Each skips its step when the output exists
(unfolds: when newer than the blocks). Logs go to `results/R<R>/systematics_logs/`.

## The band (as published)

| component | members | envelope |
|---|---|---|
| energy scale | detector pT shifted by ±pT √(((1−R_T)·0.011)² + (R_T·0.032)²); DSM thresholds ±1 ADC in the jet-to-patch match; both together | max / min of six |
| tracking | 1 % of the tracks removed from the data jets | one-sided |
| underlying event | subtraction × 0.86 / 1.18 on data, response, both | max / min of six |
| embedding statistics | 1000 Poisson replicas of the response, re-unfolded | RMS, symmetric |
| regularisation (Bayes only) | 2 and 6 iterations; prior tilted by (pT/20)^±0.5 | max / min |
| trigger efficiency (HT2 only) | ±3 % on the data / embedding factor | flat |

Luminosity (5.6 %) is not in the band.
