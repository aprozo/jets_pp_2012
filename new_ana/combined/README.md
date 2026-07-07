# combined — JPX promotion combination (matrix inversion)

Dmitry's published nominal: the JP0+JP1+JP2 **promotion combination**, unfolded
**once** by unregularized matrix inversion on a cell-filtered fine-binned
response. This is the method that covers the full 6.9–52 GeV range in one
spectrum — each pT region is measured by the trigger that is efficient there,
and the emulated turn-on cancels because the same shouldFire partition gates
the data and the response.

## The method

1. **Exclusive partition.** An event belongs to the highest jet-patch category
   it *should* fire: cat2 = should_JP2, cat1 = should_JP1 && !should_JP2,
   cat0 = should_JP0 && !should_JP1. In data the should bits are the OR of the
   per-jet **hardware** trigger_match bits; in embedding the stored
   `evt_should_JP*` (same jet-based convention). Per-category jet gates
   (identical both sides): the per-jet trigger_match veto, |det_eta|<0.5,
   NEF<0.95, and the windows cat2 ≥ 8.4, cat1 > 8.2, cat0 < 22.5 GeV.
2. **Prescale recovery.** Data enters RAW (events the prescaler recorded:
   cat2 fired_JP0|1|2, cat1 fired_JP0|1, cat0 fired_JP0) and is normalized by
   the **full** JP2 luminosity. The response's measured side instead carries
   the per-run sampling probability: w2 = 1, w1 = 1/ps0 + 1/ps1 − 1/(ps0·ps1),
   w0 = 1/ps0 (`lists/run_prescales.txt`).
3. **Trigger correction.** The summed data is divided by the MEASURED
   combination-level correction C_JPX(pt) = C_sum = eps_data/eps_emb
   (`promotion.h::JpxTrigEff`, derived by `corrections/measure_Cjpx.C`): the
   sampling-weighted category-gate probability ratio for the exact promotion
   gates, on the fired_JP0 data base vs the weighted embedding. Bin-by-bin in
   the turn-on, the [19,44) pol0 fit on the plateau. One self-contained
   in-situ measurement — the inclusive per-trigger R/T-hat tables do not see
   the exclusive-category shuffle and over-correct the turn-on by ~half.
4. **Filtered fine response.** The migration is filled on the fine 0.1 GeV
   grid covering the WHOLE analysis range (reco [5,86], truth [0,86] —
   Dmitry's [0,60] caps would truncate the 52-86 feed-down buffer column and
   swing the last quoted bin by ±8%); cells whose raw pair count in a ±1 GeV
   box is ≤ 4 are removed (with the consistent b/x subtractions,
   truth-weighted via `A_xfine`); everything is then rebinned to the square
   McBins grid and
   `M_ij = (b_i/matched_i) · A_ij/x_j`
   folds fakes (b/matched ≥ 1) and matching+trigger efficiency (A/x) into one
   matrix.
5. **Floor-restricted inversion.** The square block [8.2, 86) is inverted
   directly (`x = M⁻¹ b`) — the block extends exactly as far down as the
   MEASURED C_JPX (below 8.2 the fired_JP0 base cannot measure it and the
   6.9-8.2 bin rings). The 52-86 buffer is solved but never quoted;
   below-floor feed-up is background-scaled away by the b/matched row factor.
   Optional scale-matched second-difference Tikhonov damping
   (`promotion.h::kTikhonovLambda`, 0 = Dmitry's plain inverse; the "jpxDamp"
   systematic runs λ=0.030). A fold QA prints b_data/(M × Dmitry-truth) per
   reco row to separate data-side from solve-side residuals.

## Files

* `promotion.h` — every promotion-specific constant (windows via CatJetGate,
  fine grid, box filter, block floor, lambda, the frozen JpxTrigEff table,
  lumi-derived weights, runtime Leff). Shared physics (bins, badRuns, paths,
  Systematic presets) comes from `../config.h`.
* `response.cxx` — one pass over `merged_matching_R0.5.root` → fine
  ingredients `response_JPX_R0.5_fine.root` (A, A_entries, A_x, b, x).
* `cross_section.cpp` — data partition + sum + /C_JPX (the three
  pre-correction category histograms are cached in
  `data_JPX_R<R>_categories.root`; delete after any selection change), cell
  filter + rebin + block inversion + covariance, normalization, fold QA,
  Table III comparison. Study arguments: `cross_section("nominal", lambda,
  floor)`. Writes `xsec_JPX_R0.5.root` + `comparison_with_dmitriy_R0.5_JPX.pdf`.
* `plot_final.C` — the summary deck (spectrum, ratio with stat+syst band,
  per-trigger Bayes context, number table) → `jpx_final_R0.5.pdf`.
* `run_jpx.sh` — container driver (`run_jpx.sh [step] [systName]`).

## Run

```bash
bash run_jpx.sh                  # fine response + cross section
bash run_jpx.sh cross_section    # filter/floor/lambda studies (no tree re-read)
```

## Notes

* The quoted combination starts at the block floor (8.2 GeV — where the
  measured C_JPX begins). The 6.9–8.2 bin is covered by cat0 but the
  correction is unmeasurable there (prescale-starved base) and the bin rings;
  solving from 6.9 is a study (`cross_section("nominal", -1, 6.9)`), not the
  default.
* Quoted errors are data-statistical (X = R·diag(σ_b²)·Rᵀ). The
  embedding-statistical component is small (55M matched rows) and enters the
  systematics as the response-variation band instead.
