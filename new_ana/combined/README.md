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
3. **Trigger correction.** Each data category is divided by its measured
   C(pt) = R × T-hat (`config.h::TrigEffMeas`) — the in-situ hardware→simulator
   ruler bridge, so the hardware-gated data matches the simulator-gated
   response. cat0 (JP0) has no measured correction (C = 1).
4. **Filtered fine response.** The migration is filled on Dmitry's fine grid
   (0.1 GeV cells, detector 550×[5,60], particle 600×[0,60]); cells whose raw
   pair count in a ±1 GeV box is ≤ 4 are removed (with the consistent b/x
   subtractions, truth-weighted via `A_xfine`); everything is then rebinned to
   the square McBins grid and
   `M_ij = (b_i/matched_i) · A_ij/x_j`
   folds fakes (b/matched ≥ 1) and matching+trigger efficiency (A/x) into one
   matrix.
5. **Floor-restricted inversion.** The square block [9.7, 52) is inverted
   directly (`x = M⁻¹ b`); matched content with truth outside the block
   (buffer feed-down, below-floor feed-up) is background-scaled away by the
   b/matched row factor. Optional scale-matched second-difference Tikhonov
   damping (`promotion.h::kTikhonovLambda`, 0 = Dmitry's plain inverse).

## Files

* `promotion.h` — every promotion-specific constant (windows, fine grid, box
  filter, block floor, lambda, per-run prescales, runtime Leff). Shared
  physics (bins, C(pt), badRuns, paths) comes from `../config.h`.
* `response.cxx` — one pass over `merged_matching_R0.5.root` → fine
  ingredients `response_JPX_R0.5_fine.root` (A, A_entries, A_x, b, x).
* `cross_section.cpp` — data partition + C(pt) + sum, cell filter + rebin +
  block inversion + covariance, normalization, Table III comparison. Writes
  `xsec_JPX_R0.5.root` + `comparison_with_dmitriy_R0.5_JPX.pdf`.
* `run_jpx.sh` — container driver.

## Run

```bash
bash run_jpx.sh                  # fine response + cross section
bash run_jpx.sh cross_section    # filter/floor/lambda studies (no tree re-read)
```

## Notes

* The quoted combination starts at the block floor (9.7 GeV by default). The
  6.9–9.7 bins are covered by cat0 but are prescale-deep; lowering
  `kJpxFloor` is a study, not the default.
* Quoted errors are data-statistical (X = R·diag(σ_b²)·Rᵀ). The
  embedding-statistical component is small (55M matched rows) and enters the
  systematics as the response-variation band instead.
