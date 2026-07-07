# cross_section_inverse — matrix-inversion cross section (second solver)

Dmitry's unregularized matrix inversion, the second unfolder alongside the main
Bayes pipeline in `../` (default). Same physical inputs — the etadet Miss/Fake
response, the data divided by the measured hybrid trigger correction
`C(pt)` (`../config.h::TrigEffMeas`), and each trigger restricted to its quote
window (`../config.h::QuoteLo`). **No environment flags, no non-physical
correction.** Everything is hardcoded in `../config.h`.

## Files

* `unfold.cxx` — builds the **square** etadet response (reco axis == truth axis
  == McBins), which the inversion needs to stay well-conditioned. Writes
  `response_<T>_R0.5_square.root`.
* `cross_section.cpp` — reads the square response, divides the data by `C(pt)`,
  unfolds with `RooUnfoldInvert`, normalizes, restricts to the quote window,
  compares to Dmitry's Table III. Writes `xsec_<T>_R0.5_invert.root` and
  `comparison_with_dmitriy_R0.5_invert.pdf`.
* `run_inverse.sh` — container driver.

## Run

```bash
bash run_inverse.sh              # square response + inversion (all triggers)
bash run_inverse.sh response     # square response only
bash run_inverse.sh cross_section
```

To change triggers / paths / thread count, edit `../config.h` (shared with the
Bayes pipeline).

## Note on the turn-on

The raw unregularized inverse oscillates in the steep turn-on below ~13 GeV
(near-singular migration there). Those bins are below every quote window and are
dropped, so the quoted spectrum is the stable plateau region. Dmitry tames the
turn-on with floor-restricted cell-filtering, which is not reproduced here — use
this as a plateau cross-check of the Bayes result, which is the quoted default.
