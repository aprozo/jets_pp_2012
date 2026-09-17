# analysis: the result

One trigger at a time (MB, JP1, JP2, HT2), Bayesian unfolding with four iterations on the 2023
embedding, the systematic band of `stage2/variants.h`, at R = 0.2, 0.3, 0.4, 0.5; the ratios
between triggers; the underlying event on and off; the no-UE reference set on the 5-60 GeV bins;
the radius ratios with a band.

```bash
bash new_ana/analysis/run.sh all 0.5       # bayes, michal, radii
bash new_ana/analysis/run.sh all 0.2       # bayes, michal (radii only at R = 0.5)
```

## Steps

1. `bayes` — data blocks (jet-patch and min-bias), the 2023 responses (nominal and variants), then
   per trigger every member unfolded and combined:
   `xsec_<trigger>_R<R>_syst_bayes4.{root,txt}`, and the ratios to JP1
   `ratio_<trigger>_over_jp1_e23_R<R>_bayes4.{root,txt}`.
2. `michal` — the no-UE jet definition (raw jets at both levels): responses without the UE
   subtraction, JP1 on the 5-60 GeV bins with its band, the combined levels and the min-bias level
   for comparison (`michal.C` -> `results/Michal/xsec_R<R>_noUE_bins5-60.{root,txt}`), and the
   table of the cross section with and without the subtraction (`ue_onoff.C`).
3. `radii` — sigma(R) / sigma(0.5) for R = 0.2, 0.3, 0.4 with the band, for the nominal and the
   no-UE definition (`results/R<R>/rratio_*`). Needs steps 1-2 at every radius first.

## Quoted ranges

JP1 from 8.2 GeV, JP2 from 9.7, HT2 from 11.5, min-bias to 22.5 GeV; on the 5-60 GeV bins JP1, JP2
and HT2 from 10 GeV and min-bias to 20 GeV. The `.txt` tables carry the range in their header.

## Figures (`bash new_ana/plots.sh` -> `results/figures/analysis_*.pdf`)

`plot.C`: every trigger against Table III; the ratios to JP1 with bands; the radius ratios with
bands and PYTHIA 6; the generators against the radius ratios; UE on/off per radius; the no-UE set;
the generators against JP1 (UE on and off) at every radius; the hadronisation correction.
