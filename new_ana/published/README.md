# published: the published analysis reproduced

The construction of the published Run-12 cross section (R = 0.5, |eta| < 0.5): exclusive jet-patch
levels, per-run luminosity, the 2021 embedding as response, matrix inversion on the published bins,
and the published systematic recipe. Compared with Table III (cross section) and Table I (band).

```bash
bash new_ana/published/run.sh all 0.5      # data, response, unfold, band    -> results/R0.5/
bash new_ana/published/run.sh mb 0.5       # the min-bias level with the 2023 min-bias response
bash new_ana/published/run.sh mb2021m 0.5  # the min-bias level with the 2021 min-bias response
bash new_ana/published/run.sh ladder 0.5   # every rung of the talk (needs analysis/run.sh bayes first)
```

## Steps of `run.sh all`

1. `data` — the level spectra of the jet-patch data (`stage2/build_data.C`).
2. `response` — the 2021 responses, nominal and one per response variant, four builds in parallel
   (`stage2/build_resp.C` over the per-run trees of `output/matching2021/`).
3. `unfold` — the cross section by matrix inversion: `xsec_inversion_R0.5.root`.
4. `band` — every member unfolded and combined: `xsec_jp_R0.5_syst.{root,txt}`.

## Figures (`bash new_ana/plots.sh` -> `results/figures/published_*.pdf`)

- `plot.C`: the cross section over Table III with the published band; the ladder (inversion, Bayes
  on 2021, Bayes on 2023); the band components against Table I.
- `plot_ladder.C`: the figures of `docs/talk/ladder.tex`, one ingredient changed per figure
  (method, embedding, trigger sample, min-bias); copies go to `docs/talk/`.
- `embedding_diff.C` + `embedding_diff_plot.C`: the two embeddings at fixed particle pT (energy
  scale, resolution, matching efficiency, fragmentation, vertex, background density). The detector
  level is the same; the 2021 record has pi0, eta and Sigma0 decayed, the 2023 record keeps them,
  which makes the 2021 particle jets 0.2-0.9 % softer.
- `embedding_particles.C`: the species composition of the generated particles of one pico.
