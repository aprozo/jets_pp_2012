# Jets in pp 200 GeV, Run-12 — inclusive jet cross section

Anti-k_T R=0.5, |eta_det|<0.5 inclusive jet cross section from STAR pp200 Run-12,
validated per trigger (JP0, JP1, JP2, HT2) against Dmitry Kalinkin's Table III.
Reads `TStarJetPicoDst` trees (https://github.com/wsu-yale-rhig/TStarJetPicoMaker).
Two unfolding solvers: RooUnfoldBayes (nIter=2, default) and Dmitry-style
unregularized matrix inversion (cross-check).

The pipeline is **flag-free**: no `JETS_*` environment variables anywhere;
every physics choice is hardcoded in `new_ana/config.h` and edited there. This
is deliberate — the same pipeline must run at radii Dmitry never published
(e.g. R=0.4) as a standalone physics measurement, so nothing in it may depend
on tuning against a reference.

Everything downstream of the picos runs inside `star_star.simg`
(`/gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg`).

## The one idea

**Triggers are an analysis filter, not a production split.** One Stage-1 pass
over each pico writes every per-jet `trigger_match_JP0/JP1/JP2/HT2` bit and every
per-event `fired_JP0/JP1/JP2` bit into a single tree, at one low jet-find floor
(4.8 GeV). So there is exactly **one data production and one embedding
production**; each trigger is selected — and its analysis floor
(`config.h::TrigPtFloor`) applied — at Stage-2.

## Layout

```
src/                    RunppAna Stage-1 jet finder (hardware-first trigger match)
container.sh            Stage-1 worker: one pico -> all radii, all triggers
submit/production.xml   one condor template (type = data | embedding)
run_production.sh       Stage-1 driver: build -> submit -> wait -> merge
macros/matching_mc_reco.cxx   builds MatchedTree (embedding reco<->truth)
lists/jet_pico_dst/     data.list, embedding.list  (the latest picos)

new_ana/
  config.h              binning, floors/quote windows, C(pt)=RxThat tables,
                        Systematic variation presets — the single source of truth
  corrections/          hw_ratio.C (R), measure_T.C (That), measure_Cjpx.C (C_JPX)
  unfolding/unfold.cxx  per-trigger Miss/Fake response (decoupled det-eta gate)
  cross_section.cpp     per-trigger: raw -> /C -> Bayes -> normalize -> xsec_<T>.root
  cross_section_inverse/  per-trigger cross-check: square response + RooUnfoldInvert
  combined/             THE COMBINED RESULT: JPX promotion (JP0+JP1+JP2) unfolded
                        once by cell-filtered fine-response matrix inversion
  plot_alltriggers.C    overlay every trigger + JPX vs Dmitry's Table III
  run.sh                Stage-2 driver: response -> cross_section -> plot
systematics/            config-as-code variation driver + envelope band builder
```

## Stage-1 — produce the trees

```bash
./run_production.sh all        # build -> submit data+embedding -> wait -> merge
```

Produces `output/merged_data_R0.5.root` (ResultTree) and
`output/merged_matching_R0.5.root` (MatchedTree). `star-submit-template` runs on
the host; the build and merge run in the container. After any `src/` change,
`rm -f run12prod.zip *.package` so the scheduler ships a fresh binary.

## Stage-2 — cross section

```bash
bash new_ana/run.sh all                              # per-trigger Bayes
bash new_ana/combined/run_jpx.sh                     # JPX combination (matrix inversion)
bash new_ana/cross_section_inverse/run_inverse.sh    # per-trigger inversion cross-check
```

The full-range result is the JPX promotion combination (`new_ana/combined/`,
Dmitry's published method: exclusive JP0+JP1+JP2 shouldFire partition, one
matrix inversion on the cell-filtered fine response). The standalone triggers
(Bayes) are the per-region comparison and cross-check.

For each trigger in `config.h` it builds the Miss/Fake response, unfolds the
data (RooUnfoldBayes, nIter=2) after dividing by the measured trigger
correction C(pt)=R×That, normalizes by eta acceptance / bin width / per-run
luminosity (runtime Leff from `lumi_zilong_full.root` minus badRuns), restricts
to the trigger's quote window (`config.h::QuoteLo`), and writes
`xsec_<T>_R0.5.root` (`canonical` + `reference`). `plot_alltriggers.C` overlays
them against Dmitry's `jet_cross_section_dmitriyR0.5.root`. Outputs land in
`new_ana/`.

Quote windows (below them a standalone trigger only extrapolates its turn-on):
JP1 from 8.2 (low-pT workhorse), JP2/JP0 from 13.6, HT2 from 11.5 GeV.

## Corrections (derivation only)

The trigger correction C(pt)=R×That is frozen in `config.h::TrigEffMeas`. To
re-derive it: `bash new_ana/corrections/run.sh` (R from `hw_ratio.C`, That from
`measure_T.C`), then hand-edit the tables.

## Systematics

`bash systematics/run_systematics.sh` after the nominal run. The variation
matrix is `config.h::Systematics()` (typed presets, config-as-code); see
`systematics/README.md`.
