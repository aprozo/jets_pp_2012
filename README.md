# Inclusive jet cross section in pp at 200 GeV (STAR Run-12)

Anti-k_T inclusive jet cross section from the STAR Run-12 pp 200 GeV data set,
for radii R = 0.2-0.5 and |eta_det| < 1-R, compared bin by bin with the
published Table III (R = 0.5, |eta| < 0.5).

Stage-0 (MuDst -> `TStarJetPicoDst`) lives in the sister repository
`/gpfs01/star/pwg/prozorov/TStarJetPicoMaker`
(upstream: https://github.com/wsu-yale-rhig/TStarJetPicoMaker). This repository
is Stage-1 (picos -> per-jet trees, on the batch farm) and Stage-2 (trees ->
cross section, in `new_ana/`).

Triggers are an analysis filter, not a production split. One Stage-1 pass over
each pico writes every per-jet `trigger_match_JP0/JP1/JP2/HT2` bit and every
per-event `fired_JP0/JP1/JP2` bit into a single tree at one low jet-find floor
(4.8 GeV), so there is exactly one data production and one embedding
production; the trigger, and its analysis floor, are chosen at Stage-2. There
are no run-time switches: every physics choice is hardcoded in
`new_ana/config.h`, `container.sh` and `src/`.

## Environment

Everything downstream of the picos runs inside the analysis container (ROOT,
fastjet, RooUnfold, the TStarJetPico reader):

```bash
singularity exec -e -B /gpfs01 -B /gpfs/mnt/gpfs01 \
  /gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg bash
```

The Stage-2 drivers re-exec themselves into the container, so all commands
below are given from the host. `run_production.sh` must run on the host:
`star-submit-template` does not work inside the container.

File lists must use `/gpfs01/...` paths. The container binds `/gpfs01`; a
`/gpfs/mnt/gpfs01/...` path reads zero events without reporting an error.

## Layout

Every folder has a README.md with its steps; `new_ana/README.md` maps Stage-2.

```
src/, Makefile        Stage-1 jet finder -> bin/RunppAna
container.sh          Stage-1 worker: one pico -> jets; modes data | embedding | mbdata | mbembedding
run_production.sh     Stage-1 driver: build | submit | wait | merge (+ mb-*, mbembed-*, mbembed2021-*)
submit/production.xml, submit/production_mb2021.xml   condor templates (type, filelist, radii)
macros/matching_mc_reco.cxx     pairs embedding reco and truth jets -> MatchedTree
lib/eventStructuredAu_final/    TStarJetPico reader, bound over the container copy at run time
lists/jet_pico_dst/   data.list, dataMB.list, embedding.list (2021), embedding2023.list
lists/                bad runs, bad tower, prescales, emb2021_events_per_run.txt
                      (generated events per 2021 pico, the sample normalisation)
output/               Stage-1 trees and merges (not tracked)

new_ana/              Stage-2
  config.h            shared paths, published binning, run selection
  soft_reweight.h, vertex_reweight.h   embedding event weights
  inputs/             lumi_zilong_full.root  per-run sampled luminosity per jet-patch trigger
                      lumi_VPDMB_true.root   per-run min-bias luminosity
                      eps_chain3_R<R>.root   min-bias chain efficiency (from eps_chain/)
                      jet_cross_section_publishedR0.5.root   the published Table III
  stage2/             the engine: data blocks, responses per variant, unfolding, band combination (lib.sh)
  published/          the published analysis reproduced: 2021 embedding, inversion, band vs Table I (run.sh, plot.C)
  analysis/           the result: every trigger, Bayes on 2023, all radii, UE on/off, no-UE set, radius ratios (run.sh, plot.C)
  models/             Pythia 6/8 and Herwig 7 particle-level spectra and C_had against the data (run_models.sh)
  eps_chain/          min-bias chain-efficiency measurement (run_bits.sh)
  results/R<R>/, deck/results/R<R>/   Stage-2 outputs (not tracked)
docs/talk/            ladder.tex, preamble.tex; figures written by plot_ladder.C
```

## How to run

Site paths (container image, LCG view, maker repo) live in `site.sh`; every
driver sources it. The image has no build recipe; copy it.

```bash
# 0. Stage-0 picos: TStarJetPicoMaker (see its CLAUDE.md)                       ~7 h data, ~6 h embedding
# 1. Stage-1 input lists from the Stage-0 production folders
bash scripts/make_pico_lists.sh submit/<date>/job_pp200_towers_2012 submit/<date>/job_pp200_mb_2012 \
     submit/<date>/job_emb2023_20235003_chunks10 submit/<date>/job_emb2021_20212001
# 2. Stage-1 trees (condor; one submit stage at a time, let the queue drain)      ~2 h + merge per wave
bash run_production.sh build
bash run_production.sh submit "0.2 0.3 0.4 0.5"; bash run_production.sh wait; bash run_production.sh merge "0.2 0.3 0.4 0.5"
bash run_production.sh mb-submit "0.2 0.3 0.4 0.5"; bash run_production.sh mb-wait; bash run_production.sh mb-merge "0.2 0.3 0.4 0.5"
bash run_production.sh mbembed-submit "0.2 0.3 0.4 0.5"; bash run_production.sh mbembed-wait; bash run_production.sh mbembed-merge "0.2 0.3 0.4 0.5"
star-submit-template -template submit/production.xml -entities type=embedding,filelist=lists/jet_pico_dst/embedding.list,radii=0.5
mv output/matching/matched_e21_* output/matching2021/                          # 2021 response
# 3. Min-bias chain efficiency (only after a new data Stage-1)                   ~3 h
bash new_ana/eps_chain/run_bits.sh
# 4. The quoted cross section (inversion, 2021 embedding) and the min-bias level   ~3 h
bash new_ana/published/run.sh all 0.5
# 5. Bands, trigger ratios, no-UE reference set                                  ~18 h
for R in 0.5 0.2 0.3 0.4; do bash new_ana/analysis/run.sh all $R; done
#    -> new_ana/results/R<R>/xsec_<trigger>_R<R>_syst_bayes4.txt, ratio_*, rratio_*, results/Michal/
# 6. Generators (host, LCG view)                                                 1-3 days
bash new_ana/models/run_models.sh all
# 7. Figures                                                                    ~2 min
bash new_ana/plots.sh                                        # -> new_ana/results/figures/
```

## Reproducing the result

### 0. Picos

Produced by `TStarJetPicoMaker`: `macros/makeTStarJetPico.cxx` for data
(configs `JPHT` and `MB`) and `macros/makeTStarJetPicoEmbedding.cxx` for the
2021 and 2023 embedding requests. Their output paths are the contents of
`lists/jet_pico_dst/*.list`.

### 1. Stage-1 trees (condor)

The second argument of every stage is the space-separated radius list; the
default is `0.5`. Run one submit stage at a time and let the queue drain in
between — each submit rebuilds `run12prod.zip`, and jobs that are still
starting die at unzip.

```bash
bash run_production.sh build                              # -> bin/RunppAna

bash run_production.sh submit "0.2 0.3 0.4 0.5"           # data + 2023 embedding
bash run_production.sh wait
bash run_production.sh merge  "0.2 0.3 0.4 0.5"           # -> output/merged_{data,matching}_R<R>.root

bash run_production.sh mb-submit "0.2 0.3 0.4 0.5"        # min-bias data
bash run_production.sh mb-wait
bash run_production.sh mb-merge  "0.2 0.3 0.4 0.5"        # -> output/merged_dataMB_R<R>.root

bash run_production.sh mbembed-submit "0.2 0.3 0.4 0.5"   # min-bias pass, 2023 embedding
bash run_production.sh mbembed-wait
bash run_production.sh mbembed-merge  "0.2 0.3 0.4 0.5"   # -> output/merged_matchingMB_R<R>.root

bash run_production.sh mbembed2021-submit 0.5             # min-bias pass, 2021 embedding
bash run_production.sh mbembed2021-wait                   # -> output/matchingMB2021/
```

`mbembed2021-submit` reuses the sandbox zip that `mbembed-submit` builds, and
both `mbembed` stages refuse to start into a non-empty output directory (their
merge has no de-duplication).

The 2021 jet-patch response (`output/matching2021/`, used by
`new_ana/published`) is the same template with the 2021 pico list:

```bash
star-submit-template -template submit/production.xml \
  -entities type=embedding,filelist=lists/jet_pico_dst/embedding.list,radii=0.5
```

`production.xml` routes every `matched_*` file to `output/matching/`, so this
pass and the 2023 pass cannot be in flight together; move the 2021 trees to
`output/matching2021/` afterwards.

### 2. Min-bias chain efficiency

`new_ana/inputs/eps_chain3_R<R>.root` is already in the working tree and usable
as it stands. Rebuild it only after a new Stage-1 data production:

```bash
bash new_ana/eps_chain/run_bits.sh [nchunks] [njobs]      # -> new_ana/eps_chain/bits4/
# then, inside the container, from new_ana/eps_chain/, per radius:
root -l -b -q 'vpd_trigger_eff_jp_pass2d.C("bits4","0.5")'
```

`run_bits.sh` copies the resulting `eps_chain3_R<R>.root` into `new_ana/inputs/`.

### 3. Stage-2

```bash
bash new_ana/published/run.sh all 0.5        # the published construction: inversion, 2021 embedding, band vs Table I
bash new_ana/published/run.sh mb 0.5         # the min-bias level, inversion, 2023 response (ladder rung)
bash new_ana/published/run.sh mb2021m 0.5    # the min-bias level with the 2021 response
bash new_ana/analysis/run.sh all <R>         # per trigger, Bayes 4 it., 2023 embedding: bands, ratios, no-UE set (R = 0.5 also the radius ratios)
bash new_ana/published/run.sh ladder 0.5     # every rung of the talk (after analysis/run.sh bayes 0.5)
bash new_ana/plots.sh                        # every figure -> new_ana/results/figures/
```

The engine (`new_ana/stage2/`) is shared: `common.h` (sample table, run selection,
tree reader, jet pairing), `build_data.C` (level spectra, nominal and data variants in
one pass), `build_resp.C` (one response per variant), `unfold.C` (inversion or Bayes;
the HT2 trigger-efficiency factor), `variants.h` (the single variant table),
`systematics.C` (band combination, ratios, radius ratios), `lib.sh` (the driver functions). Results are in `new_ana/results/`
(ignored): `R<R>/`, `Michal/`, `models/`, `figures/`.
