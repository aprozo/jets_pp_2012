# jets_pp_2012 — STAR Run-12 pp 200 GeV inclusive jet cross sections

Stage-1 turns TStarJetPico trees into per-jet trees (data + PYTHIA-dijet
embedding); Stage-2 unfolds and normalises them into cross sections for five
triggers (MB = VPDMB-nobsmd, JP0/JP1/JP2, HT2) and four radii (anti-kT
R = 0.2/0.3/0.4/0.5, |eta_det| < 1-R). The reference is the published Table III
(R = 0.5 only). Stage-0 (MuDst -> pico) is the sister repo
`/gpfs01/star/pwg/prozorov/TStarJetPicoMaker/` (see its CLAUDE.md).

## Environments

```bash
# Stage-1/2 analysis (ROOT + fastjet + RooUnfold + TStarJetPico reader); the image path is in site.sh:
singularity exec -e -B /gpfs01 -B /gpfs/mnt/gpfs01 \
  /gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg bash

# MC generators (Pythia 8, Herwig 7, Rivet, LHAPDF, ROOT 6.34):
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh
```

- Run `run_production.sh` on the HOST: `star-submit-template` does not work
  inside the container. The Stage-2 drivers re-exec themselves into the
  container, so start them from the host too.
- File lists must use `/gpfs01/...` paths. The container binds `/gpfs01`; a
  `/gpfs/mnt/gpfs01/...` path reads zero events and reports no error.

## Commands

```bash
# Stage-1 (condor). Stages: build|submit|wait|merge|all, mb-*, mbembed-*,
# mbembed2021-*. The RADII are the SECOND argument, default "0.5".
bash run_production.sh build
bash run_production.sh submit "0.2 0.3 0.4 0.5"
bash run_production.sh wait
bash run_production.sh merge  "0.2 0.3 0.4 0.5"

# Stage-2, the published construction (inversion, 2021 embedding, band against Table I):
bash new_ana/published/run.sh [data|response|unfold|band|all|mb|mb2021m|ladder] [R]

# Stage-2, the analysis (every trigger, Bayes 4 it., 2023 embedding, all radii, UE on/off, no-UE set, radius ratios):
bash new_ana/analysis/run.sh [bayes|michal|radii|all] [R]

# Every figure (deck style) from the result files:
bash new_ana/plots.sh                                    # -> new_ana/results/figures/

# Generator comparison (Pythia 6/8, Herwig 7 vs the Bayes data; host, LCG view):
# every generator setting is written out at the top of its program (gen_*.cc, herwig.in); bin/ is built automatically
bash new_ana/models/run_models.sh [gen <generator> [N]|chad|spectra|compare|all|build]

# Min-bias chain efficiency (only after a new Stage-1 data production):
bash new_ana/eps_chain/run_bits.sh [nchunks] [njobs]          # -> eps_chain/bits4/
root -l -b -q 'vpd_trigger_eff_jp_pass2d.C("bits4","0.5")'    # -> eps_chain3_R0.5.root
```

## Layout

| path | role |
|---|---|
| `src/`, `Makefile`, `bin/RunppAna` | Stage-1 jet finder (validated production binary — do not break) |
| `container.sh` | Stage-1 worker: one pico -> jets; modes `data`, `embedding`, `mbdata`, `mbembedding` |
| `submit/production.xml`, `submit/production_mb2021.xml` | condor templates (entities `type`, `filelist`, `radii`) |
| `macros/matching_mc_reco.cxx` | pairs embedding reco and truth jets -> MatchedTree |
| `lib/eventStructuredAu_final/` | TStarJetPico reader, bound over the container copy at run time |
| `lists/` | pico lists (`jet_pico_dst/`), bad runs, bad tower, per-run prescales, `emb2021_events_per_run.txt` |
| `output/` | Stage-1 trees and merges: `merged_data(MB)_R<R>.root`, `merged_matching(MB)_R<R>.root` (2023 embedding), `matching2021/` (2021 embedding, read by `new_ana/published` as a per-file list; no merge), `matchingMB2021/` |
| `new_ana/` | Stage-2; `new_ana/README.md` maps the folders, the run order and the result-file names |
| `new_ana/config.h` | shared Stage-2 configuration: paths, published binning, run selection |
| `new_ana/inputs/` | luminosities, chain efficiency, the published table |
| `new_ana/stage2/` | the engine: `common.h` (shared tables, run selection, tree reader, jet pairing), `build_data.C`, `build_resp.C`, `unfold.C`, `variants.h`, `systematics.C`, `lib.sh` |
| `new_ana/published/` | the published analysis reproduced (`run.sh`: inversion, 2021 embedding, band; `plot.C`; `plot_ladder.C`; `embedding_diff*.C`) |
| `new_ana/analysis/` | the result (`run.sh`: Bayes per trigger, all radii, `michal.C`, `ue_onoff.C`; `plot.C`) |
| `new_ana/eps_chain/` | min-bias chain-efficiency measurement |
| `new_ana/models/` | particle-level Pythia6 / Pythia8 / Herwig7 spectra (embedding decays, pT-hat bins) against the data; trees in `output/models/`, results in `new_ana/results/models/` |
| `docs/talk/` | `ladder.tex` + `preamble.tex`; figures written by `plot_ladder.C` |

## Physics rules that must not drift

- Luminosity basis is the per-run sampled luminosity for EVERY trigger:
  `new_ana/inputs/lumi_zilong_full.root` for the jet-patch and high-tower
  triggers, `lumi_VPDMB_true.root` (ZDC x livetime/prescale, Leff = 0.161 pb^-1)
  for min-bias.
- NEVER take the min-bias luminosity as N_events / 25 mb. The effective
  VPDMB-nobsmd cross section is 3.28 mb; the 25 mb convention is a factor 7.6 off.
- The jet-patch triggers need NO hardware efficiency correction: kOnline data
  and embedding share the simulator decision (C = 1).
- eps_chain = P(VPDMB fired AND |vz_vpd - vz| < 6 cm | jet event), applied to the
  raw min-bias spectrum BEFORE unfolding. It is RADIUS-INDEPENDENT: both modules
  (`unfold.C::MinBiasLuminosity`) use the
  error-weighted JP0 constant 0.061. The 0.050-0.062 spread across radii and the
  ~10% JP0-JP1 sample difference are the min-bias systematic.
- The three eps_chain samples order JP0 > JP1 > JP2 (0.062 / 0.056 / 0.052 at
  R = 0.5): the harder the patch threshold, the further the event is from
  min-bias. JP2 collapses above 24 GeV on a handful of raw events — do not use it.
- Exclusive trigger levels are defined by hardware fired AND simulator should
  fire. Per-level detector windows are fixed: l0 < 22.5 GeV, l1 > 8.2 GeV,
  l2 > 9.7 GeV. Dropping the l2 window lets JP2 jets at their 8.4 GeV raw floor
  into the 8.2-9.7 GeV bin and see-saws bins 2/3 by +-4%.
- Levels 1 and 2 admit events recorded only by a lower prescaled trigger while
  being normalised by their own luminosity: this is the published convention and
  carries a +-1% normalisation systematic.
- Quoted errors are DATA STATISTICS ONLY: per-level counting noise
  B_ii = sum_l n_{l,i} / L_l^2 propagated through the inversion. They equal the
  published errors above 26 GeV and are 1.6-3.9x larger below 22.5 GeV, where the
  JP0 level (L0 = L1/48) dominates the variance.
- 2023 vs 2021 embedding: the 2023 sample has a 0.2-0.9% lower jet energy scale
  and a 1-3% lower low-pT efficiency, which is a 5% (9 GeV) to 2% (40 GeV) yield
  effect. It is NOT the run mix (restricting the 2023 response to the kept data
  runs changes the fold check by <= 0.7%) and more Bayes iterations do not remove
  it. `unfold.C` mode `jp_e23cal` calibrates the 2023 detector pT to the 2021
  energy scale per particle bin.
- MB jet-quality cuts are applied SYMMETRICALLY to data and response: NEF in
  (0.05, 0.90), zLead = ptLead/pt < 0.6 on the ANALYSIS (subtracted) jet pT
  (Stage-2), and the Stage-1 leading-track veto (`-leadsdca 0.5`).
- NEVER lower the MB reco floor below 8.2 GeV in the same unfolding pass. With a
  lower floor most MB jets sit below 8.2 GeV where data and embedding differ by
  10-28%; RooUnfoldBayes scales the fake vector by one global
  sum(data)/sum(MC measured) dominated by that soft block, and a 14 GeV truth jet
  reconstructing at 7 GeV counts as matched instead of missed. To quote 5-7 GeV,
  do it in a SEPARATE low-pT pass.
- `QuoteLo` is 13.6 GeV for every trigger in the deck. Below that the deck's
  detector-level closure is 1.11-1.22 for every trigger and JP1/JP2 are still in
  their turn-on; use the inversion (exclusive levels, 2021 embedding) there.
- Binning, floors and acceptance live in `new_ana/config.h` and `stage2/unfold.C` ONLY. MB uses the
  same bins and reco floor (8.2 GeV) as the jet-patch triggers; quoted MB bins
  above 22.5 GeV are merged.
- Hadronic correction subtracts only tracks that pass the ANALYSIS track list
  (flag, MuDst hit convention, 0.2 < pT < 200, |eta| < 2.5, DCA < 3 cm,
  pT-dependent DCA 2 -> 1 cm). It is implemented in the reader
  (`TStarJetPicoTowerCuts::HadronicCorrection`, maker repo, copied to
  `lib/eventStructuredAu_final`). Removing it costs ~4% of jets their soft tower
  constituents and doubles the bins-2/3 see-saw of the inversion.
- `build_data.C` requires a luminosity entry in the stream's own table (JP2 for the
  jet-patch mode, VPDMB for the min-bias mode). Runs before
  13047003 have none (hot tower 743 floods the MB stream with single-tower
  "jets" at 16-40 GeV); the bad-run lists alone leave 93 JP runs without
  luminosity.
- The response uses truth jets from EVERY generated event; the detector side is
  gated on a reconstructed vertex with |vz| < 60 cm within 5 cm of the thrown one.
  Each pt-hat sample is normalised by sigma x corr / N_generated
  (`lists/emb2021_events_per_run.txt`).

## Traps

1. Every `run_production.sh` submit stage does `rm -f run12prod.zip`. Rebuilding
   the sandbox while earlier condor jobs are still starting kills them at unzip:
   submit stages SEQUENTIALLY and wait for the queue to drain.
2. `mbembed-submit` refuses a non-empty `output/matchingMB`: that merge has no
   de-duplication. The JP data/MB data merges DO dedup by pico index (newest
   mtime) — keep it, repeated production rounds leave duplicates.
3. Data Stage-1 runs with `-geantnum 1`, so ResultTree (runid, eventid) =
   (adjusted `TString::Hash` of the pico BASENAME — it COLLIDES across files —,
   tree-entry index). Event-level joins must also match `runid1` (the real run).
4. The `Jets` branch (TStarJetVectorJet TClonesArray) has NO dictionary in the
   container: `SetBranchAddress` + `GetEntry` SEGFAULTS. Read kinematics through
   the streamer, e.g. `TTreeFormula("atan2(Jets.fP.fY,Jets.fP.fX)", tree)`.
5. MC constituents carry `sDCAxy = -999` sentinel, so every track-quality jet veto
   is forced OFF for `intype mcpico` (ppAnalysis MCPICO override block). Do not
   re-enable it — it silently kills half the truth jets.
6. `bg_density` and all jet-observable branches are PER-JET arrays aligned with
   `pt_corrected`. There is no per-event rho and no jet-phi array.
7. MatchedTree schema: `reco_pt` is the RAW jet pT and `reco_pt_corrected` the
   UE-subtracted one, symmetric with `mc_pt` / `mc_pt_corrected`;
   `reco_jet_area`, `reco_bg_density` and their `mc_` equivalents let the UE be
   varied on the response as well as on the data (`pt_corrected + f*area*rho`).
8. QA `vz_diff` / `vz_vpd` histograms are DOUBLE-FILLED for events passing the MB
   gate; per-event fractions must use `hEventCounter`, not histogram ratios.
9. TProfile pt-hat slices combine as sigma*N-weighted per-bin means — never
   `Scale()` + `Add()`. Count histograms scale by sigmaGen(mb) x 1e9 / nAccepted.
10. `n_constituents` includes fastjet AREA GHOSTS (~22 baseline): trends only.
11. Never ratio a merged bin against fine reference bins point-by-point (fake
    sawtooth) — integrate the reference over the wide bin.
12. Two concurrent `mb-merge` runs race on `scratch_merge_mb/` — one at a time.
13. `production.xml` routes every `matched_*` file to `output/matching/`. The 2021
    jet-patch pass and the 2023 pass therefore cannot be in flight together; move
    the 2021 trees to `output/matching2021/` before submitting the 2023 pass.

## Style

- Few files; header-only modules (no .hh/.cxx method duplication).
- No getters/setters; plain public structs; one long linear macro beats many
  small ones passing state around.
- Everything configurable sits in one `config.h` per deck; edit and rerun.
- Every study macro states its provenance and physics in a header comment.
- No Workflow-tool fanouts; delegate to targeted Opus agents only, without sub-agents.
- Documents contain physics language only — no code references.
- Change nothing without validating against the previous production on a test file.
- Statistical errors of the inversion: `canonical` carries the exact variance sum_l n_l/L_l^2 of the
  level-weighted data; `canonical_pubErr` carries the publication's estimator b^2/N_entries, which ignores the
  level weights (JP0 weight 48x JP1) and is 2.5-4.4x too small in the 8.2-22.5 GeV detector bins. With that
  estimator our bars reproduce the published ones in bins 3-12 (bins 1-2 come out
  1.4x / 1.25x larger); quote the exact ones.
- The 2021 and 2023 embeddings differ ONLY in the stored generator record: 2021 has pi0/eta/Sigma0 decayed,
  2023 keeps them undecayed (GEANT decays them, so the detector level is identical). With the 0.2 GeV constituent
  cut and the particle-level rho subtraction the 2021 truth jets are 0.2-0.9% softer, hence 2023 unfolds 2-4%
  higher. Which particle level is presented is a choice of the analysis; neither sample is altered.
  Diagnostics: new_ana/published/embedding_diff*.C, embedding_particles.C (the maker stores the PDG code of a
  generated track in its dEdx field).
- Every folder has a README.md with its steps; keep them current when a driver or an output name changes.
