# src: Stage-1, the jet finder

`bin/RunppAna` reads one TStarJetPico file and writes one tree per radius with every jet of every
event, all trigger bits and the trigger-simulator decisions, so that the trigger is chosen at
Stage-2. It is the validated production binary: change nothing without checking a test file
against the previous production.

- `RunppAna.cxx` — the program: opens the pico, loops over events, writes `ResultTree` (data) or
  `JetTree` / `JetTreeMc` inputs of the matching (embedding); the output branches are declared here.
- `ppAnalysis.cxx/.hh` — the event loop: track and tower selection, anti-kT with active area, the
  off-axis-cone underlying-event density, the jet-to-patch and jet-to-tower matching, the DSM ADCs.
- `ppParameters.hh` — the command-line parameters (radius, floors, constituent cuts).
- `JetAnalyzer.cxx/.hh`, `JetQAHistogramManager.cxx/.hh` — the fastjet wrapper and the QA histograms.

Built inside the container by `make` (`run_production.sh build`). Driven by `../container.sh`
(one pico -> jets, modes `data`, `mbdata`, `embedding`, `mbembedding`); the embedding modes then run
`../macros/matching_mc_reco.cxx`, which pairs reconstructed and truth jets into `MatchedTree`
(misses and fakes included). Condor submission: `../run_production.sh` with `../submit/production.xml`.
