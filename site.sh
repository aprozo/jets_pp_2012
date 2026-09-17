# Site configuration, sourced by every driver (run_production.sh, the Stage-2 run.sh scripts,
# scripts/*.sh). These three paths are the only thing to edit for a new account or host.
#   source site.sh
SIMG=/gpfs01/star/pwg/prozorov/jets_pp_2012/star_star.simg     # Stage-1/2 container: ROOT, fastjet, RooUnfold, TStarJetPico reader
LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt  # generators for new_ana/models (host, bash)
MAKER_REPO=/gpfs01/star/pwg/prozorov/TStarJetPicoMaker         # Stage-0 repository, source of the pico lists
