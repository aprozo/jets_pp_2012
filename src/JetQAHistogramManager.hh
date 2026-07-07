/*
   @file JetQAHistogramManager.hh
   @brief QA histogram manager for events, tracks, towers, and jets
*/

#ifndef JET_QA_HISTOGRAM_MANAGER_HH
#define JET_QA_HISTOGRAM_MANAGER_HH

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TString.h"

#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetVector.h"

#include "fastjet/PseudoJet.hh"

class JetQAHistogramManager {
public:
   void Init();

   void
   FillEvent(double vx, double vy, double vz, double vz_vpd, double refmult, int mult, int njets, float event_sum_pt);
   void FillJet(const fastjet::PseudoJet &jet, bool match_jp, bool match_ht, TList *tracksList);

   void Write(TFile *f, const TString &name = "QA_histograms");

   // Run-binned constituent QA: set the current run before FillJet so the per-run
   // profiles (track/tower vs run) get the right run id. Lets the run-period
   // structure (e.g. the jet-patch ADC / tower-energy calibration step) be seen
   // at FULL run coverage (the per-pico integrated QA only resolves dominant runs).
   void SetRun(int run) { current_run = run; }
   int current_run = 0;

   // Run-binned constituent profiles (run number axis, hadd-able)
   TProfile *run_tower_et;
   TProfile *run_track_pt;
   TProfile *run_track_dca;
   TProfile *run_track_sdca_xy;
   TProfile *run_jp_patch_adc;

   // Event
   TH1D *vx;
   TH1D *vy;
   TH1D *vz;
   TH1D *vz_vpd;
   TH1D *vz_diff;
   TH1D *refmult;
   TH1D *mult;
   TH1D *njets;
   TH1D *event_sum_pt;

   // Tracks
   TH1D *track_pt;
   TH1D *track_eta;
   TH1D *track_phi;
   TH1D *track_dca;
   TH1D *track_sdca_xy;
   TH2D *pt_sDCAxy_pos;
   TH2D *pt_sDCAxy_neg;
   TH2D *track_pt_eta;

   // Towers
   TH1D *tower_et;
   TH1D *tower_eta;
   TH1D *tower_phi;
   TH1D *tower_id;
   TH2D *tower_et_id;

   // Jets
   TH1D *jet_pt;
   TH1D *jet_eta;
   TH1D *jet_phi;
   TH1D *jet_nef;
   TH1D *jet_nconst;
   TH1D *jet_ptlead;
   TH1D *jet_ptlead_tower;
   TH2D *jet_ptlead_tower_vs_nconstituents;
   TH1D *jet_ptlead_track;
   TH2D *jet_ptlead_track_vs_nconstituents;
   TH2D *jet_pt_eta;
   TH2D *jet_pt_phi;
   TH2D *jet_nef_pt;
   TH2D *jet_eta_phi;

   // Triggered jets (JP or HT matched)
   TH1D *jet_pt_trig;
   TH1D *jet_eta_trig;
   TH1D *jet_phi_trig;
   TH1D *jet_nef_trig;
   TH1D *jet_nconst_trig;
   TH1D *jet_ptlead_trig;
   TH2D *jet_pt_eta_trig;
   TH2D *jet_pt_phi_trig;
   TH2D *jet_nef_pt_trig;
   TH2D *jet_eta_phi_trig;
};

#endif
