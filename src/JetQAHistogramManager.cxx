/*
   @file JetQAHistogramManager.cxx
   @brief QA histogram manager for events, tracks, towers, and jets
*/

#include "JetQAHistogramManager.hh"
#include "ppAnalysis.hh"

#include "TMath.h"
#include "TVector2.h"

#include "JetAnalyzer.hh"

#include <vector>

static double WrapPhi(double phi)
{
   return TVector2::Phi_mpi_pi(phi);
}

void JetQAHistogramManager::Init()
{
   // Event
   vx = new TH1D("vx", "Primary vertex x; vx, cm; N", 300, -0.5, 0.5);
   vy = new TH1D("vy", "Primary vertex y; vy, cm; N", 300, -0.5, 0.5);
   vz = new TH1D("vz", "Primary vertex z; vz, cm; N", 280, -70, 70);
   vz_vpd = new TH1D("vz_vpd", "Primary vertex z from VPD; vz, cm; N", 400, -100, 100);
   vz_diff = new TH1D("vz_diff", "vertex z VPD-TPC; vz, cm; N", 500, -30, 30);
   refmult = new TH1D("refmult", "Reference multiplicity; refmult; N", 150, 0, 150);
   mult = new TH1D("mult", "Selected multiplicity; mult; N", 150, 0, 150);
   njets = new TH1D("njets", "Jets per event; N_{jets}; N", 10, 0, 10);
   event_sum_pt = new TH1D("event_sum_pt", "Event sum p_{T}; sum p_{T} [GeV/c]; N", 100, 0, 100);

   // Run-binned constituent profiles: run number on the x axis, so they cover
   // every run and survive hadd.
   run_tower_et      = new TProfile("run_tower_et", "<tower E_{T}> vs run; run; <E_{T}> [GeV]", 35000, 13040000, 13075000);
   run_track_pt      = new TProfile("run_track_pt", "<track p_{T}> vs run; run; <p_{T}> [GeV]", 35000, 13040000, 13075000);
   run_track_dca     = new TProfile("run_track_dca", "<track DCA> vs run; run; <DCA> [cm]", 35000, 13040000, 13075000);
   run_track_sdca_xy = new TProfile("run_track_sdca_xy", "<track sDCA_{xy}> vs run; run; <sDCA_{xy}>", 35000, 13040000, 13075000);
   run_jp_patch_adc  = new TProfile("run_jp_patch_adc", "<jet-patch ADC> vs run; run; <ADC>", 35000, 13040000, 13075000);

   // Tracks
   track_pt = new TH1D("track_pt", "Track p_{T}; p_{T} [GeV/c]; N", 300, 0, 30);
   track_eta = new TH1D("track_eta", "Track #eta; #eta; N", 120, -1.2, 1.2);
   track_phi = new TH1D("track_phi", "Track #phi; #phi; N", 128, -TMath::Pi(), TMath::Pi());
   track_dca = new TH1D("track_dca", "Track DCA; DCA [cm]; N", 200, 0, 5);
   track_sdca_xy = new TH1D("track_sdca_xy", "Track sDCA_{xy}; sDCA_{xy} [cm]; N", 200, -2, 2);
   pt_sDCAxy_pos = new TH2D("pt_sDCAxy_pos", "sDCAxy positive tracks; sDCA, cm; pt, GeV/c; N", 100, -4, 4, 100, 0, 30);
   pt_sDCAxy_neg = new TH2D("pt_sDCAxy_neg", "sDCAxy negative tracks; sDCA, cm;  pt, GeV/c; N", 100, -4, 4, 100, 0, 30);
   track_pt_eta = new TH2D("track_pt_eta", "Track p_{T} vs #eta; p_{T} [GeV/c]; #eta", 300, 0, 30, 120, -1.2, 1.2);

   // Towers
   tower_et = new TH1D("tower_et", "Tower E_{T}; E_{T} [GeV]; N", 300, 0, 30);
   tower_eta = new TH1D("tower_eta", "Tower #eta; #eta; N", 120, -1.2, 1.2);
   tower_phi = new TH1D("tower_phi", "Tower #phi; #phi; N", 128, -TMath::Pi(), TMath::Pi());
   tower_id = new TH1D("tower_id", "Tower ID; ID; N", 4801, 0, 4801);
   tower_et_id = new TH2D("tower_et_id", "Tower E_{T} vs id; E_{T} [GeV]; #eta", 300, 0, 30, 4801, 0, 4801);

   // Jets
   jet_pt = new TH1D("jet_pt", "Jet p_{T}; p_{T} [GeV/c]; N", 200, 0, 100);
   jet_eta = new TH1D("jet_eta", "Jet #eta; #eta; N", 120, -1.2, 1.2);
   jet_phi = new TH1D("jet_phi", "Jet #phi; #phi; N", 128, -TMath::Pi(), TMath::Pi());
   jet_nef = new TH1D("jet_nef", "Jet NEF; NEF; N", 100, 0, 1.0);
   jet_nconst = new TH1D("jet_nconst", "Jet constituents; N_{const}; N", 200, 0, 200);
   jet_ptlead = new TH1D("jet_ptlead", "Leading constituent p_{T}; p_{T} [GeV/c]; N", 120, 0, 30);
   jet_ptlead_tower = new TH1D("jet_ptlead_tower", "Leading tower p_{T}; p_{T} [GeV/c]; N", 300, 0, 30);
   jet_ptlead_tower_vs_nconstituents =
      new TH2D("jet_ptlead_tower_vs_nconstituents", "Leading tower p_{T} vs N_{constituents}; p_{T} [GeV/c]; N", 300, 0,
               30, 30, 0, 30);
   jet_ptlead_track = new TH1D("jet_ptlead_track", "Leading track p_{T}; p_{T} [GeV/c]; N", 300, 0, 30);
   jet_ptlead_track_vs_nconstituents =
      new TH2D("jet_ptlead_track_vs_nconstituents", "Leading track p_{T} vs N_{constituents}; p_{T} [GeV/c]; N", 300, 0,
               30, 30, 0, 30);
   jet_pt_eta = new TH2D("jet_pt_eta", "Jet p_{T} vs #eta; p_{T} [GeV/c]; #eta", 200, 0, 100, 120, -1.2, 1.2);
   jet_pt_phi =
      new TH2D("jet_pt_phi", "Jet p_{T} vs #phi; p_{T} [GeV/c]; #phi", 200, 0, 100, 128, -TMath::Pi(), TMath::Pi());
   jet_nef_pt = new TH2D("jet_nef_pt", "Jet NEF vs p_{T}; p_{T} [GeV/c]; NEF", 200, 0, 100, 100, 0, 1.0);
   jet_eta_phi =
      new TH2D("jet_eta_phi", "Jet #eta vs #phi; #eta; #phi", 120, -1.2, 1.2, 128, -TMath::Pi(), TMath::Pi());

   jet_pt_trig = new TH1D("jet_pt_trig", "Triggered jet p_{T}; p_{T} [GeV/c]; N", 200, 0, 100);
   jet_eta_trig = new TH1D("jet_eta_trig", "Triggered jet #eta; #eta; N", 120, -1.2, 1.2);
   jet_phi_trig = new TH1D("jet_phi_trig", "Triggered jet #phi; #phi; N", 128, -TMath::Pi(), TMath::Pi());
   jet_nef_trig = new TH1D("jet_nef_trig", "Triggered jet NEF; NEF; N", 100, 0, 1.0);
   jet_nconst_trig = new TH1D("jet_nconst_trig", "Triggered jet constituents; N_{const}; N", 200, 0, 200);
   jet_ptlead_trig = new TH1D("jet_ptlead_trig", "Triggered leading constituent p_{T}; p_{T} [GeV/c]; N", 120, 0, 30);
   jet_pt_eta_trig =
      new TH2D("jet_pt_eta_trig", "Triggered jet p_{T} vs #eta; p_{T} [GeV/c]; #eta", 200, 0, 100, 120, -1.2, 1.2);
   jet_pt_phi_trig = new TH2D("jet_pt_phi_trig", "Triggered jet p_{T} vs #phi; p_{T} [GeV/c]; #phi", 200, 0, 100, 128,
                              -TMath::Pi(), TMath::Pi());
   jet_nef_pt_trig =
      new TH2D("jet_nef_pt_trig", "Triggered jet NEF vs p_{T}; p_{T} [GeV/c]; NEF", 200, 0, 100, 100, 0, 1.0);
   jet_eta_phi_trig = new TH2D("jet_eta_phi_trig", "Triggered jet #eta vs #phi; #eta; #phi", 120, -1.2, 1.2, 128,
                               -TMath::Pi(), TMath::Pi());
}

void JetQAHistogramManager::FillEvent(double vx_in, double vy_in, double vz_in, double vz_vpd_in, double refmult_in,
                                      int mult_in, int njets_in, float event_sum_pt_in)
{
   vx->Fill(vx_in);
   vy->Fill(vy_in);
   vz->Fill(vz_in);
   vz_vpd->Fill(vz_vpd_in);
   vz_diff->Fill(vz_vpd_in - vz_in);
   refmult->Fill(refmult_in);
   mult->Fill(mult_in);
   njets->Fill(njets_in);
   event_sum_pt->Fill(event_sum_pt_in);
}

void JetQAHistogramManager::FillJet(const fastjet::PseudoJet &jet, bool match_jp, bool match_ht, TList *tracksList)
{
   const double jet_pt_val = jet.perp();
   const int n_const = jet.constituents().size();
   double neutral_pt = 0.0;

   const std::vector<fastjet::PseudoJet> charged_constituents = fastjet::sorted_by_pt(OnlyCharged(jet.constituents()));
   const std::vector<fastjet::PseudoJet> neutral_constituents = fastjet::sorted_by_pt(OnlyNeutral(jet.constituents()));

   for (const auto &p : neutral_constituents) {
      neutral_pt += p.perp();

      const double et = p.Et();
      const int id = p.user_info<JetAnalysisUserInfo>().GetNumber();
      tower_et->Fill(et);
      tower_eta->Fill(p.eta());
      tower_phi->Fill(WrapPhi(p.phi()));
      tower_et_id->Fill(et, id);
      tower_id->Fill(id);
      run_tower_et->Fill(current_run, et);
   }
   for (const auto &p : charged_constituents) {
      const double pt = p.perp();
      track_pt->Fill(pt);
      track_eta->Fill(p.eta());
      track_phi->Fill(WrapPhi(p.phi()));
      track_pt_eta->Fill(pt, p.eta());

      int container_id = p.user_info<JetAnalysisUserInfo>().GetNumber();
      TStarJetPicoPrimaryTrack *track = (TStarJetPicoPrimaryTrack *)tracksList->At(container_id);
      track_dca->Fill(track->GetDCA());
      track_sdca_xy->Fill(track->GetsDCAxy());
      run_track_pt->Fill(current_run, pt);
      run_track_dca->Fill(current_run, track->GetDCA());
      run_track_sdca_xy->Fill(current_run, track->GetsDCAxy());
      if (track->GetCharge() > 0)
         pt_sDCAxy_pos->Fill(track->GetsDCAxy(), pt);
      else
         pt_sDCAxy_neg->Fill(track->GetsDCAxy(), pt);
   }

   // Leading constituent: is it charged or neutral? (GetQuarkCharge() is the
   // constituent charge in units of e/3.)
   const std::vector<fastjet::PseudoJet> sorted_constituents = fastjet::sorted_by_pt(jet.constituents());
   double pt_lead = sorted_constituents[0].perp();
   if (sorted_constituents[0].is_pure_ghost() ||
       !sorted_constituents[0].has_user_info<JetAnalysisUserInfo>())
      return;
   int charge = sorted_constituents[0].user_info<JetAnalysisUserInfo>().GetQuarkCharge();
   if (charge != 0) {
      jet_ptlead_track->Fill(pt_lead);
      jet_ptlead_track_vs_nconstituents->Fill(pt_lead, n_const);
   } else {
      jet_ptlead_tower->Fill(pt_lead);
      jet_ptlead_tower_vs_nconstituents->Fill(pt_lead, n_const);
   }

   const double nef = neutral_pt / jet.perp();

   jet_pt->Fill(jet.perp());
   jet_eta->Fill(jet.eta());
   jet_phi->Fill(WrapPhi(jet.phi()));
   jet_nef->Fill(nef);
   jet_nconst->Fill(n_const);
   jet_ptlead->Fill(pt_lead);

   jet_pt_eta->Fill(jet.perp(), jet.eta());
   jet_pt_phi->Fill(jet.perp(), WrapPhi(jet.phi()));
   jet_nef_pt->Fill(jet.perp(), nef);
   jet_eta_phi->Fill(jet.eta(), WrapPhi(jet.phi()));

   if (match_jp || match_ht) {
      jet_pt_trig->Fill(jet.perp());
      jet_eta_trig->Fill(jet.eta());
      jet_phi_trig->Fill(WrapPhi(jet.phi()));
      jet_nef_trig->Fill(nef);
      jet_nconst_trig->Fill(n_const);
      jet_ptlead_trig->Fill(pt_lead);
      jet_pt_eta_trig->Fill(jet.perp(), jet.eta());
      jet_pt_phi_trig->Fill(jet.perp(), WrapPhi(jet.phi()));
      jet_nef_pt_trig->Fill(jet.perp(), nef);
      jet_eta_phi_trig->Fill(jet.eta(), WrapPhi(jet.phi()));
   }
}

void JetQAHistogramManager::Write(TFile *f, const TString &name)
{
   if (!f)
      return;
   f->cd();
   f->mkdir(name);
   f->cd(name);
   vx->Write();
   vy->Write();
   vz->Write();
   vz_vpd->Write();
   vz_diff->Write();
   refmult->Write();
   mult->Write();
   njets->Write();
   event_sum_pt->Write();

   run_tower_et->Write();
   run_track_pt->Write();
   run_track_dca->Write();
   run_track_sdca_xy->Write();
   run_jp_patch_adc->Write();

   track_pt->Write();
   track_eta->Write();
   track_phi->Write();
   track_dca->Write();
   track_sdca_xy->Write();
   pt_sDCAxy_pos->Write();
   pt_sDCAxy_neg->Write();
   track_pt_eta->Write();

   tower_et->Write();
   tower_eta->Write();
   tower_phi->Write();
   tower_id->Write();
   tower_et_id->Write();

   jet_pt->Write();
   jet_eta->Write();
   jet_phi->Write();
   jet_nef->Write();
   jet_nconst->Write();
   jet_ptlead->Write();
   jet_ptlead_tower->Write();
   jet_ptlead_tower_vs_nconstituents->Write();  
   jet_ptlead_track->Write();
   jet_ptlead_track_vs_nconstituents->Write();
   jet_pt_eta->Write();
   jet_pt_phi->Write();
   jet_nef_pt->Write();
   jet_eta_phi->Write();

   jet_pt_trig->Write();
   jet_eta_trig->Write();
   jet_phi_trig->Write();
   jet_nef_trig->Write();
   jet_nconst_trig->Write();
   jet_ptlead_trig->Write();
   jet_pt_eta_trig->Write();
   jet_pt_phi_trig->Write();
   jet_nef_pt_trig->Write();
   jet_eta_phi_trig->Write();
}
