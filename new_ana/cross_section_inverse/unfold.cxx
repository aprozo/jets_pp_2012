// Response build for the matrix-inversion pipeline. Same Miss/Fake, decoupled
// detector-eta OR-gate response as ../unfolding/unfold.cxx, but on a SQUARE grid
// (reco axis == truth axis == McBins) so Dmitry's unregularized inversion
// (cross_section.cpp here) stays well-conditioned. Writes
// <workdir>/response_<T>_R<R>_square.root. No environment flags.

#include "RooUnfoldResponse.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../config.h"
#include "../soft_reweight.h"
#include "../vertex_reweight.h"

using namespace CrossSectionConfig;
AnalysisConfig cfg;

static void single_unfold(const std::string &trigger, const std::string &jetR)
{
   ROOT::EnableImplicitMT(kImtThreads);

   const TString inputFile = TString(Form("%smerged_matching_R%s.root", cfg.datapath.c_str(), jetR.c_str()));
   if (gSystem->AccessPathName(inputFile)) {
      std::cout << "Skipping " << trigger << " R=" << jetR << ": " << inputFile << " missing" << std::endl;
      return;
   }

   const char *MCPT = "mc_pt_corrected"; // UE-subtracted particle pT (paper truth axis)
   gStyle->SetOptStat(0);

   const TString outName =
      Form("%sresponse_%s_R%s_square.root", cfg.workdir.c_str(), trigger.c_str(), jetR.c_str());
   TFile *fout = new TFile(outName, "RECREATE");

   // Square response: reco axis == truth axis == McBins.
   const std::vector<double> &bins = McBins();
   std::cout << "[response] " << trigger << " SQUARE (reco==truth): " << (bins.size() - 1) << " bins"
             << std::endl;

   ROOT::RDataFrame rawJets(std::string("MatchedTree"), std::string(inputFile.Data()));

   // Decoupled detector-eta gate: keep the row if the RECO jet is in the
   // detector acceptance (|det_eta|<0.5, == the data) OR the TRUTH is in the
   // particle acceptance (|y|<0.5).
   auto rawJetsCut = rawJets.Filter("((reco_pt > -500 && std::abs(reco_det_eta) < 0.5) || "
                                    "(mc_pt > -500 && std::abs(mc_eta) < 0.5))",
                                    "decoupled |det_eta|<0.5 (reco) OR |eta|<0.5 (truth)");

   auto allJets = rawJetsCut
      .Define("vertex_weight", [](double vz)    { return VertexReweight::weight(vz); },  {"event_vz"})
      .Define("soft_weight",   [](double pthat) { return SoftReweight::weight(pthat); }, {"pthat_mid"})
      .Define("total_weight",  "mc_weight * vertex_weight * soft_weight");

   const double ptMcCut = bins.front();
   const double recoFloor = std::max(TrigPtFloor(trigger), bins.front());

   const TString cutMatched =
      Form("%s > %f && reco_pt > %f && reco_trigger_match_%s && reco_neutral_fraction < 0.95 && "
           "abs(reco_det_eta) < 0.5",
           MCPT, ptMcCut, recoFloor, trigger.c_str());
   const TString cutMissed =
      Form("%s > %f && (reco_pt < %f || !reco_trigger_match_%s || reco_neutral_fraction >= 0.95 || "
           "abs(reco_det_eta) >= 0.5)",
           MCPT, ptMcCut, recoFloor, trigger.c_str());
   const TString cutFake =
      Form("mc_pt < -500 && reco_pt > %f && reco_trigger_match_%s && reco_neutral_fraction < 0.95 && "
           "abs(reco_det_eta) < 0.5",
           recoFloor, trigger.c_str());

   auto matched = allJets.Filter(cutMatched.Data(), "matched");
   auto missed  = allJets.Filter(cutMissed.Data(),  "missed");
   auto fakes   = allJets.Filter(cutFake.Data(),    "fakes");

   TH1D measAxis("Measured", ";p_{T}^{reco};dN/dp_{T}", (int)bins.size() - 1, bins.data());
   TH1D truthAxis("Truth", ";p_{T}^{mc};dN/dp_{T}", (int)bins.size() - 1, bins.data());

   auto *response = new RooUnfoldResponse(&measAxis, &truthAxis);
   response->SetName("my_response");
   response->UseOverflow(false);

   matched.Foreach(
      [response](double rpt, double mpt, double w) { response->Fill(rpt, mpt, w); },
      {"reco_pt", MCPT, "total_weight"});
   fakes.Foreach(
      [response](double rpt, double w) { response->Fake(rpt, w); },
      {"reco_pt", "total_weight"});
   missed.Foreach(
      [response](double mpt, double w) { response->Miss(mpt, w); },
      {MCPT, "total_weight"});

   const long long nMatched = *matched.Count();
   const long long nMissed  = *missed.Count();
   const long long nFakes   = *fakes.Count();
   std::cout << "Response built (" << trigger << " R=" << jetR << ", square): "
             << "matched=" << nMatched << "  missed=" << nMissed << "  fakes=" << nFakes << std::endl;

   fout->cd();
   response->Write();
   fout->Close();
}

void unfold()
{
   DefineCustomColors();
   for (const auto &trigger : cfg.triggers)
      for (const auto &jetR : cfg.jetRs)
         single_unfold(trigger, jetR);
}
