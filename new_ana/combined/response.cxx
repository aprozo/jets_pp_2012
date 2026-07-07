// Fine-binned promotion-combined response ingredients (JPX = JP0+JP1+JP2).
//
// One pass over the Stage-1 embedding merge fills, on Dmitry's fine grid:
//   A_fine          matched reco x truth migration, weight = meas_w
//                   (total_weight x per-run promotion prescale weight)
//   A_entries_fine  the same, UNWEIGHTED (the cell filter keys off raw counts)
//   A_xfine         the same, weight = total_weight (truth-side weighting:
//                   when the filter removes a cell, its truth is routed back
//                   to Miss with the TRUTH normalization, not the prescale one)
//   b_fine(+entries) ALL accepted reco (matched + fake), weight = meas_w
//   x_fine(+entries) ALL truth in |eta|<0.5 (matched + miss), weight = total_weight
//
// The exclusive category partition and per-category reco windows are identical
// to the data side (cross_section.cpp here) — that identity is what makes the
// shouldFire turn-on cancel in the unfold. No trigger-efficiency table is
// folded in: the residual hardware-vs-simulator data-side infidelity is
// corrected on the DATA by the measured C(pt) (config.h::TrigEffMeas).
//
// Writes <workdir>response_JPX_R<R>_fine.root. The cell filter, re-coarsening
// and inversion live in cross_section.cpp so floor/lambda studies never
// re-read the tree.

#include <ROOT/RDataFrame.hxx>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TSystem.h>

#include <iostream>
#include <string>

#include "promotion.h"
#include "../soft_reweight.h"
#include "../vertex_reweight.h"

using namespace CrossSectionConfig;
using namespace PromotionConfig;
AnalysisConfig cfg;

static void build_fine(const std::string &jetR)
{
   ROOT::EnableImplicitMT(kImtThreads);

   const TString inputFile = TString(Form("%smerged_matching_R%s.root", cfg.datapath.c_str(), jetR.c_str()));
   if (gSystem->AccessPathName(inputFile)) {
      std::cout << "Skipping JPX R=" << jetR << ": " << inputFile << " missing" << std::endl;
      return;
   }
   const char *MCPT = "mc_pt_corrected"; // paper truth axis (UE-subtracted particle pT)

   const auto runPs = LoadRunPrescales(cfg);
   std::cout << "[jpx] per-run prescales loaded: " << runPs.size() << " runs "
             << "(fallback ps0=" << kPs0Avg << " ps1=" << kPs1Avg << ")" << std::endl;

   auto w0run = [runPs](int rid) {
      auto it = runPs.find(rid);
      const double p0 = (it != runPs.end()) ? it->second.first : kPs0Avg;
      return 1.0 / p0;
   };
   auto w1run = [runPs](int rid) {
      auto it = runPs.find(rid);
      const double p0 = (it != runPs.end()) ? it->second.first : kPs0Avg;
      const double p1 = (it != runPs.end()) ? it->second.second : kPs1Avg;
      return 1.0 / p0 + 1.0 / p1 - 1.0 / (p0 * p1);
   };

   ROOT::RDataFrame rawJets(std::string("MatchedTree"), std::string(inputFile.Data()));

   // Exclusive per-EVENT category from the stored embedding shouldFire bits +
   // the per-category reco-jet gate (jp_match veto, detector acceptance, NEF,
   // window). No row-level eta pre-filter: the truth acceptance |mc_eta|<0.5
   // sits in the truth selections, the reco acceptance |reco_det_eta|<0.5 in
   // reco_in — a pair with only one side in acceptance contributes exactly
   // where it belongs (background via b, or miss via x).
   auto allJets =
      rawJets
         .Define("vertex_weight", [](double vz) { return VertexReweight::weight(vz); }, {"event_vz"})
         .Define("soft_weight", [](double pthat) { return SoftReweight::weight(pthat); }, {"pthat_mid"})
         .Define("total_weight", "mc_weight * vertex_weight * soft_weight")
         .Define("rcat", "evt_should_JP2 ? 2 : (evt_should_JP1 ? 1 : (evt_should_JP0 ? 0 : -1))")
         .Define("reco_in", "reco_pt > -500 && std::abs(reco_det_eta) < 0.5 && "
                            "reco_neutral_fraction < 0.95 && "
                            "((rcat==2 && reco_trigger_match_JP2 && reco_pt >= 8.4) || "
                            "(rcat==1 && reco_trigger_match_JP1 && reco_pt > 8.2) || "
                            "(rcat==0 && reco_trigger_match_JP0 && reco_pt < 22.5))")
         .Define("w0r", w0run, {"runid"})
         .Define("w1r", w1run, {"runid"})
         .Define("wcat", "rcat==2 ? 1.0 : (rcat==1 ? w1r : (rcat==0 ? w0r : 0.0))")
         .Define("meas_w", "total_weight * wcat");

   const char *truthSel = "mc_pt > -500 && std::abs(mc_eta) < 0.5";
   auto matched  = allJets.Filter(Form("(reco_in && rcat >= 0) && %s", truthSel), "matched (A)");
   auto recoAll  = allJets.Filter("reco_in && rcat >= 0", "reco-all (b)");
   auto truthAll = allJets.Filter(truthSel, "truth-all (x)");

   auto hA_w = matched.Histo2D({"A_fine", ";reco;mc", kNRecoFine, kRecoFineLo, kRecoFineHi,
                                kNMcFine, kMcFineLo, kMcFineHi},
                               "reco_pt", MCPT, "meas_w");
   auto hA_e = matched.Histo2D({"A_entries_fine", ";reco;mc", kNRecoFine, kRecoFineLo, kRecoFineHi,
                                kNMcFine, kMcFineLo, kMcFineHi},
                               "reco_pt", MCPT);
   auto hA_x = matched.Histo2D({"A_xfine", ";reco;mc", kNRecoFine, kRecoFineLo, kRecoFineHi,
                                kNMcFine, kMcFineLo, kMcFineHi},
                               "reco_pt", MCPT, "total_weight");
   auto hB_w = recoAll.Histo1D({"b_fine", "", kNRecoFine, kRecoFineLo, kRecoFineHi}, "reco_pt", "meas_w");
   auto hB_e = recoAll.Histo1D({"b_entries_fine", "", kNRecoFine, kRecoFineLo, kRecoFineHi}, "reco_pt");
   auto hX_w = truthAll.Histo1D({"x_fine", "", kNMcFine, kMcFineLo, kMcFineHi}, MCPT, "total_weight");
   auto hX_e = truthAll.Histo1D({"x_entries_fine", "", kNMcFine, kMcFineLo, kMcFineHi}, MCPT);

   // Partition QA: exclusive category closure on the reco side.
   auto nCat2 = allJets.Filter("reco_in && rcat==2").Count();
   auto nCat1 = allJets.Filter("reco_in && rcat==1").Count();
   auto nCat0 = allJets.Filter("reco_in && rcat==0").Count();

   const TString outName = Form("%sresponse_JPX_R%s_fine.root", cfg.workdir.c_str(), jetR.c_str());
   TFile fout(outName, "RECREATE");
   hA_w->Write();
   hA_e->Write();
   hA_x->Write();
   hB_w->Write();
   hB_e->Write();
   hX_w->Write();
   hX_e->Write();
   TNamed("prescale_source", Form("run_prescales.txt (%zu runs), fallback ps0=%g ps1=%g",
                                  runPs.size(), kPs0Avg, kPs1Avg))
      .Write();
   fout.Close();

   std::cout << "[jpx] response ingredients written: " << outName << std::endl;
   std::cout << "[jpx] reco rows  cat2=" << *nCat2 << "  cat1=" << *nCat1 << "  cat0=" << *nCat0
             << std::endl;
}

void response()
{
   for (const auto &jetR : cfg.jetRs)
      build_fine(jetR);
}
