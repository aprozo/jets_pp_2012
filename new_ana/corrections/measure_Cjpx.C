// measure_Cjpx.C --------------------------------------------------------------
// Measure the JPX-combination trigger correction C_JPX(pt) = That_sum x R —
// the same hybrid philosophy as the per-trigger config.h::TrigEffMeas
// (turn-on SHAPE from data/embedding, plateau RULER from hardware/software),
// but for the EXACT gates the promotion combination uses.
//
// The combination's effective gate probability, per reco-pT bin:
//   eps(pt) = sum_cat s_cat * P(jet passes the cat gate | pt)
// with the sampling weights s = {w0, w1, 1} (the prescale-recovery weights).
//   - DATA: on the unbiased fired_JP0 base, hardware bits (rcat from the OR
//     of per-jet hardware trigger_match), per-jet category gates.
//   - EMBEDDING: on the evt_should_JP0 base (same jet-based convention),
//     simulator bits, response-side gates — WEIGHTED with the same
//     total_weight (mc x vz x soft) the response build uses. The embedding is
//     a pt-hat-stitched sample: unweighted shares are dominated by
//     over-represented high-pt-hat events and are meaningless.
//   C_sum(pt) = eps_data(pt) / eps_emb(pt)
//
// That_sum(pt) = C_sum(pt) / C_sum(plateau fit), clamped to 1 on the plateau:
// the residual data-vs-simulator turn-on SHAPE, which the response Misses
// cannot model. The un-normalized plateau level of C_sum is R x T_true;
// T_true (the physical data/embedding simulator-efficiency difference) is
// deliberately LEFT IN (Dmitry-faithful, carried as an embedding
// systematic); the plateau ruler applied instead is R (hardware/software
// within data, from hw_ratio.C) — smoothed to its plateau fit, because R is
// a smooth hardware property and per-bin noise would inject fake structure.
//
// Usage:
//   root -l -b -q 'measure_Cjpx.C+'
// Paste the printed C_JPX table into combined/promotion.h::JpxTrigEff.
// -----------------------------------------------------------------------------
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <cstdio>
#include <string>
#include <vector>

#include "../combined/promotion.h"
#include "../soft_reweight.h"
#include "../vertex_reweight.h"

using namespace PromotionConfig;
using CrossSectionConfig::AnalysisConfig;

void measure_Cjpx(
   const char *dataTree = "/gpfs01/star/pwg/prozorov/study_pp2012/fable5/jets_pp_2012/output/merged_data_R0.5.root",
   const char *matchedTree = "/gpfs01/star/pwg/prozorov/study_pp2012/fable5/jets_pp_2012/output/merged_matching_R0.5.root")
{
   ROOT::EnableImplicitMT(CrossSectionConfig::kImtThreads);
   const std::vector<double> QB = {6.9, 8.2, 9.7, 11.5, 13.6, 16.1, 19, 22.5, 26.6, 31.4, 37.2, 44, 52};
   const int nb = (int)QB.size() - 1;

   AnalysisConfig cfg;
   const PromotionWeights W = Weights(cfg);
   printf("[Cjpx] sampling weights: w0=%.6g w1=%.6g w2=1\n", W.w0, W.w1);

   // ---- DATA: fired_JP0 base; per-jet hardware category gates ---------------
   ROOT::RDataFrame d("ResultTree", dataTree);
   auto anyOf = [](const ROOT::VecOps::RVec<bool> &v) {
      for (auto b : v)
         if (b) return true;
      return false;
   };
   auto dd = d.Filter("fired_JP0")
                .Define("evt_sf0", anyOf, {"trigger_match_JP0"})
                .Define("evt_sf1", anyOf, {"trigger_match_JP1"})
                .Define("evt_sf2", anyOf, {"trigger_match_JP2"})
                .Filter("evt_sf0", "jet-based should_JP0 base")
                .Define("rcat_evt", "evt_sf2 ? 2 : (evt_sf1 ? 1 : (evt_sf0 ? 0 : -1))")
                .Define("sel", "abs(det_eta) < 0.5 && neutral_fraction <= 0.95")
                .Define("hit2", "sel && rcat_evt == 2 && " + CatJetGate(2, "", "pt_corrected"))
                .Define("hit1", "sel && rcat_evt == 1 && " + CatJetGate(1, "", "pt_corrected"))
                .Define("hit0", "sel && rcat_evt == 0 && " + CatJetGate(0, "", "pt_corrected"))
                .Define("ptAll", "pt_corrected[sel]")
                .Define("pt2", "pt_corrected[hit2]")
                .Define("pt1", "pt_corrected[hit1]")
                .Define("pt0", "pt_corrected[hit0]");
   auto hA = dd.Histo1D({"hA", "", nb, QB.data()}, "ptAll");
   auto h2 = dd.Histo1D({"h2", "", nb, QB.data()}, "pt2");
   auto h1 = dd.Histo1D({"h1", "", nb, QB.data()}, "pt1");
   auto h0 = dd.Histo1D({"h0", "", nb, QB.data()}, "pt0");

   // ---- EMBEDDING: response-side gates, response weights ---------------------
   ROOT::RDataFrame m("MatchedTree", matchedTree);
   auto mm = m.Filter("evt_should_JP0", "should_JP0 base (matches the data base)")
                .Filter("reco_pt > -500 && std::abs(reco_det_eta) < 0.5 && reco_neutral_fraction < 0.95")
                .Define("vertex_weight", [](double vz) { return VertexReweight::weight(vz); }, {"event_vz"})
                .Define("soft_weight", [](double pthat) { return SoftReweight::weight(pthat); },
                        {"pthat_mid"})
                .Define("total_weight", "mc_weight * vertex_weight * soft_weight")
                .Define("rcat", "evt_should_JP2 ? 2 : (evt_should_JP1 ? 1 : (evt_should_JP0 ? 0 : -1))")
                .Define("ehit2", "rcat == 2 && " + CatJetGate(2, "reco_", "reco_pt"))
                .Define("ehit1", "rcat == 1 && " + CatJetGate(1, "reco_", "reco_pt"))
                .Define("ehit0", "rcat == 0 && " + CatJetGate(0, "reco_", "reco_pt"));
   auto hEA = mm.Histo1D({"hEA", "", nb, QB.data()}, "reco_pt", "total_weight");
   auto hE2 = mm.Filter("ehit2").Histo1D({"hE2", "", nb, QB.data()}, "reco_pt", "total_weight");
   auto hE1 = mm.Filter("ehit1").Histo1D({"hE1", "", nb, QB.data()}, "reco_pt", "total_weight");
   auto hE0 = mm.Filter("ehit0").Histo1D({"hE0", "", nb, QB.data()}, "reco_pt", "total_weight");

   // ---- sampling-weighted sum efficiencies + C_sum ---------------------------
   gStyle->SetOptStat(0);
   TH1D *epsD = (TH1D *)h2->Clone("epsD");
   epsD->Add(h1.GetPtr(), W.w1);
   epsD->Add(h0.GetPtr(), W.w0);
   epsD->Divide(hA.GetPtr());
   TH1D *epsE = (TH1D *)hE2->Clone("epsE");
   epsE->Add(hE1.GetPtr(), W.w1);
   epsE->Add(hE0.GetPtr(), W.w0);
   epsE->Divide(hEA.GetPtr());
   TH1D *Csum = (TH1D *)epsD->Clone("Csum");
   Csum->Divide(epsE);

   // Plateau fit of C_sum (its level = R x T_true; only the SHAPE below the
   // plateau is corrected).
   TF1 fpl("fpl", "pol0", 16.1, 44.0);
   Csum->Fit(&fpl, "QR0");
   const double Cplat = fpl.GetParameter(0);
   printf("\n[Cjpx] C_sum plateau fit [16.1,44) = %.4f +- %.4f (chi2/ndf %.2f)\n", Cplat,
          fpl.GetParError(0), fpl.GetNDF() > 0 ? fpl.GetChisquare() / fpl.GetNDF() : 0.0);

   // R plateau ruler: smooth fit of the hw/sw ratio (hw_ratio.C output).
   double Rplat = 1.0;
   {
      TFile fr("R_JP2.root", "READ");
      auto *hR = fr.IsZombie() ? nullptr : (TH1D *)fr.Get("hR");
      if (hR) {
         TF1 fR("fR", "pol0", 16.1, 44.0);
         hR->Fit(&fR, "QR0");
         Rplat = fR.GetParameter(0);
         printf("[Cjpx] R(JP2) plateau fit [16.1,44) = %.4f +- %.4f\n", Rplat, fR.GetParError(0));
      } else {
         printf("[Cjpx] R_JP2.root not found here — run hw_ratio.C first (Rplat=1 for now)\n");
      }
   }

   printf("\nFrozen table for combined/promotion.h::JpxTrigEff  (C = That_sum x R_plateau):\n");
   printf("%-14s %8s %8s %8s %10s %10s\n", "bin", "eps_D", "eps_E", "C_sum", "That_sum", "C_JPX");
   for (int i = 1; i <= nb; ++i) {
      const double lo = QB[i - 1];
      const double that = (lo >= 16.1 - 1e-6) ? 1.0 : Csum->GetBinContent(i) / Cplat;
      const double cjpx = (lo >= 16.1 - 1e-6) ? Rplat : that;
      printf("[%5.1f,%5.1f) %8.4f %8.4f %8.4f %10.4f %10.4f\n", QB[i - 1], QB[i],
             epsD->GetBinContent(i), epsE->GetBinContent(i), Csum->GetBinContent(i), that, cjpx);
   }

   // ---- figure ----------------------------------------------------------------
   TCanvas c("c", "", 700, 800);
   c.Divide(1, 2, 0, 0);
   c.cd(1);
   gPad->SetGrid();
   epsD->SetTitle(";;#varepsilon = #Sigma_{cat} s_{cat} P(cat gate)");
   epsD->GetYaxis()->SetRangeUser(0, 1.09);
   epsD->SetMarkerStyle(20);
   epsD->SetMarkerColor(kBlue + 1);
   epsD->SetLineColor(kBlue + 1);
   epsD->Draw("E1");
   epsE->SetMarkerStyle(24);
   epsE->SetMarkerColor(kRed + 1);
   epsE->SetLineColor(kRed + 1);
   epsE->Draw("E1 SAME");
   TLegend lg(0.55, 0.15, 0.88, 0.35);
   lg.SetBorderSize(0);
   lg.AddEntry(epsD, "DATA (fired_JP0 base)", "lp");
   lg.AddEntry(epsE, "EMBEDDING (weighted)", "lp");
   lg.Draw();
   c.cd(2);
   gPad->SetGrid();
   Csum->SetTitle(";jet p_{T} [GeV];C_{sum} = #varepsilon_{data}/#varepsilon_{emb}");
   Csum->GetYaxis()->SetRangeUser(0.5, 1.14);
   Csum->SetMarkerStyle(21);
   Csum->Draw("E1");
   TLine one(QB.front(), 1, QB.back(), 1);
   one.SetLineStyle(2);
   one.Draw();
   TLine pl(QB.front(), Cplat, QB.back(), Cplat);
   pl.SetLineColor(kGreen + 2);
   pl.Draw();
   c.SaveAs("measure_Cjpx.pdf");
   TFile fout("measure_Cjpx.root", "RECREATE");
   epsD->Write();
   epsE->Write();
   Csum->Write();
   fout.Close();
}
