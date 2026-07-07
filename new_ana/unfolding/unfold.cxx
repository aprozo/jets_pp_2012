// Dmitry-style response build (per trigger):
//   * matched MC<->reco pairs go into RooUnfoldResponse::Fill
//   * unmatched truth ("missed" jets) go into Miss   — folds matching efficiency into M
//   * unmatched reco  ("fake"  jets) go into Fake    — folds fake rate into M
//
// The response alone is sufficient to unfold; no separate trigEff / fakeRate /
// matchEff multiplication downstream. RooUnfoldBayes with ~2 iterations
// approximates Dmitry's direct-inversion solution while staying stable.
//
// Reads the single Stage-1 embedding merge <datapath>/merged_matching_R<R>.root
// and selects the trigger via reco_trigger_match_<T>. Writes
// <workdir>/response_<T>_R<R>.root (read by cross_section.cpp) plus a train/test
// closure plot under pdf/.

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TSystem.h"

#include <fstream>
#include <map>
#include <sstream>

#include "../config.h"
#include "../soft_reweight.h"
#include "../vertex_reweight.h"

using namespace CrossSectionConfig;
AnalysisConfig cfg;

static TH1D *divideByBinWidth(TH1D *h)
{
   for (int i = 1; i <= h->GetNbinsX(); ++i)
      if (h->GetBinWidth(i) != 0) h->SetBinContent(i, h->GetBinContent(i) / h->GetBinWidth(i));
   return h;
}

static void plotIterations(TCanvas *can, TString outPdf, RooUnfoldResponse *response, TH1D *hTruth, TH1D *hMeasured)
{
   const std::vector<int> iters = {1, 2, 3, 4, 5};
   can->Divide(2, 1);
   TLegend *leg = new TLegend(0.30, 0.62, 0.44, 0.88);
   leg->AddEntry(hTruth, "Truth", "l");
   leg->AddEntry(hMeasured, "Measured", "l");
   can->cd(1);
   gPad->SetLogy();

   TH1D *hMeasN = divideByBinWidth((TH1D *)hMeasured->Clone("hMeasN"));
   TH1D *hTruthN = divideByBinWidth((TH1D *)hTruth->Clone("hTruthN"));
   hTruthN->SetMarkerStyle(20);
   hTruthN->SetLineColor(kViolet);
   hMeasN->SetMarkerStyle(21);
   hMeasN->SetLineColor(kTeal - 1);
   hMeasN->Draw("hist");
   hTruthN->Draw("hist same");

   TLine *line = new TLine();
   line->SetLineStyle(2);
   line->SetLineColor(kGray + 2);

   for (size_t k = 0; k < iters.size(); ++k) {
      can->cd(1);
      RooUnfoldBayes u(response, hMeasured, iters[k]);
      TH1D *hU = (TH1D *)u.Hreco(RooUnfold::kCovariance)->Clone(Form("hU%zu", k));
      hU->SetLineColor(2000 + (int)k);
      hU->SetMarkerColor(2000 + (int)k);
      hU->SetMarkerStyle(20);
      leg->AddEntry(hU, Form("Iter%d", iters[k]), "pel");
      TH1D *hUN = divideByBinWidth((TH1D *)hU->Clone(Form("hUN%zu", k)));
      hUN->Draw("PE same");

      can->cd(2);
      TH1D *r = (TH1D *)hU->Clone(Form("ratio%zu", k));
      r->SetLineColor(2000 + (int)k);
      r->SetMarkerColor(2000 + (int)k);
      r->SetMarkerStyle(20);
      r->Divide(hTruth);
      r->GetYaxis()->SetTitle("Unfolded/Truth");
      r->GetYaxis()->SetRangeUser(0.8, 1.2);
      r->Draw(k == 0 ? "PE" : "PE same");
      line->DrawLine(r->GetXaxis()->GetXmin(), 1.05, r->GetXaxis()->GetXmax(), 1.05);
      line->DrawLine(r->GetXaxis()->GetXmin(), 0.95, r->GetXaxis()->GetXmax(), 0.95);
   }
   can->cd();
   leg->Draw();
   can->SaveAs(outPdf);
}

static void single_unfold(const std::string &trigger, const std::string &jetR, const Systematic &syst)
{
   // IMT speeds the merged-tree read up ~3-4x. Bounded (config.h kImtThreads) —
   // unbounded IMT over gpfs on multi-GB files has crashed.
   ROOT::EnableImplicitMT(kImtThreads);

   // One Stage-1 embedding merge for every trigger; select via reco_trigger_match_<T>.
   const TString inputFile =
      TString(Form("%smerged_matching_R%s.root", cfg.datapath.c_str(), jetR.c_str()));
   if (gSystem->AccessPathName(inputFile)) {
      std::cout << "Skipping " << trigger << " R=" << jetR << ": " << inputFile << " missing" << std::endl;
      return;
   }

   // TRUTH-AXIS CONVENTION: the paper reports dsigma/d(pT - dpT_UE) at the
   // PARTICLE level; truth axis = mc_pt_corrected (UE-subtracted particle pT).
   const char *MCPT = "mc_pt_corrected";

   gStyle->SetOptStat(0);
   gStyle->SetPaintTextFormat(".2f");
   gSystem->mkdir("pdf", kTRUE);

   const TString outName = Form("%sresponse_%s_R%s%s.root", cfg.workdir.c_str(), trigger.c_str(),
                                jetR.c_str(), SystTag(syst).c_str());
   TFile *fout = new TFile(outName, "RECREATE");

   const std::vector<double> &recoBins = RecoBins();
   const std::vector<double> &mcBins = McBins();

   ROOT::RDataFrame rawJets(std::string("MatchedTree"), std::string(inputFile.Data()));

   // Row acceptance — DECOUPLED detector-eta gate (always on). The data
   // (cross_section.cpp::Raw) gates on |det_eta|<0.5. With the wide Stage-1
   // jet-find (|eta_phys|<1.0) an AND-coupled physics-eta filter drops the
   // ~+9% eta-edge reco jets the data counts (-> +8% pT-uniform overshoot);
   // swapping to AND-coupled det-eta drops migration-out misses (-> ~-10%).
   // The eta migration is nearly balanced, so keep the row if the RECO is in
   // the detector acceptance (|det_eta|<0.5, == the data) OR the TRUTH is in
   // the particle acceptance (|y|<0.5). The inner matched/miss/fake cuts below
   // classify on reco_det_eta, so edge fakes AND migration-out misses are both
   // carried.
   auto rawJetsCut = rawJets.Filter("((reco_pt > -500 && std::abs(reco_det_eta) < 0.5) || "
                                    "(mc_pt > -500 && std::abs(mc_eta) < 0.5))",
                                    "decoupled |det_eta|<0.5 (reco) OR |eta|<0.5 (truth)");

   // Vertex-z and soft-pT reweights ON.
   auto allJets = rawJetsCut
      .Define("vertex_weight", [](double vz)    { return VertexReweight::weight(vz); },  {"event_vz"})
      .Define("soft_weight",   [](double pthat) { return SoftReweight::weight(pthat); }, {"pthat_mid"})
      .Define("total_weight",  "mc_weight * vertex_weight * soft_weight")
      // Shape systematic: shift the RECONSTRUCTED energy scale used by the
      // response, per jet — flat (jesShift), EMC rt-weighted
      // tower(3.2%)/track(1.1%) scale (emcSign), track-efficiency equivalent
      // 1%*(1-rt) (trkSign), optional Gaussian smear (jerSmear). Identity for
      // the nominal. The DATA is not shifted — shifting the response reco
      // re-maps the unchanged data. Sentinel -999 (unmatched reco) is
      // preserved so the presence tests below still work.
      .DefineSlot("reco_pt_s",
                  [syst](unsigned int, double rpt, double rt) -> double {
                     if (rpt < -500.0) return rpt;
                     double f = 1.0 + syst.RecoShift(rt);
                     if (syst.jerSmear > 0.0) f += gRandom->Gaus(0.0, syst.jerSmear);
                     return rpt * f;
                  },
                  {"reco_pt", "reco_neutral_fraction"});

   // matching_mc_reco.cxx writes mc_pt = -999 / reco_pt = -999 for unmatched halves.
   const double ptMcCut = mcBins.front();

   // Per-trigger reco-pT validity floor — same TrigPtFloor() as
   // cross_section.cpp::Raw(), so data and response are consistent by
   // construction (never below the first reco bin edge here). Reco jets below
   // the floor are excluded from the measured side (matched + fake); a truth
   // jet whose only reco partner is below the floor (or not trigger-matched) is
   // a MISS, so the matching/trigger efficiency folded into the response stays
   // consistent with the data selection.
   const double recoFloor = std::max(TrigPtFloor(trigger), pt_reco_bins.front());

   // The data spectrum we unfold is filtered to jets matched to a fired trigger
   // patch (cross_section.cpp::Raw -> trigger_match_<T>). The response encodes
   // the same reco-side selection:
   //   * matched = truth+reco pair where reco jet is in a fired patch
   //   * missed  = truth jet with NO reco match, OR reco match not trigger-matched
   //   * fakes   = reco jet with no truth match, but trigger-matched
   // SYMMETRIC NEF<0.95 and |det_eta|<0.5 mirror the data selection.
   const TString cutMatched =
      Form("%s > %f && reco_pt_s > %f && reco_trigger_match_%s && reco_neutral_fraction < 0.95 && "
           "abs(reco_det_eta) < 0.5",
           MCPT, ptMcCut, recoFloor, trigger.c_str());
   const TString cutMissed =
      Form("%s > %f && (reco_pt_s < %f || !reco_trigger_match_%s || reco_neutral_fraction >= 0.95 || "
           "abs(reco_det_eta) >= 0.5)",
           MCPT, ptMcCut, recoFloor, trigger.c_str());
   const TString cutFake =
      Form("mc_pt < -500 && reco_pt_s > %f && reco_trigger_match_%s && reco_neutral_fraction < 0.95 && "
           "abs(reco_det_eta) < 0.5",
           recoFloor, trigger.c_str());

   auto matched = allJets.Filter(cutMatched.Data(), "matched");
   auto missed  = allJets.Filter(cutMissed.Data(),  "missed");
   auto fakes   = allJets.Filter(cutFake.Data(),    "fakes");

   // Train/test split on matched only — Miss/Fake fold into the response from
   // the training events so the closure test ratios are meaningful.
   const float testFraction = 0.2;
   auto matchedTrain = matched.DefineSlot("rng", [](unsigned int) -> double { return gRandom->Uniform(); }, {})
                          .Filter(Form("rng > %f", testFraction), "train (matched)");
   auto matchedTest  = matched.DefineSlot("rng", [](unsigned int) -> double { return gRandom->Uniform(); }, {})
                          .Filter(Form("rng <= %f", testFraction), "test (matched)");
   auto missedTrain  = missed.DefineSlot("rng", [](unsigned int) -> double { return gRandom->Uniform(); }, {})
                          .Filter(Form("rng > %f", testFraction), "train (miss)");
   auto missedTest   = missed.DefineSlot("rng", [](unsigned int) -> double { return gRandom->Uniform(); }, {})
                          .Filter(Form("rng <= %f", testFraction), "test (miss)");
   auto fakesTrain   = fakes.DefineSlot("rng", [](unsigned int) -> double { return gRandom->Uniform(); }, {})
                          .Filter(Form("rng > %f", testFraction), "train (fake)");
   auto fakesTest    = fakes.DefineSlot("rng", [](unsigned int) -> double { return gRandom->Uniform(); }, {})
                          .Filter(Form("rng <= %f", testFraction), "test (fake)");

   // Closure-check histograms. The Miss/Fake-aware response unfolds
   // (matched_reco + fake_reco) into (matched_truth + missed_truth), so the
   // reference TRUTH must include misses and the input MEASURED must include
   // fakes.
   auto hMeasuredMatched = matchedTest.Histo1D({"_measMatched", "", (int)recoBins.size() - 1, recoBins.data()},
                                               "reco_pt_s", "total_weight");
   auto hMeasuredFakes   = fakesTest.Histo1D({"_measFakes", "", (int)recoBins.size() - 1, recoBins.data()},
                                             "reco_pt_s", "total_weight");
   auto hTruthMatched    = matchedTest.Histo1D({"_truthMatched", "", (int)mcBins.size() - 1, mcBins.data()},
                                               MCPT, "total_weight");
   auto hTruthMissed     = missedTest.Histo1D({"_truthMissed", "", (int)mcBins.size() - 1, mcBins.data()},
                                              MCPT, "total_weight");

   TH1D measAxis("Measured", ";p_{T}^{reco};dN/dp_{T}", (int)recoBins.size() - 1, recoBins.data());
   TH1D truthAxis("Truth", ";p_{T}^{mc};dN/dp_{T}", (int)mcBins.size() - 1, mcBins.data());

   auto *response = new RooUnfoldResponse(&measAxis, &truthAxis);
   response->SetName("my_response");
   response->UseOverflow(false);

   matchedTrain.Foreach(
      [response](double rpt, double mpt, double w) { response->Fill(rpt, mpt, w); },
      {"reco_pt_s", MCPT, "total_weight"});
   fakesTrain.Foreach(
      [response](double rpt, double w) { response->Fake(rpt, w); },
      {"reco_pt_s", "total_weight"});
   missedTrain.Foreach(
      [response](double mpt, double w) { response->Miss(mpt, w); },
      {MCPT, "total_weight"});

   const long long nMatched = *matched.Count();
   const long long nMissed  = *missed.Count();
   const long long nFakes   = *fakes.Count();
   std::cout << "Response built (" << trigger << " R=" << jetR << "): "
             << "matched=" << nMatched
             << "  missed=" << nMissed
             << "  fakes=" << nFakes << std::endl;

   TH1D *hMeasuredTest =
      (TH1D *)hMeasuredMatched->Clone(Form("MeasuredTest_%s_R%s", trigger.c_str(), jetR.c_str()));
   hMeasuredTest->SetDirectory(nullptr);
   hMeasuredTest->Add(hMeasuredFakes.GetPtr());
   hMeasuredTest->SetName("MeasuredTest");

   TH1D *hTruthTest =
      (TH1D *)hTruthMatched->Clone(Form("TruthTest_%s_R%s", trigger.c_str(), jetR.c_str()));
   hTruthTest->SetDirectory(nullptr);
   hTruthTest->Add(hTruthMissed.GetPtr());
   hTruthTest->SetName("TruthTest");

   fout->cd();
   hMeasuredTest->Write();
   hTruthTest->Write();
   response->Write();
   StampProvenance(syst);

   TCanvas *can = new TCanvas("can", "Closure", 1400, 600);
   const TString outPdf =
      Form("pdf/closure_check_dmitry_%s_R%s.pdf", trigger.c_str(), jetR.c_str());
   can->SaveAs(outPdf + "[");
   plotIterations(can, outPdf, response, hTruthTest, hMeasuredTest);
   can->SaveAs(outPdf + "]");
   fout->Close();
}

// systName selects a preset from config.h::Systematics(); default "nominal" is
// identity (no reco shift), so the nominal call `unfold.cxx+` is unchanged. A
// shape systematic (jesShift/jerSmear) writes response_<T>_R<R>_<name>.root.
void unfold(const char *systName = "nominal")
{
   DefineCustomColors();
   const Systematic syst = FindSystematic(systName);
   std::cout << "[unfold] systematic = " << syst.name << " (jesShift " << syst.jesShift
             << ", jerSmear " << syst.jerSmear << ")" << std::endl;
   for (const auto &trigger : cfg.triggers)
      for (const auto &jetR : cfg.jetRs)
         single_unfold(trigger, jetR, syst);
}
