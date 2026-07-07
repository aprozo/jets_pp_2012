// Matrix-inversion cross section — Dmitry's unregularized inversion path, the
// second solver alongside the main Bayes pipeline. Same physical inputs: the
// etadet Miss/Fake response (square, built by unfold.cxx here), the data divided
// by the measured hybrid trigger correction C(pt) (config.h TrigEffMeas), and
// each trigger restricted to its quote window (config.h QuoteLo). No environment
// flags; no non-physical correction.
//
// The raw inverse oscillates in the steep turn-on below ~13 GeV (near-singular
// migration there), but those bins are below every quote window and dropped, so
// the quoted spectrum is the stable plateau region. Dmitry tames the turn-on
// with floor-restricted cell-filtering, not reproduced here.
//
// Outputs (in config.h kWorkDir):
//   xsec_<T>_R<R>_invert.root                 hist "canonical" + "reference"
//   comparison_with_dmitriy_R<R>_invert.pdf

#include "RooUnfoldInvert.h"
#include "RooUnfoldResponse.h"
#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "../config.h"
using namespace CrossSectionConfig;

AnalysisConfig cfg;

static void AppendBadRunsFromFile(const std::string &path)
{
   std::ifstream fp(path);
   if (!fp) {
      std::cerr << "[badruns] file missing: " << path << " — skipped\n";
      return;
   }
   int run, n = 0;
   while (fp >> run) {
      cfg.badRuns.push_back(run);
      ++n;
   }
   std::cout << "[badruns] appended " << n << " runs from " << path << std::endl;
}

class Spectrum {
public:
   Spectrum(std::string trig, std::string R) : trigger(std::move(trig)), jetR(std::move(R))
   {
      dataFileName = std::string(Form("%smerged_data_R%s.root", cfg.datapath.c_str(), jetR.c_str()));
      responseFile.reset(TFile::Open(
         Form((cfg.workdir + "response_%s_R%s_square.root").c_str(), trigger.c_str(), jetR.c_str()), "READ"));
      response = responseFile ? (RooUnfoldResponse *)responseFile->Get("my_response") : nullptr;
      if (!response) {
         std::cerr << "[FATAL] response_" << trigger << "_R" << jetR
                   << "_square.root missing 'my_response' — run unfold.cxx first" << std::endl;
         gSystem->Exit(2);
      }
      std::cout << "[unfold] response = response_" << trigger << "_R" << jetR << "_square.root" << std::endl;
   }

   void Build()
   {
      h = Raw(dataFileName);

      std::cout << "[unfold] solver = RooUnfoldInvert (unregularized matrix inversion, Dmitry-style)"
                << std::endl;
      RooUnfoldInvert u(response, h.get());
      TH1D *hu = (TH1D *)u.Hreco(RooUnfold::kCovariance);
      TH1D *cu = (TH1D *)hu->Clone(Form("%s_unfolded_R%s", trigger.c_str(), jetR.c_str()));
      cu->SetDirectory(0);
      h.reset(cu);

      NormBinWidthAndAcc();

      double Loverride = -1.0;
      TFile lf((cfg.workdir + "lumi_zilong_full.root").c_str(), "READ");
      auto *lumi = lf.IsZombie() ? nullptr : (TH1D *)lf.Get(Form("luminosity_%s", trigger.c_str()));
      if (lumi) {
         double sum = 0.0;
         int nKept = 0, nDropped = 0;
         for (int i = 1; i <= lumi->GetNbinsX(); ++i) {
            const char *label = lumi->GetXaxis()->GetBinLabel(i);
            if (!label || !*label) continue;
            if (std::find(cfg.badRuns.begin(), cfg.badRuns.end(), std::stoi(label)) != cfg.badRuns.end()) {
               ++nDropped;
               continue;
            }
            sum += lumi->GetBinContent(i);
            ++nKept;
         }
         Loverride = sum;
         std::cout << "[lumi] Leff(" << trigger << "): " << nKept << " kept / " << nDropped
                   << " bad runs" << std::endl;
      }
      if (Loverride <= 0) {
         std::cerr << "[FATAL] no Leff for trigger " << trigger << std::endl;
         gSystem->Exit(2);
      }
      std::cout << "[lumi] Leff(" << trigger << ") = " << Loverride << " pb^-1" << std::endl;
      h->Scale(1.0 / Loverride);

      // Restrict to the trigger's efficient quote window.
      const double qlo = QuoteLo(trigger);
      int nDrop = 0;
      for (int i = 1; i <= h->GetNbinsX(); ++i)
         if (h->GetXaxis()->GetBinCenter(i) < qlo) {
            h->SetBinContent(i, 0.0);
            h->SetBinError(i, 0.0);
            ++nDrop;
         }
      std::cout << "[quote] " << trigger << ": restricted to pT >= " << qlo << " (dropped " << nDrop
                << " low-pT bins)" << std::endl;
   }

   TH1D *Hist() const { return h.get(); }
   const std::string trigger, jetR;

private:
   std::unique_ptr<TH1D> h;
   std::string dataFileName;
   std::unique_ptr<TFile> responseFile;
   RooUnfoldResponse *response;

   std::unique_ptr<TH1D> Raw(const std::string &dataFileName) const
   {
      ROOT::RDataFrame dfRaw("ResultTree", dataFileName.c_str());
      ROOT::RDF::RNode df = dfRaw;
      const std::string ptCol = "pt_corrected";
      const std::vector<double> &recoBins = McBins(); // square: data reco axis == truth axis
      const std::string br = "trigger_match_" + trigger;

      ROOT::RDF::RNode dfn(df);
      if (trigger == "JP0" || trigger == "JP1" || trigger == "JP2") {
         std::cout << "[data] event didFire filter: fired_" << trigger << std::endl;
         dfn = dfn.Filter("fired_" + trigger);
      }

      const std::string nefCut = "neutral_fraction <= 0.95";
      const std::string floorCut = " && " + ptCol + " > " + std::to_string(TrigPtFloor(trigger));

      auto htmp = dfn.Define("m", br)
                     .Filter(
                        [](int runIndex) {
                           return runIndex >= 0 &&
                                  std::find(cfg.badRuns.begin(), cfg.badRuns.end(), runIndex) == cfg.badRuns.end();
                        },
                        {"runid1"})
                     .Define("eta_sel", ("m && abs(det_eta) < 0.5 && " + nefCut + floorCut).c_str())
                     .Define("pt_trig", (ptCol + "[eta_sel]").c_str())
                     .Filter("pt_trig.size() > 0")
                     .Histo1D({Form("%s_raw_R%s", trigger.c_str(), jetR.c_str()),
                               Form("%s raw R=%s;p_{T};counts", trigger.c_str(), jetR.c_str()),
                               (int)recoBins.size() - 1, recoBins.data()},
                              "pt_trig");

      TH1D *c = (TH1D *)htmp->Clone(Form("%s_raw_clone_R%s", trigger.c_str(), jetR.c_str()));
      c->SetDirectory(0);

      // Measured hybrid trigger correction (JP0/JP1/JP2): divide by C(pt).
      if (trigger == "JP0" || trigger == "JP1" || trigger == "JP2") {
         std::cout << "[data] dividing fired_" << trigger << " data by C(pt) (config.h TrigEffMeas)"
                   << std::endl;
         for (int i = 1; i <= c->GetNbinsX(); ++i) {
            const double p = TrigEffMeas(trigger, c->GetXaxis()->GetBinCenter(i));
            if (p > 0) {
               c->SetBinContent(i, c->GetBinContent(i) / p);
               c->SetBinError(i, c->GetBinError(i) / p);
            }
         }
      }
      return std::unique_ptr<TH1D>(c);
   }

   void NormBinWidthAndAcc()
   {
      // eta acceptance |eta_det| < kDetEtaMax (radius-independent)
      h->Scale(1.0 / EtaAcceptance());
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
         const double w = h->GetBinWidth(i);
         h->SetBinContent(i, h->GetBinContent(i) / w);
         h->SetBinError(i, h->GetBinError(i) / w);
      }
   }
};

static void CompareWithDmitriy(const std::string &jetR, const std::vector<Spectrum> &spectrums)
{
   const std::string refPath =
      std::string(Form((cfg.workdir + "jet_cross_section_dmitriyR%s.root").c_str(), jetR.c_str()));
   if (gSystem->AccessPathName(refPath.c_str())) {
      std::cerr << "[warn] " << refPath << " missing; skipping Dmitriy comparison" << std::endl;
      return;
   }
   TFile rf(refPath.c_str(), "READ");
   TH1D *ref = (TH1D *)rf.Get("crossSection_systematic")->Clone(Form("ref_R%s", jetR.c_str()));
   ref->SetDirectory(0);

   TCanvas *c = new TCanvas(Form("c_R%s", jetR.c_str()), "", 900, 900);
   c->Divide(1, 2);

   auto color_for = [](const std::string &t) -> int {
      if (t == "JP2") return kAzure - 1;
      if (t == "HT2") return kRed + 1;
      if (t == "JP1") return kOrange - 1;
      if (t == "JP0") return kGreen + 2;
      return kMagenta + 1;
   };
   auto marker_for = [](const std::string &t) -> int {
      if (t == "JP2") return 20;
      if (t == "HT2") return 21;
      if (t == "JP1") return 22;
      if (t == "JP0") return 33;
      return 34;
   };

   auto *p1 = (TPad *)c->cd(1);
   p1->SetPad(0, 0.5, 1, 1);
   p1->SetBottomMargin(0);
   p1->SetLogy();

   ref->SetMarkerStyle(29);
   ref->SetMarkerSize(2);
   ref->SetMarkerColor(kViolet);
   ref->SetLineColor(kViolet);
   ref->GetYaxis()->SetRangeUser(1.2, 1e7);
   ref->Draw("E1");

   TLegend *leg = new TLegend(0.60, 0.50, 0.88, 0.88);
   leg->AddEntry(ref, "Dmitriy Table III", "lep");

   for (const auto &s : spectrums) {
      TH1D *h = s.Hist();
      for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
         const double x = h->GetXaxis()->GetBinCenter(bin);
         if (x < ref->GetXaxis()->GetXmin() || x > ref->GetXaxis()->GetXmax()) {
            h->SetBinContent(bin, 0.0);
            h->SetBinError(bin, 0.0);
         }
      }
      h->SetMarkerStyle(marker_for(s.trigger));
      h->SetMarkerColor(color_for(s.trigger));
      h->SetLineColor(color_for(s.trigger));
      h->Draw("E1 SAME");
      leg->AddEntry(h, (s.trigger + " (invert)").c_str(), "lep");
   }
   leg->Draw();

   auto *p2 = (TPad *)c->cd(2);
   p2->SetPad(0, 0, 1, 0.5);
   p2->SetTopMargin(0);
   p2->SetBottomMargin(0.30);

   TH1D *frame = (TH1D *)ref->Clone(Form("ratio_frame_R%s", jetR.c_str()));
   frame->Reset();
   frame->SetDirectory(0);
   frame->SetTitle(Form("Ratio to Dmitriy R=%s; p_{T} [GeV/c]; Ratio", jetR.c_str()));
   frame->GetYaxis()->SetRangeUser(0.5, 1.5);
   frame->Draw();

   TH1D *refBand = (TH1D *)ref->Clone(Form("ref_unity_band_R%s", jetR.c_str()));
   refBand->Reset();
   refBand->SetDirectory(0);
   for (int b = 1; b <= refBand->GetNbinsX(); ++b) {
      const double cRef = ref->GetBinContent(b);
      const double eRef = ref->GetBinError(b);
      refBand->SetBinContent(b, 1.0);
      refBand->SetBinError(b, (cRef > 0.0) ? (eRef / cRef) : 0.0);
   }
   refBand->SetFillColorAlpha(kViolet, 0.20);
   refBand->SetFillStyle(1001);
   refBand->SetLineColor(kViolet);
   refBand->SetMarkerSize(0);
   refBand->Draw("E2 SAME");

   TLine *line = new TLine(ref->GetXaxis()->GetXmin(), 1.0, ref->GetXaxis()->GetXmax(), 1.0);
   line->SetLineColor(kViolet);
   line->SetLineWidth(2);
   line->Draw("SAME");

   for (const auto &s : spectrums) {
      TH1D *h = s.Hist();
      TH1D *r = new TH1D(Form("ratio_%s_R%s", s.trigger.c_str(), jetR.c_str()), "", ref->GetNbinsX(),
                         ref->GetXaxis()->GetXbins()->GetArray());
      r->SetDirectory(0);
      std::cout << "[ratio invert/Dmitry] " << s.trigger << " R=" << jetR << " (quote window)\n";
      for (int bin = 1; bin <= r->GetNbinsX(); ++bin) {
         const double x = r->GetXaxis()->GetBinCenter(bin);
         const double rv = ref->GetBinContent(bin);
         const double mine = h->GetBinContent(h->FindBin(x));
         r->SetBinContent(bin, rv > 0 ? mine / rv : 0.0);
         r->SetBinError(bin, rv > 0 ? h->GetBinError(h->FindBin(x)) / rv : 0.0);
         if (rv > 0 && x >= QuoteLo(s.trigger))
            printf("    pT=%5.2f  invert=%9.3g  Dmitry=%9.3g  ratio=%5.3f\n", x, mine, rv, mine / rv);
      }
      r->SetMarkerStyle(marker_for(s.trigger));
      r->SetMarkerColor(color_for(s.trigger));
      r->SetLineColor(color_for(s.trigger));
      r->Draw("E1 SAME");
   }

   c->SaveAs(Form("%scomparison_with_dmitriy_R%s_invert.pdf", cfg.workdir.c_str(), jetR.c_str()));

   for (const auto &s : spectrums) {
      TFile fout(Form("%sxsec_%s_R%s_invert.root", cfg.workdir.c_str(), s.trigger.c_str(), jetR.c_str()),
                 "RECREATE");
      s.Hist()->Write("canonical");
      ref->Write("reference");
      fout.Close();
   }
}

static void QuietRootErrors(Int_t level, Bool_t abort, const char *location, const char *msg)
{
   if (msg && std::strstr(msg, "already deleted") != nullptr) return;
   DefaultErrorHandler(level, abort, location, msg);
}

void cross_section()
{
   SetErrorHandler(QuietRootErrors);
   ROOT::EnableImplicitMT(kImtThreads);
   DefineCustomColors();
   gStyle->SetOptStat(0);
   TH1::SetDefaultSumw2();

   AppendBadRunsFromFile(cfg.workdir + "../lists/dmitry_extras.list");

   for (const auto &R : cfg.jetRs) {
      std::vector<Spectrum> spectrums;
      spectrums.reserve(cfg.triggers.size());
      for (const auto &trig : cfg.triggers) {
         const std::string dataPath = Form("%smerged_data_R%s.root", cfg.datapath.c_str(), R.c_str());
         if (gSystem->AccessPathName(dataPath.c_str())) {
            std::cerr << "[warn] " << dataPath << " missing; skipping " << trig << std::endl;
            continue;
         }
         spectrums.emplace_back(trig, R);
         spectrums.back().Build();
      }
      CompareWithDmitriy(R, spectrums);
   }
}
