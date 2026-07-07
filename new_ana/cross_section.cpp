// FINAL unifier of the Stage-2 pipeline. For each trigger it builds the
// trigger-gated raw reco-pT spectrum from the data, divides by the compiled-in
// measured trigger correction C(pt) (config.h::TrigEffMeas), unfolds with
// RooUnfoldBayes(nIter=2) against the Miss/Fake response built by unfold.cxx,
// normalizes (eta acceptance, bin width, per-run luminosity), compares to
// Dmitry's published Table III, and WRITES the result histograms.
//
// Reads the single Stage-1 data merge <datapath>/merged_data_R<R>.root and
// selects the trigger via trigger_match_<T> / fired_<T>.
//
// Outputs (in config.h kWorkDir):
//   xsec_<T>_R<R>.root                hist "canonical" (my spectrum) + "reference"
//   comparison_with_dmitriy_R<R>.pdf
//   qa.pdf                            build trace

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
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

#include "config.h"
using namespace CrossSectionConfig;

AnalysisConfig cfg;

// Append run numbers from a text file (one int per line) into cfg.badRuns, so
// spectrum and runLumi reject the same set.
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
   Spectrum(std::string trig, std::string R, const Systematic &s)
      : trigger(std::move(trig)), jetR(std::move(R)), syst(s)
   {
      // One Stage-1 data merge for every trigger; select via fired_<T> / trigger_match_<T>.
      dataFileName = std::string(Form("%smerged_data_R%s.root", cfg.datapath.c_str(), jetR.c_str()));

      // The TFile must outlive `response`; keep it as a member to anchor its
      // lifetime. A SHAPE systematic (jesShift/jerSmear) has its own rebuilt
      // response (response_<T>_R<R>_<name>.root); everything else reuses the
      // nominal one. Efficiency/fake/trigger corrections are absorbed into the
      // matrix and unfolded in one RooUnfoldBayes step.
      const std::string respTag = syst.needsResponse() ? SystTag(syst) : std::string("");
      responseFile.reset(TFile::Open(
         Form((cfg.workdir + "response_%s_R%s%s.root").c_str(), trigger.c_str(), jetR.c_str(), respTag.c_str()),
         "READ"));
      response = responseFile ? (RooUnfoldResponse *)responseFile->Get("my_response") : nullptr;
      std::cout << "[unfold] using response response_" << trigger << "_R" << jetR << respTag << ".root"
                << std::endl;
   }

   void Dump(TH1 *h, TCanvas *c, const std::string &msg = "")
   {
      c->cd();
      c->SetLogy();

      TH1D *hClone = (TH1D *)h->Clone(Form("%s_clone", h->GetName()));
      hClone->SetDirectory(nullptr);
      hClone->SetMinimum(1);
      hClone->SetMaximum(0.5e8);

      const bool hasFrame = (c->GetListOfPrimitives()->GetSize() > 0);
      hClone->DrawClone(hasFrame ? "same" : "");

      TLegend *leg = (TLegend *)c->GetListOfPrimitives()->FindObject("legend");
      if (!leg) {
         leg = new TLegend(0.60, 0.70, 0.88, 0.88);
         leg->SetName("legend");
         c->GetListOfPrimitives()->Add(leg);
      }
      leg->AddEntry(hClone, msg.c_str(), "l");
      leg->Draw();
   }

   void Build()
   {
      std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>();
      //=================== raw ===============
      h = Raw(dataFileName);
      h->SetLineColor(2001);
      Dump(h.get(), c.get(), "raw");

      //=================== unfold ===============
      Unfold();
      h->SetLineColor(2004);
      Dump(h.get(), c.get(), "after unfolding");

      //=================== normalize ===============
      NormBinWidthAndAcc();
      // Per-trigger integrated effective luminosity, computed AT RUNTIME from
      // lumi_zilong_full.root::luminosity_<T> (one labeled bin per run;
      // prescale + livetime already folded in), summed over runs NOT in
      // cfg.badRuns. badRuns edits propagate automatically.
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
         std::cout << "[lumi] Leff(" << trigger << ") from lumi_zilong_full.root: " << nKept
                   << " kept / " << nDropped << " bad runs" << std::endl;
      }
      if (Loverride <= 0) {
         std::cerr << "[FATAL] no Leff for trigger " << trigger
                   << " — luminosity_" << trigger << " missing in lumi_zilong_full.root" << std::endl;
         gSystem->Exit(2);
      }
      const double L = Loverride * syst.lumiScale;
      std::cout << "[lumi] Leff(" << trigger << ") = " << L << " pb^-1 (lumiScale " << syst.lumiScale << ")"
                << std::endl;
      h->Scale(1.0 / L);

      // Restrict to the trigger's efficient quote window: drop bins below
      // QuoteLo, where the standalone spectrum is unfolding extrapolation below
      // the turn-on rather than a measurement.
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

      c->SaveAs((cfg.workdir + "qa.pdf").c_str());
   }

   void Unfold()
   {
      std::cout << "[unfold] solver = RooUnfoldBayes (nIter=" << syst.nIter << ")" << std::endl;
      RooUnfoldBayes u(response, h.get(), syst.nIter);
      TH1D *hu = (TH1D *)u.Hreco(RooUnfold::kCovariance);

      TH1D *c = (TH1D *)hu->Clone(Form("%s_unfolded_R%s", trigger.c_str(), jetR.c_str()));
      c->SetDirectory(0);
      h.reset(c);
   }

   TH1D *Hist() const { return h.get(); }
   const Systematic &Syst() const { return syst; }
   const std::string trigger, jetR;

private:
   std::unique_ptr<TH1D> h;
   Systematic syst;
   std::string dataFileName;
   std::unique_ptr<TFile> responseFile; // anchors response lifetime
   RooUnfoldResponse *response;

   std::unique_ptr<TH1D> Raw(const std::string &dataFileName) const
   {
      ROOT::RDataFrame dfRaw("ResultTree", dataFileName.c_str());
      ROOT::RDF::RNode df = dfRaw;
      const std::string ptCol = "pt_corrected"; // UE-subtracted reco pT
      const std::vector<double> &recoBins = RecoBins();
      const std::string br = "trigger_match_" + trigger;

      // Event-level didFire requirement (hardware fired_<T>) for the JP triggers.
      // For JP1/JP0 it is the uniform-lumi prescale sample matching the
      // sampled-lumi normalization in Build(); for JP2 it defines the clean JP2
      // luminosity sample. HT2 uses only the per-jet tower match.
      const bool eventFired = (trigger == "JP0" || trigger == "JP1" || trigger == "JP2");
      ROOT::RDF::RNode dfn(df);
      if (eventFired) {
         std::cout << "[data] event didFire filter: fired_" << trigger << std::endl;
         dfn = dfn.Filter("fired_" + trigger);
      }

      // NEF<0.95 suppresses hot-tower fake jets that inflate the high-pT tail.
      const std::string nefCut = "neutral_fraction <= 0.95";

      // Per-trigger reco-pT validity floor on the UE-subtracted reco pT.
      const double trigPtFloor = TrigPtFloor(trigger);
      const std::string floorCut = " && " + ptCol + " > " + std::to_string(trigPtFloor);

      auto htmp = dfn.Define("m", br)
                     .Filter(
                        [](int runIndex) {
                           // Skip runs that are unknown or marked as bad.
                           return runIndex >= 0 &&
                                  std::find(cfg.badRuns.begin(), cfg.badRuns.end(), runIndex) == cfg.badRuns.end();
                        },
                        {"runid1"})
                     // |det_eta|<0.5 on top of the trigger match (Dmitry
                     // 00_05_00_05 = |jet.eta|<0.5 && |jet.detEta|<0.5).
                     .Define("eta_sel", ("m && abs(det_eta) < 0.5 && " + nefCut + floorCut).c_str())
                     .Define("pt_trig", (ptCol + "[eta_sel]").c_str())
                     .Filter("pt_trig.size() > 0")
                     .Histo1D({Form("%s_raw_R%s", trigger.c_str(), jetR.c_str()),
                               Form("%s raw R=%s;p_{T};counts", trigger.c_str(), jetR.c_str()),
                               (int)recoBins.size() - 1, recoBins.data()},
                              "pt_trig");

      TH1D *c = (TH1D *)htmp->Clone(Form("%s_raw_clone_R%s", trigger.c_str(), jetR.c_str()));
      c->SetDirectory(0);

      // MEASURED trigger correction (JP0/JP1/JP2): divide the recorded fired_<T>
      // data by the in-situ measured C(pt) = hybrid T-hat/R (config.h::TrigEffMeas).
      // Applied DATA-side (not folded into the response). HT2 -> C=1.
      if (trigger == "JP0" || trigger == "JP1" || trigger == "JP2") {
         std::cout << "[data] dividing fired_" << trigger
                   << " data by C(pt) (config.h TrigEffMeas), trigEffScale " << syst.trigEffScale << std::endl;
         for (int i = 1; i <= c->GetNbinsX(); ++i) {
            const double p = TrigEffMeas(trigger, c->GetXaxis()->GetBinCenter(i)) * syst.trigEffScale;
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
      // pt bin width
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
         const double w = h->GetBinWidth(i);
         h->SetBinContent(i, h->GetBinContent(i) / w);
         h->SetBinError(i, h->GetBinError(i) / w);
      }
   }
};

static void CompareWithDmitriy(const std::string &jetR, const std::vector<Spectrum> &spectrums)
{
   // Reference = PAPER Table III values (stat (+) syst as the violet band).
   const std::string refPath = std::string(
      Form((cfg.workdir + "jet_cross_section_dmitriyR%s.root").c_str(), jetR.c_str()));
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
   leg->AddEntry(ref, "Dmitriy", "lep");

   for (size_t i = 0; i < spectrums.size(); ++i) {
      TH1D *h = spectrums[i].Hist();
      for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
         const double x = h->GetXaxis()->GetBinCenter(bin);
         if (x < ref->GetXaxis()->GetXmin() || x > ref->GetXaxis()->GetXmax()) {
            h->SetBinContent(bin, 0.0);
            h->SetBinError(bin, 0.0);
         }
      }
      h->SetMarkerStyle(marker_for(spectrums[i].trigger));
      h->SetMarkerColor(color_for(spectrums[i].trigger));
      h->SetLineColor(color_for(spectrums[i].trigger));
      h->Draw("E1 SAME");
      leg->AddEntry(h, spectrums[i].trigger.c_str(), "lep");
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
      const double rel = (cRef > 0.0) ? (eRef / cRef) : 0.0;
      refBand->SetBinContent(b, 1.0);
      refBand->SetBinError(b, rel);
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

   for (size_t i = 0; i < spectrums.size(); ++i) {
      TH1D *hist = (TH1D *)spectrums[i].Hist()->Clone(Form("ratio_%s_R%s", spectrums[i].trigger.c_str(), jetR.c_str()));
      TH1D *r = new TH1D(Form("ratio_binned_%s_R%s", spectrums[i].trigger.c_str(), jetR.c_str()),
                         Form("Ratio to Dmitriy R=%s; p_{T} [GeV/c]; Ratio", jetR.c_str()), ref->GetNbinsX(),
                         ref->GetXaxis()->GetXbins()->GetArray());

      std::cout << "[ratio canonical/Dmitry] " << spectrums[i].trigger << " R=" << jetR << "\n";
      for (int bin = 1; bin <= r->GetNbinsX(); ++bin) {
         const double x = r->GetXaxis()->GetBinCenter(bin);
         const int binHist = hist->FindBin(x);
         r->SetBinContent(bin, hist->GetBinContent(binHist));
         r->SetBinError(bin, hist->GetBinError(binHist));
         const double rv = ref->GetBinContent(bin);
         if (rv > 0)
            printf("    pT=%5.2f  ours=%9.3g  ref=%9.3g  ratio=%5.3f\n", x, hist->GetBinContent(binHist), rv,
                   hist->GetBinContent(binHist) / rv);
      }
      r->SetDirectory(0);

      for (int b = 1; b <= r->GetNbinsX(); ++b) {
         const double rv = ref->GetBinContent(b);
         const double cur = r->GetBinContent(b);
         const double err = r->GetBinError(b);
         r->SetBinContent(b, rv > 0 ? cur / rv : 0.0);
         r->SetBinError(b,   rv > 0 ? err / rv : 0.0);
      }

      r->SetMarkerStyle(marker_for(spectrums[i].trigger));
      r->SetMarkerColor(color_for(spectrums[i].trigger));
      r->SetLineColor(color_for(spectrums[i].trigger));
      r->GetYaxis()->SetRangeUser(0.5, 1.5);
      r->GetYaxis()->SetTitle("Ratio to Dmitriy");
      r->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      r->Draw("E1 SAME");
   }

   c->SaveAs(Form("%scomparison_with_dmitriy_R%s.pdf", cfg.workdir.c_str(), jetR.c_str()));
}

// Persist the final canonical spectrum per trigger (+ Dmitry reference + full
// variation provenance) so plotting/systematics macros can consume it without
// rerunning. Written for every variation as xsec_<T>_R<R><SystTag>.root.
static void WriteXsec(const std::string &jetR, const std::vector<Spectrum> &spectrums)
{
   const std::string refPath =
      Form((cfg.workdir + "jet_cross_section_dmitriyR%s.root").c_str(), jetR.c_str());
   TH1D *ref = nullptr;
   if (!gSystem->AccessPathName(refPath.c_str())) {
      TFile rf(refPath.c_str(), "READ");
      TH1D *r = (TH1D *)rf.Get("crossSection_systematic");
      if (r) { ref = (TH1D *)r->Clone(Form("ref_R%s", jetR.c_str())); ref->SetDirectory(0); }
   }
   for (const auto &s : spectrums) {
      TFile fout(Form("%sxsec_%s_R%s%s.root", cfg.workdir.c_str(), s.trigger.c_str(), jetR.c_str(),
                      SystTag(s.Syst()).c_str()),
                 "RECREATE");
      s.Hist()->Write("canonical");
      if (ref) ref->Write("reference");
      StampProvenance(s.Syst());
      fout.Close();
   }
}

// ROOT 6.22 prints a flurry of "TList::Clear: ... already deleted" warnings at
// macro shutdown; filter them out via a custom error handler.
static void QuietRootErrors(Int_t level, Bool_t abort, const char *location, const char *msg)
{
   if (msg && std::strstr(msg, "already deleted") != nullptr)
      return;
   DefaultErrorHandler(level, abort, location, msg);
}

// systName selects a preset from config.h::Systematics(); default "nominal" is
// the physics result. Writes xsec_<T>_R<R><tag>.root for the variation; the
// Dmitry comparison PDF is drawn only for the nominal (a shifted spectrum vs
// Dmitry is not meaningful).
void cross_section(const char *systName = "nominal")
{
   SetErrorHandler(QuietRootErrors);
   ROOT::EnableImplicitMT(kImtThreads); // bounded IMT for the big data reads (~3-4x)
   DefineCustomColors();
   gStyle->SetOptStat(0);
   TH1::SetDefaultSumw2();

   const Systematic syst = FindSystematic(systName);
   std::cout << "[cross_section] systematic = " << syst.name << std::endl;

   // RunppAna bypasses both the bad-run mask and the good-run whitelist, so
   // ResultTree carries events from runs Dmitry rejects. Treat them as bad here
   // so spectrum and runLumi reject the same set.
   AppendBadRunsFromFile(cfg.workdir + "../lists/dmitry_extras.list");

   TCanvas c;
   c.Print((cfg.workdir + "qa.pdf[").c_str());

   for (const auto &R : cfg.jetRs) {
      std::vector<Spectrum> spectrums;
      spectrums.reserve(cfg.triggers.size());

      for (const auto &trig : cfg.triggers) {
         const std::string dataPath = Form("%smerged_data_R%s.root", cfg.datapath.c_str(), R.c_str());
         if (gSystem->AccessPathName(dataPath.c_str())) {
            std::cerr << "[warn] " << dataPath << " missing; skipping " << trig << " R=" << R << std::endl;
            continue;
         }
         spectrums.emplace_back(trig, R, syst);
         spectrums.back().Build();
      }

      WriteXsec(R, spectrums);
      if (syst.name == "nominal")
         CompareWithDmitriy(R, spectrums);
   }

   c.Print((cfg.workdir + "qa.pdf]").c_str());
}
