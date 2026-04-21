#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "/home/prozorov/dev/star/jets_pp_2012/new_ana/config.h"
using namespace CrossSectionConfig;

AnalysisConfig cfg;

static std::unique_ptr<TH1D> load(const std::string &fileName,
                                  const char *histName, const char *newName) {
  TFile f(fileName.c_str(), "READ");
  if (f.IsZombie()) {
    throw std::runtime_error("Cannot open file: " + fileName);
  }

  TH1D *src = dynamic_cast<TH1D *>(f.Get(histName));
  if (!src) {
    throw std::runtime_error("Histogram '" + std::string(histName) +
                             "' not found (or not TH1D) in file: " + fileName);
  }
  TH1D *c =
      static_cast<TH1D *>(src->Clone(newName)); // Clone allocates new hist
  c->SetDirectory(nullptr);                     // detach from file
  return std::unique_ptr<TH1D>(c);              // take ownership
}

class Spectrum {
public:
  Spectrum(std::string trig, std::string R)
      : trigger(std::move(trig)), jetR(std::move(R)) {
    // //=================== load data ===============
    dataFileName = Form("%s/merged_data_%s_R%s.root", cfg.datapath.c_str(),
                        trigger.c_str(), jetR.c_str());

    std::string triggEffName = cfg.workdir + "matching/plots/embedding_root_" +
                               trigger + "_R" + jetR + "/" + trigger +
                               "_over_all.root";

    trigEff = load(triggEffName, "h_ratio", "trigEff");
    //=================== load fake rate ===============
    std::string fakeRateName = cfg.workdir + "matching/plots/embedding_root_" +
                               trigger + "_R" + jetR +
                               "/fakereco_over_all.root";
    fakeRate = load(fakeRateName, "h_ratio", "fakeRate");
    // make 1-fakeRate
    for (int i = 1; i <= fakeRate->GetNbinsX(); ++i) {
      const double v = fakeRate->GetBinContent(i);
      fakeRate->SetBinContent(i, 1.0 - v);
    }

    //=================== load matching efficiency ===============
    std::string matchEffName = cfg.workdir + "matching/plots/embedding_root_" +
                               trigger + "_R" + jetR + "/matched_over_mc.root";
    matchEff = load(matchEffName, "h_ratio", "matchEff");
    //=================== load response matrix ===============
    TFile fileUnfolding(
        Form((cfg.workdir + "unfolding/response_%s_R%s.root").c_str(),
             trigger.c_str(), jetR.c_str()),
        "READ");
    response = (RooUnfoldResponse *)fileUnfolding.Get("my_response");
  }

  void Dump(TH1 *h, TCanvas *c, const std::string &msg = "") {
    c->cd();
    c->SetLogy();

    TH1D *hClone = (TH1D *)h->Clone(Form("%s_clone", h->GetName()));
    hClone->SetMinimum(1);
    hClone->SetMaximum(0.5e8);

    const bool hasFrame = (c->GetListOfPrimitives()->GetSize() > 0);
    // create clone of histogram

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

  void Build() {
    std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>();
    //=================== raw ===============
    h = Raw(dataFileName);
    h->SetLineColor(2001);
    Dump(h.get(), c.get(), "raw");

    //=================== apply trigger efficiency ===============
    h->Divide(trigEff.get());
    h->SetLineColor(2002);
    Dump(h.get(), c.get(), "after trig eff");

    //=================== apply fake rate ===============
    h->Multiply(fakeRate.get());
    h->SetLineColor(2003);
    Dump(h.get(), c.get(), "after fake rate");

    //=================== unfold ===============
    Unfold();
    h->SetLineColor(2004);
    Dump(h.get(), c.get(), "after unfolding");
    //=================== apply matching efficiency ===============

    h->Divide(matchEff.get());
    h->SetLineColor(2005);
    Dump(h.get(), c.get(), "after match eff");

    //=================== normalize ===============
    NormBinWidthAndAcc();
    const double LeffTot = getRunLumi()->Integral();
    h->Scale(1.0 / LeffTot);

    c->SaveAs((cfg.workdir + "qa.pdf").c_str());
  }

  void Unfold() {
    RooUnfoldBayes u(response, h.get(), cfg.nIterations);
    TH1D *hu = (TH1D *)u.Hunfold();

    TH1D *c = (TH1D *)hu->Clone(
        Form("%s_unfolded_R%s", trigger.c_str(), jetR.c_str()));
    c->SetDirectory(0);
    h.reset(c);
  }

  TH1D *Hist() const { return h.get(); }
  const std::string trigger, jetR;

private:
  std::unique_ptr<TH1D> h;
  std::string dataFileName;
  std::unique_ptr<TH1D> fakeRate;
  std::unique_ptr<TH1D> trigEff;
  std::unique_ptr<TH1D> matchEff;
  RooUnfoldResponse *response;

  std::unique_ptr<TH1D> Raw(const std::string &dataFileName) const {
    ROOT::RDataFrame df("ResultTree", dataFileName.c_str());
    const std::string br =
        "trigger_match_" + std::string(trigger == "JP2" ? "JP2" : "HT2");

    auto htmp =
        df.Define("m", br)
            .Filter(
                [](int runIndex) {
                  // Skip runs that are unknown or marked as bad
                  return runIndex >= 0 &&
                         std::find(cfg.badRuns.begin(), cfg.badRuns.end(),
                                   runIndex) == cfg.badRuns.end();
                },
                {"runid1"})
            .Define("pt_trig", "pt[m]")
            .Filter("pt_trig.size() > 0")
            .Histo1D({Form("%s_raw_R%s", trigger.c_str(), jetR.c_str()),
                      Form("%s raw R=%s;p_{T};counts", trigger.c_str(),
                           jetR.c_str()),
                      (int)pt_reco_bins.size() - 1, pt_reco_bins.data()},
                     "pt_trig");

    TH1D *c = (TH1D *)htmp->Clone(
        Form("%s_raw_clone_R%s", trigger.c_str(), jetR.c_str()));
    c->SetDirectory(0);
    return std::unique_ptr<TH1D>(c);
  }

  // std::unique_ptr<TH2D> Raw2D(const std::string &dataFileName) const {
  //   ROOT::RDataFrame df_base("ResultTree", dataFileName.c_str());
  //   const std::string br =
  //       "trigger_match_" + std::string(trigger == "JP2" ? "JP2" : "HT2");

  //   const int nRunBins = static_cast<int>(cfg.runMap.size());

  //   auto htmp =
  //       df_base
  //           .Define("runIndex",
  //                   [](int runNumber) { return cfg.runMap.at(runNumber) - 1;
  //                   },
  //                   {"runid1"})
  //           .Define("trigger", br)
  //           .Define("pt_trig", "pt[trigger]")
  //           .Filter("pt_trig.size() > 0")
  //           .Histo2D({Form("%s_pt_vs_run_trigger", trigger.c_str()),
  //                     Form("%s pt vs run;run bin;p_{T} [GeV/c];counts",
  //                          trigger.c_str()),
  //                     nRunBins, 0.0, double(nRunBins),
  //                     (int)pt_reco_bins.size() - 1, pt_reco_bins.data()},
  //                    "runIndex", "pt_trig");

  //   TH2D *c = (TH2D *)htmp->Clone(
  //       Form("%s_raw_clone_R%s", trigger.c_str(), jetR.c_str()));
  //   c->SetDirectory(0);
  //   return std::unique_ptr<TH2D>(c);
  // }

  void NormBinWidthAndAcc() {
    // scale by eta acceptance
    h->Scale(1.0 / (2.0 * (1.0 - std::stod(jetR))));
    // scale by pt bin width
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      const double w = h->GetBinWidth(i);
      h->SetBinContent(i, h->GetBinContent(i) / w);
      h->SetBinError(i, h->GetBinError(i) / w);
    }
  }

  TH1D *getRunLumi() {
    TFile lf((cfg.workdir + "lumi.root").c_str(), "READ");

    auto *maxn = (TH1D *)lf.Get(Form("nevents_%s", trigger.c_str()));
    auto *lumi = (TH1D *)lf.Get(Form("luminosity_%s", trigger.c_str()));

    TFile df(dataFileName.c_str(), "READ");
    // auto *maxn = (TH1D *)lf.Get(Form("nevents_JP2", trigger.c_str()));
    // auto *lumi = (TH1D *)lf.Get(Form("luminosity_JP2", trigger.c_str()));
    // TFile df(
    //     Form("%s/merged_data_JP2_R%s.root", cfg.datapath.c_str(),
    //     jetR.c_str()), "READ");

    auto *analyzed = (TH1D *)df.Get("hEventsRun");

    TH1D *runLumi = (TH1D *)analyzed->Clone("runLumi");

    runLumi->SetDirectory(0);
    runLumi->Divide(maxn);

    runLumi->Multiply(lumi);

    for (int i = 1; i <= runLumi->GetNbinsX(); ++i) {
      const int run = std::stoi(runLumi->GetXaxis()->GetBinLabel(i));
      if (std::find(cfg.badRuns.begin(), cfg.badRuns.end(), run) !=
          cfg.badRuns.end()) {
        runLumi->SetBinContent(i, 0.0);
        runLumi->SetBinError(i, 0.0);
      }
    }

    return runLumi;
  }
};

static void CompareWithDmitriy(const std::string &jetR,
                               const std::vector<Spectrum> &spectrums,
                               const std::string &suffix = "") {

  TFile rf(Form((cfg.workdir + "jet_cross_section_dmitriyR%s.root").c_str(),
                jetR.c_str()),
           "READ");
  TH1D *ref = (TH1D *)rf.Get("crossSection_systematic")
                  ->Clone(Form("ref_R%s", jetR.c_str()));
  ref->SetDirectory(0);

  TCanvas *c = new TCanvas(Form("c_R%s", jetR.c_str()), "", 900, 900);
  c->Divide(1, 2);

  std::vector<int> col = {kAzure - 1, kRed + 1, kOrange - 1, kGreen + 2,
                          kMagenta + 1};
  std::vector<int> mar = {20, 21, 22, 33, 34};

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
    // limit histogram to the ref range
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      const double x = h->GetXaxis()->GetBinCenter(bin);
      if (x < ref->GetXaxis()->GetXmin() || x > ref->GetXaxis()->GetXmax()) {
        h->SetBinContent(bin, 0.0);
        h->SetBinError(bin, 0.0);
      }
    }

    h->SetMarkerStyle(mar[i % mar.size()]);
    h->SetMarkerColor(col[i % col.size()]);
    h->SetLineColor(col[i % col.size()]);
    h->Draw("E1 SAME");
    leg->AddEntry(h, spectrums[i].trigger.c_str(), "lep");
  }
  leg->Draw();

  auto *p2 = (TPad *)c->cd(2);
  p2->SetPad(0, 0, 1, 0.5);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.30);

  // >>> added: draw an empty frame first (axes), then the unity band behind
  // points
  TH1D *frame = (TH1D *)ref->Clone(Form("ratio_frame_R%s", jetR.c_str()));
  frame->Reset();
  frame->SetDirectory(0);
  frame->SetTitle(
      Form("Ratio to Dmitriy R=%s; p_{T} [GeV/c]; Ratio", jetR.c_str()));
  frame->GetYaxis()->SetRangeUser(0.0, 2.5);
  frame->Draw(); // axes only

  // >>> added: unity band with ref relative uncertainties: 1 ± (eref/ref)
  TH1D *refBand = (TH1D *)ref->Clone(Form("ref_unity_band_R%s", jetR.c_str()));
  refBand->Reset();
  refBand->SetDirectory(0);

  for (int b = 1; b <= refBand->GetNbinsX(); ++b) {
    const double cRef = ref->GetBinContent(b);
    const double eRef = ref->GetBinError(b);
    const double rel =
        (cRef > 0.0) ? (eRef / cRef) : 0.0; // protect div-by-zero

    refBand->SetBinContent(b, 1.0);
    refBand->SetBinError(b, rel);
  }

  refBand->SetFillColorAlpha(kViolet, 0.20);
  refBand->SetFillStyle(1001);
  refBand->SetLineColor(kViolet);
  refBand->SetMarkerSize(0);
  refBand->Draw("E2 SAME");

  // >>> added: central unity line (draw now so it stays visible)
  TLine *line = new TLine(ref->GetXaxis()->GetXmin(), 1.0,
                          ref->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(kViolet);
  line->SetLineWidth(2);
  line->Draw("SAME");

  // --- your ratios on top ---
  for (size_t i = 0; i < spectrums.size(); ++i) {
    TH1D *hist = (TH1D *)spectrums[i].Hist()->Clone(
        Form("ratio_%s_R%s", spectrums[i].trigger.c_str(), jetR.c_str()));

    // create a new histogram with binning of ref and fill bins with values from
    // hist
    TH1D *r = new TH1D(
        Form("ratio_binned_%s_R%s", spectrums[i].trigger.c_str(), jetR.c_str()),
        Form("Ratio to Dmitriy R=%s; p_{T} [GeV/c]; Ratio", jetR.c_str()),
        ref->GetNbinsX(), ref->GetXaxis()->GetXbins()->GetArray());

    for (int bin = 1; bin <= r->GetNbinsX(); ++bin) {
      const double x = r->GetXaxis()->GetBinCenter(bin);
      const int binHist = hist->FindBin(x);
      r->SetBinContent(bin, hist->GetBinContent(binHist));
      r->SetBinError(bin, hist->GetBinError(binHist));
    }

    r->SetDirectory(0);

    r->Divide(ref);

    r->SetMarkerStyle(mar[i % mar.size()]);
    r->SetMarkerColor(col[i % col.size()]);
    r->SetLineColor(col[i % col.size()]);
    r->GetYaxis()->SetRangeUser(0.0, 2.5);
    r->GetYaxis()->SetTitle("Ratio to Dmitriy");
    r->GetXaxis()->SetTitle("p_{T} [GeV/c]");

    r->Draw("E1 SAME");
  }

  c->SaveAs(Form("%s/comparison_with_dmitriy_R%s%s.pdf", cfg.workdir.c_str(),
                 jetR.c_str(), suffix.c_str()));
}

static void compare(const Spectrum &spectrum1, const Spectrum &spectrum2) {
  TCanvas *c =
      new TCanvas(Form("c_%s_vs_%s_R%s", spectrum1.trigger.c_str(),
                       spectrum2.trigger.c_str(), spectrum1.jetR.c_str()),
                  "", 900, 900);
  c->Divide(1, 2);
  auto *p1 = (TPad *)c->cd(1);
  p1->SetPad(0, 0.5, 1, 1);
  p1->SetBottomMargin(0);
  p1->SetLogy();
  std::vector<int> col = {kAzure - 1, kRed + 1, kOrange - 1, kGreen + 2,
                          kMagenta + 1};
  std::vector<int> mar = {20, 21, 22, 33, 34};
  TH1D *ref = (TH1D *)spectrum2.Hist()->Clone(
      Form("ref_%s_R%s", spectrum2.trigger.c_str(), spectrum2.jetR.c_str()));

  ref->SetMarkerStyle(29);
  ref->SetMarkerSize(2);
  ref->SetMarkerColor(kViolet);
  ref->SetLineColor(kViolet);
  ref->GetYaxis()->SetRangeUser(1.2, 1e7);
  ref->Draw("E1");

  TLegend *leg = new TLegend(0.60, 0.50, 0.88, 0.88);
  leg->AddEntry(ref, spectrum2.trigger.c_str(), "lep");

  TH1D *h = spectrum1.Hist();
  h->SetMarkerStyle(mar[0]);
  h->SetMarkerColor(col[0]);
  h->SetLineColor(col[0]);
  h->Draw("E1 SAME");
  leg->AddEntry(h, spectrum1.trigger.c_str(), "lep");

  leg->Draw();

  auto *p2 = (TPad *)c->cd(2);
  p2->SetPad(0, 0, 1, 0.5);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.30);

  TH1D *r = (TH1D *)spectrum1.Hist()->Clone(
      Form("ratio_%s_vs_%s_R%s", spectrum1.trigger.c_str(),
           spectrum2.trigger.c_str(), spectrum1.jetR.c_str()));

  r->SetDirectory(0);

  r->Divide(ref);

  r->SetMarkerStyle(mar[0]);
  r->SetMarkerColor(col[0]);
  r->SetLineColor(col[0]);
  r->GetYaxis()->SetRangeUser(0.0, 2.5);
  r->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  r->Draw("E1");

  TLine *line = new TLine(ref->GetXaxis()->GetXmin(), 1.0,
                          ref->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(kViolet);
  line->Draw();

  c->SaveAs((cfg.workdir + "qa.pdf").c_str());
}

void cross_section() {
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  TCanvas c;
  c.Print((cfg.workdir + "qa.pdf[").c_str()); // start pdf

  for (const auto &R : cfg.jetRs) {
    std::vector<Spectrum> spectrums;
    spectrums.reserve(cfg.triggers.size());

    for (const auto &trig : cfg.triggers) {
      spectrums.emplace_back(trig, R);
      spectrums.back().Build();
    }

    CompareWithDmitriy(R, spectrums);

    if (spectrums.size() >= 2) {
      compare(spectrums[0], spectrums[1]);
    }
  }

  c.Print((cfg.workdir + "qa.pdf]").c_str()); // end pdf
}
