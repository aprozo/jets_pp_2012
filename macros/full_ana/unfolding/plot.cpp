#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TStyle.h"
#include <string>
#include <vector>

const std::vector<int> colors = {2000, 2002, 2003, 2001,
                                 2005, 2006, 2007, 2008};
const std::vector<int> markers = {20, 21, 22, 23, 33, 34};

const float max_jet_pt = 80.0;
const float min_jet_pt = 5.0;

void AddSTARLabels(TPad *pad, const TString &additionalText = "") {
  pad->cd();

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.05);

  // STAR label (if needed, uncomment)
  // latex->SetTextFont(62);
  // latex->DrawLatex(0.15, 0.895, "STAR");
  // latex->SetTextFont(42);

  // Raw label
  latex->SetTextFont(42);
  latex->SetTextSize(0.07);
  latex->DrawLatex(0.45, 0.7, "#bf{STAR} #it{Raw data}");

  // Collision system
  latex->SetTextFont(42);
  latex->SetTextSize(0.05);
  latex->SetTextAlign(11);
  latex->DrawLatex(0.45, 0.85, "p+p #sqrt{s} = 200 GeV (year 2012)");

  // Additional text if provided

  latex->DrawLatex(0.5, 0.79, additionalText.Data());
}

// Divide histogram by bin width
TH1D *divideByBinWidth(TH1D *hist) {
  TH1D *h = (TH1D *)hist->Clone();
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double binWidth = h->GetBinWidth(i);
    if (binWidth != 0) {
      h->SetBinContent(i, h->GetBinContent(i) / binWidth);
      h->SetBinError(i, h->GetBinError(i) / binWidth);
    }
  }
  return h;
}

// Plot unfolded results with comparison to truth
void PlotUnfoldedComparison(TString responseFile,
                            const std::vector<int> &iterations,
                            const std::string &trigger = "JP2",
                            const std::string &R = "0.2") {

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".2f");

  TFile *file = TFile::Open(responseFile);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << responseFile << std::endl;
    return;
  }

  RooUnfoldResponse *response = (RooUnfoldResponse *)file->Get("my_response");
  TH1D *hMeasured = (TH1D *)file->Get("MeasuredTest");
  TH1D *hTruth = (TH1D *)file->Get("TruthTest");

  if (!response || !hMeasured || !hTruth) {
    std::cerr << "Error: Cannot find required objects in file" << std::endl;
    file->Close();
    return;
  }

  TString outPdf =
      Form("pdf/unfolding_results_%s_R%s.pdf", trigger.c_str(), R.c_str());

  // ============ Plot 1: Spectra comparison ============
  TCanvas *c1 = new TCanvas("c1", "Unfolding Results", 900, 800);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetRightMargin(0.045);
  c1->SetTopMargin(0.045);
  c1->SaveAs(outPdf + "[");

  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetTopMargin(0.08);
  pad1->SetRightMargin(0.045);
  pad1->SetLogy();
  pad1->Draw();

  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.35);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.25);
  pad2->SetRightMargin(0.045);
  pad2->SetGridy();
  pad2->Draw();

  // Upper pad: Spectra
  pad1->cd();

  TH1D *hTruth_norm = divideByBinWidth(hTruth);
  TH1D *hMeasured_norm = divideByBinWidth(hMeasured);

  hTruth_norm->SetLineColor(kBlack);
  hTruth_norm->SetMarkerColor(kBlack);
  hTruth_norm->SetMarkerStyle(20);
  hTruth_norm->SetMarkerSize(1.0);
  hTruth_norm->SetLineWidth(2);

  hMeasured_norm->SetLineColor(colors[1]);
  hMeasured_norm->SetMarkerColor(colors[1]);
  hMeasured_norm->SetMarkerStyle(21);
  hMeasured_norm->SetMarkerSize(1.0);
  hMeasured_norm->SetLineWidth(2);

  hTruth_norm->GetXaxis()->SetLabelSize(0);
  hTruth_norm->GetYaxis()->SetTitle("Normalized Counts");
  hTruth_norm->GetYaxis()->SetTitleSize(0.055);
  hTruth_norm->GetYaxis()->SetLabelSize(0.055);
  hTruth_norm->GetYaxis()->SetTitleOffset(0.95);

  hMeasured_norm->GetXaxis()->SetRangeUser(min_jet_pt, max_jet_pt);
  hTruth_norm->GetXaxis()->SetRangeUser(min_jet_pt, max_jet_pt);

  hTruth_norm->GetYaxis()->SetRangeUser(hMeasured_norm->GetMinimum() * 5,
                                        hMeasured_norm->GetMaximum() * 2);
  hTruth_norm->Draw("PE ");

  hMeasured_norm->Draw("PE SAME");

  TLegend *leg = new TLegend(0.15, 0.12, 0.42, 0.42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->AddEntry(hTruth_norm, "MC Truth", "lep");
  leg->AddEntry(hMeasured_norm, "Measured", "lep");

  std::vector<TH1D *> unfolded_hists;
  std::vector<RooUnfoldBayes *> unfolders;

  for (size_t i = 0; i < iterations.size(); i++) {
    RooUnfoldBayes *unfold =
        new RooUnfoldBayes(response, hMeasured, iterations[i]);
    unfolders.push_back(unfold);

    TH1D *hUnf = (TH1D *)unfold->Hunfold();
    hUnf->SetName(Form("hUnfolded_iter%d", iterations[i]));
    TH1D *hUnf_norm = divideByBinWidth(hUnf);

    hUnf_norm->SetLineColor(colors[i + 2]);
    hUnf_norm->SetMarkerColor(colors[i + 2]);
    hUnf_norm->SetMarkerStyle(markers[i % markers.size()]);
    hUnf_norm->SetMarkerSize(1.2);
    hUnf_norm->SetLineWidth(2);

    hUnf_norm->Draw("PE SAME");
    unfolded_hists.push_back(hUnf);

    leg->AddEntry(hUnf_norm, Form("Unfolded (Iter %d)", iterations[i]), "lep");
  }

  leg->Draw();
  AddSTARLabels(pad1, Form("Anti-k_{T}, R = %s", R.c_str()));

  // Lower pad: Ratio
  pad2->cd();

  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);

  for (size_t i = 0; i < iterations.size(); i++) {
    TH1D *ratio =
        (TH1D *)unfolded_hists[i]->Clone(Form("ratio_iter%d", iterations[i]));
    ratio->Divide(hTruth);

    ratio->SetLineColor(colors[i + 2]);
    ratio->SetMarkerColor(colors[i + 2]);
    ratio->SetMarkerStyle(markers[i % markers.size()]);
    ratio->SetMarkerSize(1.2);
    ratio->SetLineWidth(2);
    ratio->SetTitle("");

    ratio->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
    ratio->GetXaxis()->SetTitleSize(0.12);
    ratio->GetXaxis()->SetLabelSize(0.10);
    ratio->GetXaxis()->SetTitleOffset(0.87);

    ratio->GetYaxis()->SetTitle("Unfolded/Truth");
    ratio->GetYaxis()->SetTitleSize(0.10);
    ratio->GetYaxis()->SetLabelSize(0.1);
    ratio->GetYaxis()->SetTitleOffset(0.5);
    ratio->GetYaxis()->SetRangeUser(0.7, 1.3);
    ratio->GetYaxis()->SetNdivisions(505);

    ratio->GetXaxis()->SetRangeUser(min_jet_pt, max_jet_pt);
    ratio->Draw(i == 0 ? "PE" : "PE SAME");

    if (i == 0) {
      line->DrawLine(ratio->GetXaxis()->GetXmin(), 1.0,
                     ratio->GetXaxis()->GetXmax(), 1.0);
      line->DrawLine(ratio->GetXaxis()->GetXmin(), 1.05,
                     ratio->GetXaxis()->GetXmax(), 1.05);
      line->DrawLine(ratio->GetXaxis()->GetXmin(), 0.95,
                     ratio->GetXaxis()->GetXmax(), 0.95);
    }
  }

  c1->SaveAs(outPdf);

  // ============ Plot 2: Response Matrix ============
  TCanvas *c2 = new TCanvas("c2", "Response Matrix", 900, 800);
  c2->SetRightMargin(0.15);
  c2->SetLogz();

  TH2D *hResponse = (TH2D *)file->Get("hResponseMatrix");
  if (hResponse) {
    hResponse->SetTitle("");
    hResponse->GetXaxis()->SetTitle("p_{T}^{reco} [GeV/#it{c}]");
    hResponse->GetYaxis()->SetTitle("p_{T}^{MC} [GeV/#it{c}]");
    hResponse->GetZaxis()->SetTitle("Counts");
    hResponse->GetXaxis()->SetTitleSize(0.045);
    hResponse->GetYaxis()->SetTitleSize(0.045);
    hResponse->GetZaxis()->SetTitleSize(0.045);
    hResponse->GetXaxis()->SetTitleOffset(1.1);
    hResponse->GetYaxis()->SetTitleOffset(1.2);

    hResponse->Draw("COLZ");

    AddSTARLabels((TPad *)c2, Form("Anti-k_{T}, R = %s", R.c_str()));

    c2->SaveAs(outPdf);
  }

  // ============ Plot 3: Detector Resolution ============
  TCanvas *c3 = new TCanvas("c3", "Detector Resolution", 900, 800);
  c3->SetRightMargin(0.15);

  TH1D *hMean = (TH1D *)file->Get("detectorResolution_1");
  TH1D *hSigma = (TH1D *)file->Get("detectorResolution_2");

  if (hMean && hSigma) {
    c3->Divide(1, 2);

    c3->cd(1);
    gPad->SetRightMargin(0.045);
    hMean->SetTitle("");
    hMean->GetXaxis()->SetTitle("p_{T}^{MC} [GeV/#it{c}]");
    hMean->GetYaxis()->SetTitle(
        "#LT p_{T}^{reco} - p_{T}^{MC} #GT [GeV/#it{c}]");
    hMean->SetLineColor(colors[0]);
    hMean->SetMarkerColor(colors[0]);
    hMean->SetMarkerStyle(20);
    hMean->SetLineWidth(2);
    hMean->Draw("PE");
    AddSTARLabels((TPad *)gPad, Form("Anti-k_{T}, R = %s", R.c_str()));

    c3->cd(2);
    gPad->SetRightMargin(0.045);
    hSigma->SetTitle("");
    hSigma->GetXaxis()->SetTitle("p_{T}^{MC} [GeV/#it{c}]");
    hSigma->GetYaxis()->SetTitle(
        "#sigma(p_{T}^{reco} - p_{T}^{MC}) [GeV/#it{c}]");
    hSigma->SetLineColor(colors[1]);
    hSigma->SetMarkerColor(colors[1]);
    hSigma->SetMarkerStyle(21);
    hSigma->SetLineWidth(2);
    hSigma->Draw("PE");

    c3->SaveAs(outPdf);
  }

  c1->SaveAs(outPdf + "]");

  file->Close();

  // Cleanup
  for (auto unfolder : unfolders)
    delete unfolder;
}

// Main plotting function
void plot() {
  std::vector<int> iterations = {1, 4};

  std::vector<std::string> triggers = {"HT2", "JP2"};
  std::vector<std::string> jetRs = {"0.2", "0.3", "0.4", "0.5", "0.6"};

  for (const auto &jetR : jetRs) {
    for (const auto &trigger : triggers) {
      TString responseFile =
          Form("response_%s_R%s.root", trigger.c_str(), jetR.c_str());
      // Plot individual triggers
      if (!responseFile)
        continue;
      cout << "Plotting trigger: " << trigger << " with R=" << jetR
           << std::endl;

      PlotUnfoldedComparison(responseFile, iterations, trigger, jetR);
    }
  }
}