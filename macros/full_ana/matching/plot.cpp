#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iostream>
#include <string>
#include <vector>

const std::vector<int> colors = {2000, 2002, 2003, 2004,
                                 2005, 2006, 2007, 2008};

const std::vector<int> markers = {20, 21, 22, 23};

const float max_jet_pt = 90.0;

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
  latex->DrawLatex(0.2, 0.2, "#bf{STAR} #it{Raw data}");

  // Collision system
  latex->SetTextFont(42);
  latex->SetTextSize(0.05);
  latex->SetTextAlign(11);
  latex->DrawLatex(0.45, 0.9, "p+p #sqrt{s} = 200 GeV (year 2012)");

  // Additional text if provided

  latex->DrawLatex(0.5, 0.85, additionalText.Data());
}
// Plot a single trigger efficiency
void PlotSingleTrigger(const TString &trigger, const TString &jetR,
                       const TString &column = "reco",
                       const TString &addition = "") {

  gStyle->SetOptStat(0);

  TString outputDir =
      Form("./plots/embedding_root_%s_R%s", trigger.Data(), jetR.Data());
  TString fileName = Form("%s/%s_over_reconstructed%s.root", outputDir.Data(),
                          trigger.Data(), addition.Data());

  // Create output directory
  gSystem->mkdir(outputDir.Data(), kTRUE);

  TFile *file = TFile::Open(fileName.Data());
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << fileName << std::endl;
    return;
  }

  TH1D *h_num = (TH1D *)file->Get("h_num");
  TH1D *h_denom = (TH1D *)file->Get("h_denom");
  TH1D *h_ratio = (TH1D *)file->Get("h_ratio");

  if (!h_num || !h_denom || !h_ratio) {
    std::cerr << "Error: Cannot find histograms in file" << std::endl;
    file->Close();
    return;
  }

  // Create canvas
  TCanvas *c =
      new TCanvas(Form("c_%s_%s", trigger.Data(), column.Data()),
                  Form("%s Trigger Efficiency", trigger.Data()), 900, 800);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetRightMargin(0.045);
  c->SetTopMargin(0.045);
  c->SetFrameBorderMode(0);

  // Create two pads
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.4, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetTopMargin(0.08);
  pad1->SetRightMargin(0.045);
  pad1->SetLogy();
  pad1->SetFrameBorderMode(0);
  pad1->Draw();

  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.4);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.25);
  pad2->SetRightMargin(0.045);
  pad2->SetGridy();
  pad2->SetFrameBorderMode(0);
  pad2->Draw();

  // ===== Upper pad: Spectra =====
  pad1->cd();

  // Style numerator
  h_num->SetStats(0);
  h_num->SetLineColor(colors[0]);
  h_num->SetMarkerColor(colors[0]);
  h_num->SetMarkerStyle(20);
  h_num->SetMarkerSize(1.0);
  h_num->SetLineWidth(2);

  // Style denominator
  h_denom->SetLineColor(colors[1]);
  h_denom->SetMarkerColor(colors[1]);
  h_denom->SetMarkerStyle(21);
  h_denom->SetMarkerSize(1.0);
  h_denom->SetLineWidth(2);

  // Adjust axes for upper pad
  h_num->GetXaxis()->SetLabelSize(0);
  h_num->GetXaxis()->SetTitleSize(0);
  h_num->GetYaxis()->SetTitleSize(0.055);
  h_num->GetYaxis()->SetLabelSize(0.055);
  h_num->GetYaxis()->SetTitleOffset(0.9);
  h_num->GetXaxis()->SetRangeUser(0, max_jet_pt);
  h_num->GetYaxis()->SetRangeUser(5e-12, 1e-3);

  h_num->GetYaxis()->SetTitle("Normalized counts");

  h_num->Draw("E1");
  h_denom->Draw("E1 SAME");

  // Legend
  TLegend *leg1 = new TLegend(0.55, 0.55, 0.85, 0.75);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.045);
  leg1->AddEntry(h_num, Form("%s matched", trigger.Data()), "lep");
  leg1->AddEntry(h_denom, "All jets", "lep");
  leg1->Draw();

  // Add STAR labels
  AddSTARLabels(pad1, Form("Anti-k_{T}, R = %s", jetR.Data()));

  // ===== Lower pad: Ratio =====
  pad2->cd();

  // Style ratio
  h_ratio->SetStats(0);
  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetMarkerStyle(20);
  h_ratio->SetMarkerSize(1.0);
  h_ratio->SetLineWidth(2);

  // Adjust axes for lower pad
  h_ratio->GetXaxis()->SetTitle(
      Form("%s jet p_{T} [GeV/#it{c}]", column.Data()));
  h_ratio->GetXaxis()->SetTitleSize(0.08);
  h_ratio->GetXaxis()->SetLabelSize(0.08);
  h_ratio->GetXaxis()->SetTitleOffset(1.1);

  h_ratio->GetYaxis()->SetTitle(Form("%s Efficiency", trigger.Data()));
  h_ratio->GetYaxis()->SetTitleSize(0.08);
  h_ratio->GetYaxis()->SetLabelSize(0.08);
  h_ratio->GetYaxis()->SetTitleOffset(0.6);
  h_ratio->GetYaxis()->SetNdivisions(505);

  h_ratio->GetXaxis()->SetRangeUser(0, max_jet_pt);
  h_ratio->GetYaxis()->SetRangeUser(0, 1.1);

  h_ratio->Draw("E1");

  // Add line at y=1
  TLine *line = new TLine(5, 1.0, max_jet_pt, 1.0);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->SetLineWidth(2);
  line->Draw("same");

  c->cd();
  c->SaveAs(Form("%s/%s_efficiency_%s.pdf", outputDir.Data(), trigger.Data(),
                 column.Data()));

  file->Close();
  delete file;
}

// Compare multiple triggers or conditions
void PlotMultipleTriggers(const TString &jetR) {

  gStyle->SetOptStat(0);

  std::vector<TString> triggers = {"JP2", "HT2"};
  std::vector<TString> columns = {"reco"};

  for (const auto &column : columns) {

    TCanvas *c = new TCanvas(Form("c_compare_%s", column.Data()),
                             "Trigger Comparison", 800, 600);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetRightMargin(0.045);
    c->SetTopMargin(0.045);
    c->SetFrameBorderMode(0);

    TLegend *leg_ratio = new TLegend(0.75, 0.28, 0.85, 0.48);
    leg_ratio->SetBorderSize(0);
    leg_ratio->SetFillStyle(0);
    leg_ratio->SetTextFont(42);
    leg_ratio->SetTextSize(0.04);

    bool firstDraw = true;

    for (size_t i = 0; i < triggers.size(); i++) {
      TString trigger = triggers[i];
      TString outputDir =
          Form("./plots/embedding_root_%s_R%s", trigger.Data(), jetR.Data());
      TString fileName = Form("%s/%s_over_reconstructed.root", outputDir.Data(),
                              trigger.Data());

      TFile *file = TFile::Open(fileName.Data());
      if (!file || file->IsZombie()) {
        std::cerr << "Warning: Cannot open file " << fileName << std::endl;
        continue;
      }

      TH1D *h_num = (TH1D *)file->Get("h_num");
      TH1D *h_denom = (TH1D *)file->Get("h_denom");
      TH1D *h_ratio = (TH1D *)file->Get("h_ratio");

      if (!h_num || !h_denom || !h_ratio) {
        std::cerr << "Warning: Cannot find histograms in file " << fileName
                  << std::endl;
        file->Close();
        continue;
      }

      // Clone to avoid issues when file closes
      h_num = (TH1D *)h_num->Clone(Form("h_num_%s", trigger.Data()));
      h_denom = (TH1D *)h_denom->Clone(Form("h_denom_%s", trigger.Data()));
      h_ratio = (TH1D *)h_ratio->Clone(Form("h_ratio_%s", trigger.Data()));
      h_num->SetDirectory(0);
      h_denom->SetDirectory(0);
      h_ratio->SetDirectory(0);

      file->Close();

      h_ratio->SetLineColor(colors[i]);
      h_ratio->SetMarkerColor(colors[i]);
      h_ratio->SetMarkerStyle(markers[i]);
      h_ratio->SetLineWidth(2);

      // Style numerator
      h_ratio->SetStats(0);
      h_ratio->SetMarkerStyle(20);
      h_ratio->SetMarkerSize(1.0);
      h_ratio->SetLineWidth(2);

      // Adjust axes for upper pad
      h_ratio->GetYaxis()->SetTitleSize(0.055);
      h_ratio->GetYaxis()->SetLabelSize(0.05);
      h_ratio->GetYaxis()->SetTitleOffset(0.9);

      h_ratio->SetStats(0);
      h_ratio->GetXaxis()->SetTitle(
          Form("%s jet p_{T} [GeV/#it{c}]", column.Data()));
      h_ratio->GetXaxis()->SetTitleSize(0.055);
      h_ratio->GetXaxis()->SetLabelSize(0.05);
      h_ratio->GetXaxis()->SetTitleOffset(1.1);

      h_ratio->GetYaxis()->SetTitle(
          Form("%s Trigger Efficiency", trigger.Data()));

      h_ratio->GetXaxis()->SetRangeUser(0, max_jet_pt);
      h_ratio->GetYaxis()->SetRangeUser(0, 1.2);

      if (i == 0) {
        h_ratio->Draw("E1");
      } else {
        h_ratio->Draw("E1 SAME");
      }

      leg_ratio->AddEntry(h_ratio, trigger.Data(), "lep");
    }

    // Draw legends and labels
    c->cd();
    leg_ratio->Draw();
    AddSTARLabels(c, Form("Anti-k_{T}, R = %s", jetR.Data()));

    // Add line at y=1
    // get max and min x from current drawing
    TLine *line = new TLine(3, 1.0, max_jet_pt + 1, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray + 2);
    line->SetLineWidth(2);
    line->Draw("same");

    c->cd();
    TString outputDir = "./plots";
    c->SaveAs(Form("%s/comparison_all_triggers_%s_R%s.pdf", outputDir.Data(),
                   column.Data(), jetR.Data()));
  }
}

void PlotMultipleJetRs(const TString &trigger, const TString &column = "reco") {

  gStyle->SetOptStat(0);

  std::vector<TString> jetRs = {"0.2", "0.3", "0.4", "0.5", "0.6"};

  TCanvas *c = new TCanvas(Form("c_compare_R_%s", trigger.Data()),
                           "Jet R Comparison", 800, 600);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetRightMargin(0.045);
  c->SetTopMargin(0.045);
  c->SetFrameBorderMode(0);

  TLegend *leg_ratio = new TLegend(0.75, 0.28, 0.85, 0.48);
  leg_ratio->SetBorderSize(0);
  leg_ratio->SetFillStyle(0);
  leg_ratio->SetTextFont(42);
  leg_ratio->SetTextSize(0.04);

  bool firstDraw = true;

  for (size_t i = 0; i < jetRs.size(); i++) {
    TString jetR = jetRs[i];
    TString outputDir =
        Form("./plots/embedding_root_%s_R%s", trigger.Data(), jetR.Data());
    TString fileName =
        Form("%s/%s_over_reconstructed.root", outputDir.Data(), trigger.Data());

    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
      std::cerr << "Warning: Cannot open file " << fileName << std::endl;
      continue;
    }

    TH1D *h_ratio = (TH1D *)file->Get("h_ratio");
    if (!h_ratio) {
      std::cerr << "Warning: Cannot find h_ratio in file " << fileName
                << std::endl;
      file->Close();
      continue;
    }

    h_ratio = (TH1D *)h_ratio->Clone(Form("h_ratio_R%s", jetR.Data()));
    h_ratio->SetDirectory(0);
    file->Close();

    h_ratio->SetLineColor(colors[i]);
    h_ratio->SetMarkerColor(colors[i]);
    h_ratio->SetMarkerStyle(markers[i % markers.size()]);
    h_ratio->SetLineWidth(2);
    h_ratio->SetMarkerSize(1.0);
    h_ratio->SetStats(0);

    h_ratio->GetXaxis()->SetTitle(
        Form("%s jet p_{T} [GeV/#it{c}]", column.Data()));
    h_ratio->GetXaxis()->SetTitleSize(0.055);
    h_ratio->GetXaxis()->SetLabelSize(0.05);
    h_ratio->GetXaxis()->SetTitleOffset(1.1);
    h_ratio->GetYaxis()->SetTitle(
        Form("%s Trigger Efficiency", trigger.Data()));
    h_ratio->GetYaxis()->SetTitleSize(0.055);
    h_ratio->GetYaxis()->SetLabelSize(0.05);
    h_ratio->GetYaxis()->SetTitleOffset(0.9);

    h_ratio->GetXaxis()->SetRangeUser(0, max_jet_pt);
    h_ratio->GetYaxis()->SetRangeUser(0, 1.2);

    if (firstDraw) {
      h_ratio->Draw("E1");
      firstDraw = false;
    } else {
      h_ratio->Draw("E1 SAME");
    }

    leg_ratio->AddEntry(h_ratio, Form("R = %s", jetR.Data()), "lep");
  }

  c->cd();
  leg_ratio->Draw();
  AddSTARLabels(c, Form("Anti-k_{T}, jets %s trigger", trigger.Data()));

  TLine *line = new TLine(3, 1.0, max_jet_pt + 1, 1.0);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->SetLineWidth(2);
  line->Draw("same");

  TString outputDir = "./plots";
  gSystem->mkdir(outputDir.Data(), kTRUE);
  c->SaveAs(Form("%s/comparison_all_Rs_%s_trigger%s.pdf", outputDir.Data(),
                 column.Data(), trigger.Data()));
}

// Main plotting function
void plot() {

  std::vector<TString> triggers = {"HT2", "JP2"};
  std::vector<TString> jetRs = {"0.2", "0.3", "0.4", "0.5", "0.6"};

  for (const auto &jetR : jetRs) {
    for (const auto &trigger : triggers) {
      // Plot individual triggers
      cout << "Plotting trigger: " << trigger << " with R=" << jetR
           << std::endl;
      PlotSingleTrigger(trigger, jetR, "reco");
    }
    PlotMultipleTriggers(jetR);
  }

  for (const auto &trigger : triggers) {
    // Plot individual triggers with addition
    cout << "Plotting trigger: " << trigger << " with R=0.6 and addition"
         << std::endl;
    PlotMultipleJetRs(trigger, "reco");
  }

  std::cout << "Done! Plots saved in ./plots/embedding_root_*/" << std::endl;
}