#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <iostream>
#include <string>
#include <vector>

const std::vector<int> colors = {2000, 2002, 2003, 2004,
                                 2005, 2006, 2007, 2008};

const std::vector<int> markers = {20, 21, 22, 23};

const float max_jet_pt = 80.0;

//// Add STAR publication style labels
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
  latex->DrawLatex(0.15, 0.2, "#bf{STAR} #it{Raw data}");

  // Collision system
  latex->SetTextFont(42);
  latex->SetTextSize(0.05);
  latex->SetTextAlign(11);
  latex->DrawLatex(0.45, 0.85, "p+p #sqrt{s} = 200 GeV (year 2012)");

  // Additional text if provided

  latex->DrawLatex(0.5, 0.79, additionalText.Data());
}

// Compare multiple triggers or conditions
void plotR(const TString &jetR) {

  gStyle->SetOptStat(0);

  std::vector<TString> triggers = {"JP2", "HT2"};

  TCanvas *c = new TCanvas(Form("c_compare_%s", jetR.Data()),
                           "Trigger Comparison", 800, 600);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetRightMargin(0.045);
  c->SetTopMargin(0.045);
  c->SetFrameBorderMode(0);
  c->SetLogy();

  TLegend *leg_ratio = new TLegend(0.75, 0.55, 0.85, 0.75);
  leg_ratio->SetBorderSize(0);
  leg_ratio->SetFillStyle(0);
  leg_ratio->SetTextFont(42);
  leg_ratio->SetTextSize(0.04);

  bool firstDraw = true;

  for (size_t i = 0; i < triggers.size(); i++) {
    TString trigger = triggers[i];
    TString fileName = Form("./data_%s_R%s.root", trigger.Data(), jetR.Data());

    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
      std::cerr << "Warning: Cannot open file " << fileName << std::endl;
      continue;
    }
    TH1D *hist = (TH1D *)file->Get("rawJetsTriggered");
    hist->SetDirectory(0);

    file->Close();

    hist->SetLineColor(colors[i]);
    hist->SetMarkerColor(colors[i]);
    hist->SetMarkerStyle(markers[i]);
    hist->SetLineWidth(2);

    // Adjust axes for upper pad
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleOffset(0.9);

    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("jet p_{T} [GeV/#c]");
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleOffset(0.9);

    hist->GetXaxis()->SetRangeUser(5, 60);

    hist->Draw(i == 0 ? "E1" : "E1 same");

    leg_ratio->AddEntry(hist, trigger.Data(), "lep");
  }

  // Draw legends and labels
  c->cd();
  leg_ratio->Draw();
  AddSTARLabels(c, Form("Anti-k_{T}, R = %s", jetR.Data()));

  c->cd();
  c->SaveAs("dataR" + jetR + ".pdf");
}

// Compare multiple jet R values for a given trigger
void plotTrigger(const TString &trigger) {
  gStyle->SetOptStat(0);
  std::vector<TString> jetRs = {"0.2", "0.3", "0.4", "0.5", "0.6"};
  TCanvas *c = new TCanvas(Form("c_compare_%s", trigger.Data()),
                           "Jet R Comparison", 800, 600);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetRightMargin(0.045);
  c->SetTopMargin(0.045);
  c->SetFrameBorderMode(0);
  c->SetLogy();

  TLegend *leg = new TLegend(0.75, 0.55, 0.85, 0.75);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  bool firstDraw = true;
  for (size_t i = 0; i < jetRs.size(); i++) {
    TString jetR = jetRs[i];
    TString fileName = Form("./data_%s_R%s.root", trigger.Data(), jetR.Data());
    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
      std::cerr << "Warning: Cannot open file " << fileName << std::endl;
      continue;
    }

    TH1D *hist = (TH1D *)file->Get("rawJetsTriggered");
    hist->SetDirectory(0);
    file->Close();

    hist->SetLineColor(colors[i]);
    hist->SetMarkerColor(colors[i]);
    hist->SetMarkerStyle(markers[i]);
    hist->SetLineWidth(2);

    // Adjust axes
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleOffset(0.9);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("jet p_{T} [GeV/#it{c}]");
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleOffset(1);
    hist->GetXaxis()->SetRangeUser(5, 60);
    hist->GetYaxis()->SetRangeUser(0.5, 2E6);

    hist->Draw(i == 0 ? "E1" : "E1 same");
    leg->AddEntry(hist, Form("R = %s", jetR.Data()), "lep");
  }

  // Draw legends and labels
  c->cd();
  leg->Draw();
  AddSTARLabels(c, Form("Anti-k_{T}, %s trigger", trigger.Data()));

  c->cd();
  c->SaveAs("data" + trigger + ".pdf");
}

// Main plotting function
void plot() {

  std::vector<TString> triggers = {"HT2", "JP2"};
  std::vector<TString> jetRs = {"0.2", "0.3", "0.4", "0.5", "0.6"};

  for (const auto &jetR : jetRs) {
    plotR(jetR);
  }

  for (const auto &trigger : triggers) {
    plotTrigger(trigger);
  }
}