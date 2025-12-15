#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iostream>
#include <string>
#include <vector>

#include "/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/cross_section/config.h"

using namespace CrossSectionConfig;

const std::vector<double> ptLead_reco_bins = {
    1.0,  2.0,  3.0,  4.1,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0,
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
    21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0};

// Plot ratio with variable binning
void PlotRatio(ROOT::RDF::RNode dfNum, ROOT::RDF::RNode dfDenom,
               const TString &column, const TString &histName,
               const TString &name, const TString &outDir,
               const TString &trigger, const TString &jetR) {

  gStyle->SetOptStat(0);
  // Get variable bins
  std::vector<double> bins = column == "mc" ? pt_mc_bins : pt_reco_bins;

  // Column to use
  TString ptColumn = column + "_pt";

  // Create histograms using RDataFrame
  auto h_num =
      dfNum.Histo1D({histName + name + "h_num",
                     Form("%s rate; p_{t} [GeV/c]; Counts", histName.Data()),
                     static_cast<int>(bins.size() - 1), bins.data()},
                    ptColumn.Data(), "mc_weight");

  auto h_denom =
      dfDenom.Histo1D({histName + name + "h_denom",
                       Form("%s rate; p_{t} [GeV/c]; Counts", histName.Data()),
                       static_cast<int>(bins.size() - 1), bins.data()},
                      ptColumn.Data(), "mc_weight");

  auto h_num_2D_vs_lead = dfNum.Histo2D(
      {histName + name + "h_num_2D_vs_lead",
       Form("%s rate vs lead; lead p_{t} [GeV/c]; p_{t} "
            "[GeV/c]; Counts",
            histName.Data()),
       static_cast<int>(ptLead_reco_bins.size() - 1), ptLead_reco_bins.data(),
       static_cast<int>(bins.size() - 1), bins.data()},
      "reco_ptLead", ptColumn.Data(), "mc_weight");

  auto h_denom_2D_vs_lead = dfDenom.Histo2D(
      {histName + name + "h_denom_2D_vs_lead",
       Form("%s rate vs lead; lead p_{t} [GeV/c]; p_{t} "
            "[GeV/c]; Counts",
            histName.Data()),
       static_cast<int>(ptLead_reco_bins.size() - 1), ptLead_reco_bins.data(),
       static_cast<int>(bins.size() - 1), bins.data()},
      "reco_ptLead", ptColumn.Data(), "mc_weight");

  TH2D *h_ratio_2D_vs_lead =
      (TH2D *)h_num_2D_vs_lead->Clone("h_ratio_2D_vs_lead");
  h_ratio_2D_vs_lead->Divide(h_denom_2D_vs_lead.GetPtr());

  TH1D *h_num_ptLead = h_num_2D_vs_lead->ProjectionX("h_num_ptLead");
  TH1D *h_denom_ptLead = h_denom_2D_vs_lead->ProjectionX("h_denom_ptLead");
  TH1D *h_ratio_ptLead = (TH1D *)h_num_ptLead->Clone("h_ratio_ptLead");
  h_ratio_ptLead->Divide(h_denom_ptLead);

  // Clone histograms for manipulation
  TH1D *h_num_clone = (TH1D *)h_num->Clone("h_num_clone");
  TH1D *h_denom_clone = (TH1D *)h_denom->Clone("h_denom_clone");

  // Create ratio histogram
  TH1D *h_ratio = (TH1D *)h_num_clone->Clone("h_ratio");
  h_ratio->Divide(h_num_clone, h_denom_clone, 1, 1, "B");

  // Create canvas with two pads
  TCanvas *canvas = new TCanvas(
      Form("canvas_%s_%s", histName.Data(), column.Data()),
      Form("canvas_%s_%s", histName.Data(), column.Data()), 800, 600);

  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.5, 1.0, 1.0);
  pad1->SetBottomMargin(0);
  pad1->SetTopMargin(0.1);
  pad1->SetLogy();
  pad1->Draw();

  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.5);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad2->SetGridy();
  pad2->Draw();

  // Upper pad - histograms
  pad1->cd();
  h_num_clone->SetLineColor(kBlue);
  h_num_clone->SetMarkerColor(kBlue);
  h_num_clone->SetMarkerStyle(20);
  h_num_clone->SetMarkerSize(1.2);
  h_num_clone->GetXaxis()->SetLabelSize(0);
  h_num_clone->GetYaxis()->SetTitleSize(0.06);
  h_num_clone->GetYaxis()->SetTitleOffset(0.8);
  h_num_clone->SetTitle(Form("%s rate", histName.Data()));
  h_num_clone->GetXaxis()->SetRangeUser(0, 100);
  h_num_clone->Draw("E1");

  h_denom_clone->SetLineColor(kRed);
  h_denom_clone->SetMarkerColor(kRed);
  h_denom_clone->SetMarkerStyle(21);
  h_denom_clone->SetMarkerSize(1.2);
  h_denom_clone->Draw("E1 SAME");

  TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend1->AddEntry(h_num_clone, histName.Data(), "lep");
  legend1->AddEntry(h_denom_clone, name.Data(), "lep");
  legend1->Draw();

  // Lower pad - ratio
  pad2->cd();
  h_ratio->SetLineColor(kBlue);
  h_ratio->SetMarkerColor(kBlue);
  h_ratio->SetMarkerStyle(20);
  h_ratio->SetMarkerSize(1.2);
  h_ratio->SetTitle("");
  h_ratio->GetXaxis()->SetTitle(Form("%s p_{t} [GeV/c]", column.Data()));
  h_ratio->GetXaxis()->SetTitleSize(0.12);
  h_ratio->GetXaxis()->SetLabelSize(0.08);
  h_ratio->GetYaxis()->SetTitle(
      Form("Ratio %s/%s", histName.Data(), name.Data()));
  h_ratio->GetYaxis()->SetTitleSize(0.08);
  h_ratio->GetYaxis()->SetTitleOffset(0.5);
  h_ratio->GetYaxis()->SetLabelSize(0.06);

  Double_t yMax = (histName == "fake")   ? 0.25
                  : (histName == "miss") ? 0.3
                                         : 1.15;
  h_ratio->GetYaxis()->SetRangeUser(0, yMax);
  h_ratio->GetXaxis()->SetRangeUser(0, 100);
  h_ratio->Draw("E1");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.08);
  latex.DrawLatex(0.15, 0.85,
                  Form("Trigger: %s, R=%s", trigger.Data(), jetR.Data()));

  // Save outputs
  canvas->SaveAs(
      Form("%s/%s_over_%s.pdf", outDir.Data(), histName.Data(), name.Data()));

  TFile *outFile = new TFile(
      Form("%s/%s_over_%s.root", outDir.Data(), histName.Data(), name.Data()),
      "RECREATE");
  h_num_clone->Write((TString) "h_num");
  h_denom_clone->Write((TString) "h_denom");
  h_ratio->Write((TString) "h_ratio");

  h_num_ptLead->Write((TString) "h_num_ptLead");
  h_denom_ptLead->Write((TString) "h_denom_ptLead");
  h_ratio_ptLead->Write((TString) "h_ratio_ptLead");

  h_num_2D_vs_lead->Write((TString) "h_num_2D_vs_lead");
  h_denom_2D_vs_lead->Write((TString) "h_denom_2D_vs_lead");
  h_ratio_2D_vs_lead->Write((TString) "h_ratio_2D_vs_lead");
  outFile->Close();

  delete outFile;
  delete legend1;
  delete h_ratio;
  delete h_denom_clone;
  delete h_num_clone;
  delete pad2;
  delete pad1;
  delete canvas;
}

// Main analysis function
void process(TString trigger = "HT2", TString jetR = "0.6") {
  // Enable implicit multi-threading for RDataFrame
  ROOT::EnableImplicitMT();
  // Configuration
  std::cout << "Processing trigger: " << trigger << std::endl;
  TString fileName = Form(
      "/home/prozorov/dev/star/jets_pp_2012/output/merged_matching_%s_R%s.root",
      trigger.Data(), jetR.Data());
  //   TString fileName =
  //   Form("/home/prozorov/dev/star/jets_pp_2012/output_working/"
  //                           "jets_embedding.root");
  TString treeName = "MatchedTree";
  TString outDir =
      Form("./plots/embedding_root_%s_R%s", trigger.Data(), jetR.Data());

  // Create output directory
  gSystem->mkdir(outDir.Data(), kTRUE);

  // Create RDataFrame from ROOT file
  std::cout << "Reading data from: " << fileName << std::endl;
  //   ROOT::RDataFrame df(treeName, fileName);
  //   auto allJets = df.Filter("isTriggerEvent", "select relevant jets");

  ROOT::RDataFrame allJets(treeName, fileName);

  auto mcJets = allJets.Filter("mc_pt != -9", "mc jets");
  auto recoJets = allJets.Filter("reco_pt != -9", "reco jets");
  auto fakeJets = recoJets.Filter("mc_pt == -9", "fake jets");

  auto matchedJets =
      allJets.Filter("mc_pt != -9 && reco_pt != -9", "matched jets");
  // Filter by trigger match

  auto trigJets =
      allJets.Filter(Form("reco_trigger_match_%s != 0", trigger.Data()),
                     Form("%s triggered", trigger.Data()));

  // 1-Miss Rate
  PlotRatio(matchedJets, allJets, "mc", "matched", "mc", outDir, trigger, jetR);

  // Trigger efficiency
  PlotRatio(trigJets, recoJets, "reco", trigger, "all", outDir, trigger, jetR);

  //   // 1- Fake Rate
  //   PlotRatio(mcJets, allJets, "reco", "reco", "all", outDir, trigger, jetR);

  PlotRatio(fakeJets, allJets, "reco", "fakereco", "all", outDir, trigger,
            jetR);
}

void matching() {

  //   std::vector<TString> triggers = {"HT2"};
  //   std::vector<TString> jetRs = {"0.6"};

  std::vector<TString> triggers = {"HT2", "JP2"};
  std::vector<TString> jetRs = {"0.5", "0.6"};
  // std::vector<TString> jetRs = {"0.2", "0.3", "0.4", "0.5", "0.6"};

  for (const auto &trigger : triggers) {
    for (const auto &jetR : jetRs) {
      process(trigger, jetR);
    }
  }
}
