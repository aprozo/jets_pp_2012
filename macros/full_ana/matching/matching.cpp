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

// Helper function to create variable bin array
std::vector<double> GetVariableBins() {
  std::vector<double> bins;
  // 30 bins from 5 to 65 (steps of 2)
  for (int i = 0; i < 30; i++) {
    bins.push_back(5.0 + i * 2.0);
  }
  bins.push_back(70.0);
  bins.push_back(80.0);
  bins.push_back(90.0);
  return bins;
}

// // Plot 1D comparison of MC and Reco pt
// void PlotComparison(ROOT::RDF::RNode df, const TString &outputDir,
//                     const TString &trigger, const TString &jetR) {
//   gStyle->SetOptStat(0);

//   TCanvas *c = new TCanvas("c_comparison", "Comparison", 1800, 1200);

//   auto h_mc = df.Histo1D({"h_mc",
//                           "Comparison of p_{t} distributions for matched "
//                           "jets;p_{t} [GeV/c];Weighted counts",
//                           400, 0, 100},
//                          "mc_pt", "mc_weight");
//   auto h_reco = df.Histo1D({"h_reco", "", 400, 0, 100}, "reco_pt",
//   "mc_weight");

//   h_mc->SetLineColor(kBlue);
//   h_reco->SetLineColor(kRed);

//   c->SetLogy();
//   h_mc->DrawClone("HIST");
//   h_reco->DrawClone("HIST SAME");

//   TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
//   leg->AddEntry(h_mc.GetPtr(), "mc", "f");
//   leg->AddEntry(h_reco.GetPtr(), "reco", "f");
//   leg->Draw();

//   TLatex latex;
//   latex.SetNDC();
//   latex.SetTextSize(0.04);
//   latex.DrawLatex(0.15, 0.85,
//                   Form("Trigger: %s, R=%s", trigger.Data(), jetR.Data()));

//   c->SaveAs((outputDir + "/matched_jets_pt.pdf").Data());

//   delete leg;
//   delete c;
// }

// // Plot 2D distributions
// void Plot2DDistributions(ROOT::RDF::RNode df, const TString &outputDir) {
//   gStyle->SetOptStat(0);
//   gStyle->SetPalette(kBlueGreenYellow);

//   // pt distribution
//   TCanvas *c_pt = new TCanvas("c_2d_pt", "2D pt", 1000, 600);
//   auto h2_pt = df.Histo2D(
//       {"h2_pt",
//        "2D Histogram of pt distributions for matched jets;mc pt;reco pt",
//        500, 0, 100, 500, 0, 100},
//       "mc_pt", "reco_pt", "mc_weight");

//   c_pt->SetLogz();
//   h2_pt->DrawClone("COLZ");
//   c_pt->SaveAs((outputDir + "/2d_pt.pdf").Data());
//   delete c_pt;

//   // eta distribution
//   TCanvas *c_eta = new TCanvas("c_2d_eta", "2D eta", 1000, 600);
//   auto h2_eta = df.Histo2D(
//       {"h2_eta",
//        "2D Histogram of eta distributions for matched jets;mc eta;reco eta",
//        500, -1, 1, 500, -1, 1},
//       "mc_eta", "reco_eta", "mc_weight");

//   c_eta->SetLogz();
//   h2_eta->DrawClone("COLZ");
//   c_eta->SaveAs((outputDir + "/2d_eta.pdf").Data());
//   delete c_eta;

//   // phi distribution
//   TCanvas *c_phi = new TCanvas("c_2d_phi", "2D phi", 1000, 600);
//   auto h2_phi = df.Histo2D(
//       {"h2_phi",
//        "2D Histogram of phi distributions for matched jets;mc phi;reco phi",
//        500, -3.14, 3.14, 500, -3.14, 3.14},
//       "mc_phi", "reco_phi", "mc_weight");

//   c_phi->SetLogz();
//   h2_phi->DrawClone("COLZ");
//   c_phi->SaveAs((outputDir + "/2d_phi.pdf").Data());
//   delete c_phi;
// }

// Plot ratio with variable binning
void PlotRatio(ROOT::RDF::RNode dfNum, ROOT::RDF::RNode dfDenom,
               const TString &histName, const TString &column,
               const TString &name, const TString &outputDir,
               const TString &trigger, const TString &jetR) {

  gStyle->SetOptStat(0);

  // Get variable bins
  std::vector<double> bins = GetVariableBins();

  // Column to use
  TString ptColumn = column + "_pt";

  // Create histograms using RDataFrame
  auto h_num = dfNum.Histo1D(
      {ptColumn + "h_num",
       Form("%s rate;p_{t} [GeV/c];#frac{1}{dp_{t}d#sigma_{pythia}} Counts",
            histName.Data()),
       static_cast<int>(bins.size() - 1), bins.data()},
      ptColumn.Data(), "mc_weight");

  auto h_denom = dfDenom.Histo1D(
      {ptColumn + "h_denom",
       Form("%s rate;p_{t} [GeV/c];#frac{1}{dp_{t}d#sigma_{pythia}} Counts",
            histName.Data()),
       static_cast<int>(bins.size() - 1), bins.data()},
      ptColumn.Data(), "mc_weight");

  // Clone histograms for manipulation
  TH1D *h_num_clone = (TH1D *)h_num->Clone("h_num_clone");
  TH1D *h_denom_clone = (TH1D *)h_denom->Clone("h_denom_clone");

  // Normalize by bin width
  for (int i = 1; i <= h_num_clone->GetNbinsX(); i++) {
    Double_t binWidth = h_num_clone->GetBinWidth(i);
    h_num_clone->SetBinContent(i, h_num_clone->GetBinContent(i) / binWidth);
    h_num_clone->SetBinError(i, h_num_clone->GetBinError(i) / binWidth);
    h_denom_clone->SetBinContent(i, h_denom_clone->GetBinContent(i) / binWidth);
    h_denom_clone->SetBinError(i, h_denom_clone->GetBinError(i) / binWidth);
  }

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
  canvas->SaveAs(Form("%s/%s_over_%s.pdf", outputDir.Data(), histName.Data(),
                      name.Data()));

  TFile *outFile = new TFile(Form("%s/%s_over_%s.root", outputDir.Data(),
                                  histName.Data(), name.Data()),
                             "RECREATE");
  h_num_clone->Write((TString) "h_num");
  h_denom_clone->Write((TString) "h_denom");
  h_ratio->Write((TString) "h_ratio");
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

  TString fileName = Form(
      "/home/prozorov/dev/star/jets_pp_2012/output/merged_matching_%s_R%s.root",
      trigger.Data(), jetR.Data());
  TString treeName = "MatchedTree";
  TString outputDir =
      Form("./plots/embedding_root_%s_R%s", trigger.Data(), jetR.Data());

  // Create output directory
  gSystem->mkdir(outputDir.Data(), kTRUE);

  // Create RDataFrame from ROOT file
  std::cout << "Reading data from: " << fileName << std::endl;
  ROOT::RDataFrame df(treeName, fileName);

  // Report initial count
  auto initialCount = df.Count();
  std::cout << "Total entries: " << *initialCount << std::endl;

  // Filter matched jets (both mc_pt and reco_pt != -9)
  auto matchedJets = df.Filter("mc_pt != -9 && reco_pt != -9", "Matched jets");

  auto matchedCount = matchedJets.Count();
  std::cout << "Matched jets: " << *matchedCount << std::endl;

  // Create comparison plots
  // std::cout << "Creating comparison plots..." << std::endl;
  // PlotComparison(matchedJets, outputDir, trigger, jetR);
  // Plot2DDistributions(matchedJets, outputDir);

  std::cout << "Processing trigger: " << trigger << std::endl;

  // Filter by trigger match

  auto triggerMatched =
      matchedJets.Filter(Form("reco_trigger_match_%s != 0", trigger.Data()),
                         Form("%s matched", trigger.Data()));

  auto trigCount = triggerMatched.Count();
  std::cout << "  " << trigger << " matched jets: " << *trigCount << std::endl;

  // Plot ratios
  PlotRatio(matchedJets, df, "matched", "mc", "mc", outputDir, trigger, jetR);
  PlotRatio(triggerMatched, matchedJets, trigger, "reco", "reconstructed",
            outputDir, trigger, jetR);

  std::cout << "Analysis complete!" << std::endl;
}

void matching() {

  //   std::vector<TString> triggers = {"HT2"};
  //   std::vector<TString> jetRs = {"0.2"};

  std::vector<TString> triggers = {"HT2", "JP2"};
  std::vector<TString> jetRs = {"0.2", "0.3", "0.4", "0.5", "0.6"};

  for (const auto &trigger : triggers) {
    for (const auto &jetR : jetRs) {
      process(trigger, jetR);
    }
  }
}

// // Comparison overlay function (like the last part of your Python script)
// void CompareMultipleTriggers(const TString &trigger,
//                              const TString &jetR) {
//   gStyle->SetOptStat(0);

//   TString outputDir =
//       Form("./plots/embedding_root_%s_R%s", trigger.Data(), jetR.Data());
//   std::vector<TString> additions = {""};
//   std::vector<TString> names = {"p_{t,lead}>0"};
//   std::vector<int> colors = {kBlue,       kRed - 2,  kGreen + 2,  kMagenta +
//   2,
//                              kOrange + 2, kCyan + 2, kViolet + 2, kGray + 2};

//   TCanvas *canvas = new TCanvas("canvas2", "", 800, 600);

//   TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.5, 1.0, 1.0);
//   pad1->SetBottomMargin(0);
//   pad1->SetTopMargin(0.1);
//   pad1->SetLogy();
//   pad1->Draw();

//   TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.5);
//   pad2->SetTopMargin(0);
//   pad2->SetBottomMargin(0.3);
//   pad2->SetGridy();
//   pad2->Draw();

//   TLegend *legend1 = new TLegend(0.5, 0.5, 0.9, 0.9);

//   int counter = 0;
//   bool firstDraw = true;

//   for (size_t i = 0; i < additions.size(); i++) {
//     TString fileName =
//         Form("%s/%s_over_reconstructed%s_reco.root", outputDir.Data(),
//              trigger.Data(), additions[i].Data());

//     TFile *rootFile = TFile::Open(fileName.Data());
//     if (!rootFile || rootFile->IsZombie()) {
//       std::cerr << "Cannot open file: " << fileName << std::endl;
//       continue;
//     }

//     TH1D *h_num = (TH1D *)rootFile->Get("h_num");
//     TH1D *h_denom = (TH1D *)rootFile->Get("h_denom");
//     TH1D *h_ratio = (TH1D *)rootFile->Get("h_ratio");

//     if (!h_num || !h_denom || !h_ratio) {
//       std::cerr << "Cannot find histograms in file: " << fileName <<
//       std::endl; rootFile->Close(); continue;
//     }

//     // Clone to avoid deletion when file closes
//     h_num = (TH1D *)h_num->Clone(
//         Form("h_num_%s_%s", trigger.Data(), additions[i].Data()));
//     h_denom = (TH1D *)h_denom->Clone(
//         Form("h_denom_%s_%s", trigger.Data(), additions[i].Data()));
//     h_ratio = (TH1D *)h_ratio->Clone(
//         Form("h_ratio_%s_%s", trigger.Data(), additions[i].Data()));
//     h_num->SetDirectory(0);
//     h_denom->SetDirectory(0);
//     h_ratio->SetDirectory(0);

//     rootFile->Close();
//     delete rootFile;

//     // Style histograms
//     h_num->SetLineColor(colors[counter]);
//     h_num->SetMarkerColor(colors[counter]);
//     h_num->SetMarkerStyle(20 + counter);
//     h_num->SetMarkerSize(1.2);

//     h_denom->SetLineColor(colors[counter] + 2);
//     h_denom->SetMarkerColor(colors[counter] + 2);
//     h_denom->SetMarkerStyle(20 + counter);
//     h_denom->SetMarkerSize(1.2);

//     h_ratio->SetLineColor(colors[counter]);
//     h_ratio->SetMarkerColor(colors[counter]);
//     h_ratio->SetMarkerStyle(20 + counter);

//     // Upper pad
//     pad1->cd();
//     h_num->GetYaxis()->SetRangeUser(1e-12, 1e-2);
//     h_num->GetXaxis()->SetRangeUser(0, 80);
//     h_num->GetXaxis()->SetLabelSize(0);
//     h_num->GetYaxis()->SetTitle("#frac{1}{dp_{t}d#sigma_{pythia}} Counts");
//     h_num->GetYaxis()->SetTitleSize(0.06);
//     h_num->GetYaxis()->SetTitleOffset(0.8);

//     if (firstDraw) {
//       h_num->Draw("E1");
//       firstDraw = false;
//     } else {
//       h_num->Draw("E1 SAME");
//     }
//     h_denom->Draw("E1 SAME");

//     legend1->AddEntry(h_num, Form("%s %s", trigger.Data(),
//     names[i].Data()),
//                       "lep");
//     legend1->AddEntry(h_denom, Form("all reco %s", names[i].Data()), "lep");

//     // Lower pad
//     pad2->cd();
//     h_ratio->GetXaxis()->SetTitle("reco p_{t} [GeV/c]");
//     h_ratio->GetXaxis()->SetTitleSize(0.12);
//     h_ratio->GetYaxis()->SetTitle(
//         Form("Ratio %s/%s", trigger.Data(), names[i].Data()));
//     h_ratio->GetYaxis()->SetTitleSize(0.08);
//     h_ratio->GetYaxis()->SetTitleOffset(0.5);
//     h_ratio->GetYaxis()->SetLabelSize(0.06);
//     h_ratio->GetXaxis()->SetRangeUser(0, 80);

//     if (counter == 0) {
//       h_ratio->Draw("E1");
//     } else {
//       h_ratio->Draw("E1 SAME");
//     }

//     counter++;
//   }

//   pad1->cd();
//   legend1->Draw();

//   canvas->SaveAs(
//       Form("%s/comparison_%s.pdf", outputDir.Data(), trigger.Data()));

//   delete legend1;
//   delete pad2;
//   delete pad1;
//   delete canvas;
// }
