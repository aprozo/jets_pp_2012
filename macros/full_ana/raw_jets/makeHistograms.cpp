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

// Main analysis function
void process(TString trigger = "HT2", TString jetR = "0.6") {
  // Enable implicit multi-threading for RDataFrame
  ROOT::EnableImplicitMT();
  // Configuration

  TString fileName = Form(
      "/home/prozorov/dev/star/jets_pp_2012/output/merged_data_%s_R%s.root",
      trigger.Data(), jetR.Data());
  TString treeName = "ResultTree";
  //   TString outputDir = Form("./data_%s_R%s", trigger.Data(), jetR.Data());

  // Create output directory
  //   gSystem->mkdir(outputDir.Data(), kTRUE);

  // Create RDataFrame from ROOT file
  std::cout << "Reading data from: " << fileName << std::endl;
  ROOT::RDataFrame df(treeName, fileName);

  // Report initial count
  auto initialCount = df.Count();
  std::cout << "Total entries: " << *initialCount << std::endl;

  // Filter matched jets (both mc_pt and reco_pt != -9)
  auto rawJets = df;

  std::cout << "Processing trigger: " << trigger << std::endl;

  // Filter by trigger match

  auto triggerJets =
      df.Define("matched_mask", Form("trigger_match_%s", trigger.Data()))
          .Define("matched_pt", "pt[matched_mask]")
          .Filter("matched_pt.size() > 0");

  //   triggerJets
  //       .Display({"njets", Form("trigger_match_%s", trigger.Data()), "pt",
  //                 "matched_pt"},
  //                20)
  //       ->Print();

  // print

  // Create comparison plots
  auto trigCount = triggerJets.Count();
  std::cout << "Matched " << trigger << " jets: " << *trigCount << std::endl;

  auto histAll =
      df.Histo1D({"rawJets", ";jet p_{t} [GeV/c]; Counts", 100, 0, 100}, "pt");

  auto histTriggered = triggerJets.Histo1D(
      {"rawJetsTriggered", ";jet p_{t} [GeV/c]; Counts", 100, 0, 100},
      "matched_pt");

  TH1D *eff = (TH1D *)histTriggered.GetPtr()->Clone("triggerEff");
  eff->Divide(histTriggered.GetPtr(), histAll.GetPtr(), 1, 1, "B");
  TString output_file = Form("./data_%s_R%s.root", trigger.Data(), jetR.Data());
  // save into root file
  TFile f(output_file, "RECREATE");
  histAll->Write();
  histTriggered->Write();
  eff->Write();
  f.Close();
}

void makeHistograms() {

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