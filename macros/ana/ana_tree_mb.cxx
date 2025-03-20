
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
using namespace std;
//! ----------------------------------------------------
struct InputTreeEntry {
  int runid;
  int runid1;
  int eventid;
  double weight;
  double refmult;
  int njets;
  double vz;
  int mult;
  float event_sum_pt;
  bool is_rejected;
  double neutral_fraction[1000];
  bool trigger_match[1000];
  double pt[1000];
  int n_constituents[1000];
  int index[1000];
};

void normalizeByRun(TH2D *pt_vs_run, TH1D *run) {
  for (int i = 1; i <= pt_vs_run->GetNbinsY(); ++i) {
    double runContent = run->GetBinContent(i);
    if (runContent <= 0)
      continue; // avoid division by zero
    for (int j = 1; j <= pt_vs_run->GetNbinsX(); ++j) {
      double binContent = pt_vs_run->GetBinContent(j, i);
      if (binContent > 0) {
        pt_vs_run->SetBinContent(j, i, binContent / runContent);
        pt_vs_run->SetBinError(j, i, pt_vs_run->GetBinError(j, i) / runContent);
      }
    }
  }
}

int ana_trees_simple(TString OutFile = "hists_mb.root") {

  // sumw2
  TH1::SetDefaultSumw2();
  float RCut = 0.6;
  float EtaCut = 1.0 - RCut;
  TFile *eventsFile = new TFile("events.root", "READ");
  TH1D *nEvents_MB = (TH1D *)eventsFile->Get("hEventsRun_MB");
  TFile *outputFile = new TFile(OutFile, "RECREATE");

  TH1D *MB_pt =
      new TH1D("MB_pt", "MB p_{T}; p_{T} (GeV/c); Entries", 100, 0, 100);

  int nRunBins = nEvents_MB->GetNbinsX();

  TH2D *MB_pt_vs_run =
      new TH2D("MB_pt_vs_run", "; p_{T} (GeV/c); Entries", MB_pt->GetNbinsX(),
               0, MB_pt->GetXaxis()->GetXmax(), nRunBins, 0, nRunBins);

  vector<TString> runBins;
  // set names of bins to run numbers
  for (unsigned int i = 1; i <= nRunBins; ++i) {
    TString runBin = nEvents_MB->GetXaxis()->GetBinLabel(i);
    runBins.push_back(runBin);
    MB_pt_vs_run->GetYaxis()->SetBinLabel(i, runBin);
  }

  TFile *mb_jets = new TFile("../output/MB/tree_jets_mb.root", "READ");
  if (!mb_jets || mb_jets->IsZombie()) {
    cout << "Error: MB tree file not found!" << endl;
    return -1;
  }
  TTree *mb_tree = (TTree *)mb_jets->Get("ResultTree");

  InputTreeEntry mb;
  mb_tree->SetBranchAddress("runid1", &mb.runid1);
  mb_tree->SetBranchAddress("vz", &mb.vz);
  mb_tree->SetBranchAddress("neutral_fraction", mb.neutral_fraction);
  mb_tree->SetBranchAddress("pt", mb.pt);
  mb_tree->SetBranchAddress("njets", &mb.njets);

  for (int i = 0; i < mb_tree->GetEntries(); ++i) {
    mb_tree->GetEntry(i);
    // show progress
    if (i % 100000 == 0)
      cout << "Processing MB entry: " << i << " / " << mb_tree->GetEntries()
           << endl;

    if (mb.vz < -30 || mb.vz > 30)
      continue; // vz cut

    int run = mb.runid1;
    if (run == 13049007 || run == 13048092 || run == 13049006 ||
        run == 13051099 || run == 13051074 || run == 13064067 ||
        run == 13070061 || run == 13069023 || run == 13048093 ||
        run == 13068060 || run == 13052061 || run == 13048019 ||
        run == 13061035) // Livetime JP2 very poor!!!!
      continue;

    int runbin = nEvents_MB->GetXaxis()->FindBin(Form("%i", mb.runid1));
    // fill the histogram with the event weight
    for (int j = 0; j < mb.njets; ++j) {
      MB_pt->Fill(mb.pt[j]);
      MB_pt_vs_run->Fill(mb.pt[j], runbin);
    }
  }

  normalizeByRun(MB_pt_vs_run, nEvents_MB);
  MB_pt->Divide(nEvents_MB);

  // Write histograms to the output file
  outputFile->cd();
  MB_pt->Write();
  MB_pt_vs_run->Write();
}