
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLine.h>
#include <TProfile.h>

#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <TStarJetPicoReader.h>
#include <TStarJetVector.h>
#include <TStarJetVectorJet.h>

#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <vector>

using namespace std;

//! ----------------------------------------------------
struct InputTreeEntry {
  InputTreeEntry() : jets(new TClonesArray("TStarJetVectorJet", 1000)) {}
  TClonesArray *jets;
  int runid;
  int runid1;
  int eventid;
  double weight;
  double refmult;
  int njets;
  double vz;
  int mult;
  float event_sum_pt;
  double neutral_fraction[1000];
  bool trigger_match[1000];
  double pt[1000];
  int n_constituents[1000];
  int index[1000];
};


int ana_trees(TString inputTreeName= "output/JP2/tree_sum0_part0.root", TString OutFile = "test.root") {
  float RCut = 0.6;
  float EtaCut = 1.0 - RCut;

  // TFile *eventsFile = new TFile("events.root", "READ");
  // TH1D *weights = (TH1D *)eventsFile->Get("hist_JP2_weight");
  // TH1D * nEvents_MB = (TH1D *)eventsFile->Get("hEventsRun_MB");
  // TH1D * nEvents_JP2 = (TH1D *)eventsFile->Get("hEventsRun_JP2");

  // TH1D * JP2_scaled = (TH1D *)nEvents_JP2->Clone("JP2_scaled");
  // JP2_scaled->Multiply(weights); // scale the JP2 events by the weights


  // =================================================================================================
  TFile *inputFile = new TFile(inputTreeName, "READ");
  TTree *MyChain = (TTree *)inputFile->Get("ResultTree");
  MyChain->BuildIndex("runid", "eventid");

  InputTreeEntry data_jet ;
  MyChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
  MyChain->SetBranchAddress("Jets", &data_jet.jets);
  MyChain->SetBranchAddress("eventid", &data_jet.eventid);
  MyChain->SetBranchAddress("runid", &data_jet.runid);
  MyChain->SetBranchAddress("weight", &data_jet.weight);
  MyChain->SetBranchAddress("njets", &data_jet.njets);
  MyChain->SetBranchAddress("vz", &data_jet.vz);
  MyChain->SetBranchAddress("mult", &data_jet.mult);
  MyChain->SetBranchAddress("event_sum_pt", &data_jet.event_sum_pt);
  MyChain->SetBranchAddress("neutral_fraction", data_jet.neutral_fraction);
  MyChain->SetBranchAddress("trigger_match_jp", data_jet.trigger_match);
  MyChain->SetBranchAddress("pt", data_jet.pt);
  MyChain->SetBranchAddress("n_constituents", data_jet.n_constituents);
  //! Output and histograms
  // =================================================================================================
  //   histograms
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  TFile *fout = new TFile(OutFile, "RECREATE");
  TH1D *hPt = new TH1D("hPt", " p_{t}; p_{t}, GeV/c", 500, 0, 50);
  

  set<float> accepted_list;
  for (Long64_t iEvent = 0; iEvent < MyChain->GetEntries();; ++iEvent) // event loop
  {
    if (!(iEvent % 50000))      cout << "Working on " << iEvent << " / " << N << endl;
    MyChain->GetEntry(iEvent);

    if (abs(data_jet.vz) > 30) // skip bad events with high weights
      continue;

    if (accepted_list.find(data_jet.event_sum_pt) != accepted_list.end())
      continue; // some events have identical total particle pT
    else
      accepted_list.insert(data_jet.event_sum_pt);

    for (int j = 0; j < data_jet.njets; ++j) {
      TStarJetVectorJet *jet = dynamic_cast<TStarJetVectorJet *>(data_jet.jets->At(j));
      if (abs(jet->Eta()) > EtaCut)        continue;
      hPt->Fill(jet->Pt());


    } // end of jet loop   

  }
 
  fout->cd();
  hEventsRun->Write();
  hPt->Write();


  fout->Close();

  return 0;
}
