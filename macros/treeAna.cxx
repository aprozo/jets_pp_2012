//! HAS to be compiled,
//! root -l macros/PrepUnfolding.cxx+
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>

#include "TStarJetVectorJet.h"

#include <iostream>
#include <vector>

using namespace std;
struct ResultStruct
{
  TStarJetVectorJet orig;
  double pT;
  int eventid;
  ResultStruct(TStarJetVectorJet orig, double pT, int eventid) : orig(orig),
                                                                 pT(pT),
                                                                 eventid(eventid) {}
};

int treeAna(TString inputFileName = "output/output_jets_MB.root")
{
  Float_t jetRadius = 0.4;
  TH1::SetDefaultSumw2(true);
  float jetEtaCut = 1.0 - jetRadius;
  float jetPtCut = 10;

  // take basename of the input file, remove folder
  TString outputName = inputFileName;
  outputName.ReplaceAll(".root", "");
  outputName += "_hists.root";

  //! Set up events
  TFile *dataFile = new TFile(inputFileName);
  TTree *jetTree = (TTree *)dataFile->Get("ResultTree");
  jetTree->BuildIndex("runid", "eventid");
  cout << "size is" << jetTree->GetEntries() << endl;

  TClonesArray *jets = new TClonesArray("TStarJetVectorJet");
  jetTree->GetBranch("Jets")->SetAutoDelete(kFALSE);
  jetTree->SetBranchAddress("Jets", &jets);

  const int numberOfPtBins = 13;
  // const char *PTBINS[numberOfPtBins]={"2_3","3_4","4_5","5_7","7_9","9_11","11_15","15_20","20_25","25_35","35_-1"};
  const float XSEC[numberOfPtBins] = {0.00230158, 0.000342755, 0.0000457002, 9.0012, .00000972535, 1.46253, 0.000000469889, 0.354566, 0.0000000269202, 0.00000000143453, 0.151622, 0.0249062, 0.00584527};
  const float NUMBEROFEVENT[numberOfPtBins] = {3000000, 3000000, 3000000, 3000000, 2000000, 3000000, 2000000, 3000000, 1000000, 1000000, 3000000, 3000000, 3000000};
  const float MAXPT[numberOfPtBins] = {30 * 0.75, 40 * 0.75, 50 * 0.75, 6 * 0.75, 70 * 0.75, 8 * 0.75, 90 * 0.75, 10 * 0.75, 110 * 0.75, 2000 * 0.75, 14 * 0.75, 18 * 0.75, 22 * 0.75};
  const vector<string> vptbins = {"1115_", "1520_", "2025_", "23_", "2535_", "34_", "3545_", "45_", "4555_", "55999_", "57_", "79_", "911_"};

  // vector<double> *pT = new vector<double>;
  // jetTree->SetBranchAddress("pT_lead0", &pT);

  int eventid;
  int runid;
  double eventWeight;
  int njets;
  jetTree->SetBranchAddress("eventid", &eventid);
  jetTree->SetBranchAddress("runid", &runid);
  jetTree->SetBranchAddress("weight", &eventWeight);
  jetTree->SetBranchAddress("njets", &njets);

  TFile *fout = new TFile(outputName, "RECREATE");

  TH1D *hJetPt = new TH1D("hJetPt", "", 80, 0, 80);
  TH1D *hJetPt_Fine = new TH1D("hJetPt_Fine", "", 800000, 0, 80);
  TH2D *hJetEtaPhi = new TH2D("hJetEtaPhi", "", 1000, -1, 1, 1000, 0, 6.3);

  //! Loop over measured level
  for (Long64_t iEvent = 0; iEvent < jetTree->GetEntries(); ++iEvent)
  {
    if (!(iEvent % 10000))
      cout << "Working on " << iEvent << " / " << jetTree->GetEntries() << endl;
    jetTree->GetEntry(iEvent);

    vector<ResultStruct> ResultStructs;
    for (int iJet = 0; iJet < njets; ++iJet)
    {
      TStarJetVectorJet *jet = (TStarJetVectorJet *)jets->At(iJet);
      // Bool_t isBadEvent = kFALSE;
      // for (int i = 0; i < vptbins.size(); ++i)
      // {
      //   if (eventWeight == XSEC[i] / NUMBEROFEVENT[i] && jet->perp() > MAXPT[i])
      //   {
      //     isBadEvent = kTRUE;
      //     break;
      //   }
      // }
      // if (isBadEvent)
      //   break;

      double pt = jet->Pt();

      if (fabs(jet->Eta()) < jetEtaCut && pt > jetPtCut)
      {
        ResultStructs.push_back(ResultStruct(*jet, pt, eventid));
      }
    }

    //! Record Measured
    for (auto resJet : ResultStructs)
    {
      double jetPt = resJet.orig.Pt();
      hJetPt->Fill(jetPt, eventWeight);
      hJetPt_Fine->Fill(jetPt, eventWeight);
      double jetEta = resJet.orig.eta();
      double jetPhi = resJet.orig.phi();
      hJetEtaPhi->Fill(jetEta, jetPhi, eventWeight);
    }

    ResultStructs.clear();
  } //! event loop

  fout->Write();

  cout << " Wrote to" << endl
       << outputName << endl;

  return 0;
}
