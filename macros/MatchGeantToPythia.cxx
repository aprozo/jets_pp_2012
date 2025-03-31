//! Code to read in geant + pythia output trees and match them
//! Youqi Song & Raghav Kunnawalkam Elayavalli & Kolja Kauder
//! contact - youqi.song@yale.edu
//! HAS to be compiled,
//! root -l MatchGeantToPythia.cxx+("pythia.root", "geant.root")

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

struct JetWithInfo {
  TStarJetVectorJet orig;
  double pt;
  double eta;
  double phi;
  double y;
  double neutral_fraction;
  bool trigger_match_JP2;
  double n_constituents;
  int event_id;
  double weight;
  double multiplicity;
  bool trigger_match_HT2;
  double vz;

  JetWithInfo(TStarJetVectorJet _orig, double _pt, double _eta, double _phi,
              double _y, double _neutral_fraction, bool _trigger_match_JP2,
              double n_constituents, int _event_id, double _weight,
              double _multiplicity, bool _trigger_match_HT2, double _vz)
      : orig(_orig), pt(_pt), eta(_eta), phi(_phi), y(_y),
        neutral_fraction(_neutral_fraction), trigger_match_JP2(_trigger_match_JP2),
        n_constituents(n_constituents), event_id(_event_id), weight(_weight),
        multiplicity(_multiplicity), trigger_match_HT2(_trigger_match_HT2), vz(_vz) {};
};

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
  bool trigger_match_HT2[1000];
  double neutral_fraction[1000];
  bool trigger_match_JP2[1000];
  double pt[1000];
  int n_constituents[1000];
  int index[1000];
};

struct OutputTreeEntry {
  OutputTreeEntry() {};
  OutputTreeEntry(JetWithInfo jet)
      : vz(jet.vz), n_constituents(jet.n_constituents), weight(jet.weight),
        multiplicity(jet.multiplicity), trigger_match_HT2(jet.trigger_match_HT2),
        pt(jet.pt), eta(jet.eta), phi(jet.phi),
        neutral_fraction(jet.neutral_fraction),
        trigger_match_JP2(jet.trigger_match_JP2) {};

  double vz;
  int n_constituents;
  double weight;
  double multiplicity;
  bool trigger_match_HT2;
  double pt;
  double eta;
  double phi;
  double neutral_fraction;
  bool trigger_match_JP2;
};

typedef pair<JetWithInfo, JetWithInfo> MatchedJetWithInfo;

int MatchGeantToPythia(TString mcTreeName, TString OutFile = "test.root") {
  float RCut = 0.6;
  float EtaCut = 1.0 - RCut;

  int MatchNumber = 0;
  int FakeNumber = 0;
  int MissNumber = 0;

  TString mcBaseName =
      mcTreeName(mcTreeName.Last('/') + 1, mcTreeName.Length());
  // if outfile conatins "test" word
  TString mcFolder = "/gpfs01/star/pwg/prozorov/jets_pp_2012/output/mc/";
  TString recoFolder = "/gpfs01/star/pwg/prozorov/jets_pp_2012/output/geant/";

  TString RecoFile = recoFolder + mcBaseName;
  TString McFile = mcFolder + mcBaseName;

  // if (OutFile.Contains("test")) {
  //   RecoFile="/gpfs01/star/pwg/prozorov/jets_pp_2012/tree_pt-hat2535_35_geant.root";
  //   McFile="/gpfs01/star/pwg/prozorov/jets_pp_2012/tree_pt-hat2535_35_mc.root";   
  // } 



  // =================================================================================================
  TFile *Mcf = new TFile(McFile, "READ");
  TH1D *hEventsRun = (TH1D *)Mcf->Get("hEventsRun");

  TTree *McChain = (TTree *)Mcf->Get("ResultTree");
  McChain->BuildIndex("runid", "eventid");

  InputTreeEntry mc;
  McChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
  McChain->SetBranchAddress("Jets", &mc.jets);
  McChain->SetBranchAddress("eventid", &mc.eventid);
  McChain->SetBranchAddress("runid", &mc.runid);
  McChain->SetBranchAddress("weight", &mc.weight);
  McChain->SetBranchAddress("njets", &mc.njets);
  McChain->SetBranchAddress("vz", &mc.vz);
  McChain->SetBranchAddress("mult", &mc.mult);
  McChain->SetBranchAddress("event_sum_pt", &mc.event_sum_pt);
  McChain->SetBranchAddress("trigger_match_HT2", mc.trigger_match_HT2);
  McChain->SetBranchAddress("neutral_fraction", mc.neutral_fraction);
  McChain->SetBranchAddress("trigger_match_JP2", mc.trigger_match_JP2);
  McChain->SetBranchAddress("pt", mc.pt);
  McChain->SetBranchAddress("n_constituents", mc.n_constituents);

  TFile *Recof = new TFile(RecoFile, "READ");
  TTree *RecoChain = (TTree *)Recof->Get("ResultTree");
  RecoChain->BuildIndex("runid", "eventid");
  RecoChain->GetBranch("Jets")->SetAutoDelete(kFALSE);

  InputTreeEntry reco;
  RecoChain->SetBranchAddress("Jets", &reco.jets);
  RecoChain->SetBranchAddress("eventid", &reco.eventid);
  RecoChain->SetBranchAddress("runid", &reco.runid);
  RecoChain->SetBranchAddress("weight", &reco.weight);
  RecoChain->SetBranchAddress("njets", &reco.njets);
  RecoChain->SetBranchAddress("vz", &reco.vz);
  RecoChain->SetBranchAddress("mult", &reco.mult);
  RecoChain->SetBranchAddress("event_sum_pt", &reco.event_sum_pt);
  RecoChain->SetBranchAddress("trigger_match_HT2", reco.trigger_match_HT2);
  RecoChain->SetBranchAddress("neutral_fraction", reco.neutral_fraction);
  RecoChain->SetBranchAddress("trigger_match_JP2", reco.trigger_match_JP2);
  RecoChain->SetBranchAddress("pt", reco.pt);
  RecoChain->SetBranchAddress("n_constituents", reco.n_constituents);

  //! Output and histograms
  // =================================================================================================
  //   histograms
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  TFile *fout = new TFile(OutFile, "RECREATE");

  TH1D *hDeltaR = new TH1D("hDeltaR", "#Delta R all; #Delta R", 350, 0, 3.5);
  TH1D *hPtMc = new TH1D("hPtMc", "Mc p_{t}; p_{t}, GeV/c", 500, 0, 50);
  TH1D *hPtReco = new TH1D("hPtReco", "Reco p_{t}; p_{t}, GeV/c", 500, 0, 50);
  TH2D *hPtMcReco = new TH2D(
      "hPtMcReco", "Mc p_{t} vs Reco p_{t}; Mc p_{t}, GeV/c; Reco p_{t}, GeV/c",
      500, 0, 50, 500, 0, 50);
  TH2D *hNJets = new TH2D("hNJets", "; Number of Mc jets; Number of Reco jets",
                          20, 0, 20, 20, 0, 20);
  TH1D *hMissedJets = new TH1D("hMissedJets", "Missed Jets; Events", 10, 0, 10);
  TH1D *hFakeJets = new TH1D("hFakeJets", "Fake Jets; Events", 20, 0, 20);
  TH1D *hMatchedJets =
      new TH1D("hMatchedJets", "Matched Jets; Events", 20, 0, 20);

  TH1D *hMissRate =
      new TH1D("hMissRate", "Miss Rate; p_{t}, GeV/c", 500, 0, 50);
  TH1D *hFakeRate =
      new TH1D("hFakeRate", "Fake Rate; p_{t}, GeV/c", 500, 0, 50);

  TH1D *stats = new TH1D("stats", "stats", 3, 0, 3);
  stats->GetXaxis()->SetBinLabel(1, "Match");
  stats->GetXaxis()->SetBinLabel(2, "Miss");
  stats->GetXaxis()->SetBinLabel(3, "Fake");

  TTree *MatchedTree = new TTree("MatchedTree", "Matched Jets");

  OutputTreeEntry reco_result;
  OutputTreeEntry mc_result;
  double deltaR = -9;
  MatchedTree->Branch("mc_pt", &mc_result.pt, "mc_pt/D");
  MatchedTree->Branch("mc_eta", &mc_result.eta, "mc_eta/D");
  MatchedTree->Branch("mc_phi", &mc_result.phi, "mc_phi/D");
  MatchedTree->Branch("mc_neutral_fraction", &mc_result.neutral_fraction,
                      "mc_neutral_fraction/D");
  MatchedTree->Branch("mc_trigger_match_JP2", &mc_result.trigger_match_JP2,
                      "mc_trigger_match_JP2/O");

  MatchedTree->Branch("mc_n_constituents", &mc_result.n_constituents,
                      "mc_n_constituents/I");
  MatchedTree->Branch("mc_event_id", &mc_result.vz, "mc_event_id/D");
  MatchedTree->Branch("mc_weight", &mc_result.weight, "mc_weight/D");
  MatchedTree->Branch("mc_multiplicity", &mc_result.multiplicity,
                      "mc_multiplicity/D");
  MatchedTree->Branch("mc_trigger_match_HT2", &mc_result.trigger_match_HT2,
                      "mc_trigger_match_HT2/O");

  MatchedTree->Branch("reco_pt", &reco_result.pt, "reco_pt/D");
  MatchedTree->Branch("reco_eta", &reco_result.eta, "reco_eta/D");
  MatchedTree->Branch("reco_phi", &reco_result.phi, "reco_phi/D");
  MatchedTree->Branch("reco_neutral_fraction", &reco_result.neutral_fraction,
                      "reco_neutral_fraction/D");
  MatchedTree->Branch("reco_trigger_match_JP2", &reco_result.trigger_match_JP2,
                      "reco_trigger_match_JP2/O");
  MatchedTree->Branch("reco_n_constituents", &reco_result.n_constituents,
                      "reco_n_constituents/I");
  MatchedTree->Branch("reco_event_id", &reco_result.vz, "reco_event_id/D");
  MatchedTree->Branch("reco_weight", &reco_result.weight, "reco_weight/D");
  MatchedTree->Branch("reco_multiplicity", &reco_result.multiplicity,
                      "reco_multiplicity/D");
  MatchedTree->Branch("reco_trigger_match_HT2", &reco_result.trigger_match_HT2,
                      "reco_trigger_match_HT2/O");

  MatchedTree->Branch("deltaR", &deltaR, "deltaR/D");

  int N = McChain->GetEntries();

  set<float> mc_accepted_list;
  set<float> reco_accepted_list;

  for (Long64_t mcEvi = 0; mcEvi < N; ++mcEvi) // event loop
  {
    if (!(mcEvi % 50000))
      cout << "Working on " << mcEvi << " / " << N << endl;
    McChain->GetEntry(mcEvi);
    if (abs(mc.vz) > 30) // skip bad events with high weights
      continue;
    if (mc_accepted_list.find(mc.event_sum_pt) != mc_accepted_list.end())
      continue; // some events have identical total particle pT
    else
      mc_accepted_list.insert(mc.event_sum_pt);

    vector<JetWithInfo> mcresult;
    for (int j = 0; j < mc.njets; ++j) {
      TStarJetVectorJet *mc_jet =
          dynamic_cast<TStarJetVectorJet *>(mc.jets->At(j));
      hPtMc->Fill(mc_jet->Pt());

      mcresult.push_back(
          JetWithInfo(*mc_jet, mc_jet->Pt(), mc_jet->Eta(), mc_jet->Phi(),
                      mc_jet->Rapidity(), mc.neutral_fraction[j],
                      mc.trigger_match_JP2[j], mc.n_constituents[j], mc.eventid,
                      mc.weight, mc.mult, mc.trigger_match_HT2[j], mc.vz));
    } // end of mcjet loop

    int recoEvent = RecoChain->GetEntryNumberWithIndex(mc.runid, mc.eventid);

    vector<JetWithInfo> recoresult;

    if (recoEvent >= 0) {
      RecoChain->GetEntry(recoEvent);
      if (abs(reco.vz) > 30)
        continue; // go to next event

      for (int j = 0; j < reco.njets; ++j) {
        TStarJetVectorJet *reco_jet = (TStarJetVectorJet *)reco.jets->At(j);
        if (abs(reco_jet->Eta()) > EtaCut)
          continue;
        hPtReco->Fill(reco_jet->Pt());
        recoresult.push_back(JetWithInfo(
            *reco_jet, reco_jet->Pt(), reco_jet->Eta(), reco_jet->Phi(),
            reco_jet->Rapidity(), reco.neutral_fraction[j],
            reco.trigger_match_JP2[j], reco.n_constituents[j], reco.eventid,
            reco.weight, reco.mult, reco.trigger_match_HT2[j], reco.vz));
      }
    } // end of recojet loop

    //======================================================================================================================================
    // Match them
    //======================================================================================================================================
    vector<MatchedJetWithInfo> MatchedResult;
    TStarJetVectorJet *dummyjet = new TStarJetVectorJet;

    if (mcresult.size() > 0) {
      for (vector<JetWithInfo>::iterator mcit = mcresult.begin();
           mcit != mcresult.end();) {
        bool isMatched = false;
        if (recoresult.size() > 0) {
          for (vector<JetWithInfo>::iterator recoit = recoresult.begin();
               recoit != recoresult.end();) {
            Double_t deta = mcit->y - recoit->y;
            Double_t dphi = TVector2::Phi_mpi_pi(mcit->phi - recoit->phi);
            Double_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
            hDeltaR->Fill(dr);

            if (dr < RCut) {
              MatchedResult.push_back(MatchedJetWithInfo(*mcit, *recoit));
              hPtMcReco->Fill(mcit->pt, recoit->pt);
              isMatched = true;
              recoit = recoresult.erase(recoit);
              ++MatchNumber;
              break;
            } else
              ++recoit;
          }
        }

        if (isMatched)
          mcit = mcresult.erase(mcit);
        else {
          JetWithInfo miss(*dummyjet, -9, -9, -9, -9, -9, 0, -9, mc.eventid,
                           mc.weight, mc.mult, 0, mc.vz);
          MatchedResult.push_back(MatchedJetWithInfo(*mcit, miss));
          // print info about missed jets
          hMissRate->Fill(mcit->pt);
          ++mcit;
          MissNumber++;
        }
      }
    }

    for (auto left_reco_jet : recoresult) {
      JetWithInfo fake(*dummyjet, -9, -9, -9, -9, -9, 0, -9, reco.eventid,
                       reco.weight, reco.mult, -9, reco.vz);
      hFakeRate->Fill(left_reco_jet.pt);
      MatchedResult.push_back(MatchedJetWithInfo(fake, left_reco_jet));
      ++FakeNumber;
    }

    for (auto matched_jet : MatchedResult) {
      mc_result = OutputTreeEntry(matched_jet.first);
      reco_result = OutputTreeEntry(matched_jet.second);
      deltaR = -9;
      if (reco_result.pt > 0 && mc_result.pt > 0) { // only calculate deltaR for matched jets with pT > 0
        Double_t deta = matched_jet.first.y - matched_jet.second.y;
        Double_t dphi = TVector2::Phi_mpi_pi(matched_jet.first.phi -
                                             matched_jet.second.phi);
        deltaR = TMath::Sqrt(deta * deta + dphi * dphi);
      }
      MatchedTree->Fill();
    }
  }

  cout << "Miss Number is " << MissNumber << endl;
  cout << "Fake Number is " << FakeNumber << endl;
  cout << "MatchNumber is " << MatchNumber << endl;

  stats->SetBinContent(1, MatchNumber);
  stats->SetBinContent(2, MissNumber);
  stats->SetBinContent(3, FakeNumber);

  fout->cd();
  hEventsRun->Write();
  hDeltaR->Write();
  hPtMc->Write();
  hPtReco->Write();
  hPtMcReco->Write();
  hNJets->Write();
  hMissedJets->Write();
  hFakeJets->Write();
  hMatchedJets->Write();
  hMissRate->Write();
  hFakeRate->Write();
  MatchedTree->Write();

  fout->Close();

  return 0;
}
