#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#include "TStarJetVectorJet.h"

#include <iostream>
#include <set>
#include <vector>

using namespace std;

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
   double ptLead[1000];
   double pt[1000];
   int n_constituents[1000];
   int index[1000];
};

struct MyJet {
   TStarJetVectorJet orig;
   double pt;
   double ptLead;
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

   float deltaR(const MyJet &other) const
   {
      if (pt < 0 || other.pt < 0) {
         return 10000; // Return 1000 to indicate invalid jets
      }
      float deta = eta - other.eta;
      float dphi = TVector2::Phi_mpi_pi(phi - other.phi);
      return sqrt(deta * deta + dphi * dphi);
   }
   MyJet()
      : pt(-9),
        ptLead(-9),
        eta(-9),
        phi(-9),
        y(-9),
        neutral_fraction(-9),
        trigger_match_JP2(false),
        n_constituents(-9),
        event_id(-9),
        weight(-9),
        multiplicity(-9),
        trigger_match_HT2(false),
        vz(-9) {};

   MyJet(TStarJetVectorJet _orig, double _pt, double _ptLead, double _eta, double _phi, double _y,
         double _neutral_fraction, bool _trigger_match_JP2, double n_constituents, int _event_id, double _weight,
         double _multiplicity, bool _trigger_match_HT2, double _vz)
      : orig(_orig),
        pt(_pt),
        ptLead(_ptLead),
        eta(_eta),
        phi(_phi),
        y(_y),
        neutral_fraction(_neutral_fraction),
        trigger_match_JP2(_trigger_match_JP2),
        n_constituents(n_constituents),
        event_id(_event_id),
        weight(_weight),
        multiplicity(_multiplicity),
        trigger_match_HT2(_trigger_match_HT2),
        vz(_vz) {};
};

typedef pair<MyJet, MyJet> MatchedJetPair;

vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets, const vector<MyJet> &RecoJets, const double &R);

int matching_mc_reco(TString mcTreeName = "/gpfs01/star/pwg/prozorov/jets_pp_2012/output/mc/"
                                          "tree_pt-hat3545_033_R0.5.root")
{

   TString mcBaseName = mcTreeName(mcTreeName.Last('/') + 1, mcTreeName.Length());
   TString OutFile = "matched_" + mcBaseName;

   TString mcFolder = "/gpfs01/star/pwg/prozorov/jets_pp_2012/output/mc/";
   TString recoFolder = "/gpfs01/star/pwg/prozorov/jets_pp_2012/output/geant/";

   TString RecoFile = recoFolder + mcBaseName;
   TString McFile = mcFolder + mcBaseName;

   float RCut = 0;
   // get RCut from input filename if it contains *Rx.x.root
   if (mcTreeName.Contains("R0.2"))
      RCut = 0.2;
   else if (mcTreeName.Contains("R0.3"))
      RCut = 0.3;
   else if (mcTreeName.Contains("R0.4"))
      RCut = 0.4;
   else if (mcTreeName.Contains("R0.5"))
      RCut = 0.5;
   else if (mcTreeName.Contains("R0.6"))
      RCut = 0.6;
   else {
      cout << "Cannot get RCut from input filename. Exiting." << endl;
      return -1;
   }
   const float EtaCut = 1.0 - RCut;
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
   // McChain->SetBranchAddress("vz", &mc.vz);
   McChain->SetBranchAddress("mult", &mc.mult);
   McChain->SetBranchAddress("event_sum_pt", &mc.event_sum_pt);
   McChain->SetBranchAddress("trigger_match_HT2", mc.trigger_match_HT2);
   McChain->SetBranchAddress("neutral_fraction", mc.neutral_fraction);
   McChain->SetBranchAddress("trigger_match_JP2", mc.trigger_match_JP2);
   McChain->SetBranchAddress("pt", mc.pt);
   McChain->SetBranchAddress("ptLead", mc.ptLead);
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
   // RecoChain->SetBranchAddress("vz", &reco.vz);
   RecoChain->SetBranchAddress("mult", &reco.mult);
   RecoChain->SetBranchAddress("event_sum_pt", &reco.event_sum_pt);
   RecoChain->SetBranchAddress("trigger_match_HT2", reco.trigger_match_HT2);
   RecoChain->SetBranchAddress("neutral_fraction", reco.neutral_fraction);
   RecoChain->SetBranchAddress("trigger_match_JP2", reco.trigger_match_JP2);
   RecoChain->SetBranchAddress("pt", reco.pt);
   RecoChain->SetBranchAddress("ptLead", reco.ptLead);
   RecoChain->SetBranchAddress("n_constituents", reco.n_constituents);

   TFile *fout = new TFile(OutFile, "RECREATE");
   TH1D *hDeltaR = new TH1D("hDeltaR", "#Delta R all; #Delta R", 350, 0, 3.5);
   TH1D *hPtMc = new TH1D("hPtMc", "Mc p_{t}; p_{t}, GeV/c", 500, 0, 50);
   TH1D *hPtReco = new TH1D("hPtReco", "Reco p_{t}; p_{t}, GeV/c", 500, 0, 50);
   TH2D *hPtMcReco =
      new TH2D("hPtMcReco", "Mc p_{t} vs Reco p_{t}; Mc p_{t}, GeV/c; Reco p_{t}, GeV/c", 500, 0, 50, 500, 0, 50);

   TH1D *hMiss = new TH1D("hMiss", "Miss Rate; p_{t}, GeV/c", 500, 0, 50);
   TH1D *hFake = new TH1D("hFake", "Fake Rate; p_{t}, GeV/c", 500, 0, 50);

   TH1D *stats = new TH1D("stats", "stats", 3, 0, 3);
   stats->GetXaxis()->SetBinLabel(1, "Match");
   stats->GetXaxis()->SetBinLabel(2, "Miss");
   stats->GetXaxis()->SetBinLabel(3, "Fake");

   TTree *MatchedTree = new TTree("MatchedTree", "Matched Jets");
   MyJet outRecoJet;
   MyJet outMcJet;
   double deltaR = -9;

   MatchedTree->Branch("mc_pt", &outMcJet.pt, "mc_pt/D");
   MatchedTree->Branch("mc_ptLead", &outMcJet.ptLead, "mc_ptLead/D");
   MatchedTree->Branch("mc_eta", &outMcJet.eta, "mc_eta/D");
   MatchedTree->Branch("mc_phi", &outMcJet.phi, "mc_phi/D");
   MatchedTree->Branch("mc_neutral_fraction", &outMcJet.neutral_fraction, "mc_neutral_fraction/D");
   MatchedTree->Branch("mc_trigger_match_JP2", &outMcJet.trigger_match_JP2, "mc_trigger_match_JP2/O");

   MatchedTree->Branch("mc_n_constituents", &outMcJet.n_constituents, "mc_n_constituents/I");
   // MatchedTree->Branch("mc_event_id", &outMcJet.vz, "mc_event_id/D");
   MatchedTree->Branch("mc_weight", &outMcJet.weight, "mc_weight/D");
   MatchedTree->Branch("mc_multiplicity", &outMcJet.multiplicity, "mc_multiplicity/D");
   MatchedTree->Branch("mc_trigger_match_HT2", &outMcJet.trigger_match_HT2, "mc_trigger_match_HT2/O");

   MatchedTree->Branch("reco_pt", &outRecoJet.pt, "reco_pt/D");
   MatchedTree->Branch("reco_ptLead", &outRecoJet.ptLead, "reco_ptLead/D");
   MatchedTree->Branch("reco_eta", &outRecoJet.eta, "reco_eta/D");
   MatchedTree->Branch("reco_phi", &outRecoJet.phi, "reco_phi/D");
   MatchedTree->Branch("reco_neutral_fraction", &outRecoJet.neutral_fraction, "reco_neutral_fraction/D");
   MatchedTree->Branch("reco_trigger_match_JP2", &outRecoJet.trigger_match_JP2, "reco_trigger_match_JP2/O");
   MatchedTree->Branch("reco_n_constituents", &outRecoJet.n_constituents, "reco_n_constituents/I");
   // MatchedTree->Branch("reco_event_id", &outRecoJet.vz, "reco_event_id/D");
   MatchedTree->Branch("reco_weight", &outRecoJet.weight, "reco_weight/D");
   MatchedTree->Branch("reco_multiplicity", &outRecoJet.multiplicity, "reco_multiplicity/D");
   MatchedTree->Branch("reco_trigger_match_HT2", &outRecoJet.trigger_match_HT2, "reco_trigger_match_HT2/O");

   MatchedTree->Branch("deltaR", &deltaR, "deltaR/D");

   int nEvents = McChain->GetEntries();

   set<float> accepted_events_list;
   int MatchNumber = 0;
   int FakeNumber = 0;
   int MissNumber = 0;

   for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) // event loop
   {
      McChain->GetEntry(iEvent);
      // if (abs(mc.vz) > 30) // skip bad events with high weights
      //   continue;
      if (accepted_events_list.count(mc.event_sum_pt) > 0)
         continue; // some events have identical total particle pT
      else
         accepted_events_list.insert(mc.event_sum_pt);

      vector<MyJet> mcJets;
      for (int j = 0; j < mc.njets; ++j) {
         TStarJetVectorJet *tempMcJet = dynamic_cast<TStarJetVectorJet *>(mc.jets->At(j));

         mcJets.push_back(MyJet(*tempMcJet, tempMcJet->Pt(), mc.ptLead[j], tempMcJet->Eta(), tempMcJet->Phi(),
                                tempMcJet->Rapidity(), mc.neutral_fraction[j], mc.trigger_match_JP2[j],
                                mc.n_constituents[j], mc.eventid, mc.weight, mc.mult, mc.trigger_match_HT2[j], mc.vz));
      }

      int recoEvent = RecoChain->GetEntryNumberWithIndex(mc.runid, mc.eventid);

      vector<MyJet> recoJets;
      if (recoEvent >= 0) {
         RecoChain->GetEntry(recoEvent);
         // if (abs(reco.vz) > 30)
         //   continue; // go to next event
         for (int j = 0; j < reco.njets; ++j) {
            TStarJetVectorJet *tempRecoJet = (TStarJetVectorJet *)reco.jets->At(j);
            if (abs(tempRecoJet->Eta()) > EtaCut)
               continue;
            recoJets.push_back(MyJet(*tempRecoJet, tempRecoJet->Pt(), reco.ptLead[j], tempRecoJet->Eta(),
                                     tempRecoJet->Phi(), tempRecoJet->Rapidity(), reco.neutral_fraction[j],
                                     reco.trigger_match_JP2[j], reco.n_constituents[j], reco.eventid, reco.weight,
                                     reco.mult, reco.trigger_match_HT2[j], reco.vz));
         }
      } // end of recojet loop
      //======================================================================================================================================
      // Match them
      //======================================================================================================================================

      vector<MatchedJetPair> MatchedJets;
      MatchedJets = MatchJetsEtaPhi(mcJets, recoJets, RCut);

      for (unsigned int j = 0; j < MatchedJets.size(); j++) {
         outMcJet = MatchedJets[j].first;
         outRecoJet = MatchedJets[j].second;
         deltaR = outMcJet.deltaR(outRecoJet);

         MatchedTree->Fill();
         // fill histograms and counters
         if (outMcJet.pt > 0 && outRecoJet.pt > 0) { // matched jets
            hDeltaR->Fill(deltaR);
            hPtMcReco->Fill(outMcJet.pt, outRecoJet.pt);
            MatchNumber++;
         } else if (outMcJet.pt > 0 && outRecoJet.pt <= 0) { // missed jets
            hMiss->Fill(outMcJet.pt);
            MissNumber++;
         } else if (outMcJet.pt <= 0 && outRecoJet.pt > 0) { // fake jets
            hFake->Fill(outRecoJet.pt);
            FakeNumber++;
         }
         hPtMc->Fill(outMcJet.pt);
         hPtReco->Fill(outRecoJet.pt);

      } // end loop over matched jets
   }

   cout << "Miss Number is " << MissNumber << endl;
   cout << "Fake Number is " << FakeNumber << endl;
   cout << "MatchNumber is " << MatchNumber << endl;

   stats->SetBinContent(1, MatchNumber);
   stats->SetBinContent(2, MissNumber);
   stats->SetBinContent(3, FakeNumber);

   fout->cd();
   MatchedTree->Write();

   stats->Write();
   hEventsRun->Write();
   hDeltaR->Write();
   hPtMc->Write();
   hPtReco->Write();
   hPtMcReco->Write();
   hMiss->Write();
   hFake->Write();

   TH1D *hMissRate = (TH1D *)hMiss->Clone("hMissRate"); // use "b" divide option to handle binomial errors
   hMissRate->Divide(hMiss, hPtMc, 1, 1, "b");
   hMissRate->Write();
   TH1D *hFakeRate = (TH1D *)hFake->Clone("hFakeRate");
   hFakeRate->Divide(hFake, hPtReco, 1, 1, "b");
   hFakeRate->Write();

   fout->Close();

   return 0;
}

vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets, const vector<MyJet> &RecoJets, const double &R)
{
   vector<MyJet> recoJetsCopy = RecoJets; // copy to avoid modifying the original
   vector<MatchedJetPair> matchedJets;
   MyJet dummy;

   for (const auto &mcJet : McJets) {
      bool isMatched = false;
      double minDeltaR = 10000; // Initialize to a large value
      auto bestMatch = recoJetsCopy.end();

      // Find the closest reco jet within the threshold
      for (auto rcit = recoJetsCopy.begin(); rcit != recoJetsCopy.end(); ++rcit) {
         MyJet recoJet = *rcit;
         double deltaR = mcJet.deltaR(recoJet);

         if (deltaR < 0.6 * R && deltaR < minDeltaR) {
            minDeltaR = deltaR;
            bestMatch = rcit;
            isMatched = true;
         }
      }
      // If a match was found, add it and remove from available jets
      if (isMatched) {
         matchedJets.push_back(make_pair(mcJet, *bestMatch));
         recoJetsCopy.erase(bestMatch);
      } else {
         // If no match was found for this MC jet, record it as unmatched
         matchedJets.push_back(make_pair(mcJet, dummy));
      }
   }
   // Add the remaining unmatched reco jets
   for (const auto &recoJet : recoJetsCopy)
      matchedJets.push_back(make_pair(dummy, recoJet));

   return matchedJets;
}