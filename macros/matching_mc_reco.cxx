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
   int mult;
   double vz;  ///< primary vertex z (cm)
   float event_sum_pt;
   bool isTriggerEvent;
   bool trigger_match_HT2[1000];
   double neutral_fraction[1000];
   bool trigger_match_JP2[1000];
   bool trigger_match_JP1[1000];
   bool trigger_match_JP0[1000];
   double ptLead[1000];
   double pt[1000];
   double pt_corrected[1000]; // off-axis-cones UE-subtracted reco pT (= pt for MC)
   double det_eta[1000];      // detector eta (BEMC-projected; reco side only)
   int n_constituents[1000];
   int index[1000];
};

struct MyJet {
   TStarJetVectorJet orig;
   double pt;
   double pt_corrected; // UE-subtracted (off-axis cones)
   double ptLead;
   double det_eta = -9; // detector eta (reco side; -9 for MC/missing)
   double eta;
   double phi;
   double y;
   double neutral_fraction;
   bool trigger_match_JP2;
   bool trigger_match_JP1;
   bool trigger_match_JP0;
   int n_constituents;
   int event_id;
   double weight;
   double multiplicity;
   bool trigger_match_HT2;

   float deltaR(const MyJet &other) const
   {
      // Sentinel for "missing" jet is pt = -999.  Real (UE-corrected) reco pt
      // can be slightly negative for low-pT jets, so compare against the
      // sentinel rather than 0.
      if (pt < -500 || other.pt < -500) {
         return 10000; // invalid
      }
      float deta = eta - other.eta;
      float dphi = TVector2::Phi_mpi_pi(phi - other.phi);
      return sqrt(deta * deta + dphi * dphi);
   }
   MyJet()
      // pt sentinels = -999 (not -9): pt_corrected can dip slightly below 0
      // for low-pT real jets after UE subtraction, so the "missing" sentinel
      // must be far below any plausible physical value.
      : pt(-999),
        pt_corrected(-999),
        ptLead(-999),
        eta(-9),
        phi(-9),
        y(-9),
        neutral_fraction(-9),
        trigger_match_JP2(false),
        trigger_match_JP1(false),
        trigger_match_JP0(false),
        n_constituents(-9),
        event_id(-9),
        weight(-9),
        multiplicity(-9),
        trigger_match_HT2(false) {};

   MyJet(TStarJetVectorJet _orig, double _pt, double _pt_corrected, double _ptLead, double _eta, double _phi,
         double _y, double _neutral_fraction, bool _trigger_match_JP2, bool _trigger_match_JP1,
         bool _trigger_match_JP0, double n_constituents, int _event_id, double _weight, double _multiplicity,
         bool _trigger_match_HT2)
      : orig(_orig),
        pt(_pt),
        pt_corrected(_pt_corrected),
        ptLead(_ptLead),
        eta(_eta),
        phi(_phi),
        y(_y),
        neutral_fraction(_neutral_fraction),
        trigger_match_JP2(_trigger_match_JP2),
        trigger_match_JP1(_trigger_match_JP1),
        trigger_match_JP0(_trigger_match_JP0),
        n_constituents(n_constituents),
        event_id(_event_id),
        weight(_weight),
        multiplicity(_multiplicity),
        trigger_match_HT2(_trigger_match_HT2) {};
};

typedef pair<MyJet, MyJet> MatchedJetPair;

vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets, const vector<MyJet> &RecoJets, const double &R);

int matching_mc_reco(TString mcTreeName = "/gpfs01/star/pwg/prozorov/jets_pp_2012/output/mc/"
                                          "tree_pt-hat3545_033_R0.5.root",
                     bool isTest = false, TString testDir = "")
{

   TString baseName = mcTreeName(mcTreeName.Last('/') + 1, mcTreeName.Length());
   // if it does not has .root extension, add it
   if (!baseName.EndsWith(".root")) {
      baseName += ".root";
   }
   // get Trigger from input filename
   TString Trigger = "MB";
   if (mcTreeName.Contains("JP2")) {
      Trigger = "JP2";
   } else if (mcTreeName.Contains("JP1")) {
      Trigger = "JP1";
   } else if (mcTreeName.Contains("JP0")) {
      Trigger = "JP0";
   } else if (mcTreeName.Contains("HT2")) {
      Trigger = "HT2";
   }

   TString dirName = "/gpfs01/star/pwg/prozorov/jets_pp_2012/output/" + Trigger + "/embedding/";
   if (isTest) {
      // testDir overrides the default test path; trailing slash auto-appended.
      if (testDir.Length() > 0) {
         dirName = testDir;
         if (!dirName.EndsWith("/")) dirName += "/";
      } else {
         dirName = "/gpfs01/star/pwg/prozorov/jets_pp_2012/";
      }
   }
   // strip the leading "mc_" prefix only — ReplaceAll would also eat any
   // "mc_" that happens to appear inside the filename (e.g. loadtest_mc).
   if (baseName.BeginsWith("mc_"))
      baseName.Remove(0, 3);
   TString geantBaseName = "geant_" + baseName;
   TString mcBaseName = "mc_" + baseName;

   // Write the matched output alongside the geant/mc inputs when in test mode
   // so all artifacts of one test run live in the same dir.
   TString OutFile = (isTest ? dirName : TString("")) + "matched_" + baseName;

   TString RecoFile = dirName + geantBaseName;
   TString McFile = dirName + mcBaseName;

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

   // Parse the Pythia pt-hat bin from the input filename so downstream
   // analyses can apply Dmitry's soft pT reweight (analysis note Eq. 6).
   // Filename convention: "..._pt_hat<lo><hi>_<batch>_R<R>.root", e.g.
   // pt_hat1115_002 → (lo=11, hi=15) → midpoint 13 GeV.
   double pthat_mid = -1;
   {
      const struct { const char *tag; double mid; } kPtHatBins[] = {
         {"pt_hat23_",    2.5},   {"pt_hat34_",    3.5},
         {"pt_hat45_",    4.5},   {"pt_hat57_",    6.0},
         {"pt_hat79_",    8.0},   {"pt_hat911_",  10.0},
         {"pt_hat1115_", 13.0},   {"pt_hat1520_", 17.5},
         {"pt_hat2025_", 22.5},   {"pt_hat2535_", 30.0},
         {"pt_hat3545_", 40.0},   {"pt_hat4555_", 50.0},
         {"pt_hat55999_", 65.0},  {"pt_hat55_",   65.0},
      };
      for (const auto &b : kPtHatBins) {
         if (mcTreeName.Contains(b.tag)) { pthat_mid = b.mid; break; }
      }
      if (pthat_mid < 0)
         cout << "WARN: could not parse pt-hat bin from " << mcTreeName
              << " (soft reweight will get pthat_mid=-1)" << endl;
   }
   const float EtaCut = 1.0 - RCut;
   // RECO-side eta cut loosened to |physics eta|<1.0 so reco jets in the
   // 0.5..(1.0-RCut) / |det_eta|<0.5 slice survive matching; the downstream unfold
   // applies the SYMMETRIC |reco_det_eta|<0.5 (Dmitry 00_05_00_05). The TRUTH/MC
   // cut stays on EtaCut (|physics eta|<0.5), consistent with the wide jet-find in
   // ppAnalysis.cxx.
   const float EtaCutReco = 1.0f;
   // =================================================================================================
   TFile *Mcf = new TFile(McFile, "READ");
   if (Mcf->IsZombie() || !Mcf->Get("ResultTree")) {
      cerr << "FATAL: cannot open mc input (or no ResultTree): " << McFile << endl;
      return 1;
   }
   TH1D *hEventsRun = (TH1D *)Mcf->Get("hEventsRun");

   TTree *McChain = (TTree *)Mcf->Get("ResultTree");
   McChain->BuildIndex("runid", "eventid");

   InputTreeEntry mc;
   // Initialise per-jet trigger flag arrays to false for trees that pre-date
   // the JP1/JP0 branch addition (RunppAna 2026-05-17). Without this, jets
   // read from old trees would inherit whatever garbage lives in mc.trigger_match_JP1[].
   std::fill(std::begin(mc.trigger_match_JP1), std::end(mc.trigger_match_JP1), false);
   std::fill(std::begin(mc.trigger_match_JP0), std::end(mc.trigger_match_JP0), false);
   McChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
   McChain->SetBranchAddress("Jets", &mc.jets);
   McChain->SetBranchAddress("eventid", &mc.eventid);
   McChain->SetBranchAddress("runid", &mc.runid);
   McChain->SetBranchAddress("weight", &mc.weight);
   McChain->SetBranchAddress("njets", &mc.njets);
   McChain->SetBranchAddress("mult", &mc.mult);
   McChain->SetBranchAddress("vz", &mc.vz);
   McChain->SetBranchAddress("event_sum_pt", &mc.event_sum_pt);
   McChain->SetBranchAddress("trigger_match_HT2", mc.trigger_match_HT2);
   McChain->SetBranchAddress("neutral_fraction", mc.neutral_fraction);
   McChain->SetBranchAddress("trigger_match_JP2", mc.trigger_match_JP2);
   if (McChain->GetBranch("trigger_match_JP1"))
      McChain->SetBranchAddress("trigger_match_JP1", mc.trigger_match_JP1);
   if (McChain->GetBranch("trigger_match_JP0"))
      McChain->SetBranchAddress("trigger_match_JP0", mc.trigger_match_JP0);
   McChain->SetBranchAddress("pt", mc.pt);
   McChain->SetBranchAddress("pt_corrected", mc.pt_corrected);
   McChain->SetBranchAddress("ptLead", mc.ptLead);
   McChain->SetBranchAddress("n_constituents", mc.n_constituents);

   TFile *Recof = new TFile(RecoFile, "READ");
   if (Recof->IsZombie() || !Recof->Get("ResultTree")) {
      cerr << "FATAL: cannot open reco input (or no ResultTree): " << RecoFile << endl;
      return 1;
   }
   TTree *RecoChain = (TTree *)Recof->Get("ResultTree");
   RecoChain->BuildIndex("runid", "eventid");
   RecoChain->GetBranch("Jets")->SetAutoDelete(kFALSE);

   InputTreeEntry reco;
   std::fill(std::begin(reco.trigger_match_JP1), std::end(reco.trigger_match_JP1), false);
   std::fill(std::begin(reco.trigger_match_JP0), std::end(reco.trigger_match_JP0), false);
   RecoChain->SetBranchAddress("Jets", &reco.jets);
   RecoChain->SetBranchAddress("eventid", &reco.eventid);
   RecoChain->SetBranchAddress("runid", &reco.runid);
   RecoChain->SetBranchAddress("weight", &reco.weight);
   RecoChain->SetBranchAddress("njets", &reco.njets);
   RecoChain->SetBranchAddress("mult", &reco.mult);
   RecoChain->SetBranchAddress("vz", &reco.vz);
   RecoChain->SetBranchAddress("event_sum_pt", &reco.event_sum_pt);
   RecoChain->SetBranchAddress("trigger_match_HT2", reco.trigger_match_HT2);
   RecoChain->SetBranchAddress("neutral_fraction", reco.neutral_fraction);
   RecoChain->SetBranchAddress("trigger_match_JP2", reco.trigger_match_JP2);
   if (RecoChain->GetBranch("trigger_match_JP1"))
      RecoChain->SetBranchAddress("trigger_match_JP1", reco.trigger_match_JP1);
   if (RecoChain->GetBranch("trigger_match_JP0"))
      RecoChain->SetBranchAddress("trigger_match_JP0", reco.trigger_match_JP0);
   RecoChain->SetBranchAddress("pt", reco.pt);
   RecoChain->SetBranchAddress("pt_corrected", reco.pt_corrected);
   RecoChain->SetBranchAddress("det_eta", reco.det_eta);
   RecoChain->SetBranchAddress("ptLead", reco.ptLead);
   RecoChain->SetBranchAddress("n_constituents", reco.n_constituents);
   RecoChain->SetBranchAddress("isTriggerEvent", &reco.isTriggerEvent);

   TFile *fout = new TFile(OutFile, "RECREATE");
   TH1D *hDeltaR = new TH1D("hDeltaR", "#Delta R all; #Delta R", 350, 0, 3.5);
   TH1D *hDeltaRMatched = new TH1D("hDeltaRMatched", "#Delta R matched; #Delta R", 350, 0, 3.5);
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
   bool isTriggerEvent = false;

   MatchedTree->Branch("mc_pt", &outMcJet.pt, "mc_pt/D");
   MatchedTree->Branch("mc_pt_corrected", &outMcJet.pt_corrected, "mc_pt_corrected/D");
   MatchedTree->Branch("mc_ptLead", &outMcJet.ptLead, "mc_ptLead/D");
   MatchedTree->Branch("mc_eta", &outMcJet.eta, "mc_eta/D");
   MatchedTree->Branch("mc_phi", &outMcJet.phi, "mc_phi/D");
   MatchedTree->Branch("mc_neutral_fraction", &outMcJet.neutral_fraction, "mc_neutral_fraction/D");
   MatchedTree->Branch("mc_trigger_match_JP2", &outMcJet.trigger_match_JP2, "mc_trigger_match_JP2/O");
   MatchedTree->Branch("mc_trigger_match_JP1", &outMcJet.trigger_match_JP1, "mc_trigger_match_JP1/O");
   MatchedTree->Branch("mc_trigger_match_JP0", &outMcJet.trigger_match_JP0, "mc_trigger_match_JP0/O");

   MatchedTree->Branch("mc_n_constituents", &outMcJet.n_constituents, "mc_n_constituents/I");
   MatchedTree->Branch("mc_weight", &outMcJet.weight, "mc_weight/D");
   MatchedTree->Branch("mc_multiplicity", &outMcJet.multiplicity, "mc_multiplicity/D");
   MatchedTree->Branch("mc_trigger_match_HT2", &outMcJet.trigger_match_HT2, "mc_trigger_match_HT2/O");

   MatchedTree->Branch("reco_pt", &outRecoJet.pt, "reco_pt/D");
   MatchedTree->Branch("reco_pt_corrected", &outRecoJet.pt_corrected, "reco_pt_corrected/D");
   MatchedTree->Branch("reco_ptLead", &outRecoJet.ptLead, "reco_ptLead/D");
   MatchedTree->Branch("reco_eta", &outRecoJet.eta, "reco_eta/D");
   // Detector eta of the reco jet (BEMC-projected, vz-dependent) — enables the
   // SYMMETRIC |det_eta|<0.5 cut in unfold.cxx, matching Dmitry's 00_05_00_05
   // selection (|eta|<0.5 && |detEta|<0.5 on BOTH data and response).
   MatchedTree->Branch("reco_det_eta", &outRecoJet.det_eta, "reco_det_eta/D");
   MatchedTree->Branch("reco_phi", &outRecoJet.phi, "reco_phi/D");
   MatchedTree->Branch("reco_neutral_fraction", &outRecoJet.neutral_fraction, "reco_neutral_fraction/D");
   MatchedTree->Branch("reco_trigger_match_JP2", &outRecoJet.trigger_match_JP2, "reco_trigger_match_JP2/O");
   MatchedTree->Branch("reco_trigger_match_JP1", &outRecoJet.trigger_match_JP1, "reco_trigger_match_JP1/O");
   MatchedTree->Branch("reco_trigger_match_JP0", &outRecoJet.trigger_match_JP0, "reco_trigger_match_JP0/O");
   MatchedTree->Branch("reco_n_constituents", &outRecoJet.n_constituents, "reco_n_constituents/I");
   MatchedTree->Branch("reco_weight", &outRecoJet.weight, "reco_weight/D");
   MatchedTree->Branch("reco_multiplicity", &outRecoJet.multiplicity, "reco_multiplicity/D");
   MatchedTree->Branch("reco_trigger_match_HT2", &outRecoJet.trigger_match_HT2, "reco_trigger_match_HT2/O");
   MatchedTree->Branch("isTriggerEvent", &isTriggerEvent, "isTriggerEvent/O");
   MatchedTree->Branch("deltaR", &deltaR, "deltaR/D");
   // Primary vertex z taken from the MC tree (always present for embedding;
   // matches the reco vz to within ~5 cm by construction).  Stored on every
   // MatchedTree row so per-jet weights can pull it via RDataFrame Define.
   double event_vz = -999;
   MatchedTree->Branch("event_vz", &event_vz, "event_vz/D");
   // Midpoint of the Pythia pt-hat bin this event came from (parsed once at
   // tree-open time from the input filename).  Lets downstream code apply
   // Dmitry's soft pT reweight per event without needing the true partonic pT
   // (which is not stored in the pico).
   MatchedTree->Branch("pthat_mid", &pthat_mid, "pthat_mid/D");
   // Event identity + per-EVENT emulated shouldFire proxies (2026-06-10, ported
   // from the JPx fork): evt_should_JPi = OR over ALL reco jets in the event of
   // trigger_match_JPi, computed BEFORE any eta cut — enables Dmitry's per-event
   // promotion partition downstream with the IDENTICAL definition on data.
   int out_runid = 0, out_eventid = 0;
   bool evt_should_JP0 = false, evt_should_JP1 = false, evt_should_JP2 = false;
   MatchedTree->Branch("runid", &out_runid, "runid/I");
   MatchedTree->Branch("eventid", &out_eventid, "eventid/I");
   MatchedTree->Branch("evt_should_JP0", &evt_should_JP0, "evt_should_JP0/O");
   MatchedTree->Branch("evt_should_JP1", &evt_should_JP1, "evt_should_JP1/O");
   MatchedTree->Branch("evt_should_JP2", &evt_should_JP2, "evt_should_JP2/O");

   int nEvents = McChain->GetEntries();

   set<float> accepted_events_list;
   int MatchNumber = 0;
   int FakeNumber = 0;
   int MissNumber = 0;

   for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) // event loop
   {
      McChain->GetEntry(iEvent);
      event_vz = mc.vz;

      if (accepted_events_list.count(mc.event_sum_pt) > 0)
         continue; // some events have identical total particle pT
      else
         accepted_events_list.insert(mc.event_sum_pt);

      vector<MyJet> mcJets;
      for (int j = 0; j < mc.njets; ++j) {
         TStarJetVectorJet *tempMcJet = dynamic_cast<TStarJetVectorJet *>(mc.jets->At(j));
         if (abs(tempMcJet->Eta()) > EtaCut)
            continue;

         mcJets.push_back(MyJet(*tempMcJet, tempMcJet->Pt(), mc.pt_corrected[j], mc.ptLead[j], tempMcJet->Eta(),
                                tempMcJet->Phi(), tempMcJet->Rapidity(), mc.neutral_fraction[j],
                                mc.trigger_match_JP2[j], mc.trigger_match_JP1[j], mc.trigger_match_JP0[j],
                                mc.n_constituents[j], mc.eventid, mc.weight, mc.mult,
                                mc.trigger_match_HT2[j]));
      }

      int recoEvent = RecoChain->GetEntryNumberWithIndex(mc.runid, mc.eventid);
      vector<MyJet> recoJets;

      out_runid = mc.runid;
      out_eventid = mc.eventid;
      evt_should_JP0 = evt_should_JP1 = evt_should_JP2 = false;

      if (recoEvent >= 0) {
         RecoChain->GetEntry(recoEvent);
         isTriggerEvent = reco.isTriggerEvent;
         // Per-event shouldFire proxy over ALL reco jets (no eta cut).
         for (int j = 0; j < reco.njets; ++j) {
            evt_should_JP0 |= reco.trigger_match_JP0[j];
            evt_should_JP1 |= reco.trigger_match_JP1[j];
            evt_should_JP2 |= reco.trigger_match_JP2[j];
         }
         for (int j = 0; j < reco.njets; ++j) {
            TStarJetVectorJet *tempRecoJet = (TStarJetVectorJet *)reco.jets->At(j);
            // RECO physics-eta cut uses EtaCutReco (=1.0, wide; unfold then
            // applies the symmetric det_eta<0.5).
            if (abs(tempRecoJet->Eta()) > EtaCutReco)
               continue;
            // Use the UE-subtracted pT for reco; kinematics (eta, phi, y) come
            // from the original 4-vector since UE subtraction is a scalar
            // pT shift, not a recompute of the jet axis.
            // reco MyJet.pt stays UE-subtracted (legacy downstream uses it as such).
            // pt_corrected duplicates pt for clarity / parity with MC side.
            MyJet rj(*tempRecoJet, reco.pt_corrected[j], reco.pt_corrected[j], reco.ptLead[j],
                     tempRecoJet->Eta(), tempRecoJet->Phi(), tempRecoJet->Rapidity(),
                     reco.neutral_fraction[j], reco.trigger_match_JP2[j],
                     reco.trigger_match_JP1[j], reco.trigger_match_JP0[j],
                     reco.n_constituents[j],
                     reco.eventid, reco.weight, reco.mult, reco.trigger_match_HT2[j]);
            rj.det_eta = reco.det_eta[j];
            recoJets.push_back(rj);
         }
      } // end of recojet loop}
      //======================================================================================================================================
      // Match them
      //======================================================================================================================================

      vector<MatchedJetPair> MatchedJets;
      MatchedJets = MatchJetsEtaPhi(mcJets, recoJets, RCut);

      for (unsigned int j = 0; j < MatchedJets.size(); j++) {
         outMcJet = MatchedJets[j].first;
         outRecoJet = MatchedJets[j].second;
         deltaR = outMcJet.deltaR(outRecoJet);
         // "missing" sentinel = -999; real (UE-corrected reco) pt can be < 0.
         const bool mcReal   = outMcJet.pt   > -500;
         const bool recoReal = outRecoJet.pt > -500;
         if (!mcReal && recoReal) {
            outMcJet.weight = outRecoJet.weight;
         }

         MatchedTree->Fill();
         // fill histograms and counters
         if (mcReal && recoReal) { // matched
            hDeltaRMatched->Fill(deltaR);
            hPtMcReco->Fill(outMcJet.pt, outRecoJet.pt);
            MatchNumber++;
         } else if (mcReal && !recoReal) { // missed
            hMiss->Fill(outMcJet.pt);
            MissNumber++;
         } else if (!mcReal && recoReal) { // fake
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
   hDeltaRMatched->Write();
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
         // find DeltaR histogram in system path and fill it
         TH1D *hDeltaR = (TH1D *)gDirectory->Get("hDeltaR");
         if (hDeltaR) {
            hDeltaR->Fill(deltaR);
         }

         // Tight matching cone (Dmitry thesis 5.5.1: ΔR < 0.2 chosen from the
         // ΔR distribution).  Independent of jet R to give apples-to-apples
         // generalized efficiency vs Dmitry.
         const double kMaxDR = 0.2;
         if (deltaR <= kMaxDR && deltaR < minDeltaR) {
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