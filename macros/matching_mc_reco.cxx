// Builds MatchedTree: pairs embedding reco jets with their truth partners
// (geant_*.root x mc_*.root -> matched_*.root), writing misses and fakes too.
//
//   root -l -b -q 'macros/matching_mc_reco.cxx+("mc_<stem>_R<R>.root", true, "<outdir>/")'
//
// Trigger, jet radius and Pythia pt-hat bin are all parsed from the input
// filename; container.sh drives it in the embedding / mbembedding modes.

#include <algorithm>
#include "TClonesArray.h"
#include "TFile.h"
#include "TSystem.h"
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
   int njets;
   int mult;
   double vz;         ///< primary vertex z (cm)
   double pthat = -1; ///< event pt-hat (MC tree only)
   double vz_vpd = 0; ///< VPD vertex z (reco side; MB vertex-efficiency study)
   double vx = 0;     ///< primary vertex x
   double vy = 0;     ///< primary vertex y
   float event_sum_pt;
   bool isTriggerEvent;
   bool trigger_match_HT2[1000];
   double neutral_fraction[1000];
   bool trigger_match_JP2[1000];
   bool trigger_match_JP1[1000];
   bool trigger_match_JP0[1000];
   double ptLead[1000];
   double pt[1000];
   double pt_corrected[1000]; // off-axis-cones UE-subtracted pT
   double jet_area[1000];     // FastJet active area
   double bg_density[1000];   // off-axis-cones UE density rho (GeV / unit area)
   double det_eta[1000];      // detector eta (BEMC-projected; reco side only)
   int n_constituents[1000];
   int jp_patch_adc[1000];    // simulator ADCs (DSM scale), reco side: per jet and event maxima, thresholds
   int ht_adc_max[1000];
   int evt_jp_adc_max = -1;
   int evt_ht_adc_max = -1;
   int jp_thr[3] = {0, 0, 0};
   int ht_thr[4] = {0, 0, 0, 0};
};

struct MyJet {
   TStarJetVectorJet orig;
   // pt = RAW jet pT, pt_corrected = off-axis-cones UE-subtracted, on BOTH sides.
   double pt;
   double pt_corrected;
   double area = -999;       // FastJet active area
   double bg_density = -999; // off-axis-cones UE density rho
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
   int jp_patch_adc = -1; // simulator ADCs (DSM scale), reco side
   int ht_adc_max = -1;

   float deltaR(const MyJet &other) const
   {
      // "Missing" sentinel is pt = -999; real UE-corrected pt can be slightly
      // negative, so test against the sentinel and not against 0.
      if (pt < -500 || other.pt < -500) {
         return 10000; // invalid
      }
      float deta = eta - other.eta;
      float dphi = TVector2::Phi_mpi_pi(phi - other.phi);
      return sqrt(deta * deta + dphi * dphi);
   }
   MyJet()
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
   if (!baseName.EndsWith(".root")) {
      baseName += ".root";
   }
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
   // Strip the LEADING "mc_" only — ReplaceAll would also eat an "mc_" that
   // happens to appear inside the filename (e.g. loadtest_mc).
   if (baseName.BeginsWith("mc_"))
      baseName.Remove(0, 3);
   TString geantBaseName = "geant_" + baseName;
   TString mcBaseName = "mc_" + baseName;

   // In test mode the matched output lands next to its geant/mc inputs.
   TString OutFile = (isTest ? dirName : TString("")) + "matched_" + baseName;

   TString RecoFile = dirName + geantBaseName;
   TString McFile = dirName + mcBaseName;

   float RCut = 0;
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

   // Pythia pt-hat bin -> bin midpoint, so downstream code can apply the soft
   // pT reweight. Filename convention "..._pt_hat<lo><hi>_<batch>_R<R>.root",
   // e.g. pt_hat1115_002 -> midpoint 13 GeV.
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
         // 2021 request (e21_<bin>_<run> picos). The trailing underscore keeps
         // pt2_3 from matching inside pt20_25.
         {"e21_pt2_3_",   2.5},   {"e21_pt3_4_",   3.5},
         {"e21_pt4_5_",   4.5},   {"e21_pt5_7_",   6.0},
         {"e21_pt7_9_",   8.0},   {"e21_pt9_11_", 10.0},
         {"e21_pt11_15_", 13.0},  {"e21_pt15_20_", 17.5},
         {"e21_pt20_25_", 22.5},  {"e21_pt25_35_", 30.0},
         {"e21_pt35_45_", 40.0},  {"e21_pt45_55_", 50.0},
         {"e21_pt55_-1_", 65.0},
      };
      for (const auto &b : kPtHatBins) {
         if (mcTreeName.Contains(b.tag)) { pthat_mid = b.mid; break; }
      }
      if (pthat_mid < 0)
         cout << "WARN: could not parse pt-hat bin from " << mcTreeName
              << " (soft reweight will get pthat_mid=-1)" << endl;
   }
   // Both sides are kept wide (|eta| < 1.0, the Stage-1 EtaJetCut) so that a jet
   // near the edge of the fiducial still finds its partner; the analysis
   // acceptance (symmetric |reco_det_eta| < 0.5) is applied downstream in unfold.
   const float EtaCut = 1.0;
   const float EtaCutReco = 1.0f;

   TFile *Mcf = new TFile(McFile, "READ");
   if (Mcf->IsZombie() || !Mcf->Get("ResultTree")) {
      cerr << "FATAL: cannot open mc input (or no ResultTree): " << McFile << endl;
      return 1;
   }
   TH1D *hEventsRun = (TH1D *)Mcf->Get("hEventsRun");

   TTree *McChain = (TTree *)Mcf->Get("ResultTree");
   McChain->BuildIndex("runid", "eventid");

   InputTreeEntry mc;
   // Trees that pre-date the JP1/JP0 branches leave these arrays unwritten;
   // clear them so their jets do not inherit garbage. (Same for reco below,
   // and for every optional branch guarded by GetBranch() != nullptr.)
   std::fill(std::begin(mc.trigger_match_JP1), std::end(mc.trigger_match_JP1), false);
   std::fill(std::begin(mc.trigger_match_JP0), std::end(mc.trigger_match_JP0), false);
   McChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
   McChain->SetBranchAddress("Jets", &mc.jets);
   McChain->SetBranchAddress("eventid", &mc.eventid);
   McChain->SetBranchAddress("runid", &mc.runid);
   McChain->SetBranchAddress("runid1", &mc.runid1);
   McChain->SetBranchAddress("weight", &mc.weight);
   McChain->SetBranchAddress("njets", &mc.njets);
   McChain->SetBranchAddress("mult", &mc.mult);
   mc.vz = 0;
   if (McChain->GetBranch("vz")) McChain->SetBranchAddress("vz", &mc.vz);
   if (McChain->GetBranch("pthat")) McChain->SetBranchAddress("pthat", &mc.pthat);
   McChain->SetBranchAddress("event_sum_pt", &mc.event_sum_pt);
   McChain->SetBranchAddress("trigger_match_HT2", mc.trigger_match_HT2);
   McChain->SetBranchAddress("neutral_fraction", mc.neutral_fraction);
   McChain->SetBranchAddress("trigger_match_JP2", mc.trigger_match_JP2);
   if (McChain->GetBranch("trigger_match_JP1"))
      McChain->SetBranchAddress("trigger_match_JP1", mc.trigger_match_JP1);
   if (McChain->GetBranch("trigger_match_JP0"))
      McChain->SetBranchAddress("trigger_match_JP0", mc.trigger_match_JP0);
   McChain->SetBranchAddress("pt", mc.pt);
   std::fill_n(mc.pt_corrected, sizeof(mc.pt_corrected)/sizeof(*mc.pt_corrected), 0.0);
   if (McChain->GetBranch("pt_corrected")) McChain->SetBranchAddress("pt_corrected", mc.pt_corrected);
   std::fill_n(mc.jet_area, sizeof(mc.jet_area)/sizeof(*mc.jet_area), 0.0);
   if (McChain->GetBranch("jet_area")) McChain->SetBranchAddress("jet_area", mc.jet_area);
   std::fill_n(mc.bg_density, sizeof(mc.bg_density)/sizeof(*mc.bg_density), 0.0);
   if (McChain->GetBranch("bg_density")) McChain->SetBranchAddress("bg_density", mc.bg_density);
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
   RecoChain->SetBranchAddress("runid1", &reco.runid1);
   RecoChain->SetBranchAddress("weight", &reco.weight);
   // Event-level should_JPn = the full trigger-simulator decision, the same
   // ruler the data gate uses. Trees without the branches fall back to the
   // per-jet trigger-match OR computed in the event loop below.
   bool reco_evt_should[3] = {false, false, false};
   const bool haveEvtShould = RecoChain->GetBranch("should_JP0") != nullptr;
   if (haveEvtShould) {
      RecoChain->SetBranchAddress("should_JP0", &reco_evt_should[0]);
      RecoChain->SetBranchAddress("should_JP1", &reco_evt_should[1]);
      RecoChain->SetBranchAddress("should_JP2", &reco_evt_should[2]);
   }
   bool reco_evt_should_ht2 = false;
   const bool haveEvtShouldHT2 = RecoChain->GetBranch("should_HT2") != nullptr;
   if (haveEvtShouldHT2) RecoChain->SetBranchAddress("should_HT2", &reco_evt_should_ht2);
   RecoChain->SetBranchAddress("njets", &reco.njets);
   RecoChain->SetBranchAddress("mult", &reco.mult);
   reco.vz = 0;
   if (RecoChain->GetBranch("vz")) RecoChain->SetBranchAddress("vz", &reco.vz);
   // Reconstructed vertex, propagated for the MB vertex-efficiency study.
   if (RecoChain->GetBranch("vz_vpd")) RecoChain->SetBranchAddress("vz_vpd", &reco.vz_vpd);
   if (RecoChain->GetBranch("vx")) RecoChain->SetBranchAddress("vx", &reco.vx);
   if (RecoChain->GetBranch("vy")) RecoChain->SetBranchAddress("vy", &reco.vy);
   RecoChain->SetBranchAddress("event_sum_pt", &reco.event_sum_pt);
   RecoChain->SetBranchAddress("trigger_match_HT2", reco.trigger_match_HT2);
   RecoChain->SetBranchAddress("neutral_fraction", reco.neutral_fraction);
   RecoChain->SetBranchAddress("trigger_match_JP2", reco.trigger_match_JP2);
   if (RecoChain->GetBranch("trigger_match_JP1"))
      RecoChain->SetBranchAddress("trigger_match_JP1", reco.trigger_match_JP1);
   if (RecoChain->GetBranch("trigger_match_JP0"))
      RecoChain->SetBranchAddress("trigger_match_JP0", reco.trigger_match_JP0);
   RecoChain->SetBranchAddress("pt", reco.pt);
   std::fill_n(reco.pt_corrected, sizeof(reco.pt_corrected)/sizeof(*reco.pt_corrected), 0.0);
   if (RecoChain->GetBranch("pt_corrected")) RecoChain->SetBranchAddress("pt_corrected", reco.pt_corrected);
   std::fill_n(reco.jet_area, sizeof(reco.jet_area)/sizeof(*reco.jet_area), 0.0);
   if (RecoChain->GetBranch("jet_area")) RecoChain->SetBranchAddress("jet_area", reco.jet_area);
   std::fill_n(reco.bg_density, sizeof(reco.bg_density)/sizeof(*reco.bg_density), 0.0);
   if (RecoChain->GetBranch("bg_density")) RecoChain->SetBranchAddress("bg_density", reco.bg_density);
   std::fill_n(reco.det_eta, sizeof(reco.det_eta)/sizeof(*reco.det_eta), 0.0);
   std::fill_n(reco.jp_patch_adc, 1000, -1);
   std::fill_n(reco.ht_adc_max, 1000, -1);
   if (RecoChain->GetBranch("jp_patch_adc")) RecoChain->SetBranchAddress("jp_patch_adc", reco.jp_patch_adc);
   if (RecoChain->GetBranch("ht_adc_max")) RecoChain->SetBranchAddress("ht_adc_max", reco.ht_adc_max);
   if (RecoChain->GetBranch("evt_jp_adc_max")) RecoChain->SetBranchAddress("evt_jp_adc_max", &reco.evt_jp_adc_max);
   if (RecoChain->GetBranch("evt_ht_adc_max")) RecoChain->SetBranchAddress("evt_ht_adc_max", &reco.evt_ht_adc_max);
   if (RecoChain->GetBranch("jp_thr")) RecoChain->SetBranchAddress("jp_thr", reco.jp_thr);
   if (RecoChain->GetBranch("ht_thr")) RecoChain->SetBranchAddress("ht_thr", reco.ht_thr);
   if (RecoChain->GetBranch("det_eta")) RecoChain->SetBranchAddress("det_eta", reco.det_eta);
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
   MatchedTree->Branch("mc_jet_area", &outMcJet.area, "mc_jet_area/D");
   MatchedTree->Branch("mc_bg_density", &outMcJet.bg_density, "mc_bg_density/D");
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

   // reco_pt is RAW, reco_pt_corrected is UE-subtracted (same as the mc_ pair).
   MatchedTree->Branch("reco_pt", &outRecoJet.pt, "reco_pt/D");
   MatchedTree->Branch("reco_pt_corrected", &outRecoJet.pt_corrected, "reco_pt_corrected/D");
   MatchedTree->Branch("reco_jet_area", &outRecoJet.area, "reco_jet_area/D");
   MatchedTree->Branch("reco_bg_density", &outRecoJet.bg_density, "reco_bg_density/D");
   MatchedTree->Branch("reco_ptLead", &outRecoJet.ptLead, "reco_ptLead/D");
   MatchedTree->Branch("reco_eta", &outRecoJet.eta, "reco_eta/D");
   // Detector eta (BEMC-projected, vz-dependent): feeds the symmetric
   // |eta|<0.5 && |det_eta|<0.5 selection unfold.cxx applies to data and response.
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
   MatchedTree->Branch("reco_jp_patch_adc", &outRecoJet.jp_patch_adc, "reco_jp_patch_adc/I");
   MatchedTree->Branch("reco_ht_adc_max", &outRecoJet.ht_adc_max, "reco_ht_adc_max/I");
   MatchedTree->Branch("isTriggerEvent", &isTriggerEvent, "isTriggerEvent/O");
   MatchedTree->Branch("deltaR", &deltaR, "deltaR/D");
   // Primary vertex z from the MC tree, repeated on every row so per-jet
   // weights can pull it via an RDataFrame Define.
   double event_vz = -999;
   MatchedTree->Branch("event_vz", &event_vz, "event_vz/D");
   // Reconstructed vertex (reco side) for the MB vertex-efficiency study.
   double reco_vz = -999, reco_vz_vpd = -999, reco_vx = -999, reco_vy = -999;
   int evt_jp_adc_max = -1, evt_ht_adc_max = -1, jp_thr[3] = {0, 0, 0}, ht_thr[4] = {0, 0, 0, 0};
   MatchedTree->Branch("evt_jp_adc_max", &evt_jp_adc_max, "evt_jp_adc_max/I");
   MatchedTree->Branch("evt_ht_adc_max", &evt_ht_adc_max, "evt_ht_adc_max/I");
   MatchedTree->Branch("jp_thr", jp_thr, "jp_thr[3]/I");
   MatchedTree->Branch("ht_thr", ht_thr, "ht_thr[4]/I");
   MatchedTree->Branch("reco_vz", &reco_vz, "reco_vz/D");
   MatchedTree->Branch("reco_vz_vpd", &reco_vz_vpd, "reco_vz_vpd/D");
   MatchedTree->Branch("reco_vx", &reco_vx, "reco_vx/D");
   MatchedTree->Branch("reco_vy", &reco_vy, "reco_vy/D");
   // Midpoint of this file's pt-hat bin, constant per file.
   MatchedTree->Branch("pthat_mid", &pthat_mid, "pthat_mid/D");
   double pthat_evt = -1; // per-event pt-hat from the MC tree; -1 for reco-only events
   MatchedTree->Branch("pthat", &pthat_evt, "pthat/D");
   // runid is the input-file hash that pairs MC and reco events (BuildIndex);
   // runid1 is the real run number from the pico header, so downstream code can
   // restrict the matched trees per run.
   int out_runid = 0, out_runid1 = 0, out_eventid = 0;
   // Per-EVENT trigger decision, computed before any eta cut so that the
   // downstream promotion partition uses an identical definition on data.
   bool evt_should_JP0 = false, evt_should_JP1 = false, evt_should_JP2 = false;
   MatchedTree->Branch("runid", &out_runid, "runid/I");
   MatchedTree->Branch("runid1", &out_runid1, "runid1/I");
   MatchedTree->Branch("eventid", &out_eventid, "eventid/I");
   MatchedTree->Branch("evt_should_JP0", &evt_should_JP0, "evt_should_JP0/O");
   MatchedTree->Branch("evt_should_JP1", &evt_should_JP1, "evt_should_JP1/O");
   MatchedTree->Branch("evt_should_JP2", &evt_should_JP2, "evt_should_JP2/O");
   bool evt_should_HT2 = false;
   MatchedTree->Branch("evt_should_HT2", &evt_should_HT2, "evt_should_HT2/O");

   int nEvents = McChain->GetEntries();

   set<float> accepted_events_list;
   std::set<int> visitedReco; // reco entries reached from the MC loop (the rest are written as fakes)
   int MatchNumber = 0;
   int FakeNumber = 0;
   int MissNumber = 0;

   for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) // event loop
   {
      McChain->GetEntry(iEvent);
      event_vz = mc.vz;
      pthat_evt = mc.pthat;
      reco_vz = reco_vz_vpd = reco_vx = reco_vy = -999; // reset; set below if reco matched
      evt_jp_adc_max = evt_ht_adc_max = -1;
      std::fill_n(jp_thr, 3, 0);
      std::fill_n(ht_thr, 4, 0);

      if (accepted_events_list.count(mc.event_sum_pt) > 0)
         continue; // some events have identical total particle pT
      else
         accepted_events_list.insert(mc.event_sum_pt);

      vector<MyJet> mcJets;
      for (int j = 0; j < mc.njets; ++j) {
         TStarJetVectorJet *tempMcJet = dynamic_cast<TStarJetVectorJet *>(mc.jets->At(j));
         if (abs(tempMcJet->Eta()) > EtaCut)
            continue;

         MyJet mj(*tempMcJet, tempMcJet->Pt(), mc.pt_corrected[j], mc.ptLead[j], tempMcJet->Eta(),
                  tempMcJet->Phi(), tempMcJet->Rapidity(), mc.neutral_fraction[j], mc.trigger_match_JP2[j],
                  mc.trigger_match_JP1[j], mc.trigger_match_JP0[j], mc.n_constituents[j], mc.eventid, mc.weight,
                  mc.mult, mc.trigger_match_HT2[j]);
         mj.area = mc.jet_area[j];
         mj.bg_density = mc.bg_density[j];
         mcJets.push_back(mj);
      }

      int recoEvent = RecoChain->GetEntryNumberWithIndex(mc.runid, mc.eventid);
      if (recoEvent >= 0) visitedReco.insert(recoEvent);
      vector<MyJet> recoJets;

      out_runid = mc.runid;
      out_runid1 = mc.runid1;
      out_eventid = mc.eventid;
      evt_should_JP0 = evt_should_JP1 = evt_should_JP2 = false;
      evt_should_HT2 = false;

      if (recoEvent >= 0) {
         RecoChain->GetEntry(recoEvent);
         isTriggerEvent = reco.isTriggerEvent;
         reco_vz = reco.vz;
         reco_vz_vpd = reco.vz_vpd;
         evt_jp_adc_max = reco.evt_jp_adc_max;
         evt_ht_adc_max = reco.evt_ht_adc_max;
         std::copy_n(reco.jp_thr, 3, jp_thr);
         std::copy_n(reco.ht_thr, 4, ht_thr);
         reco_vx = reco.vx;
         reco_vy = reco.vy;
         // Fallback: OR the per-jet trigger matches over ALL reco jets.
         for (int j = 0; j < reco.njets; ++j) {
            evt_should_JP0 |= reco.trigger_match_JP0[j];
            evt_should_JP1 |= reco.trigger_match_JP1[j];
            evt_should_JP2 |= reco.trigger_match_JP2[j];
         }
         if (haveEvtShould) { // event-level simulator decision, when available
            evt_should_JP0 = reco_evt_should[0];
            evt_should_JP1 = reco_evt_should[1];
            if (haveEvtShouldHT2) evt_should_HT2 = reco_evt_should_ht2;
            evt_should_JP2 = reco_evt_should[2];
         }
         for (int j = 0; j < reco.njets; ++j) {
            TStarJetVectorJet *tempRecoJet = (TStarJetVectorJet *)reco.jets->At(j);
            if (abs(tempRecoJet->Eta()) > EtaCutReco)
               continue;
            // Kinematics come from the original 4-vector: the UE subtraction is
            // a scalar pT shift and does not recompute the jet axis.
            MyJet rj(*tempRecoJet, tempRecoJet->Pt(), reco.pt_corrected[j], reco.ptLead[j],
                     tempRecoJet->Eta(), tempRecoJet->Phi(), tempRecoJet->Rapidity(),
                     reco.neutral_fraction[j], reco.trigger_match_JP2[j],
                     reco.trigger_match_JP1[j], reco.trigger_match_JP0[j],
                     reco.n_constituents[j],
                     reco.eventid, reco.weight, reco.mult, reco.trigger_match_HT2[j]);
            rj.det_eta = reco.det_eta[j];
            rj.area = reco.jet_area[j];
            rj.bg_density = reco.bg_density[j];
            rj.jp_patch_adc = reco.jp_patch_adc[j];
            rj.ht_adc_max = reco.ht_adc_max[j];
            recoJets.push_back(rj);
         }
      }

      // Match them
      vector<MatchedJetPair> MatchedJets;
      MatchedJets = MatchJetsEtaPhi(mcJets, recoJets, RCut);

      for (unsigned int j = 0; j < MatchedJets.size(); j++) {
         outMcJet = MatchedJets[j].first;
         outRecoJet = MatchedJets[j].second;
         deltaR = outMcJet.deltaR(outRecoJet);
         const bool mcReal   = outMcJet.pt   > -500;
         const bool recoReal = outRecoJet.pt > -500;
         if (!mcReal && recoReal) {
            outMcJet.weight = outRecoJet.weight;
         }

         MatchedTree->Fill();
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

   // Events whose MC pass produced no truth jet never enter the loop above, but
   // their reco jets are real fakes of the detector-level spectrum. Write them.
   {
      long nRecoOnly = 0;
      for (Long64_t iReco = 0; iReco < RecoChain->GetEntries(); ++iReco) {
         if (visitedReco.count((int)iReco)) continue;
         RecoChain->GetEntry(iReco);
         out_runid = reco.runid;
         out_runid1 = reco.runid1;
         out_eventid = reco.eventid;
         isTriggerEvent = reco.isTriggerEvent;
         reco_vz = reco.vz; reco_vz_vpd = reco.vz_vpd; reco_vx = reco.vx; reco_vy = reco.vy;
         evt_jp_adc_max = reco.evt_jp_adc_max; evt_ht_adc_max = reco.evt_ht_adc_max; std::copy_n(reco.jp_thr, 3, jp_thr); std::copy_n(reco.ht_thr, 4, ht_thr);
         event_vz = reco.vz; // no MC entry: use the reco vertex (the vz-difference cut is then inert)
         pthat_evt = -1;
         evt_should_JP0 = evt_should_JP1 = evt_should_JP2 = false; evt_should_HT2 = false;
         for (int j = 0; j < reco.njets; ++j) { evt_should_JP0 |= reco.trigger_match_JP0[j]; evt_should_JP1 |= reco.trigger_match_JP1[j]; evt_should_JP2 |= reco.trigger_match_JP2[j]; }
         if (haveEvtShould) { evt_should_JP0 = reco_evt_should[0]; evt_should_JP1 = reco_evt_should[1]; evt_should_JP2 = reco_evt_should[2]; if (haveEvtShouldHT2) evt_should_HT2 = reco_evt_should_ht2; }
         vector<MyJet> mcJetsNone, recoJetsOnly;
         for (int j = 0; j < reco.njets; ++j) {
            TStarJetVectorJet *tempRecoJet = (TStarJetVectorJet *)reco.jets->At(j);
            if (abs(tempRecoJet->Eta()) > EtaCutReco) continue;
            MyJet rj(*tempRecoJet, tempRecoJet->Pt(), reco.pt_corrected[j], reco.ptLead[j], tempRecoJet->Eta(), tempRecoJet->Phi(), tempRecoJet->Rapidity(),
                     reco.neutral_fraction[j], reco.trigger_match_JP2[j], reco.trigger_match_JP1[j], reco.trigger_match_JP0[j], reco.n_constituents[j],
                     reco.eventid, reco.weight, reco.mult, reco.trigger_match_HT2[j]);
            rj.det_eta = reco.det_eta[j]; rj.area = reco.jet_area[j]; rj.bg_density = reco.bg_density[j];
            rj.jp_patch_adc = reco.jp_patch_adc[j]; rj.ht_adc_max = reco.ht_adc_max[j];
            recoJetsOnly.push_back(rj);
         }
         vector<MatchedJetPair> fakes = MatchJetsEtaPhi(mcJetsNone, recoJetsOnly, RCut);
         for (unsigned int j = 0; j < fakes.size(); j++) {
            outMcJet = fakes[j].first; outRecoJet = fakes[j].second; deltaR = outMcJet.deltaR(outRecoJet);
            outMcJet.weight = outRecoJet.weight;
            MatchedTree->Fill(); hFake->Fill(outRecoJet.pt); FakeNumber++; hPtReco->Fill(outRecoJet.pt); ++nRecoOnly;
         }
      }
      cout << "reco-only events written as fakes: " << nRecoOnly << " jets" << endl;
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

   TH1D *hMissRate = (TH1D *)hMiss->Clone("hMissRate"); // "b" = binomial errors
   hMissRate->Divide(hMiss, hPtMc, 1, 1, "b");
   hMissRate->Write();
   TH1D *hFakeRate = (TH1D *)hFake->Clone("hFakeRate");
   hFakeRate->Divide(hFake, hPtReco, 1, 1, "b");
   hFakeRate->Write();

   fout->Close();

   return 0;
}

// Greedy nearest-neighbour matching in (eta, phi): each MC jet takes its closest
// free reco jet, every jet used at most once. Unmatched MC jets come back paired
// with a sentinel (a miss), leftover reco jets paired the other way (a fake).
vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets, const vector<MyJet> &RecoJets, const double &R)
{
   vector<MyJet> recoJetsCopy = RecoJets; // copy to avoid modifying the original
   vector<MatchedJetPair> matchedJets;
   MyJet dummy;

   for (const auto &mcJet : McJets) {
      bool isMatched = false;
      double minDeltaR = 10000;
      auto bestMatch = recoJetsCopy.end();

      for (auto rcit = recoJetsCopy.begin(); rcit != recoJetsCopy.end(); ++rcit) {
         MyJet recoJet = *rcit;
         double deltaR = mcJet.deltaR(recoJet);
         TH1D *hDeltaR = (TH1D *)gDirectory->Get("hDeltaR");
         if (hDeltaR) {
            hDeltaR->Fill(deltaR);
         }

         // Matching cone from the measured dR distribution; deliberately
         // independent of the jet R so the radii stay comparable.
         const double kMaxDR = 0.2;
         if (deltaR <= kMaxDR && deltaR < minDeltaR) {
            minDeltaR = deltaR;
            bestMatch = rcit;
            isMatched = true;
         }
      }
      if (isMatched) {
         matchedJets.push_back(make_pair(mcJet, *bestMatch));
         recoJetsCopy.erase(bestMatch);
      } else {
         matchedJets.push_back(make_pair(mcJet, dummy));
      }
   }
   for (const auto &recoJet : recoJetsCopy)
      matchedJets.push_back(make_pair(dummy, recoJet));

   return matchedJets;
}