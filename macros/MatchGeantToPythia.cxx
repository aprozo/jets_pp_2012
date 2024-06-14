//! HAS to be compiled,
//! root -l macros/PrepUnfolding.cxx+

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TProfile.h>
#include <TLine.h>

#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>

// #include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetVectorJet.h"

// #include "TStarJetPicoUtils.h"

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <exception>

// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBinByBin.h"
// #include "RootUnfoldBayes.h"
// #include "RootUnfoldSvd.h"
// #include "RootUnfoldTUnfold.h"

// #include "Binning.h"

using namespace std;

class ResultStruct
{
public:
  TStarJetVectorJet orig;
  double pT;
  vector<double> E;
  vector<double> phi;
  vector<double> eta;
  vector<double> pt;
  double weight;
  ResultStruct(TStarJetVectorJet orig, double pT, vector<double> E, vector<double> pt, vector<double> phi, vector<double> eta, double weight) : orig(orig),
                                                                                                                                                pT(pT),
                                                                                                                                                kT(kT),
                                                                                                                                                E(E),
                                                                                                                                                pt(pt),
                                                                                                                                                phi(phi),
                                                                                                                                                eta(eta),
                                                                                                                                                weight(weight){};

  ClassDef(ResultStruct, 1)
};

typedef pair<ResultStruct, ResultStruct> MatchedResultStruct;

int SortedIntoSamples_RooUnfold(int RADIUS = 4)
{

  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetHistLineWidth(2);

  TLegend *leg = 0;

  float RCut = (float)RADIUS / 10;

  srand(time(0));
  int MCRandom = 0;
  double RandomA = 0;
  double RandomB = 0;

  // alternative sorting method, gaussian distribution.
  // std::random_device rd;
  // std::mt19937 gen(rd());

  //! For pp Data
  TString PpLevelFile = "/gpfs01/star/pwg/robotmon/ppRun12_analysis_code/out/Geant_R04_new.root";
  TString McLevelFile = "/gpfs01/star/pwg/robotmon/ppRun12_analysis_code/out/Pythia_R04_new.root";

  float EtaCut = 1.0 - RCut;

  bool addpTCut = true;

  // Output
  // ------
  TString OutFileName = "/gpfs01/star/pwg/prozorov/jets_pp_2012/SortedSample_R04_matching_new.root";

  //! Set up Geant Chain
  TFile *Ppf = new TFile(PpLevelFile);
  TTree *PpChain = (TTree *)Ppf->Get("ResultTree");
  PpChain->BuildIndex("runid", "eventid");
  cout << "size is" << PpChain->GetEntries() << endl;

  TClonesArray *PpJets = new TClonesArray("TStarJetVectorJet");
  PpChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
  PpChain->SetBranchAddress("Jets", &PpJets);
  cout << "size is 2" << PpChain->GetEntries() << endl;

  vector<double> *pT = new vector<double>;
  vector<double> *E = new vector<double>;
  vector<vector<double>> *ppE_Vec = new vector<vector<double>>;
  vector<vector<double>> *pppt_Vec = new vector<vector<double>>;
  vector<vector<double>> *ppphi_Vec = new vector<vector<double>>;
  vector<vector<double>> *ppeta_Vec = new vector<vector<double>>;

  PpChain->SetBranchAddress("pT", &pT);
  PpChain->SetBranchAddress("ConstE", &ppE_Vec);
  PpChain->SetBranchAddress("Constphi", &ppphi_Vec);
  PpChain->SetBranchAddress("Constpt", &pppt_Vec);
  PpChain->SetBranchAddress("Consteta", &ppeta_Vec);

  int MatchNumber = 0;
  int FakeNumber = 0;
  int MissNumber = 0;
  int MissEventNumber = 0;
  int MatchedGeantEventNumber = 0;
  int TotalGeantEventNumber = 0;
  int FakeEventNumber = 0;

  int ppeventid;
  int pprunid;
  double ppweight;
  int ppnjets = 0;
  PpChain->SetBranchAddress("eventid", &ppeventid);
  PpChain->SetBranchAddress("runid", &pprunid);
  PpChain->SetBranchAddress("weight", &ppweight);
  PpChain->SetBranchAddress("njets", &ppnjets);

  //! Set up Pythia Chain
  TFile *Mcf = new TFile(McLevelFile);
  TTree *McChain = (TTree *)Mcf->Get("ResultTree");
  McChain->BuildIndex("runid", "eventid");
  cout << "size is" << McChain->GetEntries() << endl;

  TClonesArray *McJets = new TClonesArray("TStarJetVectorJet");
  McChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
  McChain->SetBranchAddress("Jets", &McJets);

  vector<double> *mcpT = new vector<double>;
  vector<double> *mcE = new vector<double>;

  vector<vector<double>> *mcE_Vec = new vector<vector<double>>;
  vector<vector<double>> *mcpt_Vec = new vector<vector<double>>;
  vector<vector<double>> *mcphi_Vec = new vector<vector<double>>;
  vector<vector<double>> *mceta_Vec = new vector<vector<double>>;

  McChain->SetBranchAddress("pT", &mcpT);

  McChain->SetBranchAddress("E", &mcE);

  McChain->SetBranchAddress("ConstE", &mcE_Vec);
  McChain->SetBranchAddress("Constphi", &mcphi_Vec);

  McChain->SetBranchAddress("Constpt", &mcpt_Vec);
  McChain->SetBranchAddress("Consteta", &mceta_Vec);

  int mceventid;
  int mcrunid;
  double mcweight;
  int mcnjets = 0;
  double mctotalpT;
  McChain->SetBranchAddress("totalpT", &mctotalpT);
  McChain->SetBranchAddress("eventid", &mceventid);
  McChain->SetBranchAddress("runid", &mcrunid);
  McChain->SetBranchAddress("weight", &mcweight);
  McChain->SetBranchAddress("njets", &mcnjets);

  // new from Isaac
  const int NUMBEROFPT = 13;
  // const char *PTBINS[NUMBEROFPT]={"2_3","3_4","4_5","5_7","7_9","9_11","11_15","15_20","20_25","25_35","35_-1"};
  const static float XSEC[NUMBEROFPT] = {0.00230158, 0.000342755, 0.0000457002, 9.0012, .00000972535, 1.46253, 0.000000469889, 0.354566, 0.0000000269202, 0.00000000143453, 0.151622, 0.0249062, 0.00584527};
  const static float NUMBEROFEVENT[NUMBEROFPT] = {3000000, 3000000, 3000000, 3000000, 2000000, 3000000, 2000000, 3000000, 1000000, 1000000, 3000000, 3000000, 3000000};
  const static float MAXPT[NUMBEROFPT] = {30, 40, 50, 6, 70, 8, 90, 10, 110, 2000, 14, 18, 22};
  const static vector<string> vptbins = {"1115_", "1520_", "2025_", "23_", "2535_", "34_", "3545_", "45_", "4555_", "55999_", "57_", "79_", "911_"};
  // Define Bounds on Delta R

  // R=0.4 bins

  // int nBins=12;
  // double bounds[13] = {exp(-3.7),exp(-3),exp(-2.5),exp(-2.2),exp(-2),exp(-1.8),exp(-1.6),exp(-1.4),exp(-1.2),exp(-1),exp(-0.8),exp(-0.6),exp(-0.2)}; //cool NEW non-log bounds that I can turn into them
  // int nBinsPtWeight=12;
  // double boundsPtWeight[13] = {20*exp(-3.7),20*exp(-3),20*exp(-2.5),20*exp(-2.2),20*exp(-2),20*exp(-1.8),20*exp(-1.6),20*exp(-1.4),20*exp(-1.2),20*exp(-1),20*exp(-0.8),20*exp(-0.6),20*exp(-0.2)};

  // Set Minimum Constituent pT
  double ConMinPt = 0.2;

  // Define Bounds on Jet pT
  // Truth and Measured Jet Bounds are the same for now

  int MeasJetBins = 10;
  double MeasJetBounds[11] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60};

  int TruthJetBins = 10;
  double TruthJetBounds[11] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60};

  // Output and histograms
  TFile *fout = new TFile(OutFileName, "RECREATE");
  TH1::SetDefaultSumw2(true);

  // TH2D* CheckWeight15to20 = new TH2D("GeantPt","",nBins,bounds

  TH1D *DoubleMatchedResolution = new TH1D("DoubleMatchedResolution", "", 20, 0, 2);
  TH1D *SingleMatchedResolution = new TH1D("SingleMatchedResolution", "", 20, 0, 2);

  TH1D *MissRate_EventNum = new TH1D("MissRate_EventNum", "", 7500, 0, 7500000);
  TH1D *TotalRate_True_EventNum = new TH1D("TotalRate_True_EventNum", "", 7500, 0, 7500000);

  TH1D *FakeRate_EventNum = new TH1D("FakeRate_EventNum", "", 7500, 0, 7500000);
  TH1D *TotalRate_Meas_EventNum = new TH1D("TotalRate_Meas_EventNum", "", 7500, 0, 7500000);

  TH1D *JetPt_True = new TH1D("JetPt_True", "", MeasJetBins, MeasJetBounds);
  TH1D *JetPt_Meas = new TH1D("JetPt_Meas", "", MeasJetBins, MeasJetBounds);
  TH1D *JetPt_TriggerMissed = new TH1D("JetPt_TriggerMissed", "", MeasJetBins, MeasJetBounds);
  TH1D *JetPt_JetMissed = new TH1D("JetPt_JetMissed", "", MeasJetBins, MeasJetBounds);
  TH1D *JetPt_Fake_Event = new TH1D("JetPt_Fake_Event", "", MeasJetBins, MeasJetBounds);
  TH1D *JetPt_Fake_Jet = new TH1D("JetPt_Fake_Jet", "", MeasJetBins, MeasJetBounds);
  TH2D *JetEtaPhi_True = new TH2D("JetEtaPhi_True", "", 100, -1, 1, 100, 0, 6.3);
  TH2D *JetEtaPhi_Meas = new TH2D("JetEtaPhi_Meas", "", 100, -1, 1, 100, 0, 6.3);
  TH2D *JetEtaPhi_TriggerMissed = new TH2D("JetEtaPhi_TriggerMissed", "", 100, -1, 1, 100, 0, 6.3);
  TH2D *JetEtaPhi_JetMissed = new TH2D("JetEtaPhi_JetMissed", "", 100, -1, 1, 100, 0, 6.3);
  TH2D *JetEtaPhi_Fake_Event = new TH2D("JetEtaPhi_Fake_Event", "", 100, -1, 1, 100, 0, 6.3);
  TH2D *JetEtaPhi_Fake_Jet = new TH2D("JetEtaPhi_Fake_Jet", "", 100, -1, 1, 100, 0, 6.3);

  TH2D *JetPt_Match = new TH2D("JetPt_Match", "", MeasJetBins, MeasJetBounds, MeasJetBins, MeasJetBounds);

  int constituentbins = 11;
  double constituentbounds[12] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};
  TH2D *NConstituentsmc = new TH2D("NConstituentsmc", "", MeasJetBins, MeasJetBounds, constituentbins, constituentbounds);
  TH2D *NConstituentspp = new TH2D("NConstituentspp", "", MeasJetBins, MeasJetBounds, constituentbins, constituentbounds);

  TH2D *NMatched = new TH2D("NMatched", "", MeasJetBins, MeasJetBounds, constituentbins, constituentbounds);
  TH2D *NMissed = new TH2D("NMissed", "", MeasJetBins, MeasJetBounds, constituentbins, constituentbounds);
  TH2D *NFake = new TH2D("NFake", "", MeasJetBins, MeasJetBounds, constituentbins, constituentbounds);

  TH2D *ConstituentpTmc = new TH2D("ConstituentpTmc", "", 11, 5, 60, 100, 0, 30);
  TH2D *ConstituentpTpp = new TH2D("ConstituentpTpp", "", 11, 5, 60, 100, 0, 30);

  int neutralenergybins = 13;
  double neutralenergybounds[14] = {-0.5, -0.1, -0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.05};
  TH2D *NeutralEnergyFractionmc = new TH2D("NeutralEnergyFractionmc", "", MeasJetBins, MeasJetBounds, neutralenergybins, neutralenergybounds);
  TH2D *NeutralEnergyFractionpp = new TH2D("NeutralEnergyFractionpp", "", MeasJetBins, MeasJetBounds, neutralenergybins, neutralenergybounds);

  TH1D *ChargeSampleCounts = new TH1D("ChargeSampleCounts", "", 2, 1, 3);

  TH1D *MatchDistance = new TH1D("MatchDistance", "", 100, 0, 0.1);
  TH2D *Match_DeltaR = new TH2D("Match_DeltaRDifference", "", 100, 0, 0.1, 100, 0, 0.1);
  TH2D *ChargeMismatch = new TH2D("ChargeMismatch", "", 50, 0, 0.8, 50, 0, 0.8);
  TH2D *Match_DeltaDeltaR_ChargeMismatch = new TH2D("Match_DeltaRDeltaRDifference_ChargeMismatch", "", 200, 0, 0.8, 200, -0.4, 0.4);
  TH2D *Match_DeltaDeltaR = new TH2D("Match_DeltaRDeltaRDifference", "", 1000, 0, 0.8, 1000, -0.2, 0.2);
  TH2D *Match_DeltaDeltaR_LikeCharge = new TH2D("Match_DeltaRDeltaRDifference_LikeCharge", "", 1000, 0, 0.8, 1000, -0.2, 0.2);
  TH2D *Match_DeltaDeltaR_OppCharge = new TH2D("Match_DeltaRDeltaRDifference_OppCharge", "", 1000, 0, 0.8, 1000, -0.2, 0.2);
  TH2D *Match_DeltaPhi = new TH2D("Match_DeltaPhiDifference", "", 500, 0, 2 * 3.14, 500, 0, 2 * 3.14);
  TH2D *Match_DeltaEta = new TH2D("Match_DeltaEtaDifference", "", 500, -1, 1, 500, -1, 1);

  // SetUp Output Trees

  TTree *MatchTree = new TTree("MatchTree", "Correlations with weights and sample assignment");

  double EventWeightmc;
  double JetPtmc;
  double DeltaRmc;
  double EventWeightpp;
  double JetPtpp;
  double DeltaRpp;
  // Charge Sample 1 = like, 2 = unlike
  int ChargeSampleMatch;
  int SampleMatch;

  MatchTree->Branch("EventWeightmc", &EventWeightmc);
  MatchTree->Branch("JetPtmc", &JetPtmc);
  MatchTree->Branch("DeltaRmc", &DeltaRmc);

  MatchTree->Branch("JetPtpp", &JetPtpp);
  MatchTree->Branch("DeltaRpp", &DeltaRpp);

  MatchTree->Branch("SampleMatch", &SampleMatch);
  MatchTree->Branch("ChargeSampleMatch", &ChargeSampleMatch);

  TTree *FakeTree = new TTree("FakeTree", "Correlations with weights and sample assignment");

  double EventWeightFake;

  double JetPtFake;
  double DeltaRFake;
  int ChargeSampleFake;
  int SampleFake;
  int TriggerFake;

  FakeTree->Branch("EventWeightFake", &EventWeightFake);
  FakeTree->Branch("JetPtFake", &JetPtFake);
  FakeTree->Branch("DeltaRFake", &DeltaRFake);

  FakeTree->Branch("ChargeSampleFake", &ChargeSampleFake);
  FakeTree->Branch("SampleFake", &SampleFake);
  FakeTree->Branch("TriggerFake", &TriggerFake);

  TTree *MissTree = new TTree("MissTree", "Correlations with weights and sample assignment");

  double EventWeightMiss;
  double JetPtMiss;
  double DeltaRMiss;
  int SampleMiss;
  int ChargeSampleMiss;
  int TriggerMissed;

  MissTree->Branch("EventWeightMiss", &EventWeightMiss);

  MissTree->Branch("JetPtMiss", &JetPtMiss);
  MissTree->Branch("DeltaRMiss", &DeltaRMiss);

  MissTree->Branch("SampleMiss", &SampleMiss);
  MissTree->Branch("ChargeSampleMiss", &ChargeSampleMiss);
  MissTree->Branch("TriggerMissed", &TriggerMissed);

  // checks on EEC Resolution

  int ResolutionBins = 10;
  double ResolutionBounds[11];

  for (int k = 0; k < ResolutionBins + 1; ++k)
  {

    ResolutionBounds[k] = k * 0.5 / ((double)ResolutionBins);
  }

  int N = McChain->GetEntries();
  int M = PpChain->GetEntries();

  vector<int> GeantEntries;
  for (Long64_t mcEvi = 0; mcEvi < (N); ++mcEvi)
  {

    // mcEvi=(mcEvil % N);

    if (!(mcEvi % 10000))
      cout << "Working on " << mcEvi << " / " << N << endl
           << M << endl;
    McChain->GetEntry(mcEvi);

    if (McJets->GetEntries() != mcnjets)
    {
      cerr << "McJets->GetEntries() != mcnjets" << endl;
      return -1;
    }

    //! Fill results in vectors for easier manipulation
    //! Also check whether there's something true in the acceptance
    bool TruthInAcceptance = false;
    vector<ResultStruct> mcresult;
    for (int j = 0; j < mcnjets; ++j)
    {

      // get the values for the current jet
      TStarJetVectorJet *mcjet = (TStarJetVectorJet *)McJets->At(j);
      double mcPT = mcpT->at(j);
      double mcE1 = mcE->at(j);
      double mcKT = mckT->at(j);
      double mcTf = mctf->at(j);
      double mcMu = mcmu->at(j);
      vector<double> mcE_Jet = mcE_Vec->at(j);
      vector<double> mcpt_Jet = mcpt_Vec->at(j);
      vector<double> mcphi_Jet = mcphi_Vec->at(j);
      vector<double> mceta_Jet = mceta_Vec->at(j);

      //! Fill in jet into pythia result

      if (fabs(mcjet->Eta()) < EtaCut)
      {

        mcresult.push_back(ResultStruct(*mcjet, mcZg, mcRg, mcPT, mcKT, mcTf, mcMu, mcE_Jet, mcpt_Jet, mcphi_Jet, mceta_Jet, mcweight));
        TruthInAcceptance = true;
      }
    }

    if (!TruthInAcceptance)
    {
      //! Skip this event, but don't count it as a loss
      continue;
    }

    // Still in MC level loop, get matching Geant Event
    Long64_t ppevi = -1;
    ppevi = PpChain->GetEntryNumberWithIndex(mcrunid, mceventid);
    if (ppevi >= 0)
    {
      GeantEntries.push_back(ppevi);
      MatchedGeantEventNumber++;
    }
    // bug in new embedding
    if (mctotalpT > 23.11003 && mctotalpT < 23.11004)
    {
      continue;
    } // 11-15
    if (mctotalpT > 33.749385 && mctotalpT < 33.749405)
    {
      continue;
    } // 15-20
    if (mctotalpT > 47.09071 && mctotalpT < 47.09072)
    {
      continue;
    } // 20-25
    if (mctotalpT > 9.62226 && mctotalpT < 9.62227)
    {
      continue;
    } // 2-3
    if (mctotalpT > 46.63831 && mctotalpT < 46.63832)
    {
      continue;
    } // 25-35
    if (mctotalpT > 6.90831 && mctotalpT < 6.90832)
    {
      continue;
    } // 3-4
    if (mctotalpT > 82.68752 && mctotalpT < 82.68753)
    {
      continue;
    } // 35-45
    if (mctotalpT > 100.25616 && mctotalpT < 100.25617)
    {
      continue;
    } // 45-55
    if (mctotalpT > 75.10883 && mctotalpT < 75.10884)
    {
      continue;
    } // 55-999
    if (mctotalpT > 3.75004 && mctotalpT < 3.75005)
    {
      continue;
    } // 5-7
    if (mctotalpT > 6.47623 && mctotalpT < 6.47624)
    {
      continue;
    } // 7-9
    if (mctotalpT > 4.22790 && mctotalpT < 4.22791)
    {
      continue;
    } // 9-11
    // Decide which sample event is sorted into
    MCRandom = rand();
    // cout << MCRandom << endl;

    // fill in regardless of match/miss to get reference and check which events are candidates for high weight and throw them out
    int isBad = 0;
    int weightBin = -1;
    int old = -1;
    for (int i = 0; i < vptbins.size(); ++i)
    {
      if (mcweight == XSEC[i] / NUMBEROFEVENT[i])
      {
        old = 0;
        weightBin = i;
      }
    }
    for (vector<ResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end(); ++mcit)
    {
      double Jetpt = mcit->orig.perp();
      double Jetpt2 = pow(Jetpt, 2);
      if (old == 0)
      {
        if (Jetpt > MAXPT[weightBin])
        {
          isBad = 1;
        }
      }
    }
    if (isBad == 1)
    {
      continue;
    }

    // matching geant event not found, fill in pythia event as a miss
    if (ppevi < 0)
    {

      for (vector<ResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end(); ++mcit)
      {

        // cout << mcit->weight;
        double Jetpt = mcit->orig.perp();
        double Jeteta = mcit->orig.eta();
        double Jetphi = mcit->orig.phi();
        // Fill in jet quantities like mass
        // roll random number at jet level
        // MCRandom=rand();

        JetPt_True->Fill(Jetpt, mcit->weight);
        JetPt_TriggerMissed->Fill(Jetpt, mcit->weight);
        if (Jetpt > 10)
        {
          JetEtaPhi_True->Fill(Jeteta, Jetphi, mcit->weight);
          JetEtaPhi_TriggerMissed->Fill(Jeteta, Jetphi, mcit->weight);
        }

        NConstituentsmc->Fill(Jetpt, mcit->E.size(), mcit->weight);
        // Fill in substructure variables such as EEC for entire missed jet
        for (unsigned i = 0; i < mcit->E.size(); ++i)
        {
          // Fill in missed constituents
          ConstituentpTmc->Fill(Jetpt, mcit->pt[i], mcit->weight);
          if (mcit->pt[i] < ConMinPt)
          {
            continue;
          }
          for (unsigned j = i + 1; j < mcit->E.size(); ++j)
          {
            if (mcit->pt[j] < ConMinPt)
            {
              continue;
            }
            double deltaPhiMiss = (double)abs(mcit->phi[i] - mcit->phi[j]);
            // calculate distance between the pair
            if (deltaPhiMiss > M_PI)
            {
              deltaPhiMiss = deltaPhiMiss - (2 * M_PI);
            }
            double deltaEtaMiss = (double)mcit->eta[i] - mcit->eta[j];
            double deltaRMiss = sqrt(pow(deltaPhiMiss, 2) + pow(deltaEtaMiss, 2));

            double Jetpt2 = pow(Jetpt, 2);

            // check weighting, will remove
            // EECWeightMiss=1;
            // Jetpt2=1;

            // check that you're filling in within bounds
            if (Jetpt > TruthJetBounds[0] && Jetpt < TruthJetBounds[TruthJetBins])
            {

              EventWeightMiss = (mcit->weight);
              JetPtMiss = (mcit->pT);
              ZgMiss = (mcit->zg);
              kTMiss = (mcit->kT);
              RgMiss = (mcit->rg);
              MuMiss = (mcit->mu);
              DeltaRMiss = (log(deltaRMiss));

              // roll random number at correlation level
              // MCRandom=rand();

              // Sort into samples
              if ((MCRandom % 2) == 0)
              {
                ++RandomA;
                SampleMiss = (1);
              }
              else
              {
                ++RandomB;
                SampleMiss = (2);
              }

              // cout << SampleMiss << endl;

              TriggerMissed = 1;

              MissTree->Fill();
              MissEventNumber++;
            }
            // end of filling in
          }
        }
      }
      continue;
    }

    // Get Geant event if its found
    PpChain->GetEntry(ppevi);
    if (ppweight != mcweight)
    {
      cout << "hey hey hey we got a problem here " << ppweight << " " << mcweight << endl;
    }
    vector<ResultStruct> ppresult;

    for (int j = 0; j < ppnjets; ++j)
    {
      TStarJetVectorJet *ppjet = (TStarJetVectorJet *)PpJets->At(j);
      double ppZg = zg->at(j);
      double ppRg = rg->at(j);
      double ppPT = pT->at(j);
      double ppE1 = E->at(j);
      double ppKT = kT->at(j);
      double ppTf = tf->at(j);
      double ppMu = mu->at(j);
      vector<double> ppE_Jet = ppE_Vec->at(j);
      vector<double> pppt_Jet = pppt_Vec->at(j);
      vector<double> ppphi_Jet = ppphi_Vec->at(j);
      vector<double> ppeta_Jet = ppeta_Vec->at(j);

      //! Ok, record
      if (fabs(ppjet->Eta()) < EtaCut)
      {
        ppresult.push_back(ResultStruct(*ppjet, ppZg, ppRg, ppPT, ppKT, ppTf, ppMu, ppE_Jet, pppt_Jet, ppphi_Jet, ppeta_Jet, ppweight));
      }
    }

    for (vector<ResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end(); ++ppit)
    {
      double Jetpt = ppit->orig.perp();
      double Jetpt2 = pow(Jetpt, 2);

      if (Jetpt > MAXPT[weightBin])
      {
        isBad = 1;
      }
    }
    if (isBad == 1)
    {
      continue;
    }

    // match and Sort Pythia and Geant jets together
    int totalJetsmc = mcresult.size();
    int totalJetspp = ppresult.size();
    int matchedJets = 0;
    int missedJets = 0;
    int fakeJets = 0;

    vector<MatchedResultStruct> MatchedResult;

    for (vector<ResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end();)
    {
      bool matched = false;

      for (vector<ResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end();)
      {

        // Check that jets are matched, haven't been matched before and geant event meets requirements to be jet
        if (mcit->orig.DeltaR(ppit->orig) < RCut && matched == false && ppit->orig.perp() > 10 && mcit->orig.perp() > 5)
        {
          SingleMatchedResolution->Fill((double)(ppit->orig.perp()) / (mcit->orig.perp()), mcweight);
          MatchedResult.push_back(MatchedResultStruct(*mcit, *ppit));
          ppit = ppresult.erase(ppit);
          matched = true;
          break;
        }
        else
        {

          ++ppit;
        }
      }
      if (matched)
      {
        if (mcit->orig.perp() > 10)
        {
          matchedJets++;
        }
        // if you find a match for a pythia jet, sort it away and remove from list
        mcit = mcresult.erase(mcit);
      }
      else
      {
        // pythia jet has been missed

        // Fill in missed jet quantities

        double ptJet = mcit->orig.perp();
        double JetMass = mcit->orig.m();
        // Fill in missed jet quantities like mass

        if (ptJet < 5)
        {
          mcit = mcresult.erase(mcit);
          continue;
        }
        if (ptJet > 10)
        {
          missedJets++;
        }
        ++mcit;
        // roll random number at jet level
        // MCRandom=rand();
      }
    }
    fakeJets = ppresult.size();
    if ((double)missedJets / (double)totalJetsmc > 0.5)
    {
      // cout << "miss ratio " << (double)missedJets/(double)totalJetsmc << endl;
      // cout << (double)fakeJets/(double)totalJetspp << endl;
      //	continue;
    }
    if ((double)fakeJets == (double)totalJetspp && missedJets == totalJetsmc)
    {
    }
    TotalRate_Meas_EventNum->Fill(mcEvi, matchedJets);
    TotalRate_True_EventNum->Fill(mcEvi, matchedJets);
    MissRate_EventNum->Fill(mcEvi, missedJets);
    TotalRate_True_EventNum->Fill(mcEvi, missedJets);

    // Loop through leftover pp Jets and fill in as fakes

    for (vector<ResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end();)
    {

      // Fill in overall jet quantities

      double ptJet = mcit->orig.perp();
      double JetMass = mcit->orig.m();
      double etaJet = mcit->orig.eta();
      double phiJet = mcit->orig.phi();

      JetPt_True->Fill(ptJet, mcit->weight);
      JetPt_JetMissed->Fill(ptJet, mcit->weight);
      if (ptJet > 10)
      {
        JetEtaPhi_True->Fill(etaJet, phiJet, mcit->weight);
        JetEtaPhi_JetMissed->Fill(etaJet, phiJet, mcit->weight);
      }
      // Fill in missed jet quantities like mass

      // roll random number at jet level
      // MCRandom=rand();

      NConstituentsmc->Fill(ptJet, mcit->E.size(), mcit->weight);

      // Loop over constituents and fill in

      for (unsigned i = 0; i < mcit->E.size(); ++i)
      {
        ConstituentpTmc->Fill(ptJet, mcit->pt[i], mcit->weight);
        if (mcit->pt[i] < ConMinPt)
        {
          continue;
        }
        for (unsigned j = i + 1; j < mcit->E.size(); ++j)
        {
          if (mcit->pt[j] < ConMinPt)
          {
            continue;
          }
          // calculate distance
          double deltaPhiMiss = (double)abs(mcit->phi[i] - mcit->phi[j]);
          if (deltaPhiMiss > M_PI)
          {
            deltaPhiMiss = deltaPhiMiss - (2 * M_PI);
          }
          double deltaEtaMiss = (double)mcit->eta[i] - mcit->eta[j];
          double deltaRMiss = sqrt(pow(deltaPhiMiss, 2) + pow(deltaEtaMiss, 2));

          double ptJet2 = pow(ptJet, 2);

          if (ptJet > TruthJetBounds[0] && ptJet < TruthJetBounds[TruthJetBins])
          {

            EventWeightMiss = (mcit->weight);
            JetPtMiss = (mcit->pT);
            ZgMiss = (mcit->zg);
            kTMiss = (mcit->kT);
            RgMiss = (mcit->rg);
            MuMiss = (mcit->mu);
            DeltaRMiss = (log(deltaRMiss));

            // roll random number at correlation level
            //  MCRandom=rand();

            // Sort into samples
            if ((MCRandom % 2) == 0)
            {
              ++RandomA;
              SampleMiss = (1);
            }
            else
            {
              ++RandomB;
              SampleMiss = (2);
            }

            // double chargeCheckFake=ppit->charge[i]*ppit->charge[j];

            /*if(chargeCheckFake>0){

        ChargeSampleFake=1;

            }

            if(chargeCheckFake<0){

        ChargeSampleFake=2;
            }

            ChargeSampleCounts->Fill(ChargeSampleFake);*/

            TriggerMissed = 2;
            MissTree->Fill();
            MissNumber++;
          }
        }
      }
      mcit = mcresult.erase(mcit);
    }
    for (vector<ResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end();)
    {

      // Fill in overall jet quantities
      double ptJet = ppit->orig.perp();
      if (ptJet < 10 || ptJet > MeasJetBounds[MeasJetBins])
      {
        ppit = ppresult.erase(ppit);
        continue;
      }
      double etaJet = ppit->orig.eta();
      double phiJet = ppit->orig.phi();
      double chargedpT = 0;
      JetPt_Meas->Fill(ptJet, ppweight);
      JetPt_Fake_Jet->Fill(ptJet, ppweight);
      JetEtaPhi_Meas->Fill(etaJet, phiJet, ppweight);
      JetEtaPhi_Fake_Jet->Fill(etaJet, phiJet, ppweight);

      NConstituentspp->Fill(ptJet, ppit->E.size(), ppit->weight);
      double ptJet2 = pow(ptJet, 2);

      TotalRate_Meas_EventNum->Fill(mcEvi);
      FakeRate_EventNum->Fill(mcEvi);

      // Loop over constituents and fill in

      for (unsigned i = 0; i < ppit->E.size(); ++i)
      {
        if (ppweight != ppit->weight)
        {
          cout << "err" << endl;
        }
        ConstituentpTpp->Fill(ptJet, ppit->pt[i], ppweight);
        if (ppit->pt[i] < ConMinPt)
        {
          continue;
        }
        for (unsigned j = i + 1; j < ppit->E.size(); ++j)
        {
          if (ppit->pt[j] < ConMinPt)
          {
            continue;
          }

          double deltaPhiFake = (double)abs(ppit->phi[i] - ppit->phi[j]);
          if (deltaPhiFake > M_PI)
          {

            deltaPhiFake = deltaPhiFake - (2 * M_PI);
          }
          double deltaEtaFake = (double)ppit->eta[i] - ppit->eta[j];
          double deltaRFake = sqrt(pow(deltaPhiFake, 2) + pow(deltaEtaFake, 2));

          EventWeightFake = (ppit->weight);
          JetPtFake = (ppit->pT);
          ZgFake = (ppit->zg);
          kTFake = (ppit->kT);
          RgFake = (ppit->rg);
          MuFake = (ppit->mu);
          DeltaRFake = (log(deltaRFake));

          // Sort into samples
          if ((MCRandom % 2) == 0)
          {
            ++RandomA;
            SampleFake = (1);
          }
          else
          {
            ++RandomB;
            SampleFake = (2);
          }

          // don't think I should be doing this, change later

          TriggerFake = 2;
          FakeTree->Fill();
          FakeNumber++;
        }
      }
      // if(chargedpT/JetPtFake > 1){cout << "something up on fake level" << endl;}
      ppit = ppresult.erase(ppit);
    }
    ppresult.clear();
    mcresult.clear();

    // Fill Matches

    for (vector<MatchedResultStruct>::iterator res = MatchedResult.begin(); res != MatchedResult.end(); ++res)
    {

      ResultStruct mc = res->first;
      ResultStruct pp = res->second;

      // Fill in Matched Jets
      double ptJetmc = mc.orig.perp();
      double ptJetpp = pp.orig.perp();
      double etaJetmc = mc.orig.eta();
      double etaJetpp = pp.orig.eta();
      double phiJetmc = mc.orig.phi();
      double phiJetpp = pp.orig.phi();

      NConstituentsmc->Fill(ptJetmc, mc.E.size(), mc.weight);
      NConstituentspp->Fill(ptJetpp, pp.E.size(), pp.weight);

      JetPt_Meas->Fill(ptJetpp, pp.weight);
      JetPt_True->Fill(ptJetmc, mc.weight);
      if (pp.orig.perp() > 10)
      {
        JetEtaPhi_Meas->Fill(etaJetpp, phiJetpp, pp.weight);
        JetEtaPhi_True->Fill(etaJetmc, phiJetmc, mc.weight);
      }
      JetPt_Match->Fill(ptJetpp, ptJetmc, mc.weight);

      double ptJetmc2 = pow(ptJetmc, 2);
      double ptJetpp2 = pow(ptJetpp, 2);
      // define vectors to help sort matches
      vector<int> mcMatchIndex;
      vector<int> ppMatchIndex;
      vector<int> mcMissIndex;
      vector<int> ppFakeIndex;

      // now begin constituent matching procedure

      for (unsigned iCont = 0; iCont < mc.E.size(); ++iCont)
      {
        ConstituentpTmc->Fill(ptJetmc, mc.pt[iCont], mcweight);
        if (mc.pt[iCont] < ConMinPt)
        {
          continue;
        }

        int matched = 0;
        double percentMatch = 0.4;
        double matchIndex = -1;

        for (unsigned jCont = 0; jCont < pp.E.size(); ++jCont)
        {
          if (pp.pt[jCont] < ConMinPt)
          {
            continue;
          }
          double percentDifference = abs((mc.pt[iCont] - pp.pt[jCont]) / mc.pt[iCont]);
          double deltaPhi = (double)abs(pp.phi[jCont] - mc.phi[iCont]);

          if (deltaPhi > M_PI)
          {
            deltaPhi = deltaPhi - (2 * M_PI);
          }
          if (deltaPhi < -3.14)
          {
            cout << deltaPhi << endl;
          }
          double deltaEta = (double)pp.eta[jCont] - mc.eta[iCont];
          double deltaR = sqrt(pow(deltaPhi, 2) + pow(deltaEta, 2));

          // check if geant-level Constituent has already been matched
          // will Double count if this is not here
          int contains = 0;

          for (unsigned check = 0; check < ppMatchIndex.size(); ++check)
          {
            if (ppMatchIndex[check] == jCont)
            {
              contains = 1;
            }
          }

          if (contains == 1)
          {
            continue;
          }

          // Match if within radius
          if (deltaR < 0.01)
          {

            // check that this match is better than previous match, if no previous match set to be within 50% of jet momentum
            if (percentDifference < percentMatch)
            {
              // mark that at least one viable match candidate has been found and store its index
              matched = 1;
              percentMatch = percentDifference;
              matchIndex = jCont;
            }
          }
        }

        // end of pp loop, decide what happens to mc particle now
        // if matched, store it with its most recent match
        if (matched == 1)
        {
          mcMatchIndex.push_back(iCont);
          ppMatchIndex.push_back(matchIndex);
        }

        // if not matched, store it as a miss
        if (matched == 0)
        {
          mcMissIndex.push_back(iCont);
        }
      }
      /*double percentDifference=abs((mc.pt[iCont]-pp.pt[jCont])/mc.pt[iCont]);
      double deltaPhi= (double)pp.phi[jCont]-mc.phi[iCont];

      if(deltaPhi>M_PI){

        deltaPhi=deltaPhi-(2*M_PI);
      }
      double deltaEta=(double)pp.eta[jCont]-mc.eta[iCont];
      double deltaR=sqrt(pow(deltaPhi,2)+pow(deltaEta,2));*/
      // check if jConstituent has already been matched
      // WILL Double count if this is not here, showing that there can be multiple consituents within same delta R, what does this mean?

      for (unsigned kCont = 0; kCont < pp.E.size(); ++kCont)
      {
        if (pp.pt[kCont] < ConMinPt)
        {
          continue;
        }
        int matched = 0;
        for (unsigned check = 0; check < ppMatchIndex.size(); ++check)
        {
          if (ppMatchIndex[check] == kCont)
          {
            matched = 1;
          }
        }
        if (matched == 0)
        {
          ppFakeIndex.push_back(kCont);
        }
      }

      // cout << "Number of matched particles is: " << ppMatchIndex.size() << " and fakes: " << ppFakeIndex.size() << "with total of " << pp.E.size() << endl;
      // Now feed in all of the Fake Correlators

      for (unsigned i = 0; i < ppFakeIndex.size(); ++i)
      {
        // fakes with matches
        for (unsigned j = 0; j < ppMatchIndex.size(); ++j)
        {
          double deltaPhiFake = (double)abs(pp.phi[ppFakeIndex[i]] - pp.phi[ppMatchIndex[j]]);
          if (deltaPhiFake > M_PI)
          {
            deltaPhiFake = deltaPhiFake - (2 * M_PI);
          }
          double deltaEtaFake = (double)pp.eta[ppFakeIndex[i]] - pp.eta[ppMatchIndex[j]];
          double deltaRFake = sqrt(pow(deltaPhiFake, 2) + pow(deltaEtaFake, 2));

          if (ptJetpp < TruthJetBounds[TruthJetBins] && ptJetpp > 10)
          {

            EventWeightFake = (pp.weight);
            JetPtFake = (ptJetpp);
            ZgFake = (Zgpp);
            kTFake = (kTpp);
            RgFake = (Rgpp);
            MuFake = (Mupp);
            DeltaRFake = (log(deltaRFake));

            // Sort into samples
            if ((MCRandom % 2) == 0)
            {
              ++RandomA;
              SampleFake = (1);
            }
            else
            {
              ++RandomB;
              SampleFake = (2);
            }

            TriggerFake = -1;
            FakeTree->Fill();
            FakeNumber++;
          }
        }

        // fakes with fakes
        for (unsigned j = i + 1; j < ppFakeIndex.size(); ++j)
        {
          double deltaPhiFake = (double)abs(pp.phi[ppFakeIndex[i]] - pp.phi[ppFakeIndex[j]]);
          if (deltaPhiFake > M_PI)
          {

            deltaPhiFake = deltaPhiFake - (2 * M_PI);
          }

          double deltaEtaFake = (double)pp.eta[ppFakeIndex[i]] - pp.eta[ppFakeIndex[j]];

          double deltaRFake = sqrt(pow(deltaPhiFake, 2) + pow(deltaEtaFake, 2));

          if (ptJetpp < TruthJetBounds[TruthJetBins] && ptJetpp > 10)
          {
            EventWeightFake = (pp.weight);
            JetPtFake = (ptJetpp);
            ZgFake = (Zgpp);
            kTFake = (kTpp);
            RgFake = (Rgpp);
            MuFake = (Mupp);
            DeltaRFake = (log(deltaRFake));

            // Sort into samples
            if ((MCRandom % 2) == 0)
            {
              ++RandomA;
              SampleFake = (1);
            }
            else
            {
              ++RandomB;
              SampleFake = (2);
            }

            // don't think I should be doing this, change later

            TriggerFake = 0;
            FakeTree->Fill();
            FakeNumber++;
          }
        }
      }

      // Fill in Missed Correlators

      for (unsigned i = 0; i < mcMissIndex.size(); ++i)
      {
        // misses with matches
        for (unsigned j = 0; j < mcMatchIndex.size(); ++j)
        {
          double deltaPhiMiss = (double)abs(mc.phi[mcMissIndex[i]] - mc.phi[mcMatchIndex[j]]);

          if (deltaPhiMiss > M_PI)
          {

            deltaPhiMiss = deltaPhiMiss - (2 * M_PI);
          }
          double deltaEtaMiss = (double)mc.eta[mcMissIndex[i]] - mc.eta[mcMatchIndex[j]];

          double deltaRMiss = sqrt(pow(deltaPhiMiss, 2) + pow(deltaEtaMiss, 2));

          double JetPtMeas = pp.orig.perp();

          // check weighting, will remove
          // EECWeightMiss=1;
          // ptJetmc2=1;

          if (ptJetmc > TruthJetBounds[0] && ptJetmc < TruthJetBounds[TruthJetBins])
          {

            EventWeightMiss = (mc.weight);
            JetPtMiss = (mc.pT);
            ZgMiss = (mc.zg);
            kTMiss = (mc.kT);
            RgMiss = (mc.rg);
            MuMiss = (mc.mu);
            DeltaRMiss = (log(deltaRMiss));

            // roll random number at correlation level
            // MCRandom=rand();

            // Sort into samples
            if ((MCRandom % 2) == 0)
            {
              ++RandomA;
              SampleMiss = (1);
            }
            else
            {
              ++RandomB;
              SampleMiss = (2);
            }

            TriggerMissed = 0;
            MissTree->Fill();
            MissNumber++;
          }
        }
        // misses with misses

        for (unsigned j = i + 1; j < mcMissIndex.size(); ++j)
        {
          double deltaPhiMiss = (double)abs(mc.phi[mcMissIndex[i]] - mc.phi[mcMissIndex[j]]);

          if (deltaPhiMiss > M_PI)
          {

            deltaPhiMiss = deltaPhiMiss - (2 * M_PI);
          }

          double deltaEtaMiss = (double)mc.eta[mcMissIndex[i]] - mc.eta[mcMissIndex[j]];
          double_t deltaRMiss = sqrt(pow(deltaPhiMiss, 2) + pow(deltaEtaMiss, 2));

          double JetPtMeas = pp.orig.perp();

          if (ptJetmc > TruthJetBounds[0] && ptJetmc < TruthJetBounds[TruthJetBins])
          {

            EventWeightMiss = (mc.weight);
            JetPtMiss = (mc.pT);
            ZgMiss = (mc.zg);
            kTMiss = (mc.kT);
            RgMiss = (mc.rg);
            MuMiss = (mc.mu);
            DeltaRMiss = (log(deltaRMiss));

            // roll random number at correlation level
            // MCRandom=rand();

            // Sort into samples
            if ((MCRandom % 2) == 0)
            {
              ++RandomA;
              SampleMiss = (1);
            }
            else
            {
              ++RandomB;
              SampleMiss = (2);
            }

            TriggerMissed = 0;
            MissTree->Fill();
            MissNumber++;
          }
        }
      }

      // roll random number at correlation level
      // MCRandom=rand();

      // Sort into samples

      // Fill in mc Matches

      for (unsigned i = 0; i < mcMatchIndex.size(); ++i)
      {

        Match_DeltaPhi->Fill(mc.phi[mcMatchIndex[i]], pp.phi[ppMatchIndex[i]]);
        Match_DeltaEta->Fill(mc.eta[mcMatchIndex[i]], pp.eta[ppMatchIndex[i]]);

        double deltaPhiMatchmcvspp = (double)abs(mc.phi[mcMatchIndex[i]] - pp.phi[ppMatchIndex[i]]);
        if (deltaPhiMatchmcvspp > M_PI)
        {
          deltaPhiMatchmcvspp = deltaPhiMatchmcvspp - (2 * M_PI);
        }
        double deltaEtaMatchmcvspp = (double)mc.eta[mcMatchIndex[i]] - pp.eta[ppMatchIndex[i]];
        double deltaRMatchmcvspp = sqrt(pow(deltaPhiMatchmcvspp, 2) + pow(deltaEtaMatchmcvspp, 2));
        MatchDistance->Fill(deltaRMatchmcvspp);

        for (unsigned j = i + 1; j < mcMatchIndex.size(); ++j)
        {

          // mcEEC
          double deltaPhiMatchmc = (double)abs(mc.phi[mcMatchIndex[i]] - mc.phi[mcMatchIndex[j]]);
          if (deltaPhiMatchmc > M_PI)
          {

            deltaPhiMatchmc = deltaPhiMatchmc - (2 * M_PI);
          }
          double deltaEtaMatchmc = (double)mc.eta[mcMatchIndex[i]] - mc.eta[mcMatchIndex[j]];
          double deltaRMatchmc = sqrt(pow(deltaPhiMatchmc, 2) + pow(deltaEtaMatchmc, 2));

          // ppEEC
          double deltaPhiMatchpp = (double)abs(pp.phi[ppMatchIndex[i]] - pp.phi[ppMatchIndex[j]]);
          if (deltaPhiMatchpp > M_PI)
          {

            deltaPhiMatchpp = deltaPhiMatchpp - (2 * M_PI);
          }
          double deltaEtaMatchpp = (double)pp.eta[ppMatchIndex[i]] - pp.eta[ppMatchIndex[j]];
          double deltaRMatchpp = sqrt(pow(deltaPhiMatchpp, 2) + pow(deltaEtaMatchpp, 2));

          if (ptJetmc > 5 && deltaRMatchpp < (RCut * 2) && ptJetpp < TruthJetBounds[TruthJetBins] && ptJetmc < TruthJetBounds[TruthJetBins])
          {

            EventWeightmc = (mc.weight);
            EventWeightpp = (pp.weight);
            JetPtpp = (ptJetpp);
            JetPtmc = (ptJetmc);

            DeltaRmc = (log(deltaRMatchmc));
            DeltaRpp = (log(deltaRMatchpp));
            Zgmc = (mc.zg);
            Zgpp = (pp.zg);
            Rgmc = (mc.rg);
            Rgpp = (pp.rg);
            Mumc = (mc.mu);
            Mupp = (pp.mu);
            kTmc = (mc.kT);
            kTpp = (pp.kT);

            // Sort into samples
            if ((MCRandom % 2) == 0)
            {
              ++RandomA;
              SampleMatch = (1);
            }
            else
            {
              ++RandomB;
              SampleMatch = (2);
            }

            // Sort into charge sample

            // Sort into charge sample

            /*double chargeCheckMatch =mc.charge[mcMatchIndex[i]]*mc.charge[mcMatchIndex[j]];

            if(chargeCheckMatch>0){

        ChargeSampleMatch=1;
        Match_DeltaDeltaR_LikeCharge->Fill(deltaRMatchmc,(deltaRMatchmc-deltaRMatchpp));

            }

            if(chargeCheckMatch<0){

        ChargeSampleMatch=2;
        Match_DeltaDeltaR_OppCharge->Fill(deltaRMatchmc,(deltaRMatchmc-deltaRMatchpp));
            }*/

            MatchTree->Fill();
            MatchNumber++;
          }
        }
      }

      NMissed->Fill(ptJetmc, mcMissIndex.size(), EventWeightmc);
      NFake->Fill(ptJetpp, ppFakeIndex.size(), EventWeightpp);
      NMatched->Fill(ptJetmc, mcMatchIndex.size(), EventWeightmc);

      if (mcMatchIndex.size() != ppMatchIndex.size())
      {
        cout << "matching error!" << endl;
      }

      ppFakeIndex.clear();
      mcMatchIndex.clear();
      ppMatchIndex.clear();
      mcMissIndex.clear();
    }
    MatchedResult.clear();
    // end of matched loop
  }
  // end of mc loop
  int FakeEvents = 0;

  sort(GeantEntries.begin(), GeantEntries.end());
  Long64_t GSize = GeantEntries.size();
  for (Long64_t ppCheck = 0; ppCheck < GSize - 1; ++ppCheck)
  {
    if (!(ppCheck % 10000))
    {
      cout << "Checking fake events: " << 100 * (double)ppCheck / (double)GSize << " percent done" << endl;
    }
    // cout << "Entry " << GeantEntries[ppCheck] << endl;
    if (GeantEntries[ppCheck] == GeantEntries[ppCheck + 1] - 1)
    {
      continue;
    }
    // if(trimming == 1){continue;}
    for (Long64_t ppEvi = GeantEntries[ppCheck] + 1; ppEvi < GeantEntries[ppCheck + 1]; ++ppEvi)
    {
      // cout << "In between " << GeantEntries[ppCheck+1]
      // cout << endl << "Got something here Event No: " << GeantEntries[ppEvi] << endl;

      // cout << "checking geant event " << ppEvi << " / " << M << " for fake jets " << endl;
      FakeEvents++;

      // reroll random sample sorting
      MCRandom = rand();

      // Fill in Geant Event as fake
      PpChain->GetEntry(ppEvi);
      // store jets
      vector<ResultStruct> ppresult;
      bool TruthInAcceptance = false;
      for (int j = 0; j < ppnjets; ++j)
      {
        TStarJetVectorJet *ppjet = (TStarJetVectorJet *)PpJets->At(j);
        vector<double> ppE_Jet = ppE_Vec->at(j);
        vector<double> pppt_Jet = pppt_Vec->at(j);
        vector<double> ppphi_Jet = ppphi_Vec->at(j);
        vector<double> ppeta_Jet = ppeta_Vec->at(j);
        double ppZg = zg->at(j);
        double ppRg = rg->at(j);
        double ppPT = pT->at(j);
        double ppE1 = E->at(j);
        double ppKT = kT->at(j);
        double ppTf = tf->at(j);
        double ppMu = mu->at(j);

        //! Ok, record
        if (fabs(ppjet->Eta()) < EtaCut)
        {
          ppresult.push_back(ResultStruct(*ppjet, ppZg, ppRg, ppPT, ppKT, ppTf, ppMu, ppE_Jet, pppt_Jet, ppphi_Jet, ppeta_Jet, ppweight));
          TruthInAcceptance = true;
        }
      }
      if (!TruthInAcceptance)
      {
        continue;
      }

      // Figure out what ptHard bin the event is from, and throw it away if a jet was produced more than twice the highest ptHard - will cut out events that have been weighted too highly.
      int isBad = 0;
      int weightBin = -1;
      for (unsigned i = 0; i < vptbins.size(); ++i)
      {
        if (ppweight == XSEC[i] / NUMBEROFEVENT[i])
        {
          weightBin = i;
        }
      }
      for (vector<ResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end(); ++ppit)
      {
        double Jetpt = ppit->orig.perp();
        double Jetpt2 = pow(Jetpt, 2);

        if (Jetpt > MAXPT[weightBin])
        {
          isBad = 1;
        }
      }
      if (isBad == 1)
      {
        cout << "bad fake event " << endl;
        continue;
      }
      TotalGeantEventNumber++;

      for (vector<ResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end();)
      {

        double Jetpt = ppit->orig.perp();
        // Only looking at jets in these ranges
        if (Jetpt > 10)
        {
          double Jeteta = ppit->orig.eta();
          double Jetphi = ppit->orig.phi();
          JetEtaPhi_Meas->Fill(Jeteta, Jetphi, ppweight);
          JetEtaPhi_Fake_Event->Fill(Jeteta, Jetphi, ppweight);

          if (ppweight != ppit->weight)
          {
            cout << "something's up " << endl;
          }
          JetPt_Meas->Fill(Jetpt, ppweight);
          JetPt_Fake_Event->Fill(Jetpt, ppweight);

          NConstituentspp->Fill(Jetpt, ppit->E.size(), ppit->weight);

          // cout << "good jet found " << endl;

          double JetMass = ppit->orig.m();

          // Fill in substructure variables such as EEC for entire fake jet
          for (unsigned i = 0; i < ppit->E.size(); ++i)
          {
            ConstituentpTpp->Fill(Jetpt, ppit->pt[i], ppweight);
            if (ppit->pt[i] < ConMinPt)
            {
              continue;
            }
            // Fill in fake constituents
            for (unsigned j = i + 1; j < ppit->E.size(); ++j)
            {
              if (ppit->pt[j] < ConMinPt)
              {
                continue;
              }
              double deltaPhiFake = (double)abs(ppit->phi[i] - ppit->phi[j]);
              if (deltaPhiFake > M_PI)
              {
                deltaPhiFake = deltaPhiFake - (2 * M_PI);
              }
              double deltaEtaFake = (double)ppit->eta[i] - ppit->eta[j];
              double deltaRFake = sqrt(pow(deltaPhiFake, 2) + pow(deltaEtaFake, 2));

              double Jetpt2 = pow(Jetpt, 2);

              // check that you're filling in within bounds
              if (Jetpt > 10 && Jetpt < TruthJetBounds[TruthJetBins])
              {

                EventWeightFake = (ppit->weight);
                JetPtFake = (Jetpt);
                DeltaRFake = (log(deltaRFake));
                ZgFake = (ppit->zg);
                kTFake = (ppit->kT);
                RgFake = (ppit->rg);
                MuFake = (ppit->mu);

                // sort into samples
                if ((MCRandom % 2) == 0)
                {
                  ++RandomA;
                  SampleFake = (1);
                }
                else
                {
                  ++RandomB;
                  SampleFake = (2);
                }

                TriggerFake = 1;
                FakeTree->Fill();

                FakeEventNumber++;
              }
              // end of filling in
            }
          }
        }
        ppit = ppresult.erase(ppit);
      }
      continue;
    }
  }

  // GeantPt->Draw("p E1");

  cout << "Sample 1 " << RandomA << endl;
  cout << "Sample 2 " << RandomB << endl;

  fout->Write();

  auto deltadeltaRPlot = new TCanvas;

  FakeRate_EventNum->Divide(TotalRate_Meas_EventNum);
  FakeRate_EventNum->Draw("p E1");

  auto deltadeltaRPlot2 = new TCanvas("MissC", "MissC");

  MissRate_EventNum->Divide(TotalRate_True_EventNum);
  MissRate_EventNum->Draw("p E1");

  auto ReweightedDist = new TCanvas("ReweightedDist", "ReweightedDist");

  // Match_DeltaDeltaR->Draw("colz");

  cout << "Miss Number is " << MissNumber << endl;
  cout << "From Missed Events is " << MissEventNumber << endl;
  cout << "Fake Number is " << FakeNumber << endl;
  cout << "From Fake Events is " << FakeEventNumber << endl;
  cout << "MatchNumber is " << MatchNumber << endl;

  return 0;
}
// end of function
