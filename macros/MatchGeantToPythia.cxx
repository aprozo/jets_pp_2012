//! pp Run12 zg, rg analysis 
//! Code to read in geant + pythia output trees and match them
//! Raghav Kunnawalkam Elayavalli & Kolja Kauder
//! contact - raghavke@wayne.edu
//! HAS to be compiled,
//! root -l macros/MatchGeantToPythia.cxx+


#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TProfile.h>
#include <TLine.h>

#include <TROOT.h>
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

#include <TStarJetVector.h>
#include <TStarJetVectorJet.h>

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <exception>

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"

#include "Binning.h"

using namespace std;

//! Load helper macro
#include "NewGeantWeightReject.hh"

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

//! ----------------------------------------------------
class RootResultStruct{
public:
  TStarJetVectorJet orig;
  double TwoSubJet_R0p1_Z;
  double TwoSubJet_R0p1_Theta;
  double TwoSubJet_R0p1_LeadpT;
  double TwoSubJet_R0p1_SubLeadpT;	  

  RootResultStruct ( TStarJetVectorJet orig,
			     double TwoSubJet_R0p1_Z, double TwoSubJet_R0p1_Theta,
			     double TwoSubJet_R0p1_LeadpT, double TwoSubJet_R0p1_SubLeadpT) :
    orig(orig),
    TwoSubJet_R0p1_Z(TwoSubJet_R0p1_Z),
    TwoSubJet_R0p1_Theta(TwoSubJet_R0p1_Theta),
    TwoSubJet_R0p1_LeadpT(TwoSubJet_R0p1_LeadpT),
    TwoSubJet_R0p1_SubLeadpT(TwoSubJet_R0p1_SubLeadpT)
  {};
  
  ClassDef(RootResultStruct,1)

};

typedef pair<RootResultStruct,RootResultStruct> MatchedRootResultStruct;

double findweight(double val = 10.0, TH1D *hpy8 = 0, TH1D *hhw7 = 0){
  // cout<<"val = "<<val<<endl;
  // hpy8->Print("base");
  // hhw7->Print("base");
  double val_r_py8 = hpy8->GetBinContent(hpy8->FindBin(val));
  double val_r_hw7 = hhw7->GetBinContent(hhw7->FindBin(val));
  return std::max(val_r_py8, val_r_hw7);
}

double findweight_max(double val = 10.0, TH1D *hpy8 = 0, TH1D *hhw7 = 0){
  double val_r_py8 = hpy8->GetBinContent(hpy8->FindBin(val));
  if(val_r_py8 > 0 && val_r_py8 < 1)
    val_r_py8 = 1./val_r_py8;
  double val_r_hw7 = hhw7->GetBinContent(hhw7->FindBin(val));
  if(val_r_hw7 > 0 && val_r_hw7 < 1)
    val_r_hw7 = 1./val_r_hw7;
  return std::max(val_r_py8, val_r_hw7);
}

double findweight_avg(double val = 10.0, TH1D *hpy8 = 0, TH1D *hhw7 = 0){
  double val_r_py8 = hpy8->GetBinContent(hpy8->FindBin(val));
  if(val_r_py8 > 0 && val_r_py8 < 1)
    val_r_py8 = 1./val_r_py8;
  double val_r_hw7 = hhw7->GetBinContent(hhw7->FindBin(val));
  if(val_r_hw7 > 0 && val_r_hw7 < 1)
    val_r_hw7 = 1./val_r_hw7;
  return (val_r_py8+val_r_hw7)/2;
}

double findweight_v2(double val = 10.0, TH1D *hweight = 0){
  return hweight->GetBinContent(hweight->FindBin(val));
}

int MatchGeantToPythia (int RADIUS = 4,
			int mode = 0
			) {
  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  gStyle->SetTitleBorderSize(0);	
  gStyle->SetHistLineWidth(2);
  TLegend* leg = 0;
  
  bool PrepClosure=true;
  bool addpTCut = true;

  TRandom rand;

  TString PpLevelFile;
  if(mode == 0)
    PpLevelFile = Form("Results/Geant_Run12_NoEff_NoBg_JP2_TSJ_wchpionMassSet_R0%d_Cleanpp12Pico_hadded_v2.root", RADIUS);
  else if(mode == 1)
    PpLevelFile = Form("Results/Geant_Run12_MIP_NoEff_NoBg_JP2_TSJ_wchpionMassSet_R0%d_Cleanpp12Pico_hadded_v2.root", RADIUS);
  else if(mode == 2)
    PpLevelFile = Form("Results/Geant_Run12_HC0p5_NoEff_NoBg_JP2_TSJ_wchpionMassSet_R0%d_Cleanpp12Pico_hadded_v2.root", RADIUS);
  else if(mode == 3)
    PpLevelFile = Form("Results/Geant_Run12_TowScaleUncer1_NoEff_NoBg_JP2_TSJ_wchpionMassSet_R0%d_Cleanpp12Pico_hadded_v2.root", RADIUS);
  else if(mode == 4)
    PpLevelFile = Form("Results/Geant_Run12_TrackEff0p96_NoEff_NoBg_JP2_TSJ_wchpionMassSet_R0%d_Cleanpp12Pico_hadded_v2.root", RADIUS);
    

  TString McLevelFile = Form("Results/Pythia6_Run12_TSJ_wparticleMassSet_R0%d_Cleanpp12Pico_hadded.root", RADIUS);

  
  TFile * fin = TFile::Open(Form("Results/PY6_pTZgRgZsjThetasj_shapeDiff_PY8HW7_R0%d.root", RADIUS));
  TH1D * hJetpT_PY6_Ratio_PY8 = (TH1D*)fin->Get("hJetpT_PY6_Ratio_PY8");
  TH1D * hJetpT_PY6_Ratio_HW7 = (TH1D*)fin->Get("hJetpT_PY6_Ratio_HW7");
  TH1D * hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[npubptbins];
  TH1D * hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[npubptbins];
  TH1D * hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[npubptbins];
  TH1D * hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[npubptbins];

  for(int ip = 0; ip < npubptbins; ++ip){
    hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[ip] = (TH1D*)fin->Get(Form("hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8_ptbin%d", ip));
    hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[ip] = (TH1D*)fin->Get(Form("hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8_ptbin%d", ip));
    hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[ip] = (TH1D*)fin->Get(Form("hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7_ptbin%d", ip));
    hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[ip] = (TH1D*)fin->Get(Form("hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7_ptbin%d", ip));
  }

  float MinJetPt = 5.0;
  bool UseMiss=true;
  bool UseFakes=true;

  bool RejectHiweights= true;
  
  float RCut = (float)RADIUS/10;
  
  float EtaCut = 1.0-RCut;

  //! Output
  TString OutFileName = "Results/";
  if(mode == 0){
    OutFileName = Form("Results/GeantToPythia_Match_ppRun12_200GeV_NoEff_NoBg_R0%d_v2.root", RADIUS);
  } else if (mode == 1)
    OutFileName = Form("Results/GeantToPythia_Match_ppRun12_200GeV_MIP_NoEff_NoBg_R0%d_v2.root", RADIUS);
  else if (mode == 2)
    OutFileName = Form("Results/GeantToPythia_Match_ppRun12_200GeV_HC0p5_NoEff_NoBg_R0%d_v2.root", RADIUS);
  else if (mode == 3)
    OutFileName = Form("Results/GeantToPythia_Match_ppRun12_200GeV_TowScaleUncer1_NoEff_NoBg_R0%d_v2.root", RADIUS);
  else if (mode == 4)
    OutFileName = Form("Results/GeantToPythia_Match_ppRun12_200GeV_TrackEff0p96_NoEff_NoBg_R0%d_v2.root", RADIUS);
  

  TFile* Mcf = new TFile( McLevelFile );
  TTree* McChain = (TTree*) Mcf->Get("ResultTree");
  McChain->BuildIndex("runid","eventid");

  TClonesArray* McJets = new TClonesArray("TStarJetVectorJet");
  McChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
  McChain->SetBranchAddress("Jets", &McJets);

  
  int mceventid;
  int mcrunid;
  double mcweight;     
  McChain->SetBranchAddress("eventid", &mceventid);
  McChain->SetBranchAddress("runid", &mcrunid);
  McChain->SetBranchAddress("weight",&mcweight );

  int mcnjets=0;
  McChain->SetBranchAddress("njets", &mcnjets );

  double mcTwoSubJet_R0p1_Z[1000];
  McChain->SetBranchAddress("TwoSubJet_R0p1_Z", mcTwoSubJet_R0p1_Z );

  double mcTwoSubJet_R0p1_Theta[1000];
  McChain->SetBranchAddress("TwoSubJet_R0p1_Theta", mcTwoSubJet_R0p1_Theta );

  double mcTwoSubJet_R0p1_LeadpT[1000];
  McChain->SetBranchAddress("TwoSubJet_R0p1_LeadpT", mcTwoSubJet_R0p1_LeadpT );

  double mcTwoSubJet_R0p1_SubLeadpT[1000];
  McChain->SetBranchAddress("TwoSubJet_R0p1_SubLeadpT", mcTwoSubJet_R0p1_SubLeadpT );
  
  TFile* Ppf = new TFile( PpLevelFile );
  TTree* PpChain = (TTree*) Ppf->Get("ResultTree");
  PpChain->BuildIndex("runid","eventid");

  TClonesArray* PpJets = new TClonesArray("TStarJetVectorJet");
  PpChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
  PpChain->SetBranchAddress("Jets", &PpJets);


  int ppeventid;
  int pprunid;
  double ppweight;     
  PpChain->SetBranchAddress("eventid", &ppeventid);
  PpChain->SetBranchAddress("runid", &pprunid);
  PpChain->SetBranchAddress("weight",&ppweight );

  int ppnjets=0;
  PpChain->SetBranchAddress("njets", &ppnjets );

  double ppTwoSubJet_R0p1_Z[1000];
  PpChain->SetBranchAddress("TwoSubJet_R0p1_Z", ppTwoSubJet_R0p1_Z );

  double ppTwoSubJet_R0p1_Theta[1000];
  PpChain->SetBranchAddress("TwoSubJet_R0p1_Theta", ppTwoSubJet_R0p1_Theta );

  double ppTwoSubJet_R0p1_LeadpT[1000];
  PpChain->SetBranchAddress("TwoSubJet_R0p1_LeadpT", ppTwoSubJet_R0p1_LeadpT );

  double ppTwoSubJet_R0p1_SubLeadpT[1000];
  PpChain->SetBranchAddress("TwoSubJet_R0p1_SubLeadpT", ppTwoSubJet_R0p1_SubLeadpT );

  
  //! Output and histograms
  TFile* fout = new TFile( OutFileName, "RECREATE");
  
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  
  //! Declare histograms for spectra
  TH1D * hRecoJetPt_QA = new TH1D("hRecoJetPt_QA","", 120, 0, 60);
  TH1D * hGenJetPt_QA = new TH1D("hGenJetPt_QA","", 120, 0, 60);  
  TH1D * hRecoJetPt = new TH1D("hRecoJetPt","", npubptbins, pubptbins);
  TH1D * hGenJetPt = new TH1D("hGenJetPt","", npubptbins, pubptbins);
  TH1D * hGenJetPt_v2 = new TH1D("hGenJetPt_v2","", npubptbins, pubptbins);
  TH1D * hGenMatchedJetPt = new TH1D("hGenMatchedJetPt","", npubptbins, pubptbins);
  TH1D * hGenUnMatchedJetPt = new TH1D("hGenUnMatchedJetPt","", npubptbins, pubptbins);
  TH1D * hRecoMatchedJetPt = new TH1D("hRecoMatchedJetPt","", npubptbins, pubptbins);
  TH1D * hGenMissJetPt = new TH1D("hGenMissJetPt","", npubptbins, pubptbins);
  TH1D * hRecoFakeJetPt = new TH1D("hRecoFakeJetPt","", npubptbins, pubptbins);

  TH1D * hJetFindingEfficiency = new TH1D("hJetFindingEfficiency", "", npubptbins, pubptbins);
  TH1D * hJetFakeRate = new TH1D("hJetFakeRate", "", npubptbins, pubptbins);
  
  TH1D * hRecoJetEta = new TH1D("hRecoJetEta","", 50, -1, 1);
  TH1D * hGenJetEta = new TH1D("hGenJetEta","", 50, -1, 1);
  TH1D * hRecoJetPhi = new TH1D("hRecoJetPhi","", 50, -3.15, 3.15);
  TH1D * hGenJetPhi = new TH1D("hGenJetPhi","", 50, -3.15, 3.15);
  TH2D * hMatrixJetPt = new TH2D("hMatrixJetPt","", npubptbins, pubptbins, npubptbins, pubptbins);

  //! histograms for JES/JER, zg, rg, M
  TH1D * hJER[npubptbins];  
  TH1D * hJER_v2[npubptbins];  
  TH1D * hJER_TwoSubJet_R0p1_Z[npubptbins];  
  TH1D * hJER_TwoSubJet_R0p1_Theta[npubptbins];  
  TH1D * hJER_v2_TwoSubJet_R0p1_Z[npubptbins];  
  TH1D * hJER_v2_TwoSubJet_R0p1_Theta[npubptbins];  

  TH1D * hGen_TwoSubJet_R0p1_Z[npubptbins];
  TH1D * hGen_TwoSubJet_R0p1_Theta[npubptbins];
  TH1D * hGen_TwoSubJet_R0p1_LeadpT[npubptbins];
  TH1D * hGen_TwoSubJet_R0p1_SubLeadpT[npubptbins];
  // TH1D * hGen_chFraction[npubptbins];  
  TH1D * hReco_TwoSubJet_R0p1_Z[npubptbins];
  TH1D * hReco_TwoSubJet_R0p1_Theta[npubptbins];
  TH1D * hReco_TwoSubJet_R0p1_LeadpT[npubptbins];
  TH1D * hReco_TwoSubJet_R0p1_SubLeadpT[npubptbins];
  // TH1D * hReco_chFraction[npubptbins];
  
  TH2D * hPtResponse = new TH2D("hPtResponse", "2D pT response matrix", nPtBinsMeas, ptminMeas, ptmaxMeas, nPtBinsTrue, ptminTrue, ptmaxTrue);
  TH2D * hPtResponse_transpose = new TH2D("hPtResponse_transpose", "2D pT response matrix", nPtBinsTrue, ptminTrue, ptmaxTrue, nPtBinsMeas, ptminMeas, ptmaxMeas);
  TH2D * hTwoSubJet_R0p1_ZResponse[npubptbins];
  TH2D * hTwoSubJet_R0p1_ThetaResponse[npubptbins];
  TH2D * hTwoSubJet_R0p1_LeadpTResponse[npubptbins];
  TH2D * hTwoSubJet_R0p1_SubLeadpTResponse[npubptbins];
  
  TH1D * hRecopTGroombypTDet[npubptbins];
  TH1D * hGenpTGroombypTDet[npubptbins];
  
  for(int i = 0; i<npubptbins; ++i){

    hTwoSubJet_R0p1_ZResponse[i] = new TH2D(Form("hTwoSubJet_R0p1_ZResponse_ptbin%d", i), Form("TwoSubJet_R0p1_Z response ptbin %d", i), nZgBinsTrue, zgminTrue, zgmaxTrue, nZgBinsTrue, zgminTrue, zgmaxTrue);
    hTwoSubJet_R0p1_ThetaResponse[i] = new TH2D(Form("hTwoSubJet_R0p1_ThetaResponse_ptbin%d", i), Form("TwoSubJet_R0p1_Theta response ptbin %d", i), nRgBinsTrue, rgminTrue, rgmaxTrue, nRgBinsTrue, rgminTrue, rgmaxTrue);
    hTwoSubJet_R0p1_LeadpTResponse[i] = new TH2D(Form("hTwoSubJet_R0p1_LeadpTResponse_ptbin%d", i), Form("TwoSubJet_R0p1_LeadpT response ptbin %d", i), 60, 0, 60, 60, 0, 60);
    hTwoSubJet_R0p1_SubLeadpTResponse[i] = new TH2D(Form("hTwoSubJet_R0p1_SubLeadpTResponse_ptbin%d", i), Form("TwoSubJet_R0p1_SubLeadpT response ptbin %d", i), 60, 0, 60, 60, 0, 60);
    
    hRecopTGroombypTDet[i]= new TH1D(Form("hRecopTGroombypTDet_ptbin%d", i), "", 100, 0, 1);
    hGenpTGroombypTDet[i]= new TH1D(Form("hGenpTGroombypTDet_ptbin%d", i), "", 100, 0, 1);
      
    hJER[i] = new TH1D(Form("hJER_ptbin_%d", i), "", 50, 0, 2);
    hJER_v2[i] = new TH1D(Form("hJER_v2_ptbin_%d", i), "", 50, 0, 2);
    hJER_TwoSubJet_R0p1_Z[i] = new TH1D(Form("hJER_TwoSubJet_R0p1_Z_ptbin_%d", i), "", 50, 0, 2);
    hJER_TwoSubJet_R0p1_Theta[i] = new TH1D(Form("hJER_TwoSubJet_R0p1_Theta_ptbin_%d", i), "", 50, 0, 2);
    hJER_v2_TwoSubJet_R0p1_Z[i] = new TH1D(Form("hJER_v2_TwoSubJet_R0p1_Z_ptbin_%d", i), "", 50, 0, 2);
    hJER_v2_TwoSubJet_R0p1_Theta[i] = new TH1D(Form("hJER_v2_TwoSubJet_R0p1_Theta_ptbin_%d", i), "", 50, 0, 2);
    
    hGen_TwoSubJet_R0p1_Z[i] = new TH1D(Form("hGen_TwoSubJet_R0p1_Z_ptbin_%d", i), "", nZgBinsTrue, zgminTrue, zgmaxTrue);
    hGen_TwoSubJet_R0p1_Theta[i] = new TH1D(Form("hGen_TwoSubJet_R0p1_Theta_ptbin_%d", i), "", nRgBinsTrue, rgminTrue, rgmaxTrue);
    hGen_TwoSubJet_R0p1_LeadpT[i] = new TH1D(Form("hGen_TwoSubJet_R0p1_LeadpT_ptbin_%d", i), "", 60, 0, 60);
    hGen_TwoSubJet_R0p1_SubLeadpT[i] = new TH1D(Form("hGen_TwoSubJet_R0p1_SubLeadpT_ptbin_%d", i), "", 60, 0, 60);
    // hGen_chFraction[i] = new TH1D(Form("hGen_chFraction_ptbin_%d", i), "", 50, 0.0, 1);    
    hReco_TwoSubJet_R0p1_Z[i] = new TH1D(Form("hReco_TwoSubJet_R0p1_Z_ptbin_%d", i), "", nZgBinsTrue, zgminTrue, zgmaxTrue);
    hReco_TwoSubJet_R0p1_Theta[i] = new TH1D(Form("hReco_TwoSubJet_R0p1_Theta_ptbin_%d", i), "", nRgBinsTrue, rgminTrue, rgmaxTrue);
    hReco_TwoSubJet_R0p1_LeadpT[i] = new TH1D(Form("hReco_TwoSubJet_R0p1_LeadpT_ptbin_%d", i), "", 60, 0, 60);
    hReco_TwoSubJet_R0p1_SubLeadpT[i] = new TH1D(Form("hReco_TwoSubJet_R0p1_SubLeadpT_ptbin_%d", i), "", 60, 0, 60);
    // hReco_chFraction[i] = new TH1D(Form("hReco_chFraction_ptbin_%d", i), "", 50, 0.0, 1);    

  }
  

  //! 2D unfolding:
  TH2D* hTrue= new TH2D ("hTrue", "Truth", nPtBinsTrue, ptminTrue, ptmaxTrue, nZgBinsTrue, zgminTrue, zgmaxTrue);
  TH2D* hMeas= new TH2D ("hMeas", "Measured", nPtBinsMeas, ptminMeas, ptmaxMeas, nZgBinsMeas, zgminMeas, zgmaxMeas);
  TH2D* hTrue_rg= new TH2D ("hTrue_rg", "Truth_rg", nPtBinsTrue, ptminTrue, ptmaxTrue, nRgBinsTrue, rgminTrue, rgmaxTrue);
  TH2D* hMeas_rg= new TH2D ("hMeas_rg", "Measured_rg", nPtBinsMeas, ptminMeas, ptmaxMeas, nRgBinsMeas, rgminMeas, rgmaxMeas);

  TH1D * hTrue_zg_1Dtest = new TH1D("hTrue_zg_1Dtest", "", nZgBinsTrue, zgminTrue, zgmaxTrue);
  TH1D * hTrue_rg_1Dtest = new TH1D("hTrue_rg_1Dtest", "", nRgBinsTrue, rgminTrue, rgmaxTrue);
  TH1D * hMeas_zg_1Dtest = new TH1D("hMeas_zg_1Dtest", "", nZgBinsMeas, zgminMeas, zgmaxMeas);
  TH1D * hMeas_rg_1Dtest = new TH1D("hMeas_rg_1Dtest", "", nRgBinsMeas, rgminMeas, rgmaxMeas);

  
  RooUnfoldResponse IncPtResponse (nPtBinsMeas, ptminMeas, ptmaxMeas, nPtBinsTrue, ptminTrue, ptmaxTrue);
  
  RooUnfoldResponse IncPtResponse_MCClosure (nPtBinsMeas, ptminMeas, ptmaxMeas, nPtBinsTrue, ptminTrue, ptmaxTrue);

  RooUnfoldResponse IncPtTwoSubJet_R0p1_ZResponse2D;
  IncPtTwoSubJet_R0p1_ZResponse2D.Setup (hMeas, hTrue );
  RooUnfoldResponse IncPtTwoSubJet_R0p1_ThetaResponse2D;
  IncPtTwoSubJet_R0p1_ThetaResponse2D.Setup (hMeas_rg, hTrue_rg );
  
  RooUnfoldResponse IncPtTwoSubJet_R0p1_ZResponse2D_MCClosure;
  IncPtTwoSubJet_R0p1_ZResponse2D_MCClosure.Setup (hMeas, hTrue);
  RooUnfoldResponse IncPtTwoSubJet_R0p1_ThetaResponse2D_MCClosure;
  IncPtTwoSubJet_R0p1_ThetaResponse2D_MCClosure.Setup (hMeas_rg, hTrue_rg);
  
  //! Various responses with truth shapes reweighted to represent a different prior
  RooUnfoldResponse IncBentPtResponse    ( nPtBinsMeas, ptminMeas, ptmaxMeas, nPtBinsTrue, ptminTrue, ptmaxTrue );
  RooUnfoldResponse IncPriorBentPtResponse( nPtBinsMeas, ptminMeas, ptmaxMeas, nPtBinsTrue, ptminTrue, ptmaxTrue );
  

  //! twosubjet bent stuff
  RooUnfoldResponse IncBentPtTwoSubJet_R0p1_ZResponse2D;
  IncBentPtTwoSubJet_R0p1_ZResponse2D.Setup (hMeas, hTrue);
  RooUnfoldResponse IncPtBentTwoSubJet_R0p1_ZResponse2D;
  IncPtBentTwoSubJet_R0p1_ZResponse2D.Setup (hMeas, hTrue);  
  RooUnfoldResponse IncBentPtBentTwoSubJet_R0p1_ZResponse2D;
  IncBentPtBentTwoSubJet_R0p1_ZResponse2D.Setup (hMeas, hTrue);

  RooUnfoldResponse IncBentPtTwoSubJet_R0p1_ThetaResponse2D;
  IncBentPtTwoSubJet_R0p1_ThetaResponse2D.Setup (hMeas_rg, hTrue_rg );
  RooUnfoldResponse IncPtBentTwoSubJet_R0p1_ThetaResponse2D;
  IncPtBentTwoSubJet_R0p1_ThetaResponse2D.Setup (hMeas_rg, hTrue_rg );  
  RooUnfoldResponse IncBentPtBentTwoSubJet_R0p1_ThetaResponse2D;
  IncBentPtBentTwoSubJet_R0p1_ThetaResponse2D.Setup (hMeas_rg, hTrue_rg );

  RooUnfoldResponse IncPriorBentPtTwoSubJet_R0p1_ZResponse2D;
  IncPriorBentPtTwoSubJet_R0p1_ZResponse2D.Setup (hMeas, hTrue );
  RooUnfoldResponse IncPtPriorBentTwoSubJet_R0p1_ZResponse2D;
  IncPtPriorBentTwoSubJet_R0p1_ZResponse2D.Setup (hMeas, hTrue );  
  RooUnfoldResponse IncPriorBentPtPriorBentTwoSubJet_R0p1_ZResponse2D;
  IncPriorBentPtPriorBentTwoSubJet_R0p1_ZResponse2D.Setup (hMeas, hTrue );

  RooUnfoldResponse IncPriorBentPtTwoSubJet_R0p1_ThetaResponse2D;
  IncPriorBentPtTwoSubJet_R0p1_ThetaResponse2D.Setup (hMeas_rg, hTrue_rg );
  RooUnfoldResponse IncPtPriorBentTwoSubJet_R0p1_ThetaResponse2D;
  IncPtPriorBentTwoSubJet_R0p1_ThetaResponse2D.Setup (hMeas_rg, hTrue_rg );  
  RooUnfoldResponse IncPriorBentPtPriorBentTwoSubJet_R0p1_ThetaResponse2D;
  IncPriorBentPtPriorBentTwoSubJet_R0p1_ThetaResponse2D.Setup (hMeas_rg, hTrue_rg );

  //! ptbin1 
  RooUnfoldResponse IncTwoSubJet_R0p1_ZResponse1D_ptbin1;
  IncTwoSubJet_R0p1_ZResponse1D_ptbin1.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin1;
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin1.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin1;
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin1.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncTwoSubJet_R0p1_ThetaResponse1D_ptbin1;
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin1.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin1;
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin1;
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);;
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin1;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin1.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin1;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin1.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin1;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin1.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin1;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin1.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin1;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin1.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin1;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin1.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);

  //! ptbin2 
  RooUnfoldResponse IncTwoSubJet_R0p1_ZResponse1D_ptbin2;
  IncTwoSubJet_R0p1_ZResponse1D_ptbin2.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin2;
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin2.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin2;
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin2.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncTwoSubJet_R0p1_ThetaResponse1D_ptbin2;
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin2.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin2;
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin2;
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin2;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin2.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin2;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin2.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin2;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin2.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin2;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin2.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin2;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin2.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin2;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin2.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);

  //! ptbin3 
  RooUnfoldResponse IncTwoSubJet_R0p1_ZResponse1D_ptbin3;
  IncTwoSubJet_R0p1_ZResponse1D_ptbin3.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin3;
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin3.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin3;
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin3.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncTwoSubJet_R0p1_ThetaResponse1D_ptbin3;
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin3.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin3;
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin3;
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin3;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin3.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin3;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin3.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin3;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin3.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin3;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin3.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin3;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin3.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin3;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin3.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);

  //! ptbin4 
  RooUnfoldResponse IncTwoSubJet_R0p1_ZResponse1D_ptbin4;
  IncTwoSubJet_R0p1_ZResponse1D_ptbin4.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin4;
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin4.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin4;
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin4.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncTwoSubJet_R0p1_ThetaResponse1D_ptbin4;
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin4.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin4;
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin4;
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin4;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin4.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin4;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin4.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin4;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin4.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin4;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin4.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin4;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin4.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin4;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin4.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);

  //! ptbin5 
  RooUnfoldResponse IncTwoSubJet_R0p1_ZResponse1D_ptbin5;
  IncTwoSubJet_R0p1_ZResponse1D_ptbin5.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin5;
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin5.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin5;
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin5.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncTwoSubJet_R0p1_ThetaResponse1D_ptbin5;
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin5.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin5;
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin5;
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin5;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin5.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin5;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin5.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin5;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin5.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin5;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin5.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin5;
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin5.Setup(hMeas_zg_1Dtest, hTrue_zg_1Dtest);
  RooUnfoldResponse IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin5;
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin5.Setup(hMeas_rg_1Dtest, hTrue_rg_1Dtest);

  

  
  
  TH2D* IncTruth2D = new TH2D( "IncTruth2D", "TRAIN z_{g}^{lead} vs. p_{T}^{lead}, Pythia6;p_{T}^{lead};z_{g}^{lead}", nPtBinsTrue, ptminTrue, ptmaxTrue, nZgBinsTrue, zgminTrue, zgmaxTrue);
  TH2D* IncMeas2D  = new TH2D( "IncMeas2D", "TRAIN z_{g}^{lead} vs. p_{T}^{lead}, Pythia6 #oplus GEANT;p_{T}^{lead};z_{g}^{lead}", nPtBinsMeas, ptminMeas, ptmaxMeas, nZgBinsMeas, zgminMeas, zgmaxMeas);  
  TH1D * IncJetpTMeasMCClosure1D = new TH1D("IncJetpTMeasMCClosure1D","", nPtBinsMeas, ptminMeas, ptmaxMeas);
  TH1D * IncJetpTMeasTestMCClosure1D = new TH1D("IncJetpTMeasTestMCClosure1D","", nPtBinsMeas, ptminMeas, ptmaxMeas);
  TH1D * IncJetpTTruthMCClosure1D = new TH1D("IncJetpTTruthMCClosure1D","", nPtBinsTrue, ptminTrue, ptmaxTrue);  
  TH2D* IncTestTruth2D = new TH2D( "IncTestTruth2D", "TEST z_{g}^{lead} vs. p_{T}^{lead}, Pythia6;p_{T}^{lead};z_{g}^{lead}", nPtBinsTrue, ptminTrue, ptmaxTrue, nZgBinsTrue, zgminTrue, zgmaxTrue);
  TH2D* IncTestMeas2D  = new TH2D( "IncTestMeas2D", "TEST z_{g}^{lead} vs. p_{T}^{lead}, Pythia6 #oplus GEANT;p_{T}^{lead};z_{g}^{lead}", nPtBinsMeas, ptminMeas, ptmaxMeas, nZgBinsMeas, zgminMeas, zgmaxMeas);

  TH2D* IncTruth2D_rg = new TH2D( "IncTruth2D_rg", "TRAIN R_{g}^{lead} vs. p_{T}^{lead}, Pythia6;p_{T}^{lead};R_{g}^{lead}", nPtBinsTrue, ptminTrue, ptmaxTrue, nRgBinsTrue, rgminTrue, rgmaxTrue);
  TH2D* IncMeas2D_rg  = new TH2D( "IncMeas2D_rg", "TRAIN R_{g}^{lead} vs. p_{T}^{lead}, Pythia6 #oplus GEANT;p_{T}^{lead};R_{g}^{lead}", nPtBinsMeas, ptminMeas, ptmaxMeas, nRgBinsMeas, rgminMeas, rgmaxMeas);
  TH2D* IncTestTruth2D_rg = new TH2D( "IncTestTruth2D_rg", "TEST R_{g}^{lead} vs. p_{T}^{lead}, Pythia6;p_{T}^{lead};R_{g}^{lead}", nPtBinsTrue, ptminTrue, ptmaxTrue, nRgBinsTrue, rgminTrue, rgmaxTrue);
  TH2D* IncTestMeas2D_rg  = new TH2D( "IncTestMeas2D_rg", "TEST R_{g}^{lead} vs. p_{T}^{lead}, Pythia6 #oplus GEANT;p_{T}^{lead};R_{g}^{lead}", nPtBinsMeas, ptminMeas, ptmaxMeas, nRgBinsMeas, rgminMeas, rgmaxMeas);


  TH2D* IncTruth2D_TwoSubJet_R0p1_Z = new TH2D( "IncTruth2D_TwoSubJet_R0p1_Z", "", nPtBinsTrue, ptminTrue, ptmaxTrue, nZgBinsTrue, zgminTrue, zgmaxTrue);
  TH2D* IncMeas2D_TwoSubJet_R0p1_Z  = new TH2D( "IncMeas2D_TwoSubJet_R0p1_Z", "", nPtBinsMeas, ptminMeas, ptmaxMeas, nZgBinsMeas, zgminMeas, zgmaxMeas);
  TH2D* IncTestTruth2D_TwoSubJet_R0p1_Z = new TH2D( "IncTestTruth2D_TwoSubJet_R0p1_Z", "", nPtBinsTrue, ptminTrue, ptmaxTrue, nZgBinsTrue, zgminTrue, zgmaxTrue);
  TH2D* IncTestMeas2D_TwoSubJet_R0p1_Z  = new TH2D( "IncTestMeas2D_TwoSubJet_R0p1_Z", "", nPtBinsMeas, ptminMeas, ptmaxMeas, nZgBinsMeas, zgminMeas, zgmaxMeas);

  TH2D* IncTruth2D_TwoSubJet_R0p1_Theta = new TH2D( "IncTruth2D_TwoSubJet_R0p1_Theta", "", nPtBinsTrue, ptminTrue, ptmaxTrue, nRgBinsTrue, rgminTrue, rgmaxTrue);
  TH2D* IncMeas2D_TwoSubJet_R0p1_Theta  = new TH2D( "IncMeas2D_TwoSubJet_R0p1_Theta", "", nPtBinsMeas, ptminMeas, ptmaxMeas, nRgBinsMeas, rgminMeas, rgmaxMeas);
  TH2D* IncTestTruth2D_TwoSubJet_R0p1_Theta = new TH2D( "IncTestTruth2D_TwoSubJet_R0p1_Theta", "", nPtBinsTrue, ptminTrue, ptmaxTrue, nRgBinsTrue, rgminTrue, rgmaxTrue);
  TH2D* IncTestMeas2D_TwoSubJet_R0p1_Theta  = new TH2D( "IncTestMeas2D_TwoSubJet_R0p1_Theta", "", nPtBinsMeas, ptminMeas, ptmaxMeas, nRgBinsMeas, rgminMeas, rgmaxMeas);

  
  TH2D* DeltaPtvsPt = new TH2D( "DeltaPtvsPt", "Delta p_{T} vs. p_{T}, Pythia6", nPtBinsTrue, ptminTrue, ptmaxTrue, 80, -5, 5);
  TH2D* DeltaPtvsTwoSubJet_R0p1_Z = new TH2D( "DeltaPtvsTwoSubJet_R0p1_Z", "Delta p_{T} vs. TwoSubJet_R0p1_Z, Pythia6", nZgBinsTrue, zgminTrue, zgmaxTrue, 80, -5, 5);
  TH2D* DeltaPtvsTwoSubJet_R0p1_Theta = new TH2D( "DeltaPtvsTwoSubJet_R0p1_Theta", "Delta p_{T} vs. TwoSubJet_R0p1_Theta, Pythia6", nRgBinsTrue, rgminTrue, rgmaxTrue, 80, -5, 5);

  TH2D* DeltaTwoSubJet_R0p1_ZvsPt = new TH2D( "DeltaTwoSubJet_R0p1_ZvsPt", "DeltaTwoSubJet_R0p1_Z vs. p_{T}, Pythia6", nPtBinsTrue, ptminTrue, ptmaxTrue, 80, -0.5, 0.5);
  TH2D* DeltaTwoSubJet_R0p1_ZvsTwoSubJet_R0p1_Z = new TH2D( "DeltaTwoSubJet_R0p1_ZvsTwoSubJet_R0p1_Z", "DeltaTwoSubJet_R0p1_Z vs. TwoSubJet_R0p1_Z, Pythia6", nZgBinsTrue, zgminTrue, zgmaxTrue, 80, -0.5, 0.5);
  TH2D* DeltaTwoSubJet_R0p1_ZvsTwoSubJet_R0p1_Theta = new TH2D( "DeltaTwoSubJet_R0p1_ZvsTwoSubJet_R0p1_Theta", "DeltaTwoSubJet_R0p1_Z vs. Mj, Pythia6", nRgBinsTrue, rgminTrue, rgmaxTrue, 80, -0.5, 0.5);
  
  TH2D* DeltaTwoSubJet_R0p1_ThetavsPt = new TH2D( "DeltaTwoSubJet_R0p1_ThetavsPt", "DeltaTwoSubJet_R0p1_Theta vs. p_{T}, Pythia6", nPtBinsTrue, ptminTrue, ptmaxTrue, 80, -0.5, 0.5);
  TH2D* DeltaTwoSubJet_R0p1_ThetavsTwoSubJet_R0p1_Z = new TH2D( "DeltaTwoSubJet_R0p1_ThetavsTwoSubJet_R0p1_Z", "DeltaTwoSubJet_R0p1_Theta vs. TwoSubJet_R0p1_Z, Pythia6", nZgBinsTrue, zgminTrue, zgmaxTrue, 80, -0.5, 0.5);
  TH2D* DeltaTwoSubJet_R0p1_ThetavsTwoSubJet_R0p1_Theta = new TH2D( "DeltaTwoSubJet_R0p1_ThetavsTwoSubJet_R0p1_Theta", "DeltaTwoSubJet_R0p1_Theta vs. TwoSubJet_R0p1_Theta, Pythia6", nRgBinsTrue, rgminTrue, rgmaxTrue, 80, -0.5, 0.5);
  

  //! Loop over particle level

  int missed=0;
  int N = McChain->GetEntries();
  for ( Long64_t mcEvi = 0; mcEvi<N  ; ++mcEvi ){
    if ( !(mcEvi%10000) ) cout << "Working on " << mcEvi << " / " << N << endl;
    McChain->GetEntry(mcEvi);

    if ( McJets->GetEntries() != mcnjets ){
      cerr << "McJets->GetEntries() != mcnjets" << endl;
      return -1;
    }
    
    if(RADIUS == 4 &&
       (mcEvi == 227376 || mcEvi == 227417 ||
	mcEvi == 107716 || mcEvi == 105889 ||
	mcEvi == 153437 || mcEvi == 226949 ||
	mcEvi == 226988 || mcEvi == 105616 ||
	mcEvi == 153089)) //! high weight events 
      continue;
    
    //! Fill results in vectors for easier manipulation
    //! Also check whether there's something true in the acceptance
    bool TruthInAcceptance=false;
    vector<RootResultStruct> mcresult;
    for (int j=0; j<mcnjets; ++j ){
      TStarJetVectorJet* mcjet = (TStarJetVectorJet*) McJets->At(j);

      //! Skip high weight outliers
      if ( RejectHiweights && NewGeantWeightReject ( mcjet->Pt(), mcweight, 2 ) )  {
	cout << "Skipping JET with pt=" << mcjet->Pt() << " due to high weight" << endl;
	continue;
      }

      //! Ok, record
      if ( fabs ( mcjet->Eta() ) < EtaCut ) {
	mcresult.push_back( RootResultStruct(*mcjet,
					     mcTwoSubJet_R0p1_Z[j], mcTwoSubJet_R0p1_Theta[j],
					     mcTwoSubJet_R0p1_LeadpT[j], mcTwoSubJet_R0p1_SubLeadpT[j]));	
	TruthInAcceptance=true;
      }
    }

    if ( !TruthInAcceptance ) {
      //! Skip this event, but don't count it as a loss
      continue;
    }


    //! Record Truth
    for ( vector<RootResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end(); ++mcit ){
      if ( !PrepClosure || mcEvi%2 == 0){
	IncJetpTTruthMCClosure1D->Fill(mcit->orig.Pt(), mcweight);
	// IncTruth2D->Fill( mcit->orig.Pt(), mcit->zg, mcweight );
	// IncTruth2D_rg->Fill( mcit->orig.Pt(), mcit->rg, mcweight );
	IncTruth2D_TwoSubJet_R0p1_Z->Fill( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight );
	IncTruth2D_TwoSubJet_R0p1_Theta->Fill( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight );
      }
      if ( !PrepClosure || mcEvi%2 == 1){
	// IncTestTruth2D->Fill( mcit->orig.Pt(), mcit->zg, mcweight );
	// IncTestTruth2D_rg->Fill( mcit->orig.Pt(), mcit->rg, mcweight );
	IncTestTruth2D_TwoSubJet_R0p1_Z->Fill( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight );
	IncTestTruth2D_TwoSubJet_R0p1_Theta->Fill( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight );
      }

      
      hGenJetPt_QA->Fill(mcit->orig.Pt(), mcweight);
      hGenJetPt->Fill(mcit->orig.Pt(), mcweight);
      hGenJetEta->Fill(mcit->orig.Eta(), mcweight);
      hGenJetPhi->Fill(mcit->orig.Phi(), mcweight);

      int ptbin=-1;
      for(int i = 0; i<npubptbins; ++i){
	if(mcit->orig.Pt() > pubptbins[i])
	  ptbin = i;
      }
      if(ptbin == -1)
	continue;

      
      // hGen_Zg[ptbin]->Fill(mcit->zg, mcweight);
      // hGen_Rg[ptbin]->Fill(mcit->rg, mcweight);
      hGen_TwoSubJet_R0p1_Z[ptbin]->Fill(mcit->TwoSubJet_R0p1_Z, mcweight);
      hGen_TwoSubJet_R0p1_Theta[ptbin]->Fill(mcit->TwoSubJet_R0p1_Theta, mcweight);
      hGen_TwoSubJet_R0p1_LeadpT[ptbin]->Fill(mcit->TwoSubJet_R0p1_LeadpT, mcweight);
      hGen_TwoSubJet_R0p1_SubLeadpT[ptbin]->Fill(mcit->TwoSubJet_R0p1_SubLeadpT, mcweight);
      
      
    }

    
  
    //! Get corresponding geant event
    Long64_t ppevi=-1;
    ppevi = PpChain->GetEntryNumberWithIndex( mcrunid, mceventid );

    if ( ppevi < 0 ){      
      //! Here is where we for the first time could file for loss
      if ( UseMiss ){
	for ( vector<RootResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end(); ++mcit ){
	  
	  IncPtResponse.Miss( mcit->orig.Pt(), mcweight );
	  IncPtTwoSubJet_R0p1_ZResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight );
	  IncPtTwoSubJet_R0p1_ThetaResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight );

	  hGenMissJetPt->Fill(mcit->orig.Pt(), mcweight);
	
	  if ( !PrepClosure || mcEvi%2 == 0){	    
	    IncPtResponse_MCClosure.Miss( mcit->orig.Pt(), mcweight );
	    IncPtTwoSubJet_R0p1_ZResponse2D_MCClosure.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight );
	    IncPtTwoSubJet_R0p1_ThetaResponse2D_MCClosure.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight );   
	  }

	  IncBentPtResponse.Miss ( mcit->orig.Pt(), mcweight );
	  IncBentPtTwoSubJet_R0p1_ZResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight );
	  IncBentPtTwoSubJet_R0p1_ThetaResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight );
	  IncPtBentTwoSubJet_R0p1_ZResponse2D.Miss(mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight);
	  IncPtBentTwoSubJet_R0p1_ThetaResponse2D.Miss(mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight);

	  int ptbin=-1;
	  for(int i = 0; i<npubptbins; ++i){
	    if(mcit->orig.Pt() > pubptbins[i])
	      ptbin = i;
	  }

	  if(ptbin!=-1){
	    //! get the prior shift from the maximum of the ratio difference for pythia-6 between pythia-8 and herwig-7
	    double trueshiftptwt = findweight(mcit->orig.Pt(), hJetpT_PY6_Ratio_PY8, hJetpT_PY6_Ratio_HW7);
	    IncPriorBentPtResponse.Miss ( mcit->orig.Pt(), mcweight*trueshiftptwt);

	    double priorBentTwoSubJet_R0p1_Zwt = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[ptbin], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[ptbin]);
	    double priorBentTwoSubJet_R0p1_Thetawt = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[ptbin], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[ptbin]);
	    IncPriorBentPtResponse.Miss ( mcit->orig.Pt(), mcweight*trueshiftptwt);
	    IncPriorBentPtTwoSubJet_R0p1_ZResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight*trueshiftptwt);
	    IncPriorBentPtTwoSubJet_R0p1_ThetaResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight*trueshiftptwt);
	    IncPtPriorBentTwoSubJet_R0p1_ZResponse2D.Miss(mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt);
	    IncPtPriorBentTwoSubJet_R0p1_ThetaResponse2D.Miss(mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt);
	    double priorBentPY8TwoSubJet_R0p1_Zwt = findweight_v2(mcit->orig.Pt(), hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[ptbin]);
	    double priorBentPY8TwoSubJet_R0p1_Thetawt = findweight_v2(mcit->orig.Pt(), hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[ptbin]);
	    double priorBentHW7TwoSubJet_R0p1_Zwt = findweight_v2(mcit->orig.Pt(), hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[ptbin]);
	    double priorBentHW7TwoSubJet_R0p1_Thetawt = findweight_v2(mcit->orig.Pt(), hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[ptbin]);
	  
	    double priorBentTwoSubJet_R0p1_Zwt_max;   
	    double priorBentTwoSubJet_R0p1_Thetawt_max;    
	    double priorBentTwoSubJet_R0p1_Zwt_max_env;
	    double priorBentTwoSubJet_R0p1_Thetawt_max_env;
	    double priorBentTwoSubJet_R0p1_Zwt_avg;
	    double priorBentTwoSubJet_R0p1_Thetawt_avg;
	  

	    if(mcit->orig.Pt() > 15 && mcit->orig.Pt() < 20){
	      IncTwoSubJet_R0p1_ZResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	      IncTwoSubJet_R0p1_ThetaResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	      IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	      IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	      IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	      IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	      priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[1]);
	      priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[1]);
	      priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[1]);
	      priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[1]);
	      priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[1]);
	      priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[1]);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);	    
	    }
	    if(mcit->orig.Pt() > 20 && mcit->orig.Pt() < 25){
	      IncTwoSubJet_R0p1_ZResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	      IncTwoSubJet_R0p1_ThetaResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	      IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	      IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	      IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	      IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	      priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[2]);
	      priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[2]);
	      priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[2]);
	      priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[2]);
	      priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[2]);
	      priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[2]);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	    }
	    if(mcit->orig.Pt() > 25 && mcit->orig.Pt() < 30){
	      IncTwoSubJet_R0p1_ZResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	      IncTwoSubJet_R0p1_ThetaResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	      IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	      IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	      IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	      IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	      priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[3]);
	      priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[3]);
	      priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[3]);
	      priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[3]);
	      priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[3]);
	      priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[3]);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	    }
	    if(mcit->orig.Pt() > 30 && mcit->orig.Pt() < 40){
	      IncTwoSubJet_R0p1_ZResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	      IncTwoSubJet_R0p1_ThetaResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	      IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	      IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	      IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	      IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	      priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[4]);
	      priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[4]);
	      priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[4]);
	      priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[4]);
	      priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[4]);
	      priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[4]);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	    }
	    if(mcit->orig.Pt() > 40 && mcit->orig.Pt() < 60){
	      IncTwoSubJet_R0p1_ZResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	      IncTwoSubJet_R0p1_ThetaResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	      IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	      IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	      IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	      IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	      priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[5]);
	      priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[5]);
	      priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[5]);
	      priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[5]);
	      priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[5]);
	      priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[5]);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	      IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	      IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	    }

	  }
	  

	}

	//! Skip this event
	missed++;
	continue;
      }
    }

    
  
    //! Get geant
    PpChain->GetEntry(ppevi);
    vector<RootResultStruct> ppresult;
    for (int j=0; j<ppnjets; ++j ){
      TStarJetVectorJet* ppjet = (TStarJetVectorJet*) PpJets->At(j);
      
      if ( RejectHiweights && NewGeantWeightReject ( ppjet->Pt(), mcweight, 12 ) )  {
	cout << "Skipping RECO JET with pt=" << ppjet->Pt() << " due to high weight" << endl;
	continue;
      }
      
      //! Ok, record
      if ( fabs ( ppjet->Eta() ) < EtaCut ) {
	ppresult.push_back( RootResultStruct(*ppjet,
					     ppTwoSubJet_R0p1_Z[j], ppTwoSubJet_R0p1_Theta[j],
					     ppTwoSubJet_R0p1_LeadpT[j], ppTwoSubJet_R0p1_SubLeadpT[j]));
      }
    }

    //! Record Truth
    for ( vector<RootResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end(); ++mcit ){      
      hGenJetPt_v2->Fill(mcit->orig.Pt(), mcweight);
    }
    
    
  
    // Record geant
    for ( vector<RootResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end(); ++ppit ){
      if ( !PrepClosure || mcEvi%2 == 0){
	IncJetpTMeasTestMCClosure1D->Fill(ppit->orig.Pt(), ppweight);
	IncMeas2D_TwoSubJet_R0p1_Z->Fill( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Z, ppweight );
	IncMeas2D_TwoSubJet_R0p1_Theta->Fill( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Theta, ppweight );
      }
      if ( !PrepClosure || mcEvi%2 == 1){
	IncJetpTMeasMCClosure1D->Fill(ppit->orig.Pt(), ppweight);
	IncTestMeas2D_TwoSubJet_R0p1_Z->Fill( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Z, ppweight );
	IncTestMeas2D_TwoSubJet_R0p1_Theta->Fill( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Theta, ppweight );
      }

      hRecoJetPt_QA->Fill(ppit->orig.Pt(), ppweight);
      hRecoJetPt->Fill(ppit->orig.Pt(), ppweight);
      hRecoJetEta->Fill(ppit->orig.Eta(), ppweight);
      hRecoJetPhi->Fill(ppit->orig.Phi(), ppweight);

      int ptbin=-1;
      for(int i = 0; i<npubptbins; ++i){
	if(ppit->orig.Pt() > pubptbins[i])
	  ptbin = i;
      }
      if(ptbin == -1)
	continue;

      hReco_TwoSubJet_R0p1_Z[ptbin]->Fill(ppit->TwoSubJet_R0p1_Z, ppweight);
      hReco_TwoSubJet_R0p1_Theta[ptbin]->Fill(ppit->TwoSubJet_R0p1_Theta, ppweight);
      hReco_TwoSubJet_R0p1_LeadpT[ptbin]->Fill(ppit->TwoSubJet_R0p1_LeadpT, ppweight);
      hReco_TwoSubJet_R0p1_SubLeadpT[ptbin]->Fill(ppit->TwoSubJet_R0p1_SubLeadpT, ppweight);
            
    }


    
    //! Sort them together
    vector<MatchedRootResultStruct> MatchedResult;
    for ( vector<RootResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end(); ){
      bool matched=false;
      for ( vector<RootResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end(); ){
	if ( mcit->orig.DeltaR( ppit->orig )< RCut ){
	  MatchedResult.push_back ( MatchedRootResultStruct ( *mcit, *ppit ) );
	  ppit = ppresult.erase( ppit );
	  matched=true;
	  break;
	} else{
	  ++ppit;
	}
      }
      if ( matched ) {
	hGenMatchedJetPt->Fill(mcit->orig.Pt(), mcweight);
	mcit = mcresult.erase( mcit );
      } else {
	hGenUnMatchedJetPt->Fill(mcit->orig.Pt(), mcweight);
	++mcit;
      }
    }
  

    // cout<<" size of matched result = "<<MatchedResult.size()<<endl;
  
    //! Fill Response
    for (vector<MatchedRootResultStruct>::iterator res = MatchedResult.begin(); res != MatchedResult.end(); ++res ){

      // cout<<"    in the matched iterator jet pTs, gen and reco = "<<res->second.orig.Pt()<<", "<<res->first.orig.Pt()<<endl;
      
      IncPtResponse.Fill( res->second.orig.Pt(), res->first.orig.Pt(), mcweight );
      IncPtTwoSubJet_R0p1_ZResponse2D.Fill( res->second.orig.Pt(), res->second.TwoSubJet_R0p1_Z,
					    res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Z, mcweight );
      IncPtTwoSubJet_R0p1_ThetaResponse2D.Fill( res->second.orig.Pt(), res->second.TwoSubJet_R0p1_Theta,
						res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Theta, mcweight );

      
      if(!PrepClosure || mcEvi%2 == 0){
	IncPtResponse_MCClosure.Fill( res->second.orig.Pt(), res->first.orig.Pt(), mcweight );
	IncPtTwoSubJet_R0p1_ZResponse2D_MCClosure.Fill( res->second.orig.Pt(), res->second.TwoSubJet_R0p1_Z,
							res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Z, mcweight );
	IncPtTwoSubJet_R0p1_ThetaResponse2D_MCClosure.Fill( res->second.orig.Pt(), res->second.TwoSubJet_R0p1_Theta,
							    res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Theta, mcweight );
      }

      
      float truept = res->first.orig.Pt();
      double randpt = res->second.orig.Pt()*(1-fabs(rand.Gaus(0, 0.05)));
      IncBentPtResponse.Fill ( randpt, res->first.orig.Pt(), mcweight );

      IncBentPtTwoSubJet_R0p1_ZResponse2D.Fill( randpt, res->second.TwoSubJet_R0p1_Z, res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Z, mcweight );
      IncBentPtTwoSubJet_R0p1_ThetaResponse2D.Fill( randpt, res->second.TwoSubJet_R0p1_Theta, res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Theta, mcweight );
     
      int ptbin=-1;
      for(int i = 0; i<npubptbins; ++i){
	if(res->first.orig.Pt() > pubptbins[i])
	  ptbin = i;
      }

      // cout<<"ptbin = "<<ptbin<<endl;

      if(ptbin!=-1){
	//! get the prior shift from the maximum of the ratio difference for pythia-6 between pythia-8 and herwig-7

	
	double trueshiftptwt = findweight(res->first.orig.Pt(), hJetpT_PY6_Ratio_PY8, hJetpT_PY6_Ratio_HW7);
	IncPriorBentPtResponse.Fill    ( res->second.orig.Pt(), res->first.orig.Pt(), mcweight*trueshiftptwt);
	IncPriorBentPtTwoSubJet_R0p1_ZResponse2D.Fill( res->second.orig.Pt(), res->second.TwoSubJet_R0p1_Z, res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Z, mcweight*trueshiftptwt);
	IncPriorBentPtTwoSubJet_R0p1_ThetaResponse2D.Fill( res->second.orig.Pt(), res->second.TwoSubJet_R0p1_Theta, res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Theta, mcweight*trueshiftptwt);
	double priorBentTwoSubJet_R0p1_Zwt = findweight(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[ptbin], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[ptbin]);
	double priorBentTwoSubJet_R0p1_Thetawt = findweight(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[ptbin], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[ptbin]);
	IncPtPriorBentTwoSubJet_R0p1_ZResponse2D.Fill( res->second.orig.Pt(), res->second.TwoSubJet_R0p1_Z, res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt);
	IncPtPriorBentTwoSubJet_R0p1_ThetaResponse2D.Fill( res->second.orig.Pt(), res->second.TwoSubJet_R0p1_Theta, res->first.orig.Pt(), res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt);
	double priorBentPY8TwoSubJet_R0p1_Zwt = findweight_v2(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[ptbin]);
	double priorBentPY8TwoSubJet_R0p1_Thetawt = findweight_v2(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[ptbin]);
	double priorBentHW7TwoSubJet_R0p1_Zwt = findweight_v2(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[ptbin]);
	double priorBentHW7TwoSubJet_R0p1_Thetawt = findweight_v2(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[ptbin]);
	
		
	double priorBentTwoSubJet_R0p1_Zwt_max;   
	double priorBentTwoSubJet_R0p1_Thetawt_max;    
	double priorBentTwoSubJet_R0p1_Zwt_max_env;
	double priorBentTwoSubJet_R0p1_Thetawt_max_env;
	double priorBentTwoSubJet_R0p1_Zwt_avg;
	double priorBentTwoSubJet_R0p1_Thetawt_avg;


	if(res->first.orig.Pt() > 15 && res->first.orig.Pt() < 20){
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin1.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin1.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight);
	  priorBentTwoSubJet_R0p1_Zwt_max = findweight(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[1]);
	  priorBentTwoSubJet_R0p1_Thetawt_max = findweight(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[1]);
	  priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[1]);
	  priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[1]);
	  priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[1]);
	  priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[1]);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin1.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin1.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin1.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin1.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin1.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin1.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin1.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin1.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	}
	if(res->first.orig.Pt() > 20 && res->first.orig.Pt() < 25){
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin2.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin2.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight);

	  priorBentTwoSubJet_R0p1_Zwt_max = findweight(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[2]);
	  priorBentTwoSubJet_R0p1_Thetawt_max = findweight(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[2]);
	  priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[2]);
	  priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[2]);
	  priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[2]);
	  priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[2]);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin2.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin2.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin2.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin2.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin2.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin2.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin2.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin2.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	}
	if(res->first.orig.Pt() > 25 && res->first.orig.Pt() < 30){
	  priorBentTwoSubJet_R0p1_Zwt_max = findweight(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[3]);
	  priorBentTwoSubJet_R0p1_Thetawt_max = findweight(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[3]);
	  priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[3]);
	  priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[3]);
	  priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[3]);
	  priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[3]);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin3.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin3.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin3.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin3.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin3.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin3.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin3.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin3.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin3.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin3.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	}
	if(res->first.orig.Pt() > 30 && res->first.orig.Pt() < 40){
	  priorBentTwoSubJet_R0p1_Zwt_max = findweight(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[4]);
	  priorBentTwoSubJet_R0p1_Thetawt_max = findweight(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[4]);
	  priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[4]);
	  priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[4]);
	  priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[4]);
	  priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[4]);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin4.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin4.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin4.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin4.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin4.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin4.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin4.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin4.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin4.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin4.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	}
	if(res->first.orig.Pt() > 40 && res->first.orig.Pt() < 60){
	  priorBentTwoSubJet_R0p1_Zwt_max = findweight(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[5]);
	  priorBentTwoSubJet_R0p1_Thetawt_max = findweight(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[5]);
	  priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[5]);
	  priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[5]);
	  priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[5]);
	  priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(res->first.TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[5]);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin5.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin5.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin5.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin5.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin5.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin5.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin5.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin5.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin5.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin5.Fill( res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Fill( res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	}
	

      }
      
      
      DeltaPtvsPt->Fill ( truept,  res->second.orig.Pt()-res->first.orig.Pt(), mcweight);
      DeltaPtvsTwoSubJet_R0p1_Z->Fill ( res->first.TwoSubJet_R0p1_Z,  res->second.orig.Pt()-res->first.orig.Pt(), mcweight);
      DeltaPtvsTwoSubJet_R0p1_Theta->Fill ( res->first.TwoSubJet_R0p1_Theta,  res->second.orig.Pt()-res->first.orig.Pt(), mcweight);

      DeltaTwoSubJet_R0p1_ZvsPt->Fill ( truept,  res->second.TwoSubJet_R0p1_Z-res->first.TwoSubJet_R0p1_Z, mcweight);
      DeltaTwoSubJet_R0p1_ZvsTwoSubJet_R0p1_Z->Fill ( res->first.TwoSubJet_R0p1_Z,  res->second.TwoSubJet_R0p1_Z-res->first.TwoSubJet_R0p1_Z, mcweight);
      DeltaTwoSubJet_R0p1_ZvsTwoSubJet_R0p1_Theta->Fill ( res->first.TwoSubJet_R0p1_Theta,  res->second.TwoSubJet_R0p1_Z-res->first.TwoSubJet_R0p1_Z, mcweight);
	
      DeltaTwoSubJet_R0p1_ThetavsPt->Fill ( truept,  res->second.TwoSubJet_R0p1_Theta-res->first.TwoSubJet_R0p1_Theta, mcweight);
      DeltaTwoSubJet_R0p1_ThetavsTwoSubJet_R0p1_Theta->Fill ( res->first.TwoSubJet_R0p1_Theta,  res->second.TwoSubJet_R0p1_Theta-res->first.TwoSubJet_R0p1_Theta, mcweight);
      DeltaTwoSubJet_R0p1_ThetavsTwoSubJet_R0p1_Z->Fill ( res->first.TwoSubJet_R0p1_Z,  res->second.TwoSubJet_R0p1_Theta-res->first.TwoSubJet_R0p1_Theta, mcweight);
      
      
      hMatrixJetPt->Fill(res->second.orig.Pt(), res->first.orig.Pt(), mcweight);

      // hGenMatchedJetPt->Fill(res->first.orig.Pt(), mcweight);
      hRecoMatchedJetPt->Fill(res->second.orig.Pt(), mcweight);
      
      ptbin=-1;
      for(int i = 0; i<npubptbins; ++i){
	if(res->second.orig.Pt() > pubptbins[i])
	  ptbin = i;
      }
      if(ptbin != -1){      
	hPtResponse->Fill(res->second.orig.Pt(), res->first.orig.Pt(), mcweight);
	hPtResponse_transpose->Fill(res->first.orig.Pt(), res->second.orig.Pt(), mcweight);

	hTwoSubJet_R0p1_ZResponse[ptbin]->Fill(res->second.TwoSubJet_R0p1_Z, res->first.TwoSubJet_R0p1_Z, mcweight);
	hTwoSubJet_R0p1_ThetaResponse[ptbin]->Fill(res->second.TwoSubJet_R0p1_Theta, res->first.TwoSubJet_R0p1_Theta, mcweight);
	
	hJER[ptbin]->Fill((float)(res->second.orig.Pt()/res->first.orig.Pt()), mcweight);

	hJER_TwoSubJet_R0p1_Z[ptbin]->Fill((float)(res->second.TwoSubJet_R0p1_Z/res->first.TwoSubJet_R0p1_Z), mcweight);
	hJER_TwoSubJet_R0p1_Theta[ptbin]->Fill((float)(res->second.TwoSubJet_R0p1_Theta/res->first.TwoSubJet_R0p1_Theta), mcweight);
      }
      

      ptbin=-1;
      for(int i = 0; i<npubptbins; ++i){
	if(res->first.orig.Pt() > pubptbins[i])
	  ptbin = i;
      }
      if(ptbin != -1){      
	hJER_v2[ptbin]->Fill((float)(res->second.orig.Pt()/res->first.orig.Pt()), mcweight);
	hJER_v2_TwoSubJet_R0p1_Z[ptbin]->Fill((float)(res->second.TwoSubJet_R0p1_Z/res->first.TwoSubJet_R0p1_Z), mcweight);
	hJER_v2_TwoSubJet_R0p1_Theta[ptbin]->Fill((float)(res->second.TwoSubJet_R0p1_Theta/res->first.TwoSubJet_R0p1_Theta), mcweight);

      }      
    }
  
    
    //! Fill misses and fakes
    if ( UseMiss ){
      for ( vector<RootResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end(); ++mcit ){
	IncPtResponse.Miss( mcit->orig.Pt(), mcweight );
	IncPtTwoSubJet_R0p1_ZResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight );
	IncPtTwoSubJet_R0p1_ThetaResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight );

	hGenMissJetPt->Fill(mcit->orig.Pt(), mcweight);
	
	if ( !PrepClosure || mcEvi%2 == 0){	    
	  IncPtResponse_MCClosure.Miss( mcit->orig.Pt(), mcweight );
	  IncPtTwoSubJet_R0p1_ZResponse2D_MCClosure.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight );
	  IncPtTwoSubJet_R0p1_ThetaResponse2D_MCClosure.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight );	    
	}

	IncBentPtResponse.Miss ( mcit->orig.Pt(), mcweight );
	IncBentPtTwoSubJet_R0p1_ZResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight );
	IncBentPtTwoSubJet_R0p1_ThetaResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight );
	IncPtBentTwoSubJet_R0p1_ZResponse2D.Miss(mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight);
	IncPtBentTwoSubJet_R0p1_ThetaResponse2D.Miss(mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight);

	int ptbin=-1;
	for(int i = 0; i<npubptbins; ++i){
	  if(mcit->orig.Pt() > pubptbins[i])
	    ptbin = i;
	}

	if(ptbin!=-1){
	  //! get the prior shift from the maximum of the ratio difference for pythia-6 between pythia-8 and herwig-7
	  double trueshiftptwt = findweight(mcit->orig.Pt(), hJetpT_PY6_Ratio_PY8, hJetpT_PY6_Ratio_HW7);
	  IncPriorBentPtResponse.Miss ( mcit->orig.Pt(), mcweight*trueshiftptwt);

	  double priorBentTwoSubJet_R0p1_Zwt = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[ptbin], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[ptbin]);
	  double priorBentTwoSubJet_R0p1_Thetawt = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[ptbin], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[ptbin]);
	  IncPriorBentPtResponse.Miss ( mcit->orig.Pt(), mcweight*trueshiftptwt);
	  IncPriorBentPtTwoSubJet_R0p1_ZResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight*trueshiftptwt);
	  IncPriorBentPtTwoSubJet_R0p1_ThetaResponse2D.Miss( mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight*trueshiftptwt);
	  IncPtPriorBentTwoSubJet_R0p1_ZResponse2D.Miss(mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt);
	  IncPtPriorBentTwoSubJet_R0p1_ThetaResponse2D.Miss(mcit->orig.Pt(), mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt);
	  double priorBentPY8TwoSubJet_R0p1_Zwt = findweight_v2(mcit->orig.Pt(), hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[ptbin]);
	  double priorBentPY8TwoSubJet_R0p1_Thetawt = findweight_v2(mcit->orig.Pt(), hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[ptbin]);
	  double priorBentHW7TwoSubJet_R0p1_Zwt = findweight_v2(mcit->orig.Pt(), hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[ptbin]);
	  double priorBentHW7TwoSubJet_R0p1_Thetawt = findweight_v2(mcit->orig.Pt(), hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[ptbin]);
	  
	  double priorBentTwoSubJet_R0p1_Zwt_max;   
	  double priorBentTwoSubJet_R0p1_Thetawt_max;    
	  double priorBentTwoSubJet_R0p1_Zwt_max_env;
	  double priorBentTwoSubJet_R0p1_Thetawt_max_env;
	  double priorBentTwoSubJet_R0p1_Zwt_avg;
	  double priorBentTwoSubJet_R0p1_Thetawt_avg;
	  

	  if(mcit->orig.Pt() > 15 && mcit->orig.Pt() < 20){
	    IncTwoSubJet_R0p1_ZResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	    IncTwoSubJet_R0p1_ThetaResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	    IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	    IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	    IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	    IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	    priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[1]);
	    priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[1]);
	    priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[1]);
	    priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[1]);
	    priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[1]);
	    priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[1], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[1]);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin1.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin1.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);	    
	  }
	  if(mcit->orig.Pt() > 20 && mcit->orig.Pt() < 25){
	    IncTwoSubJet_R0p1_ZResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	    IncTwoSubJet_R0p1_ThetaResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	    IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	    IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	    IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	    IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	    priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[2]);
	    priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[2]);
	    priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[2]);
	    priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[2]);
	    priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[2]);
	    priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[2], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[2]);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin2.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin2.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	  }
	  if(mcit->orig.Pt() > 25 && mcit->orig.Pt() < 30){
	    IncTwoSubJet_R0p1_ZResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	    IncTwoSubJet_R0p1_ThetaResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	    IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	    IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	    IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	    IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	    priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[3]);
	    priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[3]);
	    priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[3]);
	    priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[3]);
	    priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[3]);
	    priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[3], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[3]);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin3.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin3.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	  }
	  if(mcit->orig.Pt() > 30 && mcit->orig.Pt() < 40){
	    IncTwoSubJet_R0p1_ZResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	    IncTwoSubJet_R0p1_ThetaResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	    IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	    IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	    IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	    IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	    priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[4]);
	    priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[4]);
	    priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[4]);
	    priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[4]);
	    priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[4]);
	    priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[4], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[4]);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin4.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin4.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	  }
	  if(mcit->orig.Pt() > 40 && mcit->orig.Pt() < 60){
	    IncTwoSubJet_R0p1_ZResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight);
	    IncTwoSubJet_R0p1_ThetaResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight);
	    IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentPY8TwoSubJet_R0p1_Zwt);
	    IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentPY8TwoSubJet_R0p1_Thetawt);
	    IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentHW7TwoSubJet_R0p1_Zwt);
	    IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentHW7TwoSubJet_R0p1_Thetawt);
	    priorBentTwoSubJet_R0p1_Zwt_max = findweight(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[5]);
	    priorBentTwoSubJet_R0p1_Thetawt_max = findweight(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[5]);
	    priorBentTwoSubJet_R0p1_Zwt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[5]);
	    priorBentTwoSubJet_R0p1_Thetawt_max_env = findweight_max(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[5]);
	    priorBentTwoSubJet_R0p1_Zwt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Z, hJetTwoSubJet_R0p1_Z_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Z_PY6_Ratio_HW7[5]);
	    priorBentTwoSubJet_R0p1_Thetawt_avg = findweight_avg(mcit->TwoSubJet_R0p1_Theta, hJetTwoSubJet_R0p1_Theta_PY6_Ratio_PY8[5], hJetTwoSubJet_R0p1_Theta_PY6_Ratio_HW7[5]);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_max_env);
	    IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin5.Miss( mcit->TwoSubJet_R0p1_Z, mcweight*priorBentTwoSubJet_R0p1_Zwt_avg);
	    IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin5.Miss( mcit->TwoSubJet_R0p1_Theta, mcweight*priorBentTwoSubJet_R0p1_Thetawt_avg);
	  }

	}
      }
    }
  
    if ( UseFakes ){
      for ( vector<RootResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end(); ++ppit ){
	IncPtResponse.Fake( ppit->orig.Pt(), ppweight );
	IncPtTwoSubJet_R0p1_ZResponse2D.Fake( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Z, ppweight );
	IncPtTwoSubJet_R0p1_ThetaResponse2D.Fake( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Theta, ppweight );

	hRecoFakeJetPt->Fill(ppit->orig.Pt(), ppweight);		
	
	if ( !PrepClosure || mcEvi%2 == 0){	  
	  IncPtResponse_MCClosure.Fake( ppit->orig.Pt(), ppweight );
	  IncPtTwoSubJet_R0p1_ZResponse2D_MCClosure.Fake( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Z, ppweight );
	  IncPtTwoSubJet_R0p1_ThetaResponse2D_MCClosure.Fake( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Theta, ppweight );
	}
	
	IncBentPtResponse.Fake ( ppit->orig.Pt(), ppweight );

	IncBentPtTwoSubJet_R0p1_ZResponse2D.Fake( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Z, ppweight );
	IncBentPtTwoSubJet_R0p1_ThetaResponse2D.Fake( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Theta, ppweight );
	IncPtBentTwoSubJet_R0p1_ZResponse2D.Fake(ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Z, ppweight);
	IncPtBentTwoSubJet_R0p1_ThetaResponse2D.Fake(ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Theta, ppweight);
	IncPriorBentPtResponse.Fake ( ppit->orig.Pt(), ppweight );
	IncPriorBentPtTwoSubJet_R0p1_ZResponse2D.Fake( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Z, ppweight );
	IncPriorBentPtTwoSubJet_R0p1_ThetaResponse2D.Fake( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Theta, ppweight );
	IncPtPriorBentTwoSubJet_R0p1_ZResponse2D.Fake(ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Z, ppweight);
	IncPtPriorBentTwoSubJet_R0p1_ThetaResponse2D.Fake(ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Theta, ppweight);

	
	if(ppit->orig.Pt() > 15 && ppit->orig.Pt() < 20){
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin1.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin1.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin1.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin1.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin1.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin1.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin1.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin1.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin1.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin1.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	}
	if(ppit->orig.Pt() > 20 && ppit->orig.Pt() < 25){
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin2.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin2.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin2.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin2.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin2.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin2.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin2.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin2.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin2.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin2.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	}
	if(ppit->orig.Pt() > 25 && ppit->orig.Pt() < 30){
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin3.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin3.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin3.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin3.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin3.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin3.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin3.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin3.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin3.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin3.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	}
	if(ppit->orig.Pt() > 30 && ppit->orig.Pt() < 40){
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin4.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin4.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin4.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin4.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin4.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin4.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin4.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin4.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin4.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin4.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	}
	if(ppit->orig.Pt() > 40 && ppit->orig.Pt() < 60){
	  IncTwoSubJet_R0p1_ZResponse1D_ptbin5.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin5.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin5.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin5.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin5.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin5.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin5.Fake( ppit->TwoSubJet_R0p1_Z, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin5.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin5.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin5.Fake( ppit->TwoSubJet_R0p1_Theta, ppweight);
	}
	
      } 
    }
  
  }

  cout << "Misssed " << missed << endl;

  //! Done

  fout->Write();
  IncPtResponse.SetName("IncPtResponse");
  IncPtResponse.Write();

  IncPtTwoSubJet_R0p1_ZResponse2D.SetName("IncPtTwoSubJet_R0p1_ZResponse2D");
  IncPtTwoSubJet_R0p1_ZResponse2D.Write();

  IncPtTwoSubJet_R0p1_ThetaResponse2D.SetName("IncPtTwoSubJet_R0p1_ThetaResponse2D");
  IncPtTwoSubJet_R0p1_ThetaResponse2D.Write();
  
  IncPtResponse_MCClosure.SetName("IncPtResponse_MCClosure");
  IncPtResponse_MCClosure.Write();

  IncPtTwoSubJet_R0p1_ZResponse2D_MCClosure.SetName("IncPtTwoSubJet_R0p1_ZResponse2D_MCClosure");
  IncPtTwoSubJet_R0p1_ZResponse2D_MCClosure.Write();

  IncPtTwoSubJet_R0p1_ThetaResponse2D_MCClosure.SetName("IncPtTwoSubJet_R0p1_ThetaResponse2D_MCClosure");
  IncPtTwoSubJet_R0p1_ThetaResponse2D_MCClosure.Write();

  
  IncBentPtResponse.SetName("IncBentPtResponse");
  IncBentPtResponse.Write();

  IncPriorBentPtResponse.SetName("IncPriorBentPtResponse");
  IncPriorBentPtResponse.Write();  


  //! two subjet histos
  IncBentPtTwoSubJet_R0p1_ZResponse2D.SetName("IncBentPtTwoSubJet_R0p1_ZResponse2D");
  IncBentPtTwoSubJet_R0p1_ZResponse2D.Write();

  IncBentPtTwoSubJet_R0p1_ThetaResponse2D.SetName("IncBentPtTwoSubJet_R0p1_ThetaResponse2D");
  IncBentPtTwoSubJet_R0p1_ThetaResponse2D.Write();
  
  IncPtBentTwoSubJet_R0p1_ZResponse2D.SetName("IncPtBentTwoSubJet_R0p1_ZResponse2D");
  IncPtBentTwoSubJet_R0p1_ZResponse2D.Write();

  IncPtBentTwoSubJet_R0p1_ThetaResponse2D.SetName("IncPtBentTwoSubJet_R0p1_ThetaResponse2D");
  IncPtBentTwoSubJet_R0p1_ThetaResponse2D.Write();

  IncPriorBentPtTwoSubJet_R0p1_ZResponse2D.SetName("IncPriorBentPtTwoSubJet_R0p1_ZResponse2D");
  IncPriorBentPtTwoSubJet_R0p1_ZResponse2D.Write();

  IncPriorBentPtTwoSubJet_R0p1_ThetaResponse2D.SetName("IncPriorBentPtTwoSubJet_R0p1_ThetaResponse2D");
  IncPriorBentPtTwoSubJet_R0p1_ThetaResponse2D.Write();
  
  IncPtPriorBentTwoSubJet_R0p1_ZResponse2D.SetName("IncPtPriorBentTwoSubJet_R0p1_ZResponse2D");
  IncPtPriorBentTwoSubJet_R0p1_ZResponse2D.Write();

  IncPtPriorBentTwoSubJet_R0p1_ThetaResponse2D.SetName("IncPtPriorBentTwoSubJet_R0p1_ThetaResponse2D");
  IncPtPriorBentTwoSubJet_R0p1_ThetaResponse2D.Write();

  IncTwoSubJet_R0p1_ZResponse1D_ptbin1.SetName("IncTwoSubJet_R0p1_ZResponse1D_ptbin1");
  IncTwoSubJet_R0p1_ZResponse1D_ptbin1.Write();
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin1.SetName("IncTwoSubJet_R0p1_ThetaResponse1D_ptbin1");
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin1.Write();
  IncTwoSubJet_R0p1_ZResponse1D_ptbin2.SetName("IncTwoSubJet_R0p1_ZResponse1D_ptbin2");
  IncTwoSubJet_R0p1_ZResponse1D_ptbin2.Write();
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin2.SetName("IncTwoSubJet_R0p1_ThetaResponse1D_ptbin2");
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin2.Write();
  IncTwoSubJet_R0p1_ZResponse1D_ptbin3.SetName("IncTwoSubJet_R0p1_ZResponse1D_ptbin3");
  IncTwoSubJet_R0p1_ZResponse1D_ptbin3.Write();
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin3.SetName("IncTwoSubJet_R0p1_ThetaResponse1D_ptbin3");
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin3.Write();
  IncTwoSubJet_R0p1_ZResponse1D_ptbin4.SetName("IncTwoSubJet_R0p1_ZResponse1D_ptbin4");
  IncTwoSubJet_R0p1_ZResponse1D_ptbin4.Write();
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin4.SetName("IncTwoSubJet_R0p1_ThetaResponse1D_ptbin4");
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin4.Write();
  IncTwoSubJet_R0p1_ZResponse1D_ptbin5.SetName("IncTwoSubJet_R0p1_ZResponse1D_ptbin5");
  IncTwoSubJet_R0p1_ZResponse1D_ptbin5.Write();
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin5.SetName("IncTwoSubJet_R0p1_ThetaResponse1D_ptbin5");
  IncTwoSubJet_R0p1_ThetaResponse1D_ptbin5.Write();

  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin1.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin1");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin1.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin1.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin1");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin1.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin1.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin1");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin1.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin1.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin1");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin1.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin1.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin1");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin1.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin1.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin1");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin1.Write();

  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin2.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin2");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin2.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin2.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin2");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin2.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin2.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin2");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin2.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin2.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin2");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin2.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin2.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin2");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin2.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin2.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin2");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin2.Write();

  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin3.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin3");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin3.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin3.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin3");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin3.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin3.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin3");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin3.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin3.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin3");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin3.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin3.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin3");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin3.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin3.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin3");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin3.Write();

  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin4.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin4");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin4.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin4.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin4");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin4.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin4.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin4");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin4.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin4.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin4");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin4.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin4.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin4");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin4.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin4.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin4");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin4.Write();

  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin5.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin5");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_ptbin5.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin5.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin5");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_max_env_ptbin5.Write();
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin5.SetName("IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin5");
  IncPriorBentTwoSubJet_R0p1_ZResponse1D_avg_ptbin5.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin5.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin5");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_ptbin5.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin5.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin5");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_max_env_ptbin5.Write();
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin5.SetName("IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin5");
  IncPriorBentTwoSubJet_R0p1_ThetaResponse1D_avg_ptbin5.Write();
  
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin1.SetName("IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin1");
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin1.Write();
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin1.SetName("IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin1");
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Write();
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin2.SetName("IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin2");
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin2.Write();
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin2.SetName("IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin2");
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Write();
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin3.SetName("IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin3");
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin3.Write();
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin3.SetName("IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin3");
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Write();
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin4.SetName("IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin4");
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin4.Write();
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin4.SetName("IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin4");
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Write();
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin5.SetName("IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin5");
  IncPriorBentPY8TwoSubJet_R0p1_ZResponse1D_ptbin5.Write();
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin5.SetName("IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin5");
  IncPriorBentPY8TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Write();

  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin1.SetName("IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin1");
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin1.Write();
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin1.SetName("IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin1");
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin1.Write();
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin2.SetName("IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin2");
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin2.Write();
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin2.SetName("IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin2");
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin2.Write();
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin3.SetName("IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin3");
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin3.Write();
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin3.SetName("IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin3");
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin3.Write();
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin4.SetName("IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin4");
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin4.Write();
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin4.SetName("IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin4");
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin4.Write();
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin5.SetName("IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin5");
  IncPriorBentHW7TwoSubJet_R0p1_ZResponse1D_ptbin5.Write();
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin5.SetName("IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin5");
  IncPriorBentHW7TwoSubJet_R0p1_ThetaResponse1D_ptbin5.Write();
  

  cout << " Wrote to" << endl << OutFileName << endl;

  return 0;

}


