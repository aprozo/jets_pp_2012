//! pp Run12 zg, rg analysis 
//! Code to read analysis ntuples and produce histograms  
//! Raghav Kunnawalkam Elayavalli & Kolja Kauder
//! contact - raghavke@wayne.edu
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

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"

#include "Binning.h"

using namespace std;


class RootGroomingResultStruct{
public:
  TStarJetVectorJet orig;
  double TwoSubJet_R0p1_Z;
  double TwoSubJet_R0p1_Theta;
  double TwoSubJet_R0p1_LeadpT;
  double TwoSubJet_R0p1_SubLeadpT;	  
  RootGroomingResultStruct ( TStarJetVectorJet orig,
			     double TwoSubJet_R0p1_Z, double TwoSubJet_R0p1_Theta,
			     double TwoSubJet_R0p1_LeadpT, double TwoSubJet_R0p1_SubLeadpT) :
    orig(orig),
    TwoSubJet_R0p1_Z(TwoSubJet_R0p1_Z),
    TwoSubJet_R0p1_Theta(TwoSubJet_R0p1_Theta),
    TwoSubJet_R0p1_LeadpT(TwoSubJet_R0p1_LeadpT),
    TwoSubJet_R0p1_SubLeadpT(TwoSubJet_R0p1_SubLeadpT)
  {};
  
  ClassDef(RootGroomingResultStruct,1)

};

typedef pair<RootGroomingResultStruct,RootGroomingResultStruct> MatchedRootGroomingResultStruct;


int PrepUnfolding (int RADIUS=4) {

  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  gStyle->SetTitleBorderSize(0);	
  gStyle->SetHistLineWidth(2);

  TLegend* leg = 0;

  float RCut = (float)RADIUS/10;
  
  //! For pp Data
  TString PpLevelFile = Form("Results/Pp12_JP2_TSJ_NoEff_wchpionMassSet_R0%d_hadded_v2.root", RADIUS);

  float EtaCut = 1.0-RCut;

  bool addpTCut = true;


  // Output
  // ------
  TString OutFileName = "Results/ForUnfolding_TwoSubJet_";
  if ( addpTCut ) OutFileName += "WithpTCuts";
  if ( !addpTCut ) OutFileName += "NopTCuts";
  OutFileName += gSystem->BaseName(PpLevelFile);

  //! Set up pp events
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

  TH1D * hRecoJetPt = new TH1D("hRecoJetPt","", nPtBinsMeas, ptminMeas, ptmaxMeas);
  TH1D * hRecoJetEta = new TH1D("hRecoJetEta","", 50, -1, 1);
  TH1D * hRecoJetPhi = new TH1D("hRecoJetPhi","", 50, 3.15, 3.15);
  TH1D * hReco_TwoSubJet_R0p1_Z[npubptbins];
  TH1D * hReco_TwoSubJet_R0p1_Theta[npubptbins];
  TH1D * hReco_TwoSubJet_R0p1_LeadpT[npubptbins];
  TH1D * hReco_TwoSubJet_R0p1_SubLeadpT[npubptbins];

  for(int i = 0; i<npubptbins; ++i){    
    hReco_TwoSubJet_R0p1_Z[i] = new TH1D(Form("hReco_TwoSubJet_R0p1_Z_ptbin_%d", i), "", nZgBinsTrue, zgminTrue, zgmaxTrue);
    hReco_TwoSubJet_R0p1_Theta[i] = new TH1D(Form("hReco_TwoSubJet_R0p1_Theta_ptbin_%d", i), "", nRgBinsTrue, rgminTrue, rgmaxTrue);
    hReco_TwoSubJet_R0p1_LeadpT[i] = new TH1D(Form("hReco_TwoSubJet_R0p1_LeadpT_ptbin_%d", i), "", 60, 0, 60);
    hReco_TwoSubJet_R0p1_SubLeadpT[i] = new TH1D(Form("hReco_TwoSubJet_R0p1_SubLeadpT_ptbin_%d", i), "", 60, 0, 60);    
  }
    
  // Only need one of these
  TH2D* IncPtTwoSubJet_R0p1_ZMeas2D  = new TH2D( "IncPtTwoSubJet_R0p1_ZMeas2D", "Measured z_{SJ} vs. p_{T};p_{T};z_{SJ}", nPtBinsMeas, ptminMeas, ptmaxMeas, nZgBinsMeas, zgminMeas, zgmaxMeas);
  TH2D* IncPtTwoSubJet_R0p1_ThetaMeas2D  = new TH2D( "IncPtTwoSubJet_R0p1_ThetaMeas2D", "Measured #theta_{SJ} vs. p_{T};p_{T};#theta_{SJ}", nPtBinsMeas, ptminMeas, ptmaxMeas, nRgBinsMeas, rgminMeas, rgmaxMeas);

  //! Loop over measured level
  for ( Long64_t ppEvi = 0; ppEvi< PpChain->GetEntries() ; ++ppEvi ){
    if ( !(ppEvi%10000) ) cout << "Working on " << ppEvi << " / " << PpChain->GetEntries() << endl;

    PpChain->GetEntry(ppEvi);
    vector<RootGroomingResultStruct> ppresult;

    for (int j=0; j<ppnjets; ++j ){
      
      TStarJetVectorJet* ppjet = (TStarJetVectorJet*) PpJets->At(j);
	  
      if ( fabs ( ppjet->Eta() ) < EtaCut ) {
	ppresult.push_back( RootGroomingResultStruct(*ppjet, 
						     ppTwoSubJet_R0p1_Z[j], ppTwoSubJet_R0p1_Theta[j],
						     ppTwoSubJet_R0p1_LeadpT[j], ppTwoSubJet_R0p1_SubLeadpT[j]));
      }

    }


    //! Record Measured
    for (vector<RootGroomingResultStruct>::iterator ppit = ppresult.begin();
	 ppit != ppresult.end();
	 ++ppit ){

      IncPtTwoSubJet_R0p1_ZMeas2D->Fill( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Z, ppweight );
      IncPtTwoSubJet_R0p1_ThetaMeas2D->Fill( ppit->orig.Pt(), ppit->TwoSubJet_R0p1_Theta, ppweight );
      
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

    	
    }//! ppit iteration
    
  }//! event loop 

  fout->Write();

  cout << " Wrote to" << endl << OutFileName << endl;
  return 0;
}
