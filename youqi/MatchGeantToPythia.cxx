//! Code to read in geant + pythia output trees and match them
//! Youqi Song & Raghav Kunnawalkam Elayavalli & Kolja Kauder
//! contact - youqi.song@yale.edu
//! HAS to be compiled,
//! root -l MatchGeantToPythia.cxx+("pythia.root", "geant.root")

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
#include <TGraph.h>

#include <TStarJetVector.h>
#include <TStarJetVectorJet.h>
#include <TStarJetPicoReader.h>

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <exception>

using namespace std;

//! ----------------------------------------------------
class RootResultStruct
{
public:
    TStarJetVectorJet orig;
    double pt;
    double eta;
    double y;
    double phi;
    double m;
    int n;
    double q;
    double rg;
    double zg;
    double mg;
    double weight;
    int reject;
    int mcevi;
    int mult;
    double nef;
    double vz;
    RootResultStruct(TStarJetVectorJet orig, double pt, double eta, double y, double phi, double m, int n, double q, double rg, double zg, double mg, double nef, double weight, int reject, int mcevi, int mult, double vz) : orig(orig), pt(pt), eta(eta), y(y), phi(phi), m(m), n(n), q(q), rg(rg), zg(zg), mg(mg), nef(nef), weight(weight), reject(reject), mcevi(mcevi), mult(mult), vz(vz){};
    ClassDef(RootResultStruct, 1)
};

typedef pair<RootResultStruct, RootResultStruct> MatchedRootResultStruct;

int MatchGeantToPythia(string McFile, string PpFile, string OutFile = "test.root", int RADIUS = 4, int mode = 0)
{
    bool RejectHiweights = true;
    float RCut = (float)RADIUS / 10;
    float EtaCut = 1.0 - RCut;
    int nj = 0;

    TString dir = "/gpfs01/star/pwg/youqi/run12/final_0628/results/";

    TFile *Mcf = new TFile(dir + "pythia/" + McFile);
    TTree *McChain = (TTree *)Mcf->Get("ResultTree");
    McChain->BuildIndex("runid", "eventid");

    TClonesArray *McJets = new TClonesArray("TStarJetVectorJet");
    McChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
    McChain->SetBranchAddress("Jets", &McJets);

    int mcrunid;
    McChain->SetBranchAddress("runid", &mcrunid);
    int mceventid;
    McChain->SetBranchAddress("eventid", &mceventid);
    double mcweight;
    McChain->SetBranchAddress("weight", &mcweight);
    int mcnjets = 0;
    McChain->SetBranchAddress("njets", &mcnjets);
    double mcvz = 0;
    McChain->SetBranchAddress("vz", &mcvz);
    int mcmult;
    McChain->SetBranchAddress("mult", &mcmult);
    int mcreject;
    McChain->SetBranchAddress("reject", &mcreject);
    float mcpttot;
    McChain->SetBranchAddress("pttot", &mcpttot);
    double mcpt[1000];
    McChain->SetBranchAddress("pt", mcpt);
    double mcm[1000];
    McChain->SetBranchAddress("m", mcm);
    double mcq[1000];
    McChain->SetBranchAddress("q", mcq);
    double mcrg[1000];
    McChain->SetBranchAddress("rg", mcrg);
    double mczg[1000];
    McChain->SetBranchAddress("zg", mczg);
    double mcmg[1000];
    McChain->SetBranchAddress("mg", mcmg);
    int mcn[1000];
    McChain->SetBranchAddress("n", mcn);
    double mcnef[1000];
    McChain->SetBranchAddress("nef", mcnef);

    TFile *Ppf = new TFile(dir + "geant/" + PpFile);
    TTree *PpChain = (TTree *)Ppf->Get("ResultTree");
    PpChain->BuildIndex("runid", "eventid");

    TClonesArray *PpJets = new TClonesArray("TStarJetVectorJet");
    PpChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
    PpChain->SetBranchAddress("Jets", &PpJets);

    int ppeventid;
    PpChain->SetBranchAddress("eventid", &ppeventid);
    int pprunid;
    PpChain->SetBranchAddress("runid", &pprunid);
    int pprunid1;
    PpChain->SetBranchAddress("runid1", &pprunid1);    
    double ppweight;
    PpChain->SetBranchAddress("weight", &ppweight);
    int rcreject;
    PpChain->SetBranchAddress("reject", &rcreject);
    int rcmult;
    PpChain->SetBranchAddress("mult", &rcmult);
    int ppnjets = 0;
    PpChain->SetBranchAddress("njets", &ppnjets);
    double ppvz = 0;
    PpChain->SetBranchAddress("vz", &ppvz);
    float pppttot;
    PpChain->SetBranchAddress("pttot", &pppttot);
    double rcpt[1000];
    PpChain->SetBranchAddress("pt", rcpt);
    double rcm[1000];
    PpChain->SetBranchAddress("m", rcm);
    double rcq[1000];
    PpChain->SetBranchAddress("q", rcq);
    double rcrg[1000];
    PpChain->SetBranchAddress("rg", rcrg);
    double rczg[1000];
    PpChain->SetBranchAddress("zg", rczg);
    double rcmg[1000];
    PpChain->SetBranchAddress("mg", rcmg);
    int rcn[1000];
    PpChain->SetBranchAddress("n", rcn);
    double rcnef[1000];
    PpChain->SetBranchAddress("nef", rcnef);
   
    //! Output and histograms
    TFile *fout = new TFile(dir + "matched/" + OutFile, "RECREATE");
    TTree *MatchedTree = new TTree("MatchedTree", "Matched Jets");
    double pt;
    MatchedTree->Branch("pt", &pt, "pt/D");
    double pts;
    MatchedTree->Branch("pts", &pts, "pts/D");
    double m;
    MatchedTree->Branch("m", &m, "m/D");
    double ms;
    MatchedTree->Branch("ms", &ms, "ms/D");
    double q;
    MatchedTree->Branch("q", &q, "q/D");
    double qs;
    MatchedTree->Branch("qs", &qs, "qs/D");
    double rg;
    MatchedTree->Branch("rg", &rg, "rg/D");
    double rgs;
    MatchedTree->Branch("rgs", &rgs, "rgs/D");
    double zg;
    MatchedTree->Branch("zg", &zg, "zg/D");
    double zgs;
    MatchedTree->Branch("zgs", &zgs, "zgs/D");
    double mg;
    MatchedTree->Branch("mg", &mg, "mg/D");
    double mgs;
    MatchedTree->Branch("mgs", &mgs, "mgs/D");
    int n;
    MatchedTree->Branch("n", &n, "n/I");
    int ns;
    MatchedTree->Branch("ns", &ns, "ns/I");
    double w;
    MatchedTree->Branch("w", &w, "w/D");
    double ws;
    MatchedTree->Branch("ws", &ws, "ws/D");
    int reject;
    MatchedTree->Branch("reject", &reject, "reject/I");
    int rejects;
    MatchedTree->Branch("rejects", &rejects, "rejects/I");
    int mcevi;
    MatchedTree->Branch("mcevi", &mcevi, "mcevi/I");
    int mult;
    MatchedTree->Branch("mult", &mult, "mult/I");
    int mults;
    MatchedTree->Branch("mults", &mults, "mults/I");
    double nef;
    MatchedTree->Branch("nef", &nef, "nef/D");
    double nefs;
    MatchedTree->Branch("nefs", &nefs, "nefs/D");

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    int N = McChain->GetEntries();
    cout << "Number of Pythia events: " << N << endl;
    cout << "Number of Geant events:  " << PpChain->GetEntries() << endl;
    McChain->GetEntry(0); // this event contains a truth jet
    //McChain->GetEntry(2); // for small sample QA, this event contains a truth jet
    TStarJetVectorJet *dummyjet = (TStarJetVectorJet *)McJets->At(0);

    set<float> mcpt_list;
    set<float> pppt_list;

    for (Long64_t mcEvi = 0; mcEvi < N; ++mcEvi) // event loop
    {
        if (!(mcEvi % 500000))
            cout << "Working on " << mcEvi << " / " << N << endl;
	
        McChain->GetEntry(mcEvi);

        if (mcreject > 0) // skip bad events with high weights
	{
	    continue;
	}
	if (abs(mcvz) > 30)
	{
	    continue;
	}
	if (mcpt_list.find(mcpttot) != mcpt_list.end())
        {
            continue; // some events have identical total particle pT
        }
        else
        {
            mcpt_list.insert(mcpttot);
        }

	vector<RootResultStruct> mcresult;

        for (int j = 0; j < mcnjets; ++j)
        {
            TStarJetVectorJet *mcjet = (TStarJetVectorJet *)McJets->At(j);
	    mcresult.push_back(RootResultStruct(*mcjet, mcjet->Pt(), mcjet->Eta(), mcjet->Rapidity(), mcjet->Phi(), mcjet->M(), mcn[j], mcq[j], mcrg[j], mczg[j], mcmg[j], mcnef[j], mcweight, mcreject, mcEvi, mcmult, mcvz));   
        } // end of mcjet loop

        int ppevi = PpChain->GetEntryNumberWithIndex(mcrunid, mceventid);
        PpChain->GetEntry(ppevi);
	vector<RootResultStruct> ppresult;

	if (ppevi >= 0)
	{
		if (rcreject > 0)  continue; // go to next MC event
        	if (abs(ppvz) > 30)  continue;
		if (pppt_list.find(pppttot) != pppt_list.end())	continue;
		else
        	{
            		pppt_list.insert(pppttot);
        	}

        	for (int j = 0; j < ppnjets; ++j)
        	{
            		TStarJetVectorJet *ppjet = (TStarJetVectorJet *)PpJets->At(j); 
           		ppresult.push_back(RootResultStruct(*ppjet, ppjet->Pt(), ppjet->Eta(), ppjet->Rapidity(), ppjet->Phi(), ppjet->M(), rcn[j], rcq[j], rcrg[j], rczg[j], rcmg[j], rcnef[j], ppweight, rcreject, mcEvi, rcmult, ppvz));
	    	}
	} // end of ppjet loop

        //! Sort them together
        vector<MatchedRootResultStruct> MatchedResult;
        if (mcresult.size() > 0)
        {
            for (vector<RootResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end();)
            {
                bool matched = false;
                if (ppresult.size() > 0)
                {
                    for (vector<RootResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end();)
                    {
                        Double_t deta = mcit->y - ppit->y;
                        Double_t dphi = TVector2::Phi_mpi_pi(mcit->phi - ppit->phi);
                        Double_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
                        if (dr < RCut)
                        {
                            MatchedResult.push_back(MatchedRootResultStruct(*mcit, *ppit));
                            ppit = ppresult.erase(ppit);
                            matched = true;
                            //h_pt1s_pt1->Fill(mcit->pt1, ppit->pt1, ppweight);
                            break;
                        }
                        else
                        {
                            ++ppit;
                        }
                    }
                }
                if (matched)
                {
                    mcit = mcresult.erase(mcit);
                }
                else
                {
                    MatchedResult.push_back(MatchedRootResultStruct(*mcit, RootResultStruct(*dummyjet, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, ppweight, rcreject, mcEvi, rcmult, ppvz)));
                    ++mcit;
                }
            }
        }
        for (vector<RootResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end();) // fakes
        {
            MatchedResult.push_back(MatchedRootResultStruct(RootResultStruct(*dummyjet, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, mcweight, mcreject, mcEvi, mcmult, mcvz), *ppit));
            ++ppit;
        }
        
        for (vector<MatchedRootResultStruct>::iterator res = MatchedResult.begin(); res != MatchedResult.end(); ++res)
        {
            pt = res->first.pt;
            pts = res->second.pt;
            m = res->first.m;
            ms = res->second.m;
            n = res->first.n;
            ns = res->second.n;
            q = res->first.q;
            qs = res->second.q;
            rg = res->first.rg;
            rgs = res->second.rg;
            zg = res->first.zg;
            zgs = res->second.zg;
            mg = res->first.mg;
            mgs = res->second.mg;
            w = res->first.weight;
            ws = res->second.weight;
            reject = res->first.reject;
            rejects = res->second.reject;
            mult = res->first.mult;
            mults = res->second.mult;
            nef = res->first.nef;
            nefs = res->second.nef;

            if (res->second.mcevi > res->first.mcevi)
            {
                mcevi = res->second.mcevi;
            }
            else
            {
                mcevi = res->first.mcevi;
            }
            MatchedTree->Fill();
        }
    }

    fout->Write();
    return 0;

}
