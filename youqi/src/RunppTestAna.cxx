/* @file RunppTestAna.cxx
   @author Raghav Kunnawalkam Elayavalli and Youqi Song
*/

#include "ppTestParameters.hh"
#include "ppTestAnalysis.hh"
#include "TStarJetVectorJet.h"

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TFile.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TParameter.h>
#include "TString.h"
#include "TObjString.h"

#include <set>
#include <vector>
#include <algorithm>

#include <cmath>
#include <climits>

#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include <exception>

using namespace std;
using namespace fastjet;
using namespace contrib;

// Mostly for run 12
bool readinbadrunlist(vector<int> &badrun, TString csvfile);

int main(int argc, const char **argv)
{

  ppTestAnalysis *ppana = nullptr;
  try
  {
    ppana = new ppTestAnalysis(argc, argv);
  }
  catch (std::exception &e)
  {
    cerr << "Initialization failed with exception " << e.what() << endl;
    return -1;
  }

  if (ppana->InitChains() == false)
  {
    cerr << "Chain initialization failed" << endl;
    return -1;
  }

  // Get parameters we used
  // ----------------------
  const ppTestParameters pars = ppana->GetPars();

  // Explicitly choose bad tower list here
  // -------------------------------------
  // Otherwise too easy to hide somewhere and forget...
  shared_ptr<TStarJetPicoReader> pReader = ppana->GetpReader();

  if (pReader)
  {
    TStarJetPicoTowerCuts *towerCuts = pReader->GetTowerCuts();
    towerCuts->AddBadTowers("/gpfs01/star/pwg/youqi/run12/final/badtower_isaac.list");
  }

  // Explicitly add bad run list here
  // --------------------------------
  if (pReader)
  {
    TString csvfile = "/gpfs01/star/pwg/youqi/run12/final/badrun_isaac.list";
    vector<int> badruns;
    if (readinbadrunlist(badruns, csvfile) == false)
    {
      cerr << "Problems reading bad run list" << endl;
      return -1;
    }
    pReader->AddMaskedRuns(badruns);
  }

  // File
  // --------------------
  TFile *fout = new TFile(pars.OutFileName, "RECREATE");

  // Save results
  // ------------
  TTree *ResultTree = new TTree("ResultTree", "Result Jets");

  // Give each event a unique ID to compare event by event with different runs
  int runid;
  ResultTree->Branch("runid", &runid, "runid/I");
  int runid1;
  ResultTree->Branch("runid1", &runid1, "runid1/I");
  int eventid;
  ResultTree->Branch("eventid", &eventid, "eventid/I");
  double weight;
  ResultTree->Branch("weight", &weight, "weight/D");
  double refmult;
  ResultTree->Branch("refmult", &refmult, "refmult/D");
  int njets;
  ResultTree->Branch("njets", &njets, "njets/I");
  double vz;
  ResultTree->Branch("vz", &vz, "vz/D");
  int mult;
  ResultTree->Branch("mult", &mult, "mult/I");
  float pttot;
  ResultTree->Branch("pttot", &pttot, "pttot/F");
  int reject;
  ResultTree->Branch("reject", &reject, "reject/I");

  TClonesArray Jets("TStarJetVectorJet");
  ResultTree->Branch("Jets", &Jets);
  double nef[1000];
  ResultTree->Branch("nef", nef, "nef[njets]/D");
  double pt[1000];
  ResultTree->Branch("pt", pt, "pt[njets]/D");
  double m[1000];
  ResultTree->Branch("m", m, "m[njets]/D");
  double q[1000];
  ResultTree->Branch("q", q, "q[njets]/D");
  double rg[1000];
  ResultTree->Branch("rg", rg, "rg[njets]/D");
  double zg[1000];
  ResultTree->Branch("zg", zg, "zg[njets]/D");
  double mg[1000];
  ResultTree->Branch("mg", mg, "mg[njets]/D");
  int n[1000];
  ResultTree->Branch("n", n, "n[njets]/I");
  int index[1000];
  ResultTree->Branch("index", index, "index[njets]/I");
  
  // Helpers
  TStarJetVector *sv;

  // Go through events
  // -----------------
  cout << "Running analysis" << endl;
  try
  {
    bool ContinueReading = true;

    while (ContinueReading)
    {

      Jets.Clear();
      //refmult = 0;
      //runid = -(INT_MAX - 1);
      //eventid = -(INT_MAX - 1);
      
      EVENTRESULT ret = ppana->RunEvent(); // event observables reset here

      // Understand what happened in the event
      switch (ret)
      {
      case EVENTRESULT::PROBLEM:
        cerr << "Encountered a serious issue" << endl;
        return -1;
        break;
      case EVENTRESULT::ENDOFINPUT:
        cout << "End of Input" << endl;
        ContinueReading = false;
        continue;
        break;
      case EVENTRESULT::NOTACCEPTED:
        // continue;
        break;
      case EVENTRESULT::NOCONSTS:
        // cout << "Event empty." << endl;
        break;
      case EVENTRESULT::NOJETS:
        // cout << "No jets found." << endl;
        break;
      case EVENTRESULT::JETSFOUND:
        // The only way not to break out or go back to the top
        // cout << "Jets found." << endl;
        break;
      default:
        cerr << "Unknown return value." << endl;
        return -1;
        break;
      }

      // Now we can pull out details and results
      // ---------------------------------------
      
      weight = ppana->GetEventWeight();
      refmult = ppana->GetRefmult();
      runid = ppana->GetRunid();
      runid1 = ppana->GetRunid1();
      eventid = ppana->GetEventid();
      vz = ppana->GetVz();
      //highw = ppana->CheckHighW();
      mult = ppana->GetEventMult();
      pttot = ppana->GetPttot();
      reject = ppana->GetRejectCode();
      
      vector<ResultStruct> Result = ppana->GetResult();
      njets = Result.size();

      if (njets == 0)
      {
        ResultTree->Fill();
        continue;
      }

      int ijet = 0;
      for (auto &gr : Result)
      {

        TStarJetVector sv = TStarJetVector(MakeTLorentzVector(gr.orig));
        sv.SetCharge(gr.orig.user_info<JetAnalysisUserInfo>().GetQuarkCharge() / 3);
        new (Jets[ijet]) TStarJetVectorJet(sv);
        
	nef[ijet] = gr.orig.user_info<JetAnalysisUserInfo>().GetNumber();
        pt[ijet] = gr.orig.perp();
        m[ijet] = gr.orig.m();
        q[ijet] = gr.q;
        rg[ijet] = gr.sd.structure_of<fastjet::contrib::SoftDrop>().delta_R();
        zg[ijet] = gr.sd.structure_of<fastjet::contrib::SoftDrop>().symmetry();
        mg[ijet] = gr.sd.m();
        n[ijet] = gr.orig.constituents().size();
        index[ijet] = ijet;
        
        ijet++;
      }

      ResultTree->Fill();
    }
  }
  catch (std::string &s)
  {
    cerr << "RunEvent failed with string " << s << endl;
    return -1;
  }
  catch (std::exception &e)
  {
    cerr << "RunEvent failed with exception " << e.what() << endl;
    return -1;
  }

  fout->Write();

  cout << "Done." << endl;

  delete ppana;
  return 0;
}

//----------------------------------------------------------------------
bool readinbadrunlist(vector<int> &badrun, TString csvfile)
{

  // open infile
  std::string line;
  std::ifstream inFile(csvfile);

  std::cout << "Loading bad run id from " << csvfile.Data() << std::endl;
  ;

  if (!inFile.good())
  {
    std::cout << "Can't open " << csvfile.Data() << std::endl;
    return false;
  }

  while (std::getline(inFile, line))
  {
    if (line.size() == 0)
      continue; // skip empty lines
    if (line[0] == '#')
      continue; // skip comments

    std::istringstream ss(line);
    while (ss)
    {
      std::string entry;
      std::getline(ss, entry, ',');
      int ientry = atoi(entry.c_str());
      if (ientry)
      {
        badrun.push_back(ientry);
      }
    }
  }

  return true;
}
