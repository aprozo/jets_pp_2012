/* @file RunppAna.cxx
   @author Raghav Kunnawalkam Elayavalli,Youqi Song and Alexandr Prozorov
*/

#include "TStarJetVectorJet.h"
#include "ppAnalysis.hh"
#include "ppParameters.hh"

#include "TObjString.h"
#include "TString.h"
#include <TBranch.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParameter.h>
#include <TRandom.h>

#include <algorithm>
#include <set>
#include <vector>

#include <climits>
#include <cmath>

#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include <exception>

using namespace std;
using namespace fastjet;
using namespace contrib;

// Mostly for run 12
bool readinbadrunlist(vector<int> &badrun, TString csvfile);

int main(int argc, const char **argv) {

  ppAnalysis *ppana = nullptr;
  try {
    ppana = new ppAnalysis(argc, argv);
  } catch (std::exception &e) {
    cerr << "Initialization failed with exception " << e.what() << endl;
    return -1;
  }

  if (ppana->InitChains() == false) {
    cerr << "Chain initialization failed" << endl;
    return -1;
  }

  // Get parameters we used
  // ----------------------
  const ppParameters pars = ppana->GetPars();

  // Explicitly choose bad tower list here
  // -------------------------------------
  // Otherwise too easy to hide somewhere and forget...
  shared_ptr<TStarJetPicoReader> pReader = ppana->GetpReader();

  if (pReader) {
    TStarJetPicoTowerCuts *towerCuts = pReader->GetTowerCuts();
    towerCuts->AddBadTowers("lists/Combined_pp200Y12_badtower_isaac.list");

    TString csvfile = "lists/pp200Y12_badrun_isaac.list";
    vector<int> badruns;
    if (readinbadrunlist(badruns, csvfile) == false) {
      cerr << "Problems reading bad run list" << endl;
      return -1;
    }
    pReader->AddMaskedRuns(badruns);
  }

  // File
  // --------------------
  TFile *fout = new TFile(pars.OutFileName, "RECREATE");

  TH1D *hEventCounter = new TH1D("hEventCounter", "Event Counter", 6, 0, 6);
  hEventCounter->GetXaxis()->SetBinLabel(1, "ALL");
  hEventCounter->GetXaxis()->SetBinLabel(2, "AFTER_VERTEX");
  hEventCounter->GetXaxis()->SetBinLabel(3, "JETSFOUND");
  hEventCounter->GetXaxis()->SetBinLabel(4, "NOJETS");
  hEventCounter->GetXaxis()->SetBinLabel(5, "NOCONSTS");
  hEventCounter->GetXaxis()->SetBinLabel(6, "NOTACCEPTED");

  TString csvfile = "lists/good_run_list.list";
  vector<int> goodruns;
  if (readinbadrunlist(goodruns, csvfile) == false) {
    cerr << "Problems reading good run list" << endl;
    return -1;
  }

  TH1D *hEventsRun = new TH1D("hEventsRun", "Events per run", goodruns.size(),
                              0, goodruns.size());

  // set names of bins to run numbers
  for (unsigned int i = 0; i < goodruns.size(); ++i) {
    hEventsRun->GetXaxis()->SetBinLabel(i + 1, Form("%d", goodruns[i]));
  }

  TH1D *hEventsRunBeforeVertex =
      (TH1D *)hEventsRun->Clone("hEventsRunBeforeVertex");

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
  float event_sum_pt;
  ResultTree->Branch("event_sum_pt", &event_sum_pt, "event_sum_pt/F");
  bool is_rejected;
  ResultTree->Branch("is_rejected", &is_rejected, "is_rejected/O");

  TClonesArray Jets("TStarJetVectorJet");
  ResultTree->Branch("Jets", &Jets);
  double neutral_fraction[1000];
  ResultTree->Branch("neutral_fraction", neutral_fraction,
                     "neutral_fraction[njets]/D");
  bool is_leading_charged[1000];
  ResultTree->Branch("is_leading_charged", is_leading_charged,
                     "is_leading_charged[njets]/O");

  bool trigger_match_JP2[1000];
  ResultTree->Branch("trigger_match_JP2", trigger_match_JP2,
                     "trigger_match_JP2[njets]/O");
  bool trigger_match_HT2[1000];
  ResultTree->Branch("trigger_match_HT2", trigger_match_HT2,
                     "trigger_match_HT2[njets]/O");

  double pt[1000];
  ResultTree->Branch("pt", pt, "pt[njets]/D");

  double ptLead[1000];
  ResultTree->Branch("ptLead", ptLead, "ptLead[njets]/D");

  int n_constituents[1000];
  ResultTree->Branch("n_constituents", n_constituents,
                     "n_constituents[njets]/I");
  int index[1000];
  ResultTree->Branch("index", index, "index[njets]/I");

  // Helpers
  TStarJetVector *sv;

  // Go through events
  // -----------------
  cout << "Running analysis" << endl;
  try {
    bool ContinueReading = true;

    while (ContinueReading) {

      Jets.Clear();
      EVENTRESULT ret = ppana->RunEvent(); // event observables reset here

      // Understand what happened in the event
      switch (ret) {
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
      hEventCounter->Fill("ALL", 1);
      runid1 = ppana->GetRunid1();

 
      hEventsRunBeforeVertex->Fill(Form("%i", runid1), 1);
      if (fabs(vz) > 30) {
        continue;
      }
      hEventCounter->Fill("AFTER_VERTEX", 1);

      if (ret == EVENTRESULT::JETSFOUND) {
        hEventCounter->Fill("JETSFOUND", 1);
      } else if (ret == EVENTRESULT::NOCONSTS) {
        hEventCounter->Fill("NOCONSTS", 1);
      } else if (ret == EVENTRESULT::NOJETS) {
        hEventCounter->Fill("NOJETS", 1);
      } else if (ret == EVENTRESULT::NOTACCEPTED) {
        hEventCounter->Fill("NOTACCEPTED", 1);
      }

      // Now we can pull out details and results
      // ---------------------------------------

      weight = ppana->GetEventWeight();
      refmult = ppana->GetRefmult();
      runid = ppana->GetRunid();
      eventid = ppana->GetEventid();
      mult = ppana->GetEventMult();
      event_sum_pt = ppana->GetEventSumPt();
      is_rejected = ppana->GetRejectCode();

      vector<ResultStruct> Result = ppana->GetResult();
      njets = Result.size();
      // fill the bin with the same run name
      hEventsRun->Fill(Form("%i", runid1), 1);

      if (njets == 0)
        continue;

      int ijet = 0;
      for (auto &gr : Result) {
        TStarJetVector sv = TStarJetVector(MakeTLorentzVector(gr.orig));
        new (Jets[ijet]) TStarJetVectorJet(sv);
        neutral_fraction[ijet] =
            gr.orig.user_info<JetAnalysisUserInfo>().GetNumber();
        trigger_match_JP2[ijet] =
            gr.orig.user_info<JetAnalysisUserInfo>().IsMatchedJP();
        trigger_match_HT2[ijet] =
            gr.orig.user_info<JetAnalysisUserInfo>().IsMatchedHT();

        vector<PseudoJet> constituents =
            sorted_by_pt(gr.orig.constituents()); // sort by pt
        ptLead[ijet] = constituents[0].pt();
        pt[ijet] = gr.orig.perp();
        n_constituents[ijet] = gr.orig.constituents().size();
        index[ijet] = ijet;

        ijet++;
      }

      ResultTree->Fill();
    }
  } catch (std::string &s) {
    cerr << "RunEvent failed with string " << s << endl;
    return -1;
  } catch (std::exception &e) {
    cerr << "RunEvent failed with exception " << e.what() << endl;
    return -1;
  }

  // Save the output

  fout->Write();

  if (ppana->QA_pt_sDCAxy_pos) {
    ppana->QA_pt_sDCAxy_pos->Write();
  }
  if (ppana->QA_pt_sDCAxy_neg) {
    ppana->QA_pt_sDCAxy_neg->Write();
  }
  if (ppana->QA_jetpt_TowerID) {
    ppana->QA_jetpt_TowerID->Write();
  }
  if (ppana->QA_highjetpt_leadingtower_pt) {
    ppana->QA_highjetpt_leadingtower_pt->Write();
  }
  if (ppana->QA_highjetpt_leadingtrack_pt) {
    ppana->QA_highjetpt_leadingtrack_pt->Write();
  }

  if (ppana->QA_vx) {
    ppana->QA_vx->Write();
  }
  if (ppana->QA_vy) {
    ppana->QA_vy->Write();
  }
  if (ppana->QA_vz) {
    ppana->QA_vz->Write();
  }

  cout << "Done." << endl;

  delete ppana;
  return 0;
}

//----------------------------------------------------------------------
bool readinbadrunlist(vector<int> &badrun, TString csvfile) {

  // open infile
  std::string line;
  std::ifstream inFile(csvfile);

  std::cout << "Loading bad run id from " << csvfile.Data() << std::endl;
  ;

  if (!inFile.good()) {
    std::cout << "Can't open " << csvfile.Data() << std::endl;
    return false;
  }

  while (std::getline(inFile, line)) {
    if (line.size() == 0)
      continue; // skip empty lines
    if (line[0] == '#')
      continue; // skip comments

    std::istringstream ss(line);
    while (ss) {
      std::string entry;
      std::getline(ss, entry, ',');
      int ientry = atoi(entry.c_str());
      if (ientry) {
        badrun.push_back(ientry);
      }
    }
  }

  return true;
}
