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

int main(int argc, const char **argv)
{

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

   for (int i = 0; i < argc; ++i) {
      cout << "argv[" << i << "]=" << argv[i] << endl;
   }

   // Explicitly choose bad tower list here
   // -------------------------------------
   // Otherwise too easy to hide somewhere and forget...
   shared_ptr<TStarJetPicoReader> pReader = ppana->GetpReader();

   if (pReader) {
      TStarJetPicoTowerCuts *towerCuts = pReader->GetTowerCuts();
      // Full ~200-tower isaac hot/dead mask — approximates Dmitry's
      // StjTowerEnergyCutBemcStatus(1) DB mask, which he applies IN ADDITION
      // to the explicit 3407 (StjTowerEnergyCutTowerId). The 2026-06-12
      // 3407-only production PROVED the mask is load-bearing: without it,
      // hot-tower fakes (absent in embedding) produce rising tails — JP2
      // 1.41 at 40.6 / 1.95 at 48 GeV, HT2 (tower-triggered) 1.56 / 2.33.
      // ("lists/badtower_3407.list" kept for the minimal-mask systematic.)
      towerCuts->AddBadTowers("lists/Combined_pp200Y12_badtower_isaac.list");

      // Bad-run masking disabled: accept all runs in the pico list.
      // Dmitry-aligned filtering is applied downstream in
      // cross_section.cpp via dmitry_extras.list.
   }

   // File
   // --------------------
   TFile *fout = new TFile(pars.OutFileName, "RECREATE");

   TH1D *hEventCounter = new TH1D("hEventCounter", "Event Counter", 5, 0, 5);
   hEventCounter->GetXaxis()->SetBinLabel(1, "ALL");
   hEventCounter->GetXaxis()->SetBinLabel(2, "NOTACCEPTED");
   hEventCounter->GetXaxis()->SetBinLabel(3, "JETSFOUND");
   hEventCounter->GetXaxis()->SetBinLabel(4, "NOJETS");
   hEventCounter->GetXaxis()->SetBinLabel(5, "NOCONSTS");

   // Good-run list disabled: hEventsRun auto-grows as new run labels
   // appear via Fill(label, w). Bins/labels populated at runtime.
   TH1D *hEventsRun = new TH1D("hEventsRun", "Events per run", 1, 0, 1);
   hEventsRun->SetCanExtend(TH1::kAllAxes);

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
   bool isTriggerEvent;
   ResultTree->Branch("isTriggerEvent", &isTriggerEvent, "isTriggerEvent/O");
   // Per-event JP fire flags (real prescale-accepted hardware bits in data).
   // Basis of the independent-trigger analysis: filter fired_<T>, normalize by
   // the trigger's SAMPLED lumi (JP1 7.122, JP0 0.1486 pb^-1).
   bool fired_JP0, fired_JP1, fired_JP2;
   ResultTree->Branch("fired_JP0", &fired_JP0, "fired_JP0/O");
   ResultTree->Branch("fired_JP1", &fired_JP1, "fired_JP1/O");
   ResultTree->Branch("fired_JP2", &fired_JP2, "fired_JP2/O");
   double weight;
   ResultTree->Branch("weight", &weight, "weight/D");
   double refmult;
   ResultTree->Branch("refmult", &refmult, "refmult/D");
   int njets;
   ResultTree->Branch("njets", &njets, "njets/I");
   int mult;
   ResultTree->Branch("mult", &mult, "mult/I");
   // Primary vertex z (cm).  Needed downstream for vertex-z reweighting,
   // i.e. the user-side analogue of Dmitry's
   // SetVertexReweightingParams("embedding_tpc_vertex", ...).
   double vz;
   ResultTree->Branch("vz", &vz, "vz/D");
   float event_sum_pt;
   ResultTree->Branch("event_sum_pt", &event_sum_pt, "event_sum_pt/F");

   TClonesArray Jets("TStarJetVectorJet");
   ResultTree->Branch("Jets", &Jets);
   double neutral_fraction[1000];
   ResultTree->Branch("neutral_fraction", neutral_fraction, "neutral_fraction[njets]/D");
   int lead_tower_id[1000];
   ResultTree->Branch("lead_tower_id", lead_tower_id, "lead_tower_id[njets]/I");
   // Near-max JP-patch ADC in the jet's box (-1 if none). QA/provenance probe:
   // degenerate-data signature = trigger_match_JP2 && jp_patch_adc<=36 events.
   int jp_patch_adc[1000];
   ResultTree->Branch("jp_patch_adc", jp_patch_adc, "jp_patch_adc[njets]/I");

   // IsMatchedJP() returns the JP-match for whichever JP-trigger family was
   // selected via the -trig flag (JP0/JP1/JP2 — see match_jp dispatch on
   // pars.TriggerName.Contains in ppAnalysis.cxx). Expose the same data
   // under the three explicit branch names so downstream (matching.cpp,
   // unfolding/, cross_section.cpp) can filter via
   //   reco_trigger_match_<TRIG>
   // regardless of which JP-family was active. All three branches carry
   // identical content per run; the trigger map (-trig JP1 vs JP2) sets
   // WHICH JP family is the active match.
   bool trigger_match_JP2[1000];
   ResultTree->Branch("trigger_match_JP2", trigger_match_JP2, "trigger_match_JP2[njets]/O");
   bool trigger_match_JP1[1000];
   ResultTree->Branch("trigger_match_JP1", trigger_match_JP1, "trigger_match_JP1[njets]/O");
   bool trigger_match_JP0[1000];
   ResultTree->Branch("trigger_match_JP0", trigger_match_JP0, "trigger_match_JP0[njets]/O");
   bool trigger_match_HT2[1000];
   ResultTree->Branch("trigger_match_HT2", trigger_match_HT2, "trigger_match_HT2[njets]/O");
   double pt[1000];
   ResultTree->Branch("pt", pt, "pt[njets]/D");

   double ptLead[1000];
   ResultTree->Branch("ptLead", ptLead, "ptLead[njets]/D");

   int n_constituents[1000];
   ResultTree->Branch("n_constituents", n_constituents, "n_constituents[njets]/I");
   int index[1000];
   ResultTree->Branch("index", index, "index[njets]/I");

   // Jet area + UE density via off-axis cones (R= -R) + UE-subtracted pT.
   //   pt_corrected = pt - bg_density * jet_area
   double jet_area[1000];
   ResultTree->Branch("jet_area", jet_area, "jet_area[njets]/D");
   double bg_density[1000];
   ResultTree->Branch("bg_density", bg_density, "bg_density[njets]/D");
   double pt_corrected[1000];
   ResultTree->Branch("pt_corrected", pt_corrected, "pt_corrected[njets]/D");

   // Detector-frame jet eta (Pibero / Dmitry's StJetCandidate::detEta).
   // detEta is the jet eta projected from origin (0,0,0) to BEMC radius 225.405 cm.
   // Required for the detector-level |detEta| < 0.5 cut that mirrors
   // Dmitry's make_detector_level_cut. Formula: asinh(sinh(eta_phys) + vz/225.405).
   double det_eta[1000];
   ResultTree->Branch("det_eta", det_eta, "det_eta[njets]/D");

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

         if (ret == EVENTRESULT::NOTACCEPTED) {
            hEventCounter->Fill("NOTACCEPTED", 1);
         } else if (ret == EVENTRESULT::JETSFOUND) {
            hEventCounter->Fill("JETSFOUND", 1);
         } else if (ret == EVENTRESULT::NOCONSTS) {
            hEventCounter->Fill("NOCONSTS", 1);
         } else if (ret == EVENTRESULT::NOJETS) {
            hEventCounter->Fill("NOJETS", 1);
         }

         // Now we can pull out details and results
         // ---------------------------------------
         isTriggerEvent = ppana->IsTriggerEvent();
         fired_JP0 = ppana->FiredJP0();
         fired_JP1 = ppana->FiredJP1();
         fired_JP2 = ppana->FiredJP2();
         runid = ppana->GetRunid();
         // Luminosity bookkeeping: count events that were *actually analyzed for
         // the spectrum* — i.e. that pass the analysis trigger filter
         // (didFire AND shouldFire — see the `simu_fired` block in
         // ppAnalysis.cxx) and weren't otherwise rejected (Vz, ranking, ...).
         // The cross section is computed downstream as
         //     runLumi = (hEventsRun / lumi.root::nevents_<TRIG>) × lumi_<TRIG>
         // so hEventsRun must count events using the same definition that
         // nevents_<TRIG> does for STAR-recorded JP2 events. Hot-tower-only
         // triggers fire the hardware (didFire=true) but the simulator with
         // bad-tower mask disagrees (shouldFire=false); they're *recorded* in
         // nevents_<TRIG> but not analyzed for the spectrum, so they belong
         // in neither the numerator nor here.
         //
         // For MC (`intype == MCPICO`) there is no hardware trigger to filter
         // on, so we fall back to the original behaviour and fill hEventsRun
         // for every event reaching this point.
         const bool keep_for_lumi = (pars.intype == MCPICO)
                                  || (isTriggerEvent && ret != EVENTRESULT::NOTACCEPTED);
         if (keep_for_lumi) {
            hEventsRun->Fill(Form("%i", runid1), 1);
         }

         weight = ppana->GetEventWeight();
         refmult = ppana->GetRefmult();
         eventid = ppana->GetEventid();
         mult = ppana->GetEventMult();
         vz = ppana->GetVz();
         event_sum_pt = ppana->GetEventSumPt();

         // if (pars.InputName.Contains("hat") && pars.intype == INPICO &&
         //     ret == EVENTRESULT::NOTACCEPTED) { // fill events only for Geant only to account for missed jets
         //    ResultTree->Fill();
         //    continue;
         // }
         vector<ResultStruct> Result = ppana->GetResult();
         njets = Result.size();
         if (njets == 0)
            continue;

         int ijet = 0;
         for (auto &gr : Result) {
            TStarJetVector sv = TStarJetVector(MakeTLorentzVector(gr.orig));
            new (Jets[ijet]) TStarJetVectorJet(sv);
            neutral_fraction[ijet] = gr.orig.user_info<JetAnalysisUserInfo>().GetNumber();
            lead_tower_id[ijet]    = gr.orig.user_info<JetAnalysisUserInfo>().GetLeadTowerId();
            // Per-threshold JP-patch matches (isJP0/isJP1/isJP2 separately) — NOT
            // degenerate. Required for genuine per-trigger gates downstream.
            // (The old code copied one matched_jp flag into all three branches —
            // that degeneracy inflated the JP2 trigger efficiency ~2.8x and bent
            // the spectrum; see CLAUDE.md provenance gotchas.)
            const auto &ui = gr.orig.user_info<JetAnalysisUserInfo>();
            trigger_match_JP2[ijet] = ui.IsMatchedJP2();
            trigger_match_JP1[ijet] = ui.IsMatchedJP1();
            trigger_match_JP0[ijet] = ui.IsMatchedJP0();
            trigger_match_HT2[ijet] = ui.IsMatchedHT();
            jp_patch_adc[ijet]      = ui.GetJpAdc();

            vector<PseudoJet> constituents = sorted_by_pt(gr.orig.constituents()); // sort by pt
            ptLead[ijet] = constituents[0].pt();
            pt[ijet] = gr.orig.perp();
            n_constituents[ijet] = gr.orig.constituents().size();
            index[ijet] = ijet;
            jet_area[ijet]     = gr.area;
            bg_density[ijet]   = gr.bg_density;
            pt_corrected[ijet] = gr.pt_corrected;
            // detector-frame jet eta (BEMC front-face projection from origin):
            //   detEta = asinh(sinh(eta_phys) + vz / 225.405)
            // mirrors StJetCandidate::detEta(vertex) at BEMC_RADIUS = 225.405 cm.
            det_eta[ijet] = std::asinh(std::sinh(gr.orig.eta()) + vz / 225.405);
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

   ppana->GetHistogramManager().Write(fout, "QA_histograms");

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
