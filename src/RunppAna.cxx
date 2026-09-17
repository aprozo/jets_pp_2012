#include "TStarJetVectorJet.h"
#include "ppAnalysis.hh"
#include "ppParameters.hh"

#include "TString.h"
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include <cmath>
#include <exception>
#include <stdexcept>
#include <vector>

using namespace std;
using namespace fastjet;
using namespace contrib;

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

   // The bad-tower list is chosen explicitly here, not hidden in a parameter.
   shared_ptr<TStarJetPicoReader> pReader = ppana->GetpReader();

   if (pReader) {
      TStarJetPicoTowerCuts *towerCuts = pReader->GetTowerCuts();
      // Tower mask: the DB tower status is applied at Stage-0 (pico production);
      // the only explicit veto here is tower 3407 (hot in Run 12).
      towerCuts->AddBadTowers("lists/badtower_3407.list");
      std::cout << "Bad tower list: lists/badtower_3407.list" << std::endl;

      // Bad-run masking disabled: accept all runs in the pico list.
   }

   TFile *fout = new TFile(pars.OutFileName, "RECREATE");

   TH1D *hEventCounter = new TH1D("hEventCounter", "Event Counter", 5, 0, 5);
   hEventCounter->GetXaxis()->SetBinLabel(1, "ALL");
   hEventCounter->GetXaxis()->SetBinLabel(2, "NOTACCEPTED");
   hEventCounter->GetXaxis()->SetBinLabel(3, "JETSFOUND");
   hEventCounter->GetXaxis()->SetBinLabel(4, "NOJETS");
   hEventCounter->GetXaxis()->SetBinLabel(5, "NOCONSTS");

   // No good-run list: hEventsRun auto-grows as new run labels appear via
   // Fill(label, w), so bins and labels are populated at runtime.
   TH1D *hEventsRun = new TH1D("hEventsRun", "Events per run", 1, 0, 1);
   hEventsRun->SetCanExtend(TH1::kAllAxes);

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
   bool fired_JP0, fired_JP1, fired_JP2, fired_HT2, fired_MB, should_HT2;
   ResultTree->Branch("fired_JP0", &fired_JP0, "fired_JP0/O");
   ResultTree->Branch("fired_JP1", &fired_JP1, "fired_JP1/O");
   ResultTree->Branch("fired_JP2", &fired_JP2, "fired_JP2/O");
   ResultTree->Branch("fired_HT2", &fired_HT2, "fired_HT2/O");
   ResultTree->Branch("fired_MB", &fired_MB, "fired_MB/O");
   ResultTree->Branch("should_HT2", &should_HT2, "should_HT2/O");
   bool should_JP0, should_JP1, should_JP2; // offline-emulator decisions
   ResultTree->Branch("should_JP0", &should_JP0, "should_JP0/O");
   ResultTree->Branch("should_JP1", &should_JP1, "should_JP1/O");
   ResultTree->Branch("should_JP2", &should_JP2, "should_JP2/O");
   bool shouldHW_JP0, shouldHW_JP1, shouldHW_JP2; // kOnline (bit-7) decisions
   ResultTree->Branch("shouldHW_JP0", &shouldHW_JP0, "shouldHW_JP0/O");
   ResultTree->Branch("shouldHW_JP1", &shouldHW_JP1, "shouldHW_JP1/O");
   ResultTree->Branch("shouldHW_JP2", &shouldHW_JP2, "shouldHW_JP2/O");
   double weight;
   ResultTree->Branch("weight", &weight, "weight/D");
   double refmult;
   ResultTree->Branch("refmult", &refmult, "refmult/D");
   int njets;
   ResultTree->Branch("njets", &njets, "njets/I");
   int mult;
   ResultTree->Branch("mult", &mult, "mult/I");
   // Primary vertex z (cm); needed downstream for vertex-z reweighting.
   double vz;
   ResultTree->Branch("vz", &vz, "vz/D");
   double pthat = -1;
   ResultTree->Branch("pthat", &pthat, "pthat/D"); // event pt-hat from the MC pico header (-1 data)
   // VPD vertex z (the |vz_vpd - vz| < 6 cut is MB-only) and transverse vertex.
   double vx, vy, vz_vpd;
   ResultTree->Branch("vx", &vx, "vx/D");
   ResultTree->Branch("vy", &vy, "vy/D");
   ResultTree->Branch("vz_vpd", &vz_vpd, "vz_vpd/D");
   // VPD east/west hit counts: VPDMB-fired proxy is n_east>=1 && n_west>=1.
   int n_vpd_east, n_vpd_west;
   ResultTree->Branch("n_vpd_east", &n_vpd_east, "n_vpd_east/I");
   ResultTree->Branch("n_vpd_west", &n_vpd_west, "n_vpd_west/I");
   float event_sum_pt;
   ResultTree->Branch("event_sum_pt", &event_sum_pt, "event_sum_pt/F");
   // Simulator ADCs (DSM scale) and thresholds: event maxima here, per-jet below.
   int evt_jp_adc_max, evt_ht_adc_max, jp_thr[3], ht_thr[4];
   ResultTree->Branch("evt_jp_adc_max", &evt_jp_adc_max, "evt_jp_adc_max/I");
   ResultTree->Branch("evt_ht_adc_max", &evt_ht_adc_max, "evt_ht_adc_max/I");
   ResultTree->Branch("jp_thr", jp_thr, "jp_thr[3]/I");
   ResultTree->Branch("ht_thr", ht_thr, "ht_thr[4]/I");

   TClonesArray Jets("TStarJetVectorJet");
   ResultTree->Branch("Jets", &Jets);
   // Fixed capacity of the per-jet stack arrays below; njets is guarded
   // against it after GetResult() (never reached at PtJetMin=5).
   constexpr int kMaxJetsPerEvent = 1000;
   double neutral_fraction[kMaxJetsPerEvent];
   ResultTree->Branch("neutral_fraction", neutral_fraction, "neutral_fraction[njets]/D");
   int lead_tower_id[kMaxJetsPerEvent];
   ResultTree->Branch("lead_tower_id", lead_tower_id, "lead_tower_id[njets]/I");
   // Near-max JP-patch ADC in the jet's box (-1 if none), a QA probe.
   int jp_patch_adc[kMaxJetsPerEvent];
   ResultTree->Branch("jp_patch_adc", jp_patch_adc, "jp_patch_adc[njets]/I");
   int ht_adc_max[kMaxJetsPerEvent];
   ResultTree->Branch("ht_adc_max", ht_adc_max, "ht_adc_max[njets]/I");

   bool trigger_match_JP2[kMaxJetsPerEvent];
   ResultTree->Branch("trigger_match_JP2", trigger_match_JP2, "trigger_match_JP2[njets]/O");
   bool trigger_match_JP1[kMaxJetsPerEvent];
   ResultTree->Branch("trigger_match_JP1", trigger_match_JP1, "trigger_match_JP1[njets]/O");
   bool trigger_match_JP0[kMaxJetsPerEvent];
   ResultTree->Branch("trigger_match_JP0", trigger_match_JP0, "trigger_match_JP0[njets]/O");
   bool trigger_match_HT2[kMaxJetsPerEvent];
   ResultTree->Branch("trigger_match_HT2", trigger_match_HT2, "trigger_match_HT2[njets]/O");
   double pt[kMaxJetsPerEvent];
   ResultTree->Branch("pt", pt, "pt[njets]/D");

   double ptLead[kMaxJetsPerEvent];
   ResultTree->Branch("ptLead", ptLead, "ptLead[njets]/D");

   int n_constituents[kMaxJetsPerEvent];
   ResultTree->Branch("n_constituents", n_constituents, "n_constituents[njets]/I");
   int index[kMaxJetsPerEvent];
   ResultTree->Branch("index", index, "index[njets]/I");

   // Jet area, off-axis-cones UE density, and pt_corrected = pt - bg_density * jet_area
   double jet_area[kMaxJetsPerEvent];
   ResultTree->Branch("jet_area", jet_area, "jet_area[njets]/D");
   double bg_density[kMaxJetsPerEvent];
   ResultTree->Branch("bg_density", bg_density, "bg_density[njets]/D");
   double pt_corrected[kMaxJetsPerEvent];
   ResultTree->Branch("pt_corrected", pt_corrected, "pt_corrected[njets]/D");

   // Detector-frame jet eta: the jet eta projected from the origin onto the BEMC
   // front face, asinh(sinh(eta_phys) + vz/225.405) with R_BEMC = 225.405 cm.
   // The detector-level |det_eta| < 0.5 acceptance cut is applied downstream.
   double det_eta[kMaxJetsPerEvent];
   ResultTree->Branch("det_eta", det_eta, "det_eta[njets]/D");

   cout << "Running analysis" << endl;
   try {
      bool ContinueReading = true;

      while (ContinueReading) {

         Jets.Clear();
         EVENTRESULT ret = ppana->RunEvent(); // event observables reset here

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
            break;
         case EVENTRESULT::NOCONSTS:
            break;
         case EVENTRESULT::NOJETS:
            break;
         case EVENTRESULT::JETSFOUND:
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

         isTriggerEvent = ppana->IsTriggerEvent();
         fired_JP0 = ppana->FiredJP0();
         fired_JP1 = ppana->FiredJP1();
         fired_JP2 = ppana->FiredJP2();
         fired_HT2 = ppana->FiredHT2();
         fired_MB = ppana->FiredMB();
         should_HT2 = ppana->ShouldHT2();
         should_JP0 = ppana->ShouldJP0();
         should_JP1 = ppana->ShouldJP1();
         should_JP2 = ppana->ShouldJP2();
         shouldHW_JP0 = ppana->ShouldHwJP0();
         shouldHW_JP1 = ppana->ShouldHwJP1();
         shouldHW_JP2 = ppana->ShouldHwJP2();
         runid = ppana->GetRunid();
         // Luminosity bookkeeping: hEventsRun counts events actually analyzed for
         // the spectrum — those passing the trigger filter (didFire AND shouldFire,
         // see the `simu_fired` block in ppAnalysis.cxx) and not otherwise rejected
         // (Vz, ranking, ...). The cross section is computed downstream as
         //     runLumi = (hEventsRun / lumi.root::nevents_<TRIG>) × lumi_<TRIG>
         // so hEventsRun must use the same event definition as nevents_<TRIG>.
         // In MC there is no hardware trigger to filter on, so every event
         // reaching this point is counted.
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
         pthat = ppana->GetPthat();
         vx = ppana->GetVx();
         vy = ppana->GetVy();
         vz_vpd = ppana->GetVpdVz();
         n_vpd_east = ppana->GetNVpdEast();
         n_vpd_west = ppana->GetNVpdWest();
         event_sum_pt = ppana->GetEventSumPt();
         evt_jp_adc_max = ppana->GetJpAdcMax();
         evt_ht_adc_max = ppana->GetHtAdcMax();
         for (int i = 0; i < 3; ++i) jp_thr[i] = ppana->GetJpThr(i);
         for (int i = 0; i < 4; ++i) ht_thr[i] = ppana->GetHtThr(i);

         vector<ResultStruct> Result = ppana->GetResult();
         njets = Result.size();
         if (njets > kMaxJetsPerEvent)
            throw std::runtime_error("njets exceeds the kMaxJetsPerEvent per-jet array capacity");
         if (njets == 0)
            continue;

         int ijet = 0;
         for (auto &gr : Result) {
            TStarJetVector sv = TStarJetVector(MakeTLorentzVector(gr.orig));
            new (Jets[ijet]) TStarJetVectorJet(sv);
            neutral_fraction[ijet] = gr.orig.user_info<JetAnalysisUserInfo>().GetNumber();
            lead_tower_id[ijet]    = gr.orig.user_info<JetAnalysisUserInfo>().GetLeadTowerId();
            // Per-threshold JP-patch matches (JP0/JP1/JP2 independently), for
            // per-trigger gates downstream.
            const auto &ui = gr.orig.user_info<JetAnalysisUserInfo>();
            trigger_match_JP2[ijet] = ui.IsMatchedJP2();
            trigger_match_JP1[ijet] = ui.IsMatchedJP1();
            trigger_match_JP0[ijet] = ui.IsMatchedJP0();
            trigger_match_HT2[ijet] = ui.IsMatchedHT();
            jp_patch_adc[ijet]      = ui.GetJpAdc();
            ht_adc_max[ijet]        = ui.GetHtAdc();

            vector<PseudoJet> constituents = sorted_by_pt(gr.orig.constituents()); // sort by pt
            ptLead[ijet] = constituents[0].pt();
            pt[ijet] = gr.orig.perp();
            n_constituents[ijet] = gr.orig.constituents().size();
            index[ijet] = ijet;
            jet_area[ijet]     = gr.area;
            bg_density[ijet]   = gr.bg_density;
            pt_corrected[ijet] = gr.pt_corrected;
            // detector-frame jet eta, BEMC front-face projection from the origin
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

   fout->Write();

   ppana->GetHistogramManager().Write(fout, "QA_histograms");

   cout << "Done." << endl;

   delete ppana;
   return 0;
}
