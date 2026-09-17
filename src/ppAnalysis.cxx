#include <map>

#include "ppAnalysis.hh"
#include <climits>
#include <iostream>
#include <stdlib.h> // for atof, atoi
#include <string>

using std::cerr;
using std::cout;
using std::endl;

bool match_jp(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R, const TString &which = "JP2", double vz = 0.0);
bool match_ht(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R);
int ht_adc_max_in_jet(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers);
int jp_patch_adc_near(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R, double vz);
void setTriggerBitMap(TStarJetPicoTriggerInfo *trig, TStarJetPicoEventHeader *header);
bool getBarrelJetPatchEtaPhi(int jetPatch, float &eta, float &phi);

// Off-axis-cones UE density (GeV per unit (η,φ) area): two cones at
// (η_jet, φ_jet ± π/2) of radius R; ρ = avg(ΣpT) / (π R²).
double off_axis_cones_density(const PseudoJet &jet, const vector<PseudoJet> &particles, double R);

double getPythiaWeight(TString filename);

ppAnalysis::ppAnalysis(const int argc, const char **const argv)
{
   vector<string> arguments(argv + 1, argv + argc);
   bool argsokay = true;
   NEvents = -1;
   for (auto parg = arguments.begin(); parg != arguments.end(); ++parg) {
      string arg = *parg;
      if (arg == "-R") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.R = atof(parg->data());
      } else if (arg == "-lja") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.LargeJetAlgorithm = AlgoFromString(*parg);
      } else if (arg == "-pj") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.PtJetMin = atof((parg)->data());
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.PtJetMax = atof((parg)->data());
      } else if (arg == "-ec") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.EtaConsCut = atof((parg)->data());
      } else if (arg == "-pc") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.PtConsMin = atof((parg)->data());
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.PtConsMax = atof((parg)->data());
      } else if (arg == "-hadcorr") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.HadronicCorr = atof(parg->data());
      } else if (arg == "-o") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.OutFileName = *parg;
      } else if (arg == "-i") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.InputName = *parg;
      } else if (arg == "-c") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.ChainName = *parg;
      } else if (arg == "-trig") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.TriggerName = *parg;
      } else if (arg == "-intype") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         if (*parg == "pico") {
            pars.intype = INPICO;
            continue;
         }
         if (*parg == "mcpico") {
            pars.intype = MCPICO;
            continue;
         }
         argsokay = false;
         break;
      } else if (arg == "-N") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         NEvents = atoi(parg->data());
      } else if (arg == "-sdca") {
         // |signed dca_xy| flat cap (cm): removes low-pT tracks reconstructed as
         // high-pT (large transverse impact parameter) that fake high-pT jets.
         if (++parg == arguments.end()) { argsokay = false; break; }
         pars.sDCAxyCut = atof(parg->data());
         cout << "Setting |sDCAxy| cut to " << pars.sDCAxyCut << " cm" << endl;
      } else if (arg == "-dca") {
         // flat 3-D DCA cap (cm), applied by the reader (SetDCACut).
         if (++parg == arguments.end()) { argsokay = false; break; }
         pars.DcaCut = atof(parg->data());
         cout << "Setting |DCA| cut to " << pars.DcaCut << " cm" << endl;
      } else if (arg == "-leadsdca") {
         // Veto a jet whose LEADING charged constituent has |sDCAxy| above this
         // (cm) — single-fake-track jets. Default off (99999).
         if (++parg == arguments.end()) { argsokay = false; break; }
         pars.LeadTrackSdcaCut = atof(parg->data());
         cout << "Setting leading-track |sDCAxy| veto to " << pars.LeadTrackSdcaCut << " cm" << endl;
      } else if (arg == "-fakeeff") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.FakeEff = atof(parg->data());
         cout << "Setting fake efficiency to " << pars.FakeEff << endl;
         if (pars.FakeEff < 0 || pars.FakeEff > 1) {
            argsokay = false;
            break;
         }
      } else if (arg == "-towunc") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.IntTowScale = atoi(parg->data());
         pars.fTowScale = 1.0 + pars.IntTowScale * pars.fTowUnc;
         cout << "Setting tower scale to " << pars.fTowScale << endl;
         if (pars.IntTowScale < -1 || pars.IntTowScale > 1) {
            argsokay = false;
            break;
         }
      } else if (arg == "-geantnum") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         if (*parg != "0" && *parg != "1") {
            argsokay = false;
            break;
         }
         pars.UseGeantNumbering = bool(atoi(parg->data()));
      } else if (arg == "-jetnef") {
         if (++parg == arguments.end()) {
            argsokay = false;
            break;
         }
         pars.MaxJetNEF = atof(parg->data());
         cout << "Setting Max Jet NEF to " << pars.MaxJetNEF << endl;
         if (pars.MaxJetNEF < 0 || pars.MaxJetNEF > 1) {
            argsokay = false;
            break;
         }
      } else {
         argsokay = false;
         break;
      }
   }

   if (!argsokay) {
      cerr << "usage: " << argv[0] << endl
           << " [-o OutFileName]" << endl
           << " [-N Nevents (<0 for all)]" << endl
           << " [-R radius]" << endl
           << " [-lja LargeJetAlgorithm]" << endl
           << " [-i infilepattern]" << endl
           << " [-c chainname]" << endl
           << " [-intype pico|mcpico(embedding)]" << endl
           << " [-trig trigger name (e.g. HT)]" << endl
           << " [-pj PtJetMin PtJetMax]" << endl
           << " [-ec EtaConsCut]" << endl
           << " [-pc PtConsMin PtConsMax]" << endl
           << " [-hadcorr HadronicCorrection]  -- Set to a negative value for "
              "MIP correction."
           << endl
           << " [-psc PtSubConsMin PtSubConsMax]" << endl
           << " [-fakeeff (0..1)] -- enable fake efficiency for systematics. "
              "0.95 is a reasonable example."
           << endl
           << " [-towunc -1|0|1 ] -- Shift tower energy by this times " << pars.fTowUnc << endl
           << " [-geantnum true|false] (Force Geant run event id hack)" << endl
           << endl
           << endl
           << "NOTE: Wildcarded file patterns should be in single quotes." << endl
           << endl;
      throw std::runtime_error("Not a valid list of options");
   }

   if (pars.PtJetMin <= 0) {
      throw std::runtime_error("PtJetMin needs to be positive (0.001 will work).");
   }

   if (pars.intype == MCPICO) {
      // This refers to the JetTreeMc branch
      if (pars.ChainName != "JetTreeMc") {
         throw std::runtime_error("Unsuitable chain name for pars.intype==MCPICO");
      }
   }
   if (pars.ChainName == "JetTreeMc") {
      if (pars.intype != MCPICO) {
         throw std::runtime_error("Unsuitable chain name for pars.intype==MCPICO");
      }
   }

   // Jet-finding eta window, deliberately wider than the final acceptance: a jet
   // with |det_eta| < 0.5 at |vz| < 60 has |physics eta| < ~0.72, so 1.0 is safely
   // inclusive. The final |det_eta| < 0.5 detector cut is applied downstream on the
   // det_eta branch. The CONSTITUENT cut is pars.EtaConsCut (particle loop).
   EtaJetCut = 1.0;
   EtaGhostCut = EtaJetCut + 2.0 * pars.R;

   select_jet_eta = SelectorAbsEtaMax(EtaJetCut);
   select_jet_pt = SelectorPtRange(pars.PtJetMin, pars.PtJetMax);
   select_jet = select_jet_eta * select_jet_pt;

   pars.Recursive = pars.InputName.Contains("Pythia") && false;

   JetDef = JetDefinition(pars.LargeJetAlgorithm, pars.R);

   cout << " R = " << pars.R << endl;
   cout << " Original jet algorithm : " << pars.LargeJetAlgorithm << endl;
   cout << " PtJetMin = " << pars.PtJetMin << endl;
   cout << " PtJetMax = " << pars.PtJetMax << endl;
   cout << " PtConsMin = " << pars.PtConsMin << endl;
   cout << " PtConsMax = " << pars.PtConsMax << endl;
   cout << " Constituent eta cut = " << pars.EtaConsCut << endl;
   cout << " Jet eta cut = " << EtaJetCut << endl;
   cout << " Ghosts out to eta = " << EtaGhostCut << endl;
   cout << " Reading tree named \"" << pars.ChainName << "\" from " << pars.InputName << endl;
   cout << " intype = " << pars.intype << endl;
   cout << " Writing to " << pars.OutFileName << endl;
   cout << " ----------------------------" << endl;
}
//----------------------------------------------------------------------
ppAnalysis::~ppAnalysis()
{
   if (pJA) {
      delete pJA;
      pJA = 0;
   }
}

//----------------------------------------------------------------------
bool ppAnalysis::InitChains()
{
   Events = new TChain(pars.ChainName);
   Events->Add(pars.InputName);
   if (NEvents < 0)
      NEvents = INT_MAX;

   pReader = SetupReader(Events, pars);

   InitializeReader(pReader, pars.InputName, NEvents, PicoDebugLevel, pars.HadronicCorr);
   if (pars.intype == MCPICO) {
      TurnOffCuts(pReader);
      pars.MaxJetNEF = 1.0;
      pars.sDCAxyCut = 99999;
      pars.PtConsMin = 0.0;
      pars.PtConsMax = 99999.0;
      // No DCA / flag concept for particle-level MC tracks.
      pars.ApplyTdcaPtDep = false;
      pars.FlagMin = INT_MIN;
      // MC constituents carry the sDCAxy sentinel -999, so the leading-track
      // |sDCAxy| jet veto is meaningless at particle level and must be forced
      // off — otherwise fabs(-999) > cut vetoes every jet with a charged
      // constituent among the two leading ones.
      pars.LeadTrackSdcaCut = 99999;
   }

   QA_hist.Init();

   cout << "N = " << NEvents << endl;

   cout << "Done initializing chains. " << endl;
   return true;
}
//----------------------------------------------------------------------
// Main routine for one event.
EVENTRESULT ppAnalysis::RunEvent()
{
   // Reset results from the last event
   Result.clear();
   weight = 1.;
   mult = 0;
   refmult = 0;
   runid = -(INT_MAX - 1);
   runid1 = -(INT_MAX - 1);
   eventid = -(INT_MAX - 1);
   event_sum_pt = 0.;
   njets = 0;
   particles.clear();

   if (!pReader->NextEvent()) {
      pReader->PrintStatus();
      return EVENTRESULT::ENDOFINPUT;
   }
   pReader->PrintStatus(10);

   pFullEvent = pReader->GetOutputContainer()->GetArray();
   TStarJetPicoEventHeader *header = pReader->GetEvent()->GetHeader();

   set<int> event_triggers;
   for (int i = 0; i < header->GetNOfTriggerIds(); i++) {
      event_triggers.insert(header->GetTriggerId(i));
   }

   // pp 2012 200 GeV trigger ids
   map<TString, set<int>> trigger_map_2012;
   trigger_map_2012["JP2"] = {370621};
   trigger_map_2012["JP1"] = {370611};
   trigger_map_2012["JP0"] = {370601};
   trigger_map_2012["HT2"] = {370531, 500205}; // 500205: legacy id in some embedding trees
   trigger_map_2012["MB"] = {370011};

   // Per-event JP fire flags from the header trigger ids: in data these are the
   // real prescale-accepted hardware bits.
   firedJP0 = event_triggers.count(370601) != 0;
   firedJP1 = event_triggers.count(370611) != 0;
   firedJP2 = event_triggers.count(370621) != 0;
   firedHT2 = event_triggers.count(370531) != 0;
   firedMB  = event_triggers.count(370011) != 0 || event_triggers.count(370001) != 0;

   TString current_trigger = "";
   if (pars.TriggerName.Contains("JP2"))
      current_trigger = "JP2";
   else if (pars.TriggerName.Contains("JP1"))
      current_trigger = "JP1";
   else if (pars.TriggerName.Contains("JP0"))
      current_trigger = "JP0";
   else if (pars.TriggerName.Contains("HT2"))
      current_trigger = "HT2";
   else if (pars.TriggerName.Contains("MB"))
      current_trigger = "MB";

   isTriggerEvent = false;

   if (pars.intype == INPICO) {
      for (auto trig_id : trigger_map_2012[current_trigger]) {
         if (event_triggers.count(trig_id) != 0) {
            isTriggerEvent = true;
            break;
         }
      }
   }

   vector<TStarJetPicoTriggerInfo *> triggers;
   TStarJetPicoTowerCuts *towerCuts = pReader->GetTowerCuts();
   vector<int> HT2_trigger_ids;
   int count_bad_tower_triggers = 0;

   // Build the simu-side trigger-object list; count BHT2 triggers fired by
   // bad towers (for the hot-tower-only event veto below).
   shouldHT2 = false;
   haveSimuDecision = false;
   bool simuJP0 = false, simuJP1 = false, simuJP2 = false;
   for (int i = 0; i < header->GetNOfTrigObjs(); ++i) {
      auto trig = pReader->GetEvent()->GetTrigObj(i);
      // bit-8 objects carry the FULL StTriggerSimuMaker isTrigger() decision —
      // every detector leg, not just barrel patches (id = trigger id, ADC = 0/1).
      // They are decisions, not patches: record them and keep them out of the
      // patch-matching list.
      if (trig->GetBit(8)) {
         haveSimuDecision = true;
         const bool yes = trig->GetADC() == 1;
         if (trig->GetId() == 370601) simuJP0 = yes;
         if (trig->GetId() == 370611) simuJP1 = yes;
         if (trig->GetId() == 370621) simuJP2 = yes;
         if (trig->GetId() == 370531) shouldHT2 = yes;
         continue;
      }
      if (trig->GetTriggerFlag() == 5 && trig->GetADC() == 0)
         continue; // skip triggers with BBC decision
      setTriggerBitMap(trig, header);

      if (trig->isBHT2())
         HT2_trigger_ids.push_back(trig->GetId());

      // For BHT triggers the trigger id IS the firing tower id: if that tower is
      // in the bad-tower mask the trigger is "hot tower" contamination — count it,
      // but keep the trigger object so jp_match / ht_match can still decide.
      for (Int_t ntower = 0; ntower < header->GetNOfTowers(); ntower++) {
         TStarJetPicoTower *ptower = pReader->GetEvent()->GetTower(ntower);
         if (ptower->GetId() == trig->GetId() && !towerCuts->IsTowerOK(ptower, pReader->GetEvent())) {
            count_bad_tower_triggers++;
            break;
         }
      }
      triggers.push_back(trig);
   }

   // Drop events whose every BHT2 trigger came from a bad tower.
   if (current_trigger == "HT2" && HT2_trigger_ids.size() > 0 && count_bad_tower_triggers == (int)HT2_trigger_ids.size()) {
      return EVENTRESULT::NOTACCEPTED;
   }

   // Per-event shouldFire from the trigger objects:
   shouldJP0 = shouldJP1 = shouldJP2 = false;
   shouldHwJP0 = shouldHwJP1 = shouldHwJP2 = false;
   for (auto t : triggers) {
      if (t->GetBit(7)) { // bit 7 = kOnline (hardware-equivalent) copies, whose
         if (t->isJP0()) shouldHwJP0 = true; // family bits come from the ONLINE
         if (t->isJP1()) shouldHwJP1 = true; // thresholds
         if (t->isJP2()) shouldHwJP2 = true;
         continue;
      }
      if (t->isJP0()) shouldJP0 = true;
      if (t->isJP1()) shouldJP1 = true;
      if (t->isJP2()) shouldJP2 = true;
   }
   // When the pico stores the full-simulator decision it supersedes the
   // barrel-patch derivation above.
   if (haveSimuDecision) {
      shouldJP0 = simuJP0;
      shouldJP1 = simuJP1;
      shouldJP2 = simuJP2;
   }
   // Event maxima of the simulator ADCs (hardware-first, as jp_patch_adc_near) and the
   // thresholds, so a +-1 ADC threshold shift can be replayed downstream.
   jp_adc_max = ht_adc_max = -1;
   {
      bool haveHW = false;
      for (auto t : triggers)
         if (t->GetBit(7)) { haveHW = true; break; }
      for (auto t : triggers) {
         if ((t->isJP0() || t->isJP1() || t->isJP2() || t->GetBit(9)) && t->GetBit(7) == haveHW && t->GetADC() > jp_adc_max) jp_adc_max = t->GetADC(); // bit 9: patch below the JP0 threshold (ADC only)
         if ((t->isBHT1() || t->isBHT2() || t->isBHT3()) && t->GetADC() > ht_adc_max) ht_adc_max = t->GetADC();
      }
      for (int i = 0; i < 3; ++i) jp_thr[i] = header->GetJetPatchThreshold(i);
      for (int i = 0; i < 4; ++i) ht_thr[i] = header->GetHighTowerThreshold(i);
   }

   // isTriggerEvent = didFire AND shouldFire.
   //   didFire    = trigger id present in the header trigger-id list (hardware);
   //   shouldFire = the simulator agrees (stored decision, or a trigger object
   //                with the matching isJP0/isJP1/isJP2/isBHT2 bit).
   // Hot-tower-only triggers fire the hardware but not the simulator, so they
   // are excluded here.
   if (pars.intype == INPICO && current_trigger.Length() > 0) {
      bool simu_fired = false;
      if (haveSimuDecision) {
         if (current_trigger == "JP2") simu_fired = shouldJP2;
         else if (current_trigger == "JP1") simu_fired = shouldJP1;
         else if (current_trigger == "JP0") simu_fired = shouldJP0;
         else if (current_trigger == "HT2") simu_fired = shouldHT2;
      } else
      for (auto t : triggers) {
         if (current_trigger == "JP2" && t->isJP2()) {
            simu_fired = true;
            break;
         } else if (current_trigger == "JP1" && t->isJP1()) {
            simu_fired = true;
            break;
         } else if (current_trigger == "JP0" && t->isJP0()) {
            simu_fired = true;
            break;
         } else if (current_trigger == "HT2" && t->isBHT2()) {
            simu_fired = true;
            break;
         }
      }
      // MB has no simu equivalent here — keep `isTriggerEvent` as-is.
      if (current_trigger != "MB")
         isTriggerEvent = isTriggerEvent && simu_fired;
   }

   for (auto trig : triggers) {

      if (trig->isJP2() || trig->isJP1() || trig->isJP0()) {
         // don't do anything if eta and phi are not 0
         if (trig->GetEta() != 0 && trig->GetPhi() != 0)
            continue;
         float eta, phi;
         if (getBarrelJetPatchEtaPhi(trig->GetId(), eta, phi)) {
            trig->SetEta(eta);
            trig->SetPhi(phi);
         }
      }
   }

   refmult = header->GetProperReferenceMultiplicity();
   eventid = header->GetEventId();
   runid1 = header->GetRunId();
   // 2021 embedding: every event header reports the same runId, so the real run
   // survives only in the pico name (e21_<bin>_<run>.root, one run per file).
   // Parse it from there to keep runid1 and the per-run bookkeeping meaningful.
   if (pars.InputName.Contains("e21_pt")) {
      TString base = gSystem->BaseName(pReader->GetInputChain()->GetCurrentFile()->GetName());
      base.ReplaceAll(".root", "");
      TString runtok = base(base.Last('_') + 1, base.Length()); // trailing _<run>
      if (runtok.IsDigit()) runid1 = runtok.Atoi();
   }
   QA_hist.SetRun(runid1); // run-binned constituent QA (track/tower vs run)
   vz = header->GetPrimaryVertexZ();
   pthat = (pars.intype == MCPICO) ? header->GetReferenceCentralityWeight() : -1.0; // MC header pt-hat, -1 for data
   vy = header->GetPrimaryVertexY();
   vx = header->GetPrimaryVertexX();
   vz_vpd = header->GetVpdVz();
   n_vpd_east = header->GetnumberOfVpdEastHits();
   n_vpd_west = header->GetnumberOfVpdWestHits();

   QA_hist.vx->Fill(vx);
   QA_hist.vy->Fill(vy);
   QA_hist.vz->Fill(vz);
   QA_hist.vz_vpd->Fill(vz_vpd);
   QA_hist.vz_diff->Fill(vz_vpd - vz);

   // if vpd is available, use it to discard events in MB trigger only
   if (pars.TriggerName.Contains("MB") && abs(vz_vpd - vz) > 6) {
      return EVENTRESULT::NOTACCEPTED;
   }

   // For GEANT: Need to devise a runid that's unique but also
   // reproducible to match Geant and GeantMc data.
   if (pars.UseGeantNumbering) {
      TString cname = gSystem->BaseName(pReader->GetInputChain()->GetCurrentFile()->GetName());
      UInt_t filehash = cname.Hash();
      while (filehash > INT_MAX - 100000)
         filehash -= INT_MAX / 4;
      if (filehash < 1000000)
         filehash += 1000001;
      runid = filehash;
      eventid = pReader->GetNOfCurrentEvent();
   }

   TList *tracksList = pReader->GetListOfSelectedTracks();

   // Fill the particle container. Note the event-level max-track-pT veto is OFF
   // (pars.MaxEventPtCut = 1000); the equivalent cut is applied per jet via
   // pars.MaxJetTrackPt in the jet loop below.
   for (int i = 0; i < pFullEvent->GetEntries(); ++i) {
      TStarJetVector *sv = (TStarJetVector *)pFullEvent->At(i);
      int trackid = sv->GetTrackID();
      int container_id = sv->GetTowerID();
      float track_sdcaxy = -999; // carried into the constituent user_info (for the leading-track veto)

      if (trackid < 0 && pars.intype == INPICO) { // it means it is a track -  not tower
         TStarJetPicoPrimaryTrack *track = (TStarJetPicoPrimaryTrack *)tracksList->At(i);
         trackid = i;
         float sDCAxy = track->GetsDCAxy();
         track_sdcaxy = sDCAxy;
         float pt = track->GetPt();
         // flat |sDCAxy| cap (kept for legacy / systematics; default disabled)
         if (fabs(sDCAxy) > pars.sDCAxyCut)
            continue;
         // Pico maker keeps flag >= 0, so we drop flag == 0 here.
         if (track->GetFlag() < pars.FlagMin)
            continue;
         // Hit counts in the MuDst convention: nHits = nHitsFit - 1 (the vertex point is not a
         // hit) and nHitsPoss = TPC nHitsPoss + 1. Require nHits > 12 and nHits/nHitsPoss >= 0.51.
         {
            const int nh = track->GetNOfFittedHits() - 1, np = track->GetNOfPossHits() + 1;
            if (nh <= 12 || (np > 0 && (double)nh / np < 0.51))
               continue;
         }
         // pT-dependent cap on the FULL 3-D DCA (dcaGlobal().mag() = GetDCA(),
         // WITH z — not the transverse dca_xy), on top of the flat |DCA| < 3 cm
         // the reader already applied. Thresholds: see ppParameters.hh.
         if (pars.ApplyTdcaPtDep) {
            const double a = track->GetDCA();
            double dcaMax;
            if (pt < pars.TdcaPt1) {
               dcaMax = pars.TdcaDcaMax1;
            } else if (pt < pars.TdcaPt2) {
               dcaMax = pars.TdcaDcaMax1
                        + (pars.TdcaDcaMax2 - pars.TdcaDcaMax1) /
                              (pars.TdcaPt2 - pars.TdcaPt1) * (pt - pars.TdcaPt1);
            } else {
               dcaMax = pars.TdcaDcaMax2;
            }
            if (a > fabs(dcaMax))
               continue;
         }
      }

      // The constituent pT window is a TRACK cut only: towers enter with any pT,
      // after the Stage-0 Et > 0.2 GeV threshold and the hadronic correction.
      if (sv->GetCharge() != 0 && (sv->Pt() < pars.PtConsMin || sv->Pt() > pars.PtConsMax))
         continue;
      if (fabs(sv->Eta()) > pars.EtaConsCut)
         continue;

      //! setting pid mass for particles from pythia
      if (pars.intype == MCPICO) {
         TParticlePDG *pid = (TParticlePDG *)PDGdb.GetParticle((Int_t)sv->mc_pdg_pid());
         double sv_mass = pid->Mass();
         double E = TMath::Sqrt(sv->P() * sv->P() + sv_mass * sv_mass);
         sv->SetPxPyPzE(sv->Px(), sv->Py(), sv->Pz(), E);
         event_sum_pt += sv->Pt();
      }

      //! setting pion mass for charged particles and zero mass for towers
      if (pars.intype == INPICO) {
         TParticlePDG *pid = (TParticlePDG *)PDGdb.GetParticle(211);
         double sv_mass;
         if (sv->GetCharge() != 0)
            sv_mass = pid->Mass(); //! tracks get pion mass
         else
            sv_mass = 0.0; //! towers get photon mass ~ 0

         double E = TMath::Sqrt(sv->P() * sv->P() + sv_mass * sv_mass);

         sv->SetPxPyPzE(sv->Px(), sv->Py(), sv->Pz(), E);
         event_sum_pt += sv->Pt();
      }

      // Tracks: efficiency systematic (randomly drop a fraction of tracks).
      if (sv->GetCharge() != 0) {
         Double_t mran = gRandom->Uniform(0, 1);
         if (mran > pars.FakeEff) {
            continue;
         }
      }

      // Towers: gain systematic.
      if (!sv->GetCharge()) {
         (*sv) *= pars.fTowScale;
      }

      particles.push_back(PseudoJet(*sv));
      int id = sv->GetCharge() != 0 ? trackid : container_id;
      auto *cons_ui = new JetAnalysisUserInfo(sv->GetCharge(), sv->mc_pdg_pid(), "", id);
      cons_ui->SetsDCAxy(track_sdcaxy);
      particles.back().set_user_info(cons_ui);
   }

   // For pythia, use the cross section as event weight
   if (pars.InputName.Contains("hat") || pars.InputName.Contains("e21_pt")) {
      TString currentfile = pReader->GetInputChain()->GetCurrentFile()->GetName();
      weight = getPythiaWeight(currentfile);
      if (fabs(weight - 1) < 1e-4) {
         throw std::runtime_error("mcweight unchanged!");
      }
   }

   if (pJA) {
      delete pJA;
      pJA = 0;
   }
   // Active area with explicit ghosts (ghost area 0.04, one repeat). The ghost
   // random sequence gives a ~10% area jitter per jet.
   fastjet::GhostedAreaSpec gspec(EtaGhostCut, 1, 0.04);
   fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, gspec);
   pJA = new JetAnalyzer(particles, JetDef, area_def);

   JetAnalyzer &JA = *pJA;
   vector<PseudoJet> JAResult = sorted_by_pt(select_jet(JA.inclusive_jets()));
   if (JAResult.size() == 0) {
      QA_hist.FillEvent(vx, vy, vz, vz_vpd, refmult, mult, 0, event_sum_pt);
      return EVENTRESULT::NOJETS;
   }

   // Leading-jet vs pt-hat outlier veto, effectively DISABLED (pthat_mult = 1e6):
   // the pt-hat samples are combined with the per-event weight instead.
   double pthat_mult = 1e6;

   if ((pars.InputName.Contains("hat23_") && JAResult[0].perp() > 3.0 * pthat_mult) ||
       (pars.InputName.Contains("hat34_") && JAResult[0].perp() > 4.0 * pthat_mult) ||
       (pars.InputName.Contains("hat45_") && JAResult[0].perp() > 5.0 * pthat_mult) ||
       (pars.InputName.Contains("hat57_") && JAResult[0].perp() > 7.0 * pthat_mult) ||
       (pars.InputName.Contains("hat79_") && JAResult[0].perp() > 9.0 * pthat_mult) ||
       (pars.InputName.Contains("hat911_") && JAResult[0].perp() > 11.0 * pthat_mult) ||
       (pars.InputName.Contains("hat1115_") && JAResult[0].perp() > 15.0 * pthat_mult) ||
       (pars.InputName.Contains("hat1520_") && JAResult[0].perp() > 20.0 * pthat_mult) ||
       (pars.InputName.Contains("hat2025_") && JAResult[0].perp() > 25.0 * pthat_mult) ||
       (pars.InputName.Contains("hat2535_") && JAResult[0].perp() > 35.0 * pthat_mult) ||
       (pars.InputName.Contains("hat3545_") && JAResult[0].perp() > 45.0 * pthat_mult) ||
       (pars.InputName.Contains("hat4555_") && JAResult[0].perp() > 55.0 * pthat_mult) ||
       (pars.InputName.Contains("hat55999_") && JAResult[0].perp() > 1000.0)) {
      return EVENTRESULT::NOTACCEPTED;
   }

   // same veto, 2021-embedding naming (e21_<bin>_<run> picos)
   if ((pars.InputName.Contains("e21_pt2_3_") && JAResult[0].perp() > 3.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt3_4_") && JAResult[0].perp() > 4.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt4_5_") && JAResult[0].perp() > 5.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt5_7_") && JAResult[0].perp() > 7.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt7_9_") && JAResult[0].perp() > 9.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt9_11_") && JAResult[0].perp() > 11.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt11_15_") && JAResult[0].perp() > 15.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt15_20_") && JAResult[0].perp() > 20.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt20_25_") && JAResult[0].perp() > 25.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt25_35_") && JAResult[0].perp() > 35.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt35_45_") && JAResult[0].perp() > 45.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt45_55_") && JAResult[0].perp() > 55.0 * pthat_mult) ||
       (pars.InputName.Contains("e21_pt55_-1_") && JAResult[0].perp() > 1000.0)) {
      return EVENTRESULT::NOTACCEPTED;
   }

   for (unsigned ijet = 0; ijet < JAResult.size(); ijet++) {
      PseudoJet &CurrentJet = JAResult[ijet];

      // Leading-track fake veto: a single fake high-pT track (large |sDCAxy|)
      // fakes a high-pT jet, typically with 1-2 constituents. If either of the
      // jet's TWO leading constituents is a HARD (pT > 2 GeV) charged track with
      // |sDCAxy| above pars.LeadTrackSdcaCut, drop the whole jet. The pT > 2
      // threshold is essential: soft tracks have genuinely wide sDCAxy (the
      // pT-dependent DCA cut allows up to 2 cm at low pT), so vetoing on them
      // would throw away real jets. Default off (99999).
      if (pars.LeadTrackSdcaCut < 9999) {
         const vector<PseudoJet> lead_sorted = sorted_by_pt(CurrentJet.constituents());
         bool fake_lead = false;
         const size_t ncheck = std::min<size_t>(2, lead_sorted.size());
         for (size_t ic = 0; ic < ncheck; ++ic) {
            if (lead_sorted[ic].perp() < 2.0) continue; // soft tracks: normal DCA spread, exempt
            if (!lead_sorted[ic].has_user_info<JetAnalysisUserInfo>()) continue;
            const auto &lui = lead_sorted[ic].user_info<JetAnalysisUserInfo>();
            if (lui.GetQuarkCharge() != 0 && fabs(lui.GetsDCAxy()) > pars.LeadTrackSdcaCut) {
               fake_lead = true;
               break;
            }
         }
         if (fake_lead) continue; // fake-track jet: skip
      }

      vector<PseudoJet> charged_constituents = sorted_by_pt(OnlyCharged(CurrentJet.constituents()));
      vector<PseudoJet> neutral_constituents = sorted_by_pt(OnlyNeutral(CurrentJet.constituents()));

      PseudoJet NeutralPart = join(OnlyNeutral(CurrentJet.constituents()));

      bool is_matched_jp = match_jp(CurrentJet, triggers, pars.R, pars.TriggerName, vz);
      // Per-threshold JP matches (isJP0/isJP1/isJP2 separately): each trigger
      // gates on its OWN patch match downstream.
      bool is_matched_jp0 = match_jp(CurrentJet, triggers, pars.R, "JP0", vz);
      bool is_matched_jp1 = match_jp(CurrentJet, triggers, pars.R, "JP1", vz);
      bool is_matched_jp2 = match_jp(CurrentJet, triggers, pars.R, "JP2", vz);
      bool is_matched_ht = match_ht(CurrentJet, triggers, pars.R);

      double jetptne = 0.0;
      for (PseudoJet &n : NeutralPart.constituents()) {
         jetptne += n.perp();
      }

      // Neutral energy fraction = sum(tower pT) / (sum(tower pT) + sum(track pT)),
      // SCALAR sums (the standard STAR jet-tree definition, not the jet's vector pT).
      double jetptch = 0.0;
      for (const PseudoJet &c : charged_constituents)
         jetptch += c.perp();
      const double jetnef = (jetptne + jetptch) > 0 ? jetptne / (jetptne + jetptch) : 0.0;
      JetAnalysisUserInfo *userinfo = new JetAnalysisUserInfo();
      // The multi-purpose "number" field carries the neutral energy fraction.
      userinfo->SetNumber(jetnef);
      userinfo->SetMatchJP(is_matched_jp);
      userinfo->SetMatchJP0(is_matched_jp0);
      userinfo->SetMatchJP1(is_matched_jp1);
      userinfo->SetMatchJP2(is_matched_jp2);
      userinfo->SetMatchHT(is_matched_ht);
      userinfo->SetHtAdc(ht_adc_max_in_jet(CurrentJet, triggers));
      // Leading-neutral tower BEMC id (-1 if no neutral constituents); for a
      // tower the constituent user_info "number" field holds its id.
      int leadTowId = -1;
      if (!neutral_constituents.empty()) {
         leadTowId = (int)neutral_constituents.front().user_info<JetAnalysisUserInfo>().GetNumber();
      }
      userinfo->SetLeadTowerId(leadTowId);
      // QA probe: near-max JP-patch ADC matched to this jet (-1 if no JP0+ patch
      // in the box).
      {
         const int jpadc = jp_patch_adc_near(CurrentJet, triggers, pars.R, vz);
         userinfo->SetJpAdc(jpadc);
         if (jpadc >= 0) QA_hist.run_jp_patch_adc->Fill(runid1, jpadc);
      }

      if (pars.MaxJetNEF < 1.0 && jetnef > pars.MaxJetNEF)
         continue;

      // Per-jet max-track-pT veto, INPICO only: MC truth jets keep their high-pT
      // particles.
      if (pars.intype == INPICO && pars.MaxJetTrackPt > 0) {
         bool drop_jet = false;
         for (const PseudoJet &c : charged_constituents) {
            if (c.perp() > pars.MaxJetTrackPt) {
               drop_jet = true;
               break;
            }
         }
         if (drop_jet)
            continue;
      }

      CurrentJet.set_user_info(userinfo);

      QA_hist.FillJet(CurrentJet, is_matched_jp, is_matched_ht, tracksList);
      // Jet area + UE density (off-axis cones).
      double area = 0.0, bg_rho = 0.0, pt_corr = CurrentJet.perp();
      try {
         area = CurrentJet.area(); // ghosted-active area
         bg_rho = off_axis_cones_density(CurrentJet, particles, pars.R);
         pt_corr = CurrentJet.perp() - bg_rho * area;
      } catch (...) {
         // A cluster sequence without an area definition throws; we always make
         // one above, so this is only a guard.
      }
      ResultStruct rs(CurrentJet);
      rs.area = area;
      rs.bg_density = bg_rho;
      rs.pt_corrected = pt_corr;
      Result.push_back(rs);
   }
   // By default, sort for original jet pt
   sort(Result.begin(), Result.end(), ResultStruct::origptgreater);

   QA_hist.FillEvent(vx, vy, vz, vz_vpd, refmult, mult, Result.size(), event_sum_pt);

   return EVENTRESULT::JETSFOUND;
}
//----------------------------------------------------------------------
void InitializeReader(std::shared_ptr<TStarJetPicoReader> pReader, const TString InputName, const Long64_t NEvents,
                      const int PicoDebugLevel, const double HadronicCorr)
{

   TStarJetPicoReader &reader = *pReader;

   if (HadronicCorr < 0) {
      reader.SetApplyFractionHadronicCorrection(kFALSE);
      reader.SetApplyMIPCorrection(kTRUE);
      reader.SetRejectTowerElectrons(kTRUE);
   } else {
      reader.SetApplyFractionHadronicCorrection(kTRUE);
      reader.SetFractionHadronicCorrection(HadronicCorr);
      reader.SetApplyMIPCorrection(kFALSE);
      reader.SetRejectTowerElectrons(kFALSE);
   }

   reader.Init(NEvents);
   TStarJetPicoDefinitions::SetDebugLevel(PicoDebugLevel);
}
//----------------------------------------------------------------------
shared_ptr<TStarJetPicoReader> SetupReader(TChain *chain, const ppParameters &pars)
{
   TStarJetPicoDefinitions::SetDebugLevel(0); // 10 for more output

   shared_ptr<TStarJetPicoReader> pReader = make_shared<TStarJetPicoReader>();
   TStarJetPicoReader &reader = *pReader;
   reader.SetInputChain(chain);

   // Event and track selection
   // -------------------------
   TStarJetPicoEventCuts *evCuts = reader.GetEventCuts();
   evCuts->SetVertexZCut(pars.VzCut);
   evCuts->SetPVRankingCut(0.0);
   evCuts->SetRefMultCut(pars.RefMultCut);
   evCuts->SetMaxEventPtCut(pars.MaxEventPtCut);
   evCuts->SetMaxEventEtCut(pars.MaxEventEtCut);

   std::cout << "Using these event cuts:" << std::endl;
   std::cout << " Vz: " << evCuts->GetVertexZCut() << std::endl;
   std::cout << " Refmult: " << evCuts->GetRefMultCutMin() << " -- " << evCuts->GetRefMultCutMax() << std::endl;
   std::cout << " MaxEventPt:  " << evCuts->GetMaxEventPtCut() << std::endl;
   std::cout << " MaxEventEt:  " << evCuts->GetMaxEventEtCut() << std::endl;

   // Tracks cuts
   TStarJetPicoTrackCuts *trackCuts = reader.GetTrackCuts();
   trackCuts->SetDCACut(pars.DcaCut);
   trackCuts->SetMinNFitPointsCut(pars.NMinFit);
   trackCuts->SetFitOverMaxPointsCut(pars.FitOverMaxPointsCut);
   trackCuts->SetMaxPtCut(pars.MaxTrackPt);

   std::cout << "Using these track cuts:" << std::endl;
   std::cout << " dca : " << trackCuts->GetDCACut() << std::endl;
   std::cout << " nfit : " << trackCuts->GetMinNFitPointsCut() << std::endl;
   std::cout << " nfitratio : " << trackCuts->GetFitOverMaxPointsCut() << std::endl;
   std::cout << " maxpt : " << trackCuts->GetMaxPtCut() << std::endl;

   // Towers
   TStarJetPicoTowerCuts *towerCuts = reader.GetTowerCuts();
   towerCuts->SetMaxEtCut(pars.MaxEtCut);

   std::cout << "Using these tower cuts:" << std::endl;
   std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
   std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;

   reader.SetProcessV0s(false);

   return pReader;
}

//----------------------------------------------------------------------
void TurnOffCuts(std::shared_ptr<TStarJetPicoReader> pReader)
{
   pReader->SetProcessTowers(false);
   TStarJetPicoEventCuts *evCuts = pReader->GetEventCuts();
   evCuts->SetTriggerSelection("All"); // All, MB, HT, pp, ppHT, ppJP
   evCuts->SetVertexZCut(999);
   evCuts->SetRefMultCut(0);
   evCuts->SetVertexZDiffCut(999999);

   evCuts->SetMaxEventPtCut(99999);
   evCuts->SetMaxEventEtCut(99999);
   evCuts->SetMinEventEtCut(-1);

   evCuts->SetPVRankingCutOff(); // vertex ranking cut off

   // Tracks cuts
   TStarJetPicoTrackCuts *trackCuts = pReader->GetTrackCuts();
   trackCuts->SetDCACut(99999);
   trackCuts->SetMinNFitPointsCut(-1);
   trackCuts->SetFitOverMaxPointsCut(-1);
   trackCuts->SetMaxPtCut(99999);

   // Towers: there should be no tower in MC — charged and neutral alike arrive
   // as tracks.
   TStarJetPicoTowerCuts *towerCuts = pReader->GetTowerCuts();
   towerCuts->SetMaxEtCut(99999);

   cout << " TURNED OFF ALL CUTS" << endl;
}

//----------------------------------------------------------------------
double getPythiaWeight(TString filename)
{
   const vector<double> cross_section_mb = {9.00176,     1.46259,     0.354407,    0.151627,    0.0249102,
                                            0.00584656,  0.0023021,   0.000342608, 4.56842e-05, 9.71569e-06,
                                            4.69593e-07, 2.69062e-08, 1.43197e-09};

   const vector<double> n_events = {3614773,  3706843, 3709985, 3563592, 3637343, 17337984, 17233020,
                                    16422119, 3547865, 2415179, 2525739, 1203188, 1264931};

   const static vector<string> pt_hat_bins = {"hat23_",   "hat34_",   "hat45_",   "hat57_",   "hat79_",
                                              "hat911_",  "hat1115_", "hat1520_", "hat2025_", "hat2535_",
                                              "hat3545_", "hat4555_", "hat55999_"};

   for (int i = 0; i < pt_hat_bins.size(); ++i) {
      if (filename.Contains(pt_hat_bins.at(i).data())) {
         double weight = cross_section_mb[i] / n_events[i];
         return weight;
      }
   }

   // 2021 embedding: raw Pythia cross sections from the request document, same
   // convention as the table above. n_events = GENERATED events per bin over all
   // runs (totals of lists/emb2021_events_per_run.txt). The tags keep their
   // trailing underscore so that "pt2_3" cannot match inside "pt20_25".
   const vector<double> e21_xsec_mb = {9.0012,      1.46253,     0.354566,    0.151622,    0.0249062,
                                       0.00584527,  0.00230158,  0.000342755, 4.57002e-05, 9.72535e-06,
                                       4.69889e-07, 2.69202e-08, 1.43453e-09};
   const vector<double> e21_n_events = {3685750, 3688505, 3688505, 3688505, 3688505, 3688505, 3688505,
                                        3687820, 3688505, 2458949, 2458949, 1229393, 1229393};
   const static vector<string> e21_bins = {"e21_pt2_3_",   "e21_pt3_4_",   "e21_pt4_5_",   "e21_pt5_7_",
                                           "e21_pt7_9_",   "e21_pt9_11_",  "e21_pt11_15_", "e21_pt15_20_",
                                           "e21_pt20_25_", "e21_pt25_35_", "e21_pt35_45_", "e21_pt45_55_",
                                           "e21_pt55_-1_"};
   for (int i = 0; i < e21_bins.size(); ++i) {
      if (filename.Contains(e21_bins.at(i).data())) {
         return e21_xsec_mb[i] / e21_n_events[i];
      }
   }

   throw std::runtime_error(std::string("No matching pythia pt hat bin found in filename: ") + filename.Data());
}
bool match_jp(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R, const TString &which, double vz)
{
   // which: "JP0/1/2" family or "JPany"; Contains() so prefixed names match.
   // Geometry: strict patch containment for R<0.5, |d_eta|,|d_phi|<0.6 box for R>=0.5.
   const bool any = (which.Length() == 0) || which.Contains("JPany");
   const bool jp0 = any || which.Contains("JP0");
   const bool jp1 = any || which.Contains("JP1");
   const bool jp2 = any || which.Contains("JP2");
   const bool useStrict = R < 0.5f;
   // match on DETECTOR eta: detEta = asinh(sinh(eta) + vz/R_BEMC).
   const double jetDetEta = std::asinh(std::sinh(jet.eta()) + vz / 225.405);
   // bit-7 (hardware-equivalent kOnline) patches take precedence when present.
   bool haveHW = false;
   for (auto t : triggers)
      if (t->GetBit(7)) { haveHW = true; break; }
   for (auto trigger : triggers) {
      if (trigger->GetBit(7) != haveHW)
         continue;
      const bool fired =
         (jp2 && trigger->isJP2()) || (jp1 && trigger->isJP1()) || (jp0 && trigger->isJP0());
      if (!fired)
         continue;
      const bool match = useStrict
         ? isInsideJetPatch(trigger->GetId(), jetDetEta, jet.phi())
         : isInsideJetPatchBox(trigger->GetId(), jetDetEta, jet.phi());
      if (match)
         return true;
   }
   return false;
}

// Max ADC over JP patches geometrically matched to the jet (SAME box/detEta as
// match_jp), regardless of which JP threshold fired. The maker injects a trig obj
// only for patches with ADC>JP0(=20), so this is the near-max ADC over JP0+ patches
// in the jet's box; -1 => no such patch.
int jp_patch_adc_near(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R, double vz)
{
   const bool useStrict = R < 0.5f;
   const double jetDetEta = std::asinh(std::sinh(jet.eta()) + vz / 225.405);
   // hardware-first, same convention as match_jp.
   bool haveHW = false;
   for (auto t : triggers)
      if (t->GetBit(7)) { haveHW = true; break; }
   int maxAdc = -1;
   for (auto trigger : triggers) {
      if (trigger->GetBit(7) != haveHW)
         continue;
      const bool match = useStrict
         ? isInsideJetPatch(trigger->GetId(), jetDetEta, jet.phi())
         : isInsideJetPatchBox(trigger->GetId(), jetDetEta, jet.phi());
      if (match && trigger->GetADC() > maxAdc)
         maxAdc = trigger->GetADC();
   }
   return maxAdc;
}

bool match_ht(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R)
{
   for (auto trigger : triggers) {
      if (!trigger->isBHT2())
         continue;
      int trigger_towerid = trigger->GetId();
      for (PseudoJet &part : jet.constituents()) {
         if (part.is_pure_ghost())
            continue; // ghosts have no JetAnalysisUserInfo
         if (!part.has_user_info<JetAnalysisUserInfo>())
            continue;
         if (part.user_info<JetAnalysisUserInfo>().GetNumber() == trigger_towerid) {
            return true;
         }
      }
   }
   return false;
}

// Max DSM high-tower ADC among the BHT trigger objects whose tower is a jet constituent (-1 if none).
int ht_adc_max_in_jet(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers)
{
   int maxAdc = -1;
   for (auto trigger : triggers) {
      if (!(trigger->isBHT1() || trigger->isBHT2() || trigger->isBHT3()) || trigger->GetADC() <= maxAdc)
         continue;
      for (PseudoJet &part : jet.constituents()) {
         if (part.is_pure_ghost() || !part.has_user_info<JetAnalysisUserInfo>())
            continue;
         if (part.user_info<JetAnalysisUserInfo>().GetQuarkCharge() == 0 &&
             part.user_info<JetAnalysisUserInfo>().GetNumber() == trigger->GetId()) {
            maxAdc = trigger->GetADC();
            break;
         }
      }
   }
   return maxAdc;
}

void setTriggerBitMap(TStarJetPicoTriggerInfo *trig, TStarJetPicoEventHeader *header) // needed for real data
{
   std::bitset<32> original_bitmap = trig->GetBitMap();
   Int_t trigMap = original_bitmap.to_ulong();
   if (trigMap != 0) {
      return; // bitmap already set
   }
   // bitmap layout :
   // bit 1: barrel high tower 1
   // bit 2: barrel high tower 2
   // bit 3: barrel high tower 3
   // bit 4: jet patch 0
   // bit 5: jet patch 1
   // bit 6: jet patch 2
   // bit 7-31: open
   // valid only for pp12 data:
   header->SetJetPatchThreshold(0, 20);  // jp0
   header->SetJetPatchThreshold(1, 28);  // jp1
   header->SetJetPatchThreshold(2, 36);  // jp2
   header->SetHighTowerThreshold(0, 11); // bht0
   header->SetHighTowerThreshold(1, 15); // bht1
   header->SetHighTowerThreshold(2, 18); // bht2
   header->SetHighTowerThreshold(3, 8);  // bht3

   // Jet patches sit at detector eta 0.5, -0.1 or -0.5; anything else is a tower.
   Float_t eta = trig->GetEta();
   bool jp_eta_flag = false;
   if (eta > -0.10000001 && eta < -0.09999999999)
      jp_eta_flag = true;
   else if (eta == 0.5 || eta == -0.5)
      jp_eta_flag = true;

   if (trig->GetId() <= 17 && jp_eta_flag) {

      Int_t jpAdc = trig->GetADC();
      UInt_t jp0 = header->GetJetPatchThreshold(0);
      UInt_t jp1 = header->GetJetPatchThreshold(1);
      UInt_t jp2 = header->GetJetPatchThreshold(2);

      if (jp0 > 0 && jpAdc > jp0)
         trigMap |= 1 << 4;
      if (jp1 > 0 && jpAdc > jp1)
         trigMap |= 1 << 5;
      if (jp2 > 0 && jpAdc > jp2)
         trigMap |= 1 << 6;

      trig->SetBitMap(trigMap);
   } else // it means bht
   {
      Int_t bhtAdc = trig->GetADC();
      UInt_t bht1 = header->GetHighTowerThreshold(1);
      UInt_t bht2 = header->GetHighTowerThreshold(2);
      UInt_t bht3 = header->GetHighTowerThreshold(3);

      if (bht1 > 0 && bhtAdc > bht1)
         trigMap |= 1 << 1;
      if (bht2 > 0 && bhtAdc > bht2)
         trigMap |= 1 << 2;
      if (bht3 > 0 && bhtAdc > bht3)
         trigMap |= 1 << 3;

      trig->SetBitMap(trigMap);
   }
}

float getJetPatchPhi(int jetPatch)
{
   return TVector2::Phi_mpi_pi((150 - (jetPatch % 6) * 60) * TMath::DegToRad());
}

bool getBarrelJetPatchEtaPhi(int jetPatch, float &eta, float &phi)
{
   if (jetPatch < 0 || jetPatch >= 18)
      return false;
   // Patch numbering (BEMC, http://drupal.star.bnl.gov/STAR/system/files/BEMC_y2004.pdf):
   // JP0 is centred at phi = 150 deg seen from the West into the interaction
   // region and the numbering runs clockwise in 60 deg steps. Patches 0-5 sit at
   // eta = 0.5, 6-11 at eta = -0.5, 12-17 at eta = -0.1; patches 6 apart share a
   // phi (JP0 and JP6, JP1 and JP7, ...).
   if (jetPatch >= 0 && jetPatch < 6)
      eta = 0.5f;
   if (jetPatch >= 6 && jetPatch < 12)
      eta = -0.5f;
   if (jetPatch >= 12 && jetPatch < 18)
      eta = -0.1f;
   phi = getJetPatchPhi(jetPatch);
   return true;
}

bool isInsideJetPatch(const int &jetPatch, const float &jetEta, const float &jetPhi)
{
   float eta_center, phi_center;
   if (!getBarrelJetPatchEtaPhi(jetPatch, eta_center, phi_center))
      return false;
   // Jet patch size is 1.0 in eta and 60 degrees in phi
   float etaMin = eta_center - 0.5f;
   float etaMax = eta_center + 0.5f;
   float phiMin = TVector2::Phi_mpi_pi(phi_center - 30 * TMath::DegToRad());
   float phiMax = TVector2::Phi_mpi_pi(phi_center + 30 * TMath::DegToRad());
   float phi = TVector2::Phi_mpi_pi(jetPhi);

   if (jetEta < etaMin || jetEta > etaMax)
      return false;

   if (phiMin <= phiMax) {
      return (phiMin <= phi && phi < phiMax);
   } else {
      // Wrapped interval, e.g. [2.8, -2.8) in radians
      return (phi >= phiMin || phi < phiMax);
   }

   return false;
}
bool isInsideJetPatchBox(const int &jetPatch, const float &jetEta, const float &jetPhi)
{
   float eta_center, phi_center;
   if (!getBarrelJetPatchEtaPhi(jetPatch, eta_center, phi_center))
      return false;
   float deta = jetEta - eta_center;
   float dphi = TVector2::Phi_mpi_pi(jetPhi - phi_center);
   return (std::fabs(deta) < 0.6f) && (std::fabs(dphi) < 0.6f);
}

double off_axis_cones_density(const PseudoJet &jet, const vector<PseudoJet> &particles, double R)
{
   // `particles` is the raw input list: not attached to any cluster sequence, so
   // calling is_pure_ghost() on it would throw. No ghost filter is needed anyway —
   // the inputs are real tracks/towers by construction.
   const double eta_j = jet.eta();
   const double phi_j = jet.phi();
   const double phi1 = phi_j + M_PI / 2;
   const double phi2 = phi_j - M_PI / 2;
   const double R2 = R * R;
   double sum1 = 0.0, sum2 = 0.0;
   for (const auto &p : particles) {
      const double deta = p.eta() - eta_j;
      const double deta2 = deta * deta;
      const double dphi1 = TVector2::Phi_mpi_pi(p.phi() - phi1);
      if (deta2 + dphi1 * dphi1 < R2)
         sum1 += p.perp();
      const double dphi2 = TVector2::Phi_mpi_pi(p.phi() - phi2);
      if (deta2 + dphi2 * dphi2 < R2)
         sum2 += p.perp();
   }
   const double cone_area = M_PI * R2;
   return 0.5 * (sum1 + sum2) / cone_area;
}
