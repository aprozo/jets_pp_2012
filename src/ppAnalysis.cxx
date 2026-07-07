
#include "ppAnalysis.hh"
#include <climits>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h> // for atof, atoi
#include <string>

using std::cerr;
using std::cout;
using std::endl;

bool match_jp(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R);
bool match_jp(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R, const TString &which = "JP2", double vz = 0.0);
bool match_ht(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R);
int jp_patch_adc_near(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R, double vz);
void setTriggerBitMap(TStarJetPicoTriggerInfo *trig, TStarJetPicoEventHeader *header);
bool getBarrelJetPatchEtaPhi(int jetPatch, float &eta, float &phi);

// Off-axis-cones UE density estimator (GeV per unit (η,φ) area).
// Two cones at (η_jet, φ_jet ± π/2) of radius R; ρ = avg(ΣpT) / (π R²).
double off_axis_cones_density(const PseudoJet &jet, const vector<PseudoJet> &particles, double R);

double getPythiaWeight(TString filename);
// Standard ctor
ppAnalysis::ppAnalysis(const int argc, const char **const argv)
{
   // Parse arguments
   // ---------------
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

   // Consistency checks
   // ------------------

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

   // Derived rapidity cuts
   // ---------------------
   // Jet PHYSICS-eta acceptance = Dmitry's ETAFID = 1.0, so that jets with
   // |det_eta|<0.5 survive jet-finding instead of being pre-cut at |physics|<0.5
   // (a jet with |det_eta|<0.5 at |vz|<60 has |physics eta| < ~0.72, so 1.0 is
   // safely inclusive). The FINAL |det_eta|<0.5 detector acceptance is imposed
   // downstream (det_eta is written by RunppAna.cxx and cut in cross_section
   // Raw()), reproducing Dmitry's |physics|<1.0 then |det_eta|<0.5 acceptance.
   // The CONSTITUENT cut (pars.EtaConsCut, applied in the particle loop) is
   // unchanged.
   EtaJetCut = 1.0;
   EtaGhostCut = EtaJetCut + 2.0 * pars.R;

   // Jet candidate selectors
   // -----------------------
   select_jet_eta = SelectorAbsEtaMax(EtaJetCut);
   select_jet_pt = SelectorPtRange(pars.PtJetMin, pars.PtJetMax);

   // if (pars.intype == MCPICO)
   //   select_jet = SelectorPtRange(2, pars.PtJetMax);
   // else
   select_jet = select_jet_eta * select_jet_pt;

   // Repeat on subjets?
   // ------------------
   pars.Recursive = pars.InputName.Contains("Pythia") && false;

   // Initialize jet finding
   // ----------------------

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

   // For trees of TStarJetVector
   // (like a previous result)
   // ------------------------
   Events = new TChain(pars.ChainName);
   Events->Add(pars.InputName);
   if (NEvents < 0)
      NEvents = INT_MAX;

   // For picoDSTs
   // -------------

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
   }

   // initialize qa histograms
   QA_hist.Init();

   cout << "N = " << NEvents << endl;

   cout << "Done initializing chains. " << endl;
   return true;
}
//----------------------------------------------------------------------
// Main routine for one event.
EVENTRESULT ppAnalysis::RunEvent()
{
   // Reset results (from last event)
   // -------------------------------
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
   // pp 2012 200 GeV
   // JP2 event trigger - 370621
   // BHT2 trigger:  370531

   set<int> event_triggers;
   for (int i = 0; i < header->GetNOfTriggerIds(); i++) {
      event_triggers.insert(header->GetTriggerId(i));
   }

   // only if not mcpico
   map<TString, set<int>> trigger_map_2012;
   trigger_map_2012["JP2"] = {370621};
   trigger_map_2012["JP1"] = {370611};
   trigger_map_2012["JP0"] = {370601};
   trigger_map_2012["HT2"] = {370531, 500205}; // 500205 - leftover in Youqi embedding trees
   trigger_map_2012["MB"] = {370011};

   // Per-event JP fire flags (header trigger IDs). For data these are the real
   // prescale-accepted hardware bits — measured P(fired_JP0|fired_JP2)=0.0083
   // ~ 1/ps0, P(fired_JP1|fired_JP2)=0.40 ~ 1/ps1 (2026-06-11).
   firedJP0 = event_triggers.count(370601) != 0;
   firedJP1 = event_triggers.count(370611) != 0;
   firedJP2 = event_triggers.count(370621) != 0;

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
      // if (!has_trigger)
      //    return EVENTRESULT::NOTACCEPTED;
   }

   vector<TStarJetPicoTriggerInfo *> triggers;
   TStarJetPicoTowerCuts *towerCuts = pReader->GetTowerCuts();
   vector<int> HT2_trigger_ids;
   int count_bad_tower_triggers = 0;

   // Build the simu-side trigger-object list and, in parallel, count BHT2
   // triggers caused by towers in the analysis bad-tower list. This implements
   // the "didFire AND shouldFire" event filter (mirrors Dmitry's
   // make_jet_plots.cxx:336-365): a trigger object exists only if the trigger
   // simulator agreed that the patch/tower fired (shouldFire>0), and the
   // hardware-fired list (header->GetTriggerId) provides didFire.
   for (int i = 0; i < header->GetNOfTrigObjs(); ++i) {
      auto trig = pReader->GetEvent()->GetTrigObj(i);
      // https://github.com/wsu-yale-rhig/TStarJetPicoMaker/blob/82e051867037038001ea1218256ef48e3dfca9a0/StRoot/JetPicoMaker/StMuJetAnalysisTreeMaker.cxx#L772
      // check if triggerBitMap are not set to zero
      if (trig->GetTriggerFlag() == 5 && trig->GetADC() == 0)
         continue; // skip triggers with BBC decision
      setTriggerBitMap(trig, header);

      if (trig->isBHT2())
         HT2_trigger_ids.push_back(trig->GetId());

      // For BHT triggers (trig ID = firing tower ID): is the firing tower in
      // the analysis bad-tower mask? If yes, this BHT trigger is "hot tower"
      // contamination — count it but still keep the trigger object (so
      // jp_match / ht_match can decide).
      bool is_bad_tower = false;
      for (Int_t ntower = 0; ntower < header->GetNOfTowers(); ntower++) {
         TStarJetPicoTower *ptower = pReader->GetEvent()->GetTower(ntower);
         if (ptower->GetId() == trig->GetId() && !towerCuts->IsTowerOK(ptower, pReader->GetEvent())) {
            is_bad_tower = true;
            count_bad_tower_triggers++;
            break;
         }
      }
      triggers.push_back(trig);
   }

   // Hot-tower-only HT2 events: every BHT2 trigger came from a bad tower → drop.
   // (Mirrors Dmitry's analysis-time filter; biases the cross section LOW
   // otherwise because the event lands in the luminosity sum but contains
   // no real jets after the bad-tower mask removes the firing tower.)
   if (HT2_trigger_ids.size() > 0 && count_bad_tower_triggers == (int)HT2_trigger_ids.size()) {
      return EVENTRESULT::NOTACCEPTED;
   }

   // didFire AND shouldFire (the simu-vs-hardware AND), per trigger.
   // didFire    = trigger ID present in header->fTriggerIdArray (hardware).
   // shouldFire = at least one TStarJetPicoTriggerInfo has the matching
   //              isJP0/isJP1/isJP2/isBHT2 bit set (simu agrees).
   // For events where only didFire is true (hot-tower-only triggers), the
   // simu disagrees and the event is excluded from `isTriggerEvent`.
   if (pars.intype == INPICO && current_trigger.Length() > 0) {
      bool simu_fired = false;
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

   // //  ADC values:
   // // fHighTowerThreshold[4] = 11 , 15 , 18 , 8  - bht0, bht1, bht2, bht3
   // // fJetPatchThreshold[3]  = 20 , 28 , 36  -     jp0, jp1, jp2

   refmult = header->GetProperReferenceMultiplicity();
   eventid = header->GetEventId();
   runid1 = header->GetRunId();
   QA_hist.SetRun(runid1); // run-binned constituent QA (track/tower vs run)
   // Promote vz to a class member so RunppAna can expose it on ResultTree
   // for downstream vertex-z reweighting (Dmitry-equivalent of
   // SetVertexReweightingParams).
   vz = header->GetPrimaryVertexZ();
   double vy = header->GetPrimaryVertexY();
   double vx = header->GetPrimaryVertexX();
   double vz_vpd = header->GetVpdVz();

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
         filehash -= INT_MAX / 4; // some random large number
      if (filehash < 1000000)
         filehash += 1000001;
      runid = filehash;
      eventid = pReader->GetNOfCurrentEvent();
   }

   TList *tracksList = pReader->GetListOfSelectedTracks();
   TList *towersList = pReader->GetListOfSelectedTowers();
   // print content of tracks and towers

   // Event-level max-track-pT veto is applied via
   // TStarJetPicoEventCuts::SetMaxEventPtCut(pars.MaxEventPtCut) — see
   // BuildEventAndJetCuts() below. With MaxEventPtCut=30 we match
   // Dmitry's make_max_track_pt_cut(30).

   // Fill particle container
   // -----------------------
   for (int i = 0; i < pFullEvent->GetEntries(); ++i) {
      TStarJetVector *sv = (TStarJetVector *)pFullEvent->At(i);
      int trackid = sv->GetTrackID();
      int container_id = sv->GetTowerID();

      if (trackid < 0 && pars.intype == INPICO) { // it means it is a track -  not tower
         TStarJetPicoPrimaryTrack *track = (TStarJetPicoPrimaryTrack *)tracksList->At(i);
         trackid = i;
         float sDCAxy = track->GetsDCAxy();
         int charge = track->GetCharge();
         float pt = track->GetPt();
         // flat |sDCAxy| cap (kept for legacy / systematics; default disabled)
         if (fabs(sDCAxy) > pars.sDCAxyCut)
            continue;
         // Dmitry's StjTrackCutFlag(0): reject flag <= 0 (require flag > 0).
         // Pico maker keeps flag >= 0, so we drop flag == 0 here.
         if (track->GetFlag() < pars.FlagMin)
            continue;
         // Two-part DCA selection, Dmitry's Run12 alignment (Table 2 of the
         // analysis note). BOTH parts cap the FULL 3-D DCA = dcaGlobal().mag()
         // = track->GetDCA() (WITH z), NOT the transverse component:
         //   (1) flat |DCA| < 3 cm, applied by the reader via
         //       SetDCACut(pars.DcaCut) on track->GetDCA() (see SetupReader); and
         //   (2) the pT-dependent cap below — Dmitry's StjTrackCutTdcaPtDependent,
         //       which despite the "T" name cuts Tdca = dcaGlobal().mag() (the
         //       3-D magnitude); see star-jet/.../mudst/StjTPCMuDst.cxx:100:
         //         DCA < 2 cm                 for pt < 0.5 GeV
         //             < 2.5 cm - (1/GeV)*pt  for 0.5 <= pt < 1.5 GeV  (slope -1)
         //             < 1 cm                 for pt >= 1.5 GeV
         //       i.e. pt1=0.5/dca1=2.0, pt2=1.5/dca2=1.0. (Since dca1=2 < 3, this
         //       pT-dependent cap is always at least as tight as the flat 3 cm.)
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

      // Ensure kinematic similarity
      if (sv->Pt() < pars.PtConsMin || sv->Pt() > pars.PtConsMax)
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

      // TRACKS
      // ------
      if (sv->GetCharge() != 0) {
         // EFFICIENCY uncertainty
         // ----------------------
         Double_t mran = gRandom->Uniform(0, 1);
         if (mran > pars.FakeEff) {
            continue;
         }
      }

      // TOWERS
      // ------
      // Shift gain
      if (!sv->GetCharge()) {
         (*sv) *= pars.fTowScale; // for systematics
      }

      particles.push_back(PseudoJet(*sv));
      int id = sv->GetCharge() != 0 ? trackid : container_id;
      particles.back().set_user_info(new JetAnalysisUserInfo(sv->GetCharge(), sv->mc_pdg_pid(), "", id));
   }

   // mult = particles.size();

   // For pythia, use cross section as weight
   // ---------------------------------------

   if (pars.InputName.Contains("hat")) {
      TString currentfile = pReader->GetInputChain()->GetCurrentFile()->GetName();
      weight = getPythiaWeight(currentfile);
      if (fabs(weight - 1) < 1e-4) {
         throw std::runtime_error("mcweight unchanged!");
      }
   }

   // Run analysis
   // ------------
   if (pJA) {
      delete pJA;
      pJA = 0;
   }
   // Area-aware clustering so jet.area() works. Ghost area 0.04 matches Dmitry's
   // StFastJetAreaPars in run12_200GeVJetCode/RunJetFinder2012UePro.C:162.
   fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,
                                    fastjet::GhostedAreaSpec(EtaGhostCut, 1, 0.04));
   pJA = new JetAnalyzer(particles, JetDef, area_def);

   JetAnalyzer &JA = *pJA;
   vector<PseudoJet> JAResult = sorted_by_pt(select_jet(JA.inclusive_jets()));
   if (JAResult.size() == 0) {
      QA_hist.FillEvent(vx, vy, vz, vz_vpd, refmult, mult, 0, event_sum_pt);
      return EVENTRESULT::NOJETS;
   }


   // check if the event has high weight or large |vz|
   double pthat_mult = 2;

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

   int njets = JAResult.size();
   // cout << "-----------------------" << endl;

   for (unsigned ijet = 0; ijet < JAResult.size(); ijet++) {
      PseudoJet &CurrentJet = JAResult[ijet];
      vector<PseudoJet> charged_constituents = sorted_by_pt(OnlyCharged(CurrentJet.constituents()));
      vector<PseudoJet> neutral_constituents = sorted_by_pt(OnlyNeutral(CurrentJet.constituents()));

      PseudoJet NeutralPart = join(OnlyNeutral(CurrentJet.constituents()));

      // bool is_matched_jp = match_jp(CurrentJet, triggers, pars.R);
      bool is_matched_jp = match_jp(CurrentJet, triggers, pars.R, pars.TriggerName, vz);
      // Per-threshold JP matches (isJP0/isJP1/isJP2 separately) — each trigger's
      // analysis gates on its OWN patch match. NOT degenerate (the old code
      // copied one flag into all three downstream branches).
      bool is_matched_jp0 = match_jp(CurrentJet, triggers, pars.R, "JP0", vz);
      bool is_matched_jp1 = match_jp(CurrentJet, triggers, pars.R, "JP1", vz);
      bool is_matched_jp2 = match_jp(CurrentJet, triggers, pars.R, "JP2", vz);
      bool is_matched_ht = match_ht(CurrentJet, triggers, pars.R);

      // Analysis-level JP cuts (matches Dmitry's make_detector_level_cut for jp=0/1/2):
      //   - jet must be inside a fired JP-patch (jp_match)
      //   - jp1 requires pT >= 6.0 GeV, jp2 requires pT >= 8.4 GeV (jp0 has no extra threshold)
      // (See Dmitry/star-jet/StJetPlots/StJetCut.h:182-208.)

      double jetpttot = CurrentJet.perp();

      double jetptne = 0.0;
      for (PseudoJet &n : NeutralPart.constituents()) {
         jetptne += n.perp();
      }

      JetAnalysisUserInfo *userinfo = new JetAnalysisUserInfo();
      // Save neutral energy fraction in multi-purpose field
      userinfo->SetNumber(jetptne / jetpttot);
      userinfo->SetMatchJP(is_matched_jp);
      userinfo->SetMatchJP0(is_matched_jp0);
      userinfo->SetMatchJP1(is_matched_jp1);
      userinfo->SetMatchJP2(is_matched_jp2);
      userinfo->SetMatchHT(is_matched_ht);
      // Leading-neutral tower BEMC id (-1 if no neutral constituents).
      // Constituent user_info stores tower id via container_id (see L601).
      int leadTowId = -1;
      if (!neutral_constituents.empty()) {
         leadTowId = (int)neutral_constituents.front().user_info<JetAnalysisUserInfo>().GetNumber();
      }
      userinfo->SetLeadTowerId(leadTowId);
      // Near-max JP-patch ADC matched to this jet (-1 if no JP0+ patch in box).
      // Provenance/QA probe: degenerate-data signature check uses
      // trigger_match_JP2 && jp_patch_adc<=36 ~ 0.
      {
         const int jpadc = jp_patch_adc_near(CurrentJet, triggers, pars.R, vz);
         userinfo->SetJpAdc(jpadc);
         if (jpadc >= 0) QA_hist.run_jp_patch_adc->Fill(runid1, jpadc); // per-run ADC, full coverage
      }

      if (pars.MaxJetNEF < 1.0 && (jetptne / jetpttot) > pars.MaxJetNEF)
         continue;

      // Per-jet max-track-pT cut (mirrors Dmitry's make_max_track_pt_cut).
      // INPICO only — MC truth jets keep their high-pT particles.
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

      // auto leadingTrack = charged_constituents.at(0);
      // if (jetpttot > 22)
      //    QA_hist.highjetpt_leadingtrack_pt->Fill(leadingTrack.perp());

      // if (charged_constituents.size() == 1) {
      //    int container_id = leadingTrack.user_info<JetAnalysisUserInfo>().GetNumber();
      //    TStarJetPicoPrimaryTrack *track = (TStarJetPicoPrimaryTrack *)tracksList->At(container_id);
      //    QA_hist.FillTrack(track);
      // }

      // if (neutral_constituents.size() > 0) {
      //    auto leadingTower = neutral_constituents.at(0);
      //    int towerId = leadingTower.user_info<JetAnalysisUserInfo>().GetNumber();
      //    QA_hist.jetpt_TowerID->Fill(leadingTower.perp(), towerId);
      //    if (jetpttot > 22.0)
      //       // Fill QA histograms only for high pt jets
      //       QA_hist.highjetpt_leadingtower_pt->Fill(leadingTower.perp());
      // }

      CurrentJet.set_user_info(userinfo);

      QA_hist.FillJet(CurrentJet, is_matched_jp, is_matched_ht, tracksList);
      // Jet area + UE density (off-axis cones).
      double area = 0.0, bg_rho = 0.0, pt_corr = CurrentJet.perp();
      try {
         area = CurrentJet.area(); // ghosted-active area
         bg_rho = off_axis_cones_density(CurrentJet, particles, pars.R);
         pt_corr = CurrentJet.perp() - bg_rho * area;
      } catch (...) {
         // CS without area definition would throw; we created one above so this is just a guard.
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
// Helper to deal with repetitive stuff
shared_ptr<TStarJetPicoReader> SetupReader(TChain *chain, const ppParameters &pars)
{
   TStarJetPicoDefinitions::SetDebugLevel(0); // 10 for more output

   shared_ptr<TStarJetPicoReader> pReader = make_shared<TStarJetPicoReader>();
   TStarJetPicoReader &reader = *pReader;
   reader.SetInputChain(chain);

   // Event and track selection
   // -------------------------
   TStarJetPicoEventCuts *evCuts = reader.GetEventCuts();
   // evCuts->SetTriggerSelection(pars.TriggerName); // All, MB, HT, pp, ppHT, ppJP

   // Additional cuts
   evCuts->SetVertexZCut(pars.VzCut);
   evCuts->SetPVRankingCut(0.0);
   evCuts->SetRefMultCut(pars.RefMultCut);
   // evCuts->SetVertexZDiffCut(pars.VzDiffCut);
   evCuts->SetMaxEventPtCut(pars.MaxEventPtCut);
   evCuts->SetMaxEventEtCut(pars.MaxEventEtCut);

   // evCuts->SetMinEventEtCut(pars.MinEventEtCut);

   std::cout << "Using these event cuts:" << std::endl;
   std::cout << " Vz: " << evCuts->GetVertexZCut() << std::endl;
   std::cout << " Refmult: " << evCuts->GetRefMultCutMin() << " -- " << evCuts->GetRefMultCutMax() << std::endl;
   // std::cout << " Delta Vz:  " << evCuts->GetVertexZDiffCut() << std::endl;
   std::cout << " MaxEventPt:  " << evCuts->GetMaxEventPtCut() << std::endl;
   std::cout << " MaxEventEt:  " << evCuts->GetMaxEventEtCut() << std::endl;

   // This method does NOT WORK for GEANT MC trees because everything is in the
   // tracks... Do it by hand later on, using pars.ManualHtCut; Also doesn't
   // work for general trees, but there it can't be fixed

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

   // V0s: Turn off
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
   // evCuts->SetMinEventPtCut (-1);
   evCuts->SetMinEventEtCut(-1);

   evCuts->SetPVRankingCutOff(); //  Use SetPVRankingCutOff() to turn off
                                 //  vertex ranking cut.  default is OFF

   // Tracks cuts
   TStarJetPicoTrackCuts *trackCuts = pReader->GetTrackCuts();
   trackCuts->SetDCACut(99999);
   trackCuts->SetMinNFitPointsCut(-1);
   trackCuts->SetFitOverMaxPointsCut(-1);
   trackCuts->SetMaxPtCut(99999);

   // Towers: should be no tower in MC. All (charged or neutral) are handled in
   // track
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
   // if no match found, throw error
   throw std::runtime_error(std::string("No matching pythia pt hat bin found in filename: ") + filename.Data());
}
// bool match_jp(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R)
// {
//    for (auto trigger : triggers) {
//       if (!trigger->isJP2())
//          continue;
//       float eta = trigger->GetEta();
//       float phi = trigger->GetPhi();
//       float deta = jet.eta() - eta;
//       float dphi = TVector2::Phi_mpi_pi(jet.phi() - phi);
//       if (sqrt(deta * deta + dphi * dphi) < R)
//          return true;
//    }
//    return false;
// }

bool match_jp(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R, const TString &which, double vz)
{
   // `which` selects which JP-trigger family to require a fired patch from.
   // "JP2" -> isJP2, "JP1" -> isJP1, "JP0" -> isJP0, "JPany" or empty -> any of the three.
   // Use Contains() so prefixed names from container.sh ("data_JP2", "geant_JP2",
   // "mc_JP2") are recognized — exact equality fails for those.
   //
   // Matching geometry depends on R:
   //   R <  0.5: use strict patch containment (jet axis inside the 1.0×60°
   //             JP-patch box). Appropriate when the cone fits inside the
   //             patch — fewer edge effects.
   //   R >= 0.5: use Dmitry-relaxed match (|Δη|<0.6 AND |Δφ|<0.6 from patch
   //             centre) since the cone overlaps the patch edges and the
   //             strict cut would drop genuinely-triggered jets.
   const bool any = (which.Length() == 0) || which.Contains("JPany");
   const bool jp0 = any || which.Contains("JP0");
   const bool jp1 = any || which.Contains("JP1");
   const bool jp2 = any || which.Contains("JP2");
   const bool useStrict = (R < 0.5f);
   // AUDIT FIX #4: match the JP patch on the jet DETECTOR eta (Dmitry jp_match.h:17),
   // not the physics eta. detEta = asinh(sinh(eta_phys) + vz/225.405) (BEMC radius).
   const double jetDetEta = std::asinh(std::sinh(jet.eta()) + vz / 225.405);
   // HARDWARE-first match: the data picos carry the kOnline (hardware-equivalent,
   // verified == Dmitry skim 513/513) JP patches as bit-7 trigger objects. When
   // present, match EXCLUSIVELY against them with their native family bits
   // (= the L0 register decision). Embedding and legacy picos have no bit-7
   // objects, so they fall through to the offline-emulator isJP*() gate (the same
   // ruler the response uses, since sim events have no hardware).
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
         : isInsideJetPatchDmitry(trigger->GetId(), jetDetEta, jet.phi());
      if (match)
         return true;
   }
   return false;
}

bool match_jp(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R)
{
   for (auto trigger : triggers) {
      if (trigger->isJP2() && isInsideJetPatch(trigger->GetId(), jet.eta(), jet.phi()))
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
   const bool useStrict = (R < 0.5f);
   const double jetDetEta = std::asinh(std::sinh(jet.eta()) + vz / 225.405);
   // 2026-07-02: hardware-first, same convention as match_jp — on new data
   // picos the reported patch ADC is the kOnline (hardware) value.
   bool haveHW = false;
   for (auto t : triggers)
      if (t->GetBit(7)) { haveHW = true; break; }
   int maxAdc = -1;
   for (auto trigger : triggers) {
      if (trigger->GetBit(7) != haveHW)
         continue;
      const bool match = useStrict
         ? isInsideJetPatch(trigger->GetId(), jetDetEta, jet.phi())
         : isInsideJetPatchDmitry(trigger->GetId(), jetDetEta, jet.phi());
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
      // double eta = trigger->GetEta();
      // double phi = trigger->GetPhi();
      // double deta = jet.eta() - eta;
      // double dphi = TVector2::Phi_mpi_pi(jet.phi() - phi);
      // if (sqrt(deta * deta + dphi * dphi) < R) {
      //    return true;
      // }
   }
   return false;
}

void setTriggerBitMap(TStarJetPicoTriggerInfo *trig, TStarJetPicoEventHeader *header) // needed for real data
{
   // check if trigmap is not 0
   std::bitset<32> original_bitmap = trig->GetBitMap();
   Int_t trigMap = original_bitmap.to_ulong(); // get the original bitmap
   if (trigMap != 0) {
      return; // if the bitmap is already set, no need to set it again
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
   // //  ADC values:
   // // fHighTowerThreshold[4] = 11 , 15 , 18 , 8  - bht0, bht1, bht2, bht3
   // // fJetPatchThreshold[3]  = 20 , 28 , 36  -     jp0, jp1, jp2

   Float_t eta = trig->GetEta();
   // compare eta to -0.100000
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
   // Sanity check
   if (jetPatch < 0 || jetPatch >= 18)
      return false;
   // The jet patches are numbered starting with JP0 centered at 150 degrees
   // looking from the West into the IR (intersection region) and increasing
   // clockwise, i.e. JP1 at 90 degrees, JP2 at 30 degrees, etc. On the East
   // side the numbering picks up at JP6 centered again at 150 degrees and
   // increasing clockwise (again as seen from the *West* into the IR). Thus
   // JP0 and JP6 are in the same phi location in the STAR coordinate system.
   // So are JP1 and JP7, etc.
   // JP locations:
   // Jet Patch# Eta   Phi,degrees(rad)   Quadrant
   // 0          0.5   150 (2.618)           10'
   // 1          0.5   90 (1.5708)           12'
   // 2          0.5   30 (0.5236)            2'
   // 3          0.5  -30 (-0.5236)           4'
   // 4          0.5  -90 (-1.5708)           6'
   // 5          0.5  -150 (-2.618)            8'
   // 6         -0.5   150 (2.618)           10'
   // 7         -0.5   90 (1.5708)           12'
   // 8         -0.5   30 (0.5236)            2'
   // 9         -0.5  -30 (-0.5236)           4'
   // 10        -0.5  -90 (-1.5708)           6'
   // 11        -0.5  -150 (-2.618)            8'
   // 12        -0.1   150 (2.618)           10'
   // 13        -0.1   90 (1.5708)           12'
   // 14        -0.1   30 (0.5236)            2'
   // 15        -0.1  -30 (-0.5236)           4'
   // 16        -0.1  -90 (-1.5708)           6'
   // 17        -0.1  -150 (-2.618)            8'

   // http://drupal.star.bnl.gov/STAR/system/files/BEMC_y2004.pdf

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
bool isInsideJetPatchDmitry(const int &jetPatch, const float &jetEta, const float &jetPhi)
{
   // Match Dmitry's `jp_match`: |Δη| < 0.6 AND |Δφ| < 0.6 from patch center.
   // (See Dmitry/star-jet/src/common/jp_match.h.)
   float eta_center, phi_center;
   if (!getBarrelJetPatchEtaPhi(jetPatch, eta_center, phi_center))
      return false;
   float deta = jetEta - eta_center;
   float dphi = TVector2::Phi_mpi_pi(jetPhi - phi_center);
   return (std::fabs(deta) < 0.6f) && (std::fabs(dphi) < 0.6f);
}

double off_axis_cones_density(const PseudoJet &jet, const vector<PseudoJet> &particles, double R)
{
   // Off-axis cones UE density: two cones at (η_jet, φ_jet ± π/2) of radius R.
   // ρ = ((ΣpT)_cone1 + (ΣpT)_cone2) / (2 · π R²).
   //
   // `particles` is the raw input list — these are not associated with any
   // cluster sequence, so `is_pure_ghost()` would trigger
   // "Trying to access the structure of a PseudoJet which has no associated structure".
   // Inputs are real tracks/towers by construction, so no ghost filter needed.
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