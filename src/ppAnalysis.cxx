
#include "ppAnalysis.hh"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h> // for getenv, atof, atoi
#include <string>

using std::cerr;
using std::cout;
using std::endl;

bool getBarrelJetPatchEtaPhi(int jetPatch, float &eta, float &phi);

bool match_jp(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers,
              float R);
bool match_ht(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R);

void setTriggerBitMap(TStarJetPicoTriggerInfo *trig,
                      TStarJetPicoEventHeader *header);
// Standard ctor
ppAnalysis::ppAnalysis(const int argc, const char **const argv) {
  // Parse arguments
  // ---------------
  vector<string> arguments(argv + 1, argv + argc);
  bool argsokay = true;
  bool forcedgeantnum = false;
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
      if (pars.IntTowScale < -1 || pars.FakeEff > 1) {
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
      forcedgeantnum = true;
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
         << " [-towunc -1|0|1 ] -- Shift tower energy by this times "
         << pars.fTowUnc << endl
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
    throw std::runtime_error(
        "PtJetMin needs to be positive (0.001 will work).");
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
  EtaJetCut = pars.EtaConsCut - pars.R;
  EtaGhostCut = EtaJetCut + 2.0 * pars.R;

  // Jet candidate selectors
  // -----------------------
  select_jet_eta = SelectorAbsEtaMax(EtaJetCut);
  select_jet_pt = SelectorPtRange(pars.PtJetMin, pars.PtJetMax);
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
  cout << " Reading tree named \"" << pars.ChainName << "\" from "
       << pars.InputName << endl;
  cout << " intype = " << pars.intype << endl;
  cout << " Writing to " << pars.OutFileName << endl;
  cout << " ----------------------------" << endl;
}
//----------------------------------------------------------------------
ppAnalysis::~ppAnalysis() {
  if (pJA) {
    delete pJA;
    pJA = 0;
  }
}

//----------------------------------------------------------------------
bool ppAnalysis::InitChains() {

  // For trees of TStarJetVector
  // (like a previous result)
  // ------------------------
  Events = new TChain(pars.ChainName);
  Events->Add(pars.InputName);
  if (NEvents < 0)
    NEvents = INT_MAX;

  // For picoDSTs
  // -------------
  if (pars.intype == INPICO || pars.intype == MCPICO) {
    pReader = SetupReader(Events, pars);

    InitializeReader(pReader, pars.InputName, NEvents, PicoDebugLevel,
                     pars.HadronicCorr);
    if (pars.intype == MCPICO) {
      TurnOffCuts(pReader);
    }

    cout << "Don't forget to set tower cuts for pReader!!" << endl;
  }

  cout << "N = " << NEvents << endl;

  cout << "Done initializing chains. " << endl;
  return true;
}
//----------------------------------------------------------------------
// Main routine for one event.
EVENTRESULT ppAnalysis::RunEvent() {
  TStarJetVector *sv;

  TStarJetPicoEventHeader *header = 0;
  vector<TStarJetPicoTriggerInfo *> triggers;
  UInt_t filehash = 0;
  TString cname = "";

  // Reset results (from last event)
  // -------------------------------
  Result.clear();
  weight = 1.;
  is_rejected = false;
  mult = 0;
  refmult = 0;
  runid = -(INT_MAX - 1);
  runid1 = -(INT_MAX - 1);
  eventid = -(INT_MAX - 1);
  event_sum_pt = 0.;
  njets = 0;
  vz = 999;
  particles.clear();

  // cout <<"============================================" << endl;
  // cout <<"===================New event================" << endl;
  // cout <<"============================================" << endl;


  TString trigger_array = "";
  switch (pars.intype) {
    // =====================================================
  case INPICO:
  case MCPICO:
    if (!pReader->NextEvent()) {
      pReader->PrintStatus();
      return EVENTRESULT::ENDOFINPUT;
      break;
    }
    pReader->PrintStatus(10);

    pFullEvent = pReader->GetOutputContainer()->GetArray();
    header = pReader->GetEvent()->GetHeader();
    // cout << "============================================ " << endl;

    for (int i = 0; i < header->GetNOfTrigObjs(); ++i) {
      auto trig = pReader->GetEvent()->GetTrigObj(i);
      // https://github.com/wsu-yale-rhig/TStarJetPicoMaker/blob/82e051867037038001ea1218256ef48e3dfca9a0/StRoot/JetPicoMaker/StMuJetAnalysisTreeMaker.cxx#L772
      // check if triggerBitMap are not set to zero
      if (trig->GetTriggerFlag() == 5)
        continue; // skip triggers with BBC decision
      // check if triggermap is 0 and fill it

      setTriggerBitMap(trig, header);
      
      triggers.push_back(trig);
    }
    // fill remaining triggers positions for embedding files where it is not
    // filled
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

    // for (Int_t i = 0; i < header->GetNOfTriggerIds(); i++) {
    //   trigger_array += Form("%7d ", header->GetTriggerId(i));
    // }
    // cout << "event triggers : " << trigger_array << endl;

    // for (auto trig : triggers) {
    //   //  Bitmap
    //   //  last 7 bits :
    //   // |jp2|jp1|jp0|bht3|bht2|bht1|0
    //   // if (!trig->isJP2())
    //   //   continue;
    //   // cout << " BHT2";
    //   cout << "Trigger: " << trig->GetId();
    //   cout << " ADC: " << trig->GetADC();
    //   cout << " Eta: " << trig->GetEta();
    //   cout << " Phi: " << trig->GetPhi();
    //   cout << " BitMap: " << trig->GetBitMap() << endl;
    // }

    // //  ADC values:
    // // fHighTowerThreshold[4] = 11 , 15 , 18 , 8  - bht0, bht1, bht2, bht3
    // // fJetPatchThreshold[3]  = 20 , 28 , 36  -     jp0, jp1, jp2

    // // (500205) BHT2 trigger Id
    // // (500215) BHT2 trigger Id
    // // (500401) JP2 trigger Id
    // // (500411) JP2 trigger Id

    //    fTrigSel.Contains("ppHT"))
    //    mTrigId==370541 || mTrigId==370542 || mTrigId==370351)   //
    //    BHT0*BBCMB*TOF0, BHT0*BBCMB*TOF0, MTD*BHT3
    //   fTrigSel.Contains("ppJP2")) {
    //   mTrigId==370621)		//  JP2
    //   mTrigId==370601 || mTrigId==370611

    refmult = header->GetProperReferenceMultiplicity();
    eventid = header->GetEventId();
    runid1 = header->GetRunId();
    vz = header->GetPrimaryVertexZ();

    // For GEANT: Need to devise a runid that's unique but also
    // reproducible to match Geant and GeantMc data.
    if (pars.UseGeantNumbering) {
      TString cname = gSystem->BaseName(
          pReader->GetInputChain()->GetCurrentFile()->GetName());
      UInt_t filehash = cname.Hash();
      while (filehash > INT_MAX - 100000)
        filehash -= INT_MAX / 4; // some random large number
      if (filehash < 1000000)
        filehash += 1000001;
      runid = filehash;
      eventid = pReader->GetNOfCurrentEvent();
    }
    break;
  default:
    cerr << "Unknown intype " << pars.intype << endl;
    return EVENTRESULT::PROBLEM;
  }

  // Fill particle container
  // -----------------------
  for (int i = 0; i < pFullEvent->GetEntries(); ++i) {
    sv = (TStarJetVector *)pFullEvent->At(i);

    // skip towers to get charged jets

    if (sv->GetCharge() == 0 ) {
      continue;
    }

    // Ensure kinematic similarity
    if (sv->Pt() < pars.PtConsMin || sv->Pt() > pars.PtConsMax)
      continue;
    if (fabs(sv->Eta()) > pars.EtaConsCut)
      continue;

    //! setting pid mass for particles from pythia
    if (pars.intype == MCPICO) {
      TParticlePDG *pid =
          (TParticlePDG *)PDGdb.GetParticle((Int_t)sv->mc_pdg_pid());
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
    particles.back().set_user_info(new JetAnalysisUserInfo(
        3 * sv->GetCharge(), sv->mc_pdg_pid(), "", sv->GetTowerID()));
  }

  mult = particles.size();
  if (particles.size() == 0)
    return EVENTRESULT::NOCONSTS;

  // For pythia, use cross section as weight
  // ---------------------------------------

  if (pars.InputName.Contains("pt-hat")) {
    TString currentfile = pReader->GetInputChain()->GetCurrentFile()->GetName();
    weight = LookupRun12Xsec(currentfile);
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
  pJA = new JetAnalyzer(particles, JetDef);

  JetAnalyzer &JA = *pJA;
  vector<PseudoJet> JAResult = sorted_by_pt(select_jet(JA.inclusive_jets()));
  if (JAResult.size() == 0) {
    return EVENTRESULT::NOJETS;
  }

  // check if the event has high weight or large |vz|
  double pthat_mult=1.5; 

  if ((pars.InputName.Contains("hat23_") && JAResult[0].perp() > 3.0*pthat_mult) ||
      (pars.InputName.Contains("hat34_") && JAResult[0].perp() > 4.0*pthat_mult)||
      (pars.InputName.Contains("hat45_") && JAResult[0].perp() > 5.0*pthat_mult)||
      (pars.InputName.Contains("hat57_") && JAResult[0].perp() > 7.0*pthat_mult)||
      (pars.InputName.Contains("hat79_") && JAResult[0].perp() > 9.0*pthat_mult)||
      (pars.InputName.Contains("hat911_") && JAResult[0].perp() > 11.0*pthat_mult)||
      (pars.InputName.Contains("hat1115_") && JAResult[0].perp() > 15.0*pthat_mult)||
      (pars.InputName.Contains("hat1520_") && JAResult[0].perp() > 20.0*pthat_mult)||
      (pars.InputName.Contains("hat2025_") && JAResult[0].perp() > 25.0*pthat_mult)||
      (pars.InputName.Contains("hat2535_") && JAResult[0].perp() > 35.0*pthat_mult)||
      (pars.InputName.Contains("hat3545_") && JAResult[0].perp() > 45.0*pthat_mult)||
      (pars.InputName.Contains("hat4555_") && JAResult[0].perp() > 55.0*pthat_mult)||
      (pars.InputName.Contains("hat55999_") && JAResult[0].perp() > 1000.0)) {
    is_rejected = true;
    return EVENTRESULT::NOTACCEPTED;
  }

  int njets = JAResult.size();
  // cout << "-----------------------" << endl;

  for (unsigned ijet = 0; ijet < JAResult.size(); ijet++) {
    PseudoJet &CurrentJet = JAResult[ijet];
    PseudoJet NeutralPart = join(OnlyNeutral(CurrentJet.constituents()));
    PseudoJet ChargedPart = join(OnlyCharged(CurrentJet.constituents()));

    bool is_matched_jp = match_jp(CurrentJet, triggers, pars.R);
    bool is_matched_ht = match_ht(CurrentJet, triggers, pars.R);

    double jetptne = 0.0;
    double jetpttot = 0.0;

    for (PseudoJet &n : NeutralPart.constituents()) {
      jetptne += n.perp();
    }

    for (PseudoJet &part : CurrentJet.constituents()) {
      jetpttot += part.perp();
    }

    JetAnalysisUserInfo *userinfo = new JetAnalysisUserInfo();
    // Save neutral energy fraction in multi-purpose field
    userinfo->SetNumber(jetptne / jetpttot);
    userinfo->SetMatchJP(is_matched_jp);
    userinfo->SetMatchHT(is_matched_ht);

    CurrentJet.set_user_info(userinfo);

    if (pars.MaxJetNEF < 1.0 && jetptne / jetpttot > pars.MaxJetNEF)
      continue;

    // vector<PseudoJet> constituents =
    // sorted_by_pt(CurrentJet.constituents()); int nparticles =
    // CurrentJet.constituents().size(); if (nparticles == 0)
    //   continue;
    // float pTlead = constituents[0].pt();
    // double pT_lead0 = 0;
    // double pT_lead3 = 0;
    // double pT_lead5 = 0;
    // double pT_lead7 = 0;
    // if (pTlead > 0)
    // {
    //   pT_lead0 = CurrentJet.pt();
    // }
    // if (pTlead > 3)
    // {
    //   pT_lead3 = CurrentJet.pt();
    // }
    // if (pTlead > 5)
    // {
    //   pT_lead5 = CurrentJet.pt();
    // }
    // if (pTlead > 7)
    // {
    //   pT_lead7 = CurrentJet.pt();
    // }
    Result.push_back(ResultStruct(CurrentJet));
    
  }
  // By default, sort for original jet pt
  sort(Result.begin(), Result.end(), ResultStruct::origptgreater);

  return EVENTRESULT::JETSFOUND;
}
//----------------------------------------------------------------------
void InitializeReader(std::shared_ptr<TStarJetPicoReader> pReader,
                      const TString InputName, const Long64_t NEvents,
                      const int PicoDebugLevel, const double HadronicCorr) {

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
shared_ptr<TStarJetPicoReader> SetupReader(TChain *chain,
                                           const ppParameters &pars) {
  TStarJetPicoDefinitions::SetDebugLevel(0); // 10 for more output

  shared_ptr<TStarJetPicoReader> pReader = make_shared<TStarJetPicoReader>();
  TStarJetPicoReader &reader = *pReader;
  reader.SetInputChain(chain);

  // Event and track selection
  // -------------------------
  TStarJetPicoEventCuts *evCuts = reader.GetEventCuts();
  evCuts->SetTriggerSelection("All"); // All, MB, HT, pp, ppHT, ppJP
  // Additional cuts
  evCuts->SetVertexZCut(pars.VzCut);
  evCuts->SetRefMultCut(pars.RefMultCut);
  evCuts->SetVertexZDiffCut(pars.VzDiffCut);
  evCuts->SetMaxEventPtCut(pars.MaxEventPtCut);
  evCuts->SetMaxEventEtCut(pars.MaxEventEtCut);

  evCuts->SetMinEventEtCut(pars.MinEventEtCut);

  std::cout << "Using these event cuts:" << std::endl;
  std::cout << " Vz: " << evCuts->GetVertexZCut() << std::endl;
  std::cout << " Refmult: " << evCuts->GetRefMultCutMin() << " -- "
            << evCuts->GetRefMultCutMax() << std::endl;
  std::cout << " Delta Vz:  " << evCuts->GetVertexZDiffCut() << std::endl;
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
  std::cout << " nfitratio : " << trackCuts->GetFitOverMaxPointsCut()
            << std::endl;
  std::cout << " maxpt : " << trackCuts->GetMaxPtCut() << std::endl;

  // Towers
  TStarJetPicoTowerCuts *towerCuts = reader.GetTowerCuts();
  towerCuts->SetMaxEtCut(pars.MaxEtCut);

  std::cout << "Using these tower cuts:" << std::endl;
  std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
  std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut()
            << std::endl;

  // V0s: Turn off
  reader.SetProcessV0s(false);

  return pReader;
}

//----------------------------------------------------------------------
void TurnOffCuts(std::shared_ptr<TStarJetPicoReader> pReader) {

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
double LookupRun12Xsec(TString filename) {
  //   (pthatrange)  nevents  weighted-Xsection
  //   2, 3,       2409849,           9.00176
  //   3, 4,       3706843,           1.46259
  //   4, 5,       3709985,          0.354407
  //   5, 7,       3563592,          0.151627
  //   7, 9,       3637343,         0.0249102
  //   9,11,      17337984,        0.00584656
  //  11,15,      17233020,         0.0023021
  //  15,20,      16422119,       0.000342608
  //  20,25,       3547865,       4.56842e-05
  //  25,35,       2415179,       9.71569e-06
  //  35,45,       2525739,       4.69593e-07
  //  45,55,       1203188,       2.69062e-08
  //  55,99,       1264931,       1.43197e-09

  const int NUMBEROFPT = 13;
  // const char
  // *PTBINS[NUMBEROFPT]={"2_3","3_4","4_5","5_7","7_9","9_11","11_15","15_20","20_25","25_35","35_-1"};
  const static float XSEC[NUMBEROFPT] = {
      9.0012,      1.46253,     0.354566,    0.151622,    0.0249062,
      0.00584527,  0.00230158,  0.000342755, 4.57002e-05, 9.72535e-06,
      4.69889e-07, 2.69202e-08, 1.43453e-09};
  const static float NUMBEROFEVENT[NUMBEROFPT] = {
      3e6, 3e6, 3e6, 3e6, 3e6, 3e6, 3e6, 3e6, 3e6, 2e6, 2e6, 1e6, 1e6};
  const static vector<string> vptbins = {
      "pt-hat23_",   "pt-hat34_",   "pt-hat45_",   "pt-hat57_",   "pt-hat79_",
      "pt-hat911_",  "pt-hat1115_", "pt-hat1520_", "pt-hat2025_", "pt-hat2535_",
      "pt-hat3545_", "pt-hat4555_", "pt-hat55999_"};
  for (int i = 0; i < vptbins.size(); ++i) {
    if (filename.Contains(vptbins.at(i).data()))
      return XSEC[i] / NUMBEROFEVENT[i];
  }

  return -1;
}

void print_assert(const char *exp, const char *file, int line) {
  std::cerr << "Assertion '" << exp << "' failed in '" << file << "' at line '"
            << line << "'" << std::endl;
  std::terminate();
}

#define always_assert(exp)                                                     \
  ((exp) ? (void)0 : print_assert(#exp, __FILE__, __LINE__))

float getJetPatchPhi(int jetPatch) {
  return TVector2::Phi_mpi_pi((150 - (jetPatch % 6) * 60) * TMath::DegToRad());
}

bool getBarrelJetPatchEtaPhi(int jetPatch, float &eta, float &phi) {
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
  // 5          0.5  -150 (2.618)            8'
  // 6         -0.5   150 (2.618)           10'
  // 7         -0.5   90 (1.5708)           12'
  // 8         -0.5   30 (0.5236)            2'
  // 9         -0.5  -30 (-0.5236)           4'
  // 10        -0.5  -90 (-1.5708)           6'
  // 11        -0.5  -150 (2.618)            8'
  // 12        -0.1   150 (2.618)           10'
  // 13        -0.1   90 (1.5708)           12'
  // 14        -0.1   30 (0.5236)            2'
  // 15        -0.1  -30 (-0.5236)           4'
  // 16        -0.1  -90 (-1.5708)           6'
  // 17        -0.1  -150 (2.618)            8'

  // http://drupal.star.bnl.gov/STAR/system/files/BEMC_y2004.pdf

  if (jetPatch >= 0 && jetPatch < 6)
    eta = 0.5;
  if (jetPatch >= 6 && jetPatch < 12)
    eta = -0.5;
  if (jetPatch >= 12 && jetPatch < 18)
    eta = -0.1;
  phi = getJetPatchPhi(jetPatch);
  return true;
}

bool match_jp(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers,
              float R) {
  for (auto trigger : triggers) {
    if (trigger->isJP2()) {
      float eta, phi;
      eta = trigger->GetEta();
      phi = trigger->GetPhi();
      double deta = jet.eta() - eta;
      double dphi = TVector2::Phi_mpi_pi(jet.phi() - phi);
      if ((fabs(deta) < R) && (fabs(dphi) < R)) {
        return true;
      }
    }
  }
  return false;
}

bool match_ht(PseudoJet &jet, vector<TStarJetPicoTriggerInfo *> triggers, float R) {

  // print jet towers:
  PseudoJet NeutralPart = join(OnlyNeutral(jet.constituents()));
    for (auto trigger : triggers) {
  
      if (trigger->isBHT2()) {
        int trigger_towerid = trigger->GetId();        
        for (PseudoJet &part : NeutralPart.constituents()) {
          if (part.user_info<JetAnalysisUserInfo>().GetNumber() ==
              trigger_towerid) {

            return true;
          }
        }

        double eta = trigger->GetEta();
        double phi = trigger->GetPhi();
        double deta = jet.eta() - eta;
        double dphi = TVector2::Phi_mpi_pi(jet.phi() - phi);
        if ((fabs(deta) < R) && (fabs(dphi) <R))
          return true; 
      }
    }
    return false;
  }

void setTriggerBitMap(TStarJetPicoTriggerInfo *trig,
                      TStarJetPicoEventHeader *header) {

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

  Float_t eta= trig->GetEta();
  // compare eta to -0.100000
  bool jp_eta_flag=false;
  if (eta > -0.10000001 && eta <-0.09999999999)
    jp_eta_flag=true;
  else if (eta==0.5||eta==-0.5)
    jp_eta_flag=true;

  
  if (trig->GetId() <= 17 && jp_eta_flag)
  {

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
  } 
  else // it means bht
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
