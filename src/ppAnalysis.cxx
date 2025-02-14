
#include "ppAnalysis.hh"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h> // for getenv, atof, atoi
#include <string>

using std::cerr;
using std::cout;
using std::endl;

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
  TStarJetPicoTriggerInfo *triggerInfo = 0;

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
    triggerInfo = pReader->GetEvent()->GetTrigObj(1);
    triggerInfo->PrintInfo();

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
  if ((pars.InputName.Contains("hat23_") && JAResult[0].perp() > 6.0) ||
      (pars.InputName.Contains("hat34_") && JAResult[0].perp() > 8.0) ||
      (pars.InputName.Contains("hat45_") && JAResult[0].perp() > 10.0) ||
      (pars.InputName.Contains("hat57_") && JAResult[0].perp() > 14.0) ||
      (pars.InputName.Contains("hat79_") && JAResult[0].perp() > 18.0) ||
      (pars.InputName.Contains("hat911_") && JAResult[0].perp() > 22.0) ||
      (pars.InputName.Contains("hat1115_") && JAResult[0].perp() > 30.0) ||
      (pars.InputName.Contains("hat1520_") && JAResult[0].perp() > 40.0) ||
      (pars.InputName.Contains("hat2025_") && JAResult[0].perp() > 50.0) ||
      (pars.InputName.Contains("hat2535_") && JAResult[0].perp() > 70.0) ||
      (pars.InputName.Contains("hat3545_") && JAResult[0].perp() > 90.0) ||
      (pars.InputName.Contains("hat4555_") && JAResult[0].perp() > 110.0) ||
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

    double jetptne = 0.0;
    double jetpttot = 0.0;
    double q = 0;
    for (PseudoJet &c : ChargedPart.constituents()) {
      q += c.user_info<JetAnalysisUserInfo>().GetQuarkCharge() / 3.0;
    }

    for (PseudoJet &n : NeutralPart.constituents()) {
      jetptne += n.perp();
    }

    for (PseudoJet &part : CurrentJet.constituents()) {
      jetpttot += part.perp();
    }

    JetAnalysisUserInfo *userinfo = new JetAnalysisUserInfo();
    // Save neutral energy fraction in multi-purpose field
    userinfo->SetNumber(jetptne / jetpttot);
    CurrentJet.set_user_info(userinfo);

    if (pars.MaxJetNEF < 1.0 && jetptne / jetpttot > pars.MaxJetNEF)
      continue;

    // vector<PseudoJet> constituents = sorted_by_pt(CurrentJet.constituents());
    // int nparticles = CurrentJet.constituents().size();
    // if (nparticles == 0)
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
  evCuts->SetTriggerSelection(pars.TriggerName); // All, MB, HT, pp, ppHT, ppJP
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
  // tracks... Do it by hand later on, using pars.ManualHtCut; Also doesn't work
  // for general trees, but there it can't be fixed

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

  evCuts->SetPVRankingCutOff(); //  Use SetPVRankingCutOff() to turn off vertex
                                //  ranking cut.  default is OFF

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

//----------------------------------------------------------------------
ostream &operator<<(ostream &ostr, const PseudoJet &jet) {
  if (jet == 0) {
    ostr << " 0 ";
  } else {
    ostr << " pt = " << jet.pt() << " m = " << jet.m() << " y = " << jet.rap()
         << " phi = " << jet.phi()
         << " ClusSeq = " << (jet.has_associated_cs() ? "yes" : "no");
  }
  return ostr;
}
