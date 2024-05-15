/* @file ppTestAnalysis.cxx
    @author Raghav Kunnawalkam Elayavalli
    @version Revision 1.0
    @brief test pp analysis for Run12 data and embedding 
    @details Uses JetAnalyzer objects
    @date March 16, 2022
*/

#include "ppTestAnalysis.hh"
#include <stdlib.h>     // for getenv, atof, atoi


using std::cout;
using std::cerr;
using std::endl;

// Standard ctor
ppTestAnalysis::ppTestAnalysis ( const int argc, const char** const argv )
{
  // Parse arguments
  // ---------------
  vector<string> arguments(argv + 1, argv + argc);
  bool argsokay=true;
  bool forcedgeantnum=false;
  NEvents=-1;
  for ( auto parg = arguments.begin() ; parg!=arguments.end() ; ++parg){
    string arg=*parg;
    if ( arg == "-R" ){      
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.R = atof( parg->data());
    } else if ( arg == "-lja" ){     
      if (++parg==arguments.end() ){ argsokay=false; break; }
      pars.LargeJetAlgorithm = AlgoFromString ( *parg);
    } else if ( arg == "-pj" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.PtJetMin = atof((parg)->data());
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.PtJetMax = atof((parg)->data());
    } else if ( arg == "-ec" ){      
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.EtaConsCut = atof((parg)->data());      
    } else if ( arg == "-pc" ){      
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.PtConsMin = atof((parg)->data());      
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.PtConsMax = atof((parg)->data());
    } else if ( arg == "-hadcorr" ){
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.HadronicCorr = atof( parg->data());
    } else if ( arg == "-o" ){     
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.OutFileName=*parg;
    } else if ( arg == "-i" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.InputName=*parg;
    } else if ( arg == "-c" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.ChainName=*parg;
    } else if ( arg == "-trig" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      pars.TriggerName=*parg;
    } else if ( arg == "-intype" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      if ( *parg == "pico" ){
	pars.intype = INPICO;
	continue;
      }
      if ( *parg == "mcpico" ){
	pars.intype = MCPICO;
	continue;
      }
      if ( *parg == "mctree" ){
	pars.intype = MCTREE;
	continue;
      }      
      argsokay=false;
      break;
    } else if ( arg == "-N" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      NEvents=atoi(parg->data());
    } else if ( arg == "-tracksmear" ){
      if (++parg==arguments.end() ){ argsokay=false; break; }
      switch ( atoi(parg->data()) ){
      case 0:
	cout << "Track Smearing disabled." << endl;
	break;
      case 1:
	cout << "Track Smearing set to pp primaries. DeltapT/pT = 0.01+0.005 pT " << endl;
	SigmaPt = new TF1("SigmaPt","[0] + [1]*x",0,100);
	SigmaPt->FixParameter ( 0, 0.01);
	SigmaPt->FixParameter ( 1, 0.005);
	break;
      case 2:
	cout << "Track Smearing set to AA primaries. DeltapT/pT = 0.005+0.0025 pT " << endl;
	SigmaPt = new TF1("SigmaPt","[0] + [1]*x",0,100);
	SigmaPt->FixParameter ( 0, 0.005);
	SigmaPt->FixParameter ( 1, 0.0025);

	break;
      case 3:
	cout << "Track Smearing set to AA globals not implemented. DeltapT/pT = 0.01 * pT^2" << endl;
	SigmaPt = new TF1("SigmaPt","[0]*x*x",0,100);
	SigmaPt->FixParameter ( 0, 0.01);
	break;
      case 4:
	cout << "Unrecognized pT smearing option. " << endl;
	argsokay=false;
	break;	
      }
    } else if ( arg == "-fakeeff" ){
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.FakeEff = atof( parg->data());
      cout << "Setting fake efficiency to " << pars.FakeEff << endl;
      if ( pars.FakeEff<0 || pars.FakeEff>1 ){ argsokay=false; break; }
    } else if ( arg == "-towunc" ){
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.IntTowScale = atoi( parg->data());
      pars.fTowScale = 1.0 + pars.IntTowScale*pars.fTowUnc;
      cout << "Setting tower scale to " << pars.fTowScale << endl;
      if ( pars.IntTowScale<-1 || pars.FakeEff>1 ){ argsokay=false; break; }
    } else if ( arg == "-geantnum" ){      
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      if ( *parg != "0" && *parg != "1" ){ argsokay=false; break; }
      forcedgeantnum=true;
      pars.UseGeantNumbering = bool ( atoi(parg->data()));      
    } else if ( arg == "-jetnef" ){
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      pars.MaxJetNEF = atof( parg->data());
      cout << "Setting Max Jet NEF to " << pars.MaxJetNEF << endl;
      if ( pars.MaxJetNEF<0 || pars.MaxJetNEF>1 ){ argsokay=false; break; }      
    } else {
      argsokay=false;
      break;
    }
  }
  
  if ( !argsokay ) {
    cerr << "usage: " << argv[0] << endl
	 << " [-o OutFileName]"  << endl
      	 << " [-N Nevents (<0 for all)]" << endl
	 << " [-R radius]" << endl
      	 << " [-lja LargeJetAlgorithm]" << endl
      	 << " [-i infilepattern]" << endl
      	 << " [-c chainname]" << endl
      	 << " [-intype pico|mcpico|tree|mctree|herwigtree]" << endl
	 << " [-trig trigger name (e.g. HT)]" << endl
	 << " [-pj PtJetMin PtJetMax]" << endl
	 << " [-ec EtaConsCut]" << endl
	 << " [-pc PtConsMin PtConsMax]" << endl
	 << " [-hadcorr HadronicCorrection]  -- Set to a negative value for MIP correction." << endl
      	 << " [-psc PtSubConsMin PtSubConsMax]" << endl
	 << " [-tracksmear number] -- enable track pT smearing. " << endl
	 << "                      -- 1: pp primaries, 2: AuAu primaries, 3: AuAu (maybe also pp) globals." << endl
      	 << " [-fakeeff (0..1)] -- enable fake efficiency for systematics. 0.95 is a reasonable example." << endl
      	 << " [-towunc -1|0|1 ] -- Shift tower energy by this times " << pars.fTowUnc << endl
      	 << " [-geantnum true|false] (Force Geant run event id hack)" << endl      
	 << endl << endl
	 << "NOTE: Wildcarded file patterns should be in single quotes." << endl
	 << endl;      
    throw std::runtime_error("Not a valid list of options");
  }


  // Consistency checks
  // ------------------

  if ( !forcedgeantnum ){
    if ( pars.InputName.Contains( "Geant" ) || pars.InputName.Contains( "EmbedPythiaRun12" ) ){
      pars.UseGeantNumbering=true;
    }
  }
    
  if ( pars.PtJetMin<=0 ){
    throw std::runtime_error("PtJetMin needs to be positive (0.001 will work).");
  }

  if ( pars.intype==MCPICO ){
    // This refers to the JetTreeMc branch
    if ( pars.ChainName!= "JetTreeMc" ){
      throw std::runtime_error("Unsuitable chain name for pars.intype==MCPICO");
    }
  }          
  if ( pars.ChainName == "JetTreeMc" ){
    if ( pars.intype !=MCPICO ){
      throw std::runtime_error("Unsuitable chain name for pars.intype==MCPICO");
    }
  }

  // Derived rapidity cuts
  // ---------------------
  EtaJetCut     = pars.EtaConsCut - pars.R;
  EtaGhostCut   = EtaJetCut  + 2.0*pars.R;
  
  // Jet candidate selectors
  // -----------------------
  // select_jet_eta     = SelectorAbsRapMax( EtaJetCut );
  select_jet_eta     = SelectorAbsEtaMax( EtaJetCut );
  select_jet_pt      = SelectorPtRange( pars.PtJetMin, pars.PtJetMax );
  select_jet         = select_jet_eta * select_jet_pt;     
  
  // Repeat on subjets?
  // ------------------
  pars.Recursive = pars.InputName.Contains("Pythia") && false;

  // Initialize jet finding
  // ----------------------
  AreaSpec = GhostedAreaSpec ( EtaGhostCut, pars.GhostRepeat, pars.GhostArea );
  AreaDef = AreaDefinition (fastjet::active_area_explicit_ghosts, AreaSpec);
  JetDef    = JetDefinition( pars.LargeJetAlgorithm, pars.R );

  SelectClose = fastjet::SelectorCircle( pars.R );

  cout << " R = " << pars.R << endl;
  cout << " Original jet algorithm : "<< pars.LargeJetAlgorithm << endl;
  cout << " PtJetMin = " << pars.PtJetMin << endl;
  cout << " PtJetMax = " << pars.PtJetMax << endl;
  cout << " PtConsMin = " << pars.PtConsMin << endl;
  cout << " PtConsMax = " << pars.PtConsMax << endl;
  cout << " Constituent eta cut = " << pars.EtaConsCut << endl;
  cout << " Jet eta cut = " << EtaJetCut << endl;
  cout << " Ghosts out to eta = " << EtaGhostCut << endl;
  cout << " Reading tree named \""<< pars.ChainName << "\" from " << pars.InputName << endl;
  cout << " intype = " << pars.intype << endl;
  cout << " Writing to " << pars.OutFileName << endl;
  cout << " ----------------------------" << endl;

  // Provide a gaussian for track pt smearing
  // ----------------------------------------
  SmearPt = new TF1( "SmearPt","gaus(0)",-1,1);
}
//----------------------------------------------------------------------
ppTestAnalysis::~ppTestAnalysis(){
  if (pJA){
    delete pJA; pJA=0;
  }  
}

//----------------------------------------------------------------------
bool ppTestAnalysis::InitChains(){

  // For trees of TStarJetVector
  // (like a previous result)
  // ------------------------
  Events = new TChain(pars.ChainName);
  Events->Add(pars.InputName);
  if ( NEvents<0 ) NEvents = INT_MAX;

  if ( pars.intype==INTREE || pars.intype==MCTREE ){
    assert ( Events->GetEntries()>0 && "Something went wrong loading events.");
    NEvents=min(NEvents,Events->GetEntries() );
    
    pFullEvent = new TClonesArray("TStarJetVector");
    Events->GetBranch("PythiaParticles")->SetAutoDelete(kFALSE);
    Events->SetBranchAddress("PythiaParticles", &pFullEvent);

  }
  
  // For picoDSTs
  // -------------
  if ( pars.intype==INPICO || pars.intype==MCPICO ){
    pReader = SetupReader( Events, pars );

    InitializeReader(  pReader, pars.InputName, NEvents, PicoDebugLevel, pars.HadronicCorr );
    if ( pReader && pars.InputName.Contains( "picoDst_4_5" ) ) pReader->SetUseRejectAnyway( true );

    if ( pars.intype==MCPICO )      TurnOffCuts ( pReader );

    cout << "Don't forget to set tower cuts for pReader!!" << endl;
    
  }
  
  cout << "N = " << NEvents << endl;
    
  cout << "Done initializing chains. " << endl;
  return true;
}
//----------------------------------------------------------------------
// Main routine for one event.
EVENTRESULT ppTestAnalysis::RunEvent (){
  // cout << "-----------------------" << endl;
  // cout << "Entering PpZgAnalysis::RunEvent " << endl;
  TStarJetVector* sv;

  TStarJetPicoEventHeader* header=0;

  UInt_t filehash = 0;
  TString cname = "";
  
  // Reset results (from last event)
  // -------------------------------
  Result.clear();
  weight=1;
  particles.clear();

  switch (pars.intype) {
    // =====================================================
  case INPICO :      
  case MCPICO :      
    if ( !pReader->NextEvent() ) {
      // cout << "Can't find a next event" << endl;
      // done=true;
      pReader->PrintStatus();
      return EVENTRESULT::ENDOFINPUT;
      break;
    }
    pReader->PrintStatus(10);
    
    // cout << pReader->GetOutputContainer()->GetEntries() << endl; 
    pFullEvent = pReader->GetOutputContainer()->GetArray();
    
    header = pReader->GetEvent()->GetHeader();
    refmult=header->GetProperReferenceMultiplicity();
      
    eventid = header->GetEventId();
    runid   = header->GetRunId();

    // For GEANT: Need to devise a runid that's unique but also
    // reproducible to match Geant and GeantMc data.
    if ( pars.UseGeantNumbering ){
      TString cname = gSystem->BaseName( pReader->GetInputChain()->GetCurrentFile()->GetName() );
      UInt_t filehash = cname.Hash();
      while ( filehash > INT_MAX - 100000 ) filehash -= INT_MAX / 4; // some random large number
      if ( filehash < 1000000 ) filehash += 1000001;
      runid = filehash;
      // Sigh. Apparently also need to uniquefy the event id
      // since some are the same in the same file. Grr.
      // The following isn't great, because it uses the number in the chain, but it should get the job done
      eventid = pReader->GetNOfCurrentEvent();
    }  

    break;      
    // =====================================================
  case INTREE :
  case MCTREE :
    if ( evi>= NEvents ) {
      return EVENTRESULT::ENDOFINPUT;
      break;
    }
    
    if ( !(evi%200) ) cout << "Working on " << evi << " / " << NEvents << endl;
    Events->GetEntry(evi);

    cname = Events->GetCurrentFile()->GetName();
    eventid = Events->GetLeaf("eventid")->GetValue();
    runid = Events->GetLeaf("runid")->GetValue();

    if (pars.intype == MCTREE ){
      filehash = cname.Hash();
      while ( filehash > INT_MAX - 100000 ) filehash /= 10;
      if ( filehash < 1000000 ) filehash += 1000001;
      runid += filehash;
    }    

    ++evi;
    break;
    // =====================================================
  default:
    cerr << "Unknown intype " << pars.intype << endl;
    return EVENTRESULT::PROBLEM;
  }
    
    
  // Fill particle container
  // -----------------------
  for ( int i=0 ; i<pFullEvent->GetEntries() ; ++i ){
    sv = (TStarJetVector*) pFullEvent->At(i);

    // Ensure kinematic similarity
    if ( sv->Pt()< pars.PtConsMin )             continue;
    if ( fabs( sv->Eta() )>pars.EtaConsCut )    continue;


    //! setting pid mass for particles from pythia 
    if( pars.intype==MCPICO ){
      TParticlePDG * pid = (TParticlePDG*)PDGdb.GetParticle((Int_t)sv->mc_pdg_pid());
      double sv_mass = pid->Mass();
      double E= TMath::Sqrt(sv->P()*sv->P() + sv_mass*sv_mass);
      sv->SetPxPyPzE(sv->Px(), sv->Py(), sv->Pz(), E);
    }

    //! setting pion mass for charged particles and zero mass for towers 
    if ( pars.intype==INPICO ) {
      TParticlePDG * pid = (TParticlePDG*)PDGdb.GetParticle(111);
      double sv_mass;
      if(sv->GetCharge() !=0 )
	sv_mass = pid->Mass(); //! tracks get pion mass 
      else
	sv_mass = 0.0; //! towers get photon mass ~ 0 
      double E= TMath::Sqrt(sv->P()*sv->P() + sv_mass*sv_mass);

      sv->SetPxPyPzE(sv->Px(), sv->Py(), sv->Pz(), E);  
    }
        
    // TRACKS
    // ------
    if ( sv->GetCharge()!=0 ){
      // EFFICIENCY uncertainty
      // ----------------------
      Double_t mran=gRandom->Uniform(0,1);
      if ( mran>pars.FakeEff )  {
	continue;
      }
      if ( SigmaPt ){
	// cout << sv->Pt() << "  " << SigmaPt->Eval( sv->Pt() ) << endl;
	SmearPt->SetParameters ( 1,0,SigmaPt->Eval( sv->Pt() ));
	double ptscaler=1 + SmearPt->GetRandom();
	(*sv) *= ptscaler;
      }      
    }

    // TOWERS
    // ------
    // Shift gain
    if ( !sv->GetCharge()) {
      (*sv) *= pars.fTowScale; // for systematics
    }
 
    particles.push_back( PseudoJet (*sv ) );
    particles.back().set_user_info ( new JetAnalysisUserInfo( 3*sv->GetCharge(),"",sv->GetTowerID() ) );


  }
  // cout << pFullEvent->GetEntries() << "  " <<  particles.size() << endl;
  
  if ( particles.size()==0 ) return EVENTRESULT::NOCONSTS;

  // For pythia, use cross section as weight
  // ---------------------------------------
  if ( TParameter<double>* sigmaGen=(TParameter<double>*) Events->GetCurrentFile()->Get("sigmaGen") ){
    weight=sigmaGen->GetVal();
  }
  
  // For GEANT data
  // --------------
  if ( pars.InputName.Contains("Geant") ){
    TString currentfile = pReader->GetInputChain()->GetCurrentFile()->GetName();
    weight=LookupXsec ( currentfile );
    if ( fabs(weight-1)<1e-4){
      throw std::runtime_error("mcweight unchanged!");
    }
  }

  if ( pars.InputName.Contains("EmbedPythiaRun12pp200") ){
    TString currentfile = pReader->GetInputChain()->GetCurrentFile()->GetName();
    weight=LookupRun12Xsec ( currentfile );
    if ( fabs(weight-1)<1e-4){
      throw std::runtime_error("mcweight unchanged!");
    }
  }

  
  // Run analysis
  // ------------
  if (pJA){
    delete pJA; pJA=0;
  }
  pJA = new JetAnalyzer( particles, JetDef ) ;

  JetAnalyzer& JA = *pJA;      
  vector<PseudoJet> JAResult = sorted_by_pt( select_jet ( JA.inclusive_jets() ) );
  
  if ( JAResult.size()==0 ) {
    // cout << "Nothing found" << endl;
    return EVENTRESULT::NOJETS;
  } 

  int njets = JAResult.size();
  // cout << "-----------------------" << endl;

  for (unsigned ijet = 0; ijet < JAResult.size(); ijet++) {
    PseudoJet& CurrentJet = JAResult[ijet];
    PseudoJet NeutralPart = join ( OnlyNeutral( CurrentJet.constituents() ) );
    PseudoJet ChargedPart = join ( OnlyCharged( CurrentJet.constituents() ) );
    double q = 0;
    for ( PseudoJet& c : ChargedPart.constituents() ){
      q+= c.user_info<JetAnalysisUserInfo>().GetQuarkCharge();
    }
    
    JetAnalysisUserInfo* userinfo = new JetAnalysisUserInfo( q );
    // Save neutral energy fraction in multi-purpose field
    userinfo->SetNumber(NeutralPart.pt()  / CurrentJet.pt());
    CurrentJet.set_user_info ( userinfo );

    if ( pars.MaxJetNEF<1.0 &&  NeutralPart.pt()/CurrentJet.pt() > pars.MaxJetNEF ) continue;

      vector<PseudoJet> constituents = sorted_by_pt(CurrentJet.constituents());
      int nparticles = CurrentJet.constituents().size();
      if (nparticles == 0) continue;
      float pTlead = constituents[0].pt();
      double pT_lead0;
      double pT_lead3;
      double pT_lead5;
      double pT_lead7;
      if(pTlead > 0) {
          pT_lead0 = CurrentJet.pt();
      }
      if(pTlead > 3) {
          pT_lead3 = CurrentJet.pt();
      }
      if(pTlead > 5) {
          pT_lead5 = CurrentJet.pt();
      }
      if(pTlead > 7) {
          pT_lead7 = CurrentJet.pt();
      }
      Result.push_back ( ResultStruct ( CurrentJet, pT_lead0, pT_lead3, pT_lead5, pT_lead7) );

  }  
  // By default, sort for original jet pt
  sort ( Result.begin(), Result.end(), ResultStruct::origptgreater);
  
  return EVENTRESULT::JETSFOUND;
}
//----------------------------------------------------------------------
void InitializeReader(  std::shared_ptr<TStarJetPicoReader> pReader, const TString InputName, const Long64_t NEvents,
			const int PicoDebugLevel, const double HadronicCorr ){

  TStarJetPicoReader& reader = *pReader;

  if ( HadronicCorr < 0 ){
    reader.SetApplyFractionHadronicCorrection(kFALSE);
    reader.SetApplyMIPCorrection(kTRUE);
    reader.SetRejectTowerElectrons( kTRUE );
  } else {
    reader.SetApplyFractionHadronicCorrection(kTRUE);
    reader.SetFractionHadronicCorrection( HadronicCorr );
    reader.SetApplyMIPCorrection(kFALSE);
    reader.SetRejectTowerElectrons( kFALSE );
  }
    
  // Run 11: Use centrality cut
  if ( InputName.Contains("NPE") ){
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetReferenceCentralityCut (  6, 8 ); // 6,8 for 0-20%
  }
        
  reader.Init(NEvents);
  TStarJetPicoDefinitions::SetDebugLevel(PicoDebugLevel);
}
//----------------------------------------------------------------------
// Helper to deal with repetitive stuff
shared_ptr<TStarJetPicoReader> SetupReader ( TChain* chain, const ppTestParameters& pars ){
  TStarJetPicoDefinitions::SetDebugLevel(0); // 10 for more output

  shared_ptr<TStarJetPicoReader> pReader = make_shared<TStarJetPicoReader>();
  TStarJetPicoReader& reader = *pReader;
  reader.SetInputChain (chain);

  // Event and track selection
  // -------------------------
  TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
  evCuts->SetTriggerSelection( pars.TriggerName ); //All, MB, HT, pp, ppHT, ppJP
  // Additional cuts 
  evCuts->SetVertexZCut (pars.VzCut);
  evCuts->SetRefMultCut ( pars.RefMultCut );
  evCuts->SetVertexZDiffCut( pars.VzDiffCut );
  evCuts->SetMaxEventPtCut ( pars.MaxEventPtCut );
  evCuts->SetMaxEventEtCut ( pars.MaxEventEtCut );

  evCuts->SetMinEventEtCut ( pars.MinEventEtCut );

  std::cout << "Using these event cuts:" << std::endl;
  std::cout << " Vz: " << evCuts->GetVertexZCut() << std::endl;
  std::cout << " Refmult: " << evCuts->GetRefMultCutMin() << " -- " << evCuts->GetRefMultCutMax() << std::endl;
  std::cout << " Delta Vz:  " << evCuts->GetVertexZDiffCut() << std::endl;
  std::cout << " MaxEventPt:  " << evCuts->GetMaxEventPtCut() << std::endl;
  std::cout << " MaxEventEt:  " << evCuts->GetMaxEventEtCut() << std::endl;
  
  // This method does NOT WORK for GEANT MC trees because everything is in the tracks...
  // Do it by hand later on, using pars.ManualHtCut;
  // Also doesn't work for general trees, but there it can't be fixed
    
  // // TESTING ONLY:
  // evCuts->SetMaxEventPtCut ( 20000000. );
  // evCuts->SetMaxEventEtCut ( 20000000. );
  // evCuts->SetReferenceCentralityCut (  6, 8 ); // 6,8 for 0-20%
  // evCuts->SetMinEventEtCut ( -1.0 );
  // evCuts->SetMinEventEtCut ( 6.0 );


  // Tracks cuts
  TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
  trackCuts->SetDCACut( pars.DcaCut );
  trackCuts->SetMinNFitPointsCut( pars.NMinFit );
  trackCuts->SetFitOverMaxPointsCut( pars.FitOverMaxPointsCut );
  trackCuts->SetMaxPtCut ( pars.MaxTrackPt );

  std::cout << "Using these track cuts:" << std::endl;
  std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
  std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
  std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
  std::cout << " maxpt : " << trackCuts->GetMaxPtCut (  ) << std::endl;
  
  // Towers
  TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
  towerCuts->SetMaxEtCut(pars.MaxEtCut);

  std::cout << "Using these tower cuts:" << std::endl;
  std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
  std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
  
  // V0s: Turn off
  reader.SetProcessV0s(false);

  // cout << "towerCuts=" << towerCuts << endl;
  // cout << "trackCuts=" << trackCuts << endl;
  return pReader;
}
//----------------------------------------------------------------------
void TurnOffCuts ( std::shared_ptr<TStarJetPicoReader> pReader ){


  pReader->SetProcessTowers ( false );
  TStarJetPicoEventCuts* evCuts = pReader->GetEventCuts();
  evCuts->SetTriggerSelection( "All" ); //All, MB, HT, pp, ppHT, ppJP
  // No Additional cuts -- but keep the Vz cut
  // evCuts->SetVertexZCut ( 30 );
  // evCuts->SetVertexZCut (99999);
  evCuts->SetRefMultCut (0);
  evCuts->SetVertexZDiffCut(999999);
  
  evCuts->SetMaxEventPtCut (99999);
  evCuts->SetMaxEventEtCut (99999);
  evCuts->SetPVRankingCutOff();		//  Use SetPVRankingCutOff() to turn off vertex ranking cut.  default is OFF
  
  // Tracks cuts
  TStarJetPicoTrackCuts* trackCuts = pReader->GetTrackCuts();
  trackCuts->SetDCACut(99999);
  trackCuts->SetMinNFitPointsCut(-1);
  trackCuts->SetFitOverMaxPointsCut(-1);
  trackCuts->SetMaxPtCut (99999);
  
  // Towers: should be no tower in MC. All (charged or neutral) are handled in track
  TStarJetPicoTowerCuts* towerCuts = pReader->GetTowerCuts();
  towerCuts->SetMaxEtCut(99999);

  cout << " TURNED OFF ALL CUTS" << endl;
}
  
  
//----------------------------------------------------------------------
double LookupXsec( TString filename ){
  // Some data for geant
  // -------------------
  //cross-sections for simulated GEANT data sample
  // From Renee.
  // also available via
  // http://people.physics.tamu.edu/sakuma/star/jets/c101121_event_selection/s0150_mclist_001/web.php
  // Double_t MinbXsec=28.12;
  // Double_t Xsec[12];
  // Xsec[0]=28.11;//2
  // Xsec[1]=1.287;//3
  // Xsec[2]=0.3117;//4
  // Xsec[3]=0.1360;//5
  // Xsec[4]=0.02305;//7
  // Xsec[5]=0.005494;//9
  // Xsec[6]=0.002228;//11
  // Xsec[7]=0.0003895;//15
  // Xsec[8]=0.00001016;//25
  // Xsec[9]=0.0000005010;//35
  // Xsec[10]=0.0000000283;//45
  // Xsec[11]=0.000000001443;//55
  
  static const Double_t MinbXsec=28.12;
  // static const Double_t Xsec[12] = {
  //   28.11,		// 2-3
  //   1.287,		// 3-4
  //   0.3117,		// 4-5
  //   0.1360,		// 5-7
  //   0.02305,		// 7-9
  //   0.005494,		// 9-11
  //   0.002228,		// 11-15
  //   0.0003895,		// 15-25
  //   0.00001016,		// 25-35
  //   0.0000005010,	// 35-45
  //   0.0000000283,	// 45-55
  //   0.000000001443	// 55-65
  // };

  static const Double_t Xsec[12] = {
    1.0,        // Placeholder for 2-3
    1.30E+09,	// 3-4
    3.15E+08,	// 4-5
    1.37E+08,	// 5-7
    2.30E+07,	// 7-9
    5.53E+06,	// 9-11
    2.22E+06,	// 11-15
    3.90E+05,	// 15-25
    1.02E+04,	// 25-35
    5.01E+02,	// 35-45
    2.86E+01,	// 45-55
    1.46E+00	// 55-65
  };

  static const Double_t Nmc[12] = {
    1,			// 2-3
    672518,		// 3-4
    672447,		// 4-5
    393498,		// 5-7
    417659,		// 7-9
    412652,		// 9-11
    419030,		// 11-15
    396744,		// 15-25
    399919,		// 25-35
    119995,		// 35-45
    117999,		// 45-55
    119999		// 55-65
  };

  Double_t w[12];
  for ( int i=0; i<12 ; ++i ){
    w[i] = Xsec[i] / Nmc[i];
    // w[i] = Nmc[i] / Xsec[i] ;
  }

  // static const Double_t w[12] = {
  //   1,			// Placeholder
  //   1.90E+03,
  //   6.30E+02,
  //   3.43E+02,
  //   5.49E+01,
  //   1.33E+01,
  //   5.30E+00,
  //   9.81E-01,
  //   2.56E-02,
  //   4.56E-03,
  //   2.43E-04,
  //   1.20E-05
  // };
    
  if ( filename.Contains("picoDst_3_4") ) return w[1];
  if ( filename.Contains("picoDst_4_5") ) return w[2];
  if ( filename.Contains("picoDst_5_7") ) return w[3];
  if ( filename.Contains("picoDst_7_9") ) return w[4];
  if ( filename.Contains("picoDst_9_11") ) return w[5];
  if ( filename.Contains("picoDst_11_15") ) return w[6];
  if ( filename.Contains("picoDst_15_25") ) return w[7];
  if ( filename.Contains("picoDst_25_35") ) return w[8];
  if ( filename.Contains("picoDst_35_45") ) return w[9];
  if ( filename.Contains("picoDst_45_55") ) return w[10];
  if ( filename.Contains("picoDst_55_65") ) return w[11];

  return 1;

}

//----------------------------------------------------------------------
double LookupRun12Xsec( TString filename ){
  
  const int NUMBEROFPT = 11;
  // const char *PTBINS[NUMBEROFPT]={"2_3","3_4","4_5","5_7","7_9","9_11","11_15","15_20","20_25","25_35","35_-1"};
  const static float XSEC[NUMBEROFPT] = {9.00581646, 1.461908221, 0.3544350863, 0.1513760388, 0.02488645725, 0.005845846143, 0.002304880181, 0.000342661835, 4.562988397e-05, 9.738041626e-06, 5.019978175e-07};
  const static float NUMBEROFEVENT[NUMBEROFPT] = {2100295, 600300, 600300, 300289, 300289, 300289, 160295, 100302, 80293, 76303, 23307};

  const static vector<string> vptbins={"pp12Pico_pt2_3","pp12Pico_pt3_4","pp12Pico_pt4_5","pp12Pico_pt5_7","pp12Pico_pt7_9","pp12Pico_pt9_11","pp12Pico_pt11_15","pp12Pico_pt15_20","pp12Pico_pt20_25","pp12Pico_pt25_35","35_-1"};
  for ( int i=0; i<vptbins.size(); ++i ){
    if ( filename.Contains(vptbins.at(i).data())) return XSEC[i] / NUMBEROFEVENT[i];
  }

  throw std::runtime_error("Not a valid filename");
  return -1;
  
}

//----------------------------------------------------------------------
/*
    convenient output
*/
ostream & operator<<(ostream & ostr, const PseudoJet& jet) {
  if (jet == 0) {
    ostr << " 0 ";
  } else {
    ostr << " pt = " << jet.pt()
         << " m = " << jet.m()
         << " y = " << jet.rap()
         << " phi = " << jet.phi()
         << " ClusSeq = " << (jet.has_associated_cs() ? "yes" : "no");
  }
  return ostr;
}
//----------------------------------------------------------------------

// ------------------------------------------------------------
// DELETE ME, ONLY FOR EXPERIMENTING

vector<PseudoJet> reshuffle ( const vector<PseudoJet> orig ){

  // vector< pair<double,double> > phieta;
  // vector<PseudoJet> shuffled;
  // for ( auto j : orig ){
  //   phieta.push_back ( pair<double,double> (j.phi(), j.eta()) );
  //   shuffled.push_back ( j );
  // }

  // std::random_device rd;
  // std::mt19937 g(rd());
  // std::shuffle ( phieta.begin(), phieta.end(), g );
  // // std::random_shuffle ( phieta.begin(), phieta.end() );
    
  // for ( int i=0 ; i<shuffled.size() ; ++i ){
  //   PseudoJet& j = shuffled.at(i);
  //   j.reset_PtYPhiM ( j.pt(), phieta.at(i).second, phieta.at(i).first, j.m() );
  // }

  // vector< double > pt;
  // vector<PseudoJet> shuffled;
  // for ( auto j : orig ){
  //   pt.push_back ( j.pt() );
  //   shuffled.push_back ( j );
  // }

  // std::random_device rd;
  // std::mt19937 g(rd());
  // std::shuffle ( pt.begin(), pt.end(), g );

  // for ( int i=0 ; i<shuffled.size() ; ++i ){
  //   PseudoJet& j = shuffled.at(i);
  //   j.reset_PtYPhiM ( pt.at(i), j.eta(), j.phi(), j.m() );
  // }

  vector< double > pt;
  for ( auto j : orig ){
    pt.push_back ( j.pt() );
  }

  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle ( pt.begin(), pt.end(), g );

  vector<PseudoJet> shuffled;

  for ( int i=0 ; i<orig.size() ; ++i ){
    PseudoJet j = orig.at(i);
    j.reset_PtYPhiM ( pt.at(i), j.eta(), j.phi(), j.m() );
    shuffled.push_back( j );
  }

  return shuffled;
}


