/* @file ppTestAnalysis.hh
    @author Raghav Kunnawalkam Elayavalli 
    @version Revision 1.0
    @brief test analysis class 
    @details Uses JetAnalyzer objects
    @date March 16, 2022
*/


#ifndef __PPTESTANALYSIS_HH
#define __PPTESTANALYSIS_HH

#include "ppTestParameters.hh"
#include "JetAnalyzer.hh"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TRandom.h"
#include "TChain.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TFile.h"
#include "TSystem.h"
#include "TParameter.h"
#include "TClonesArray.h"

#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

// Not needed for analysis per se
#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"

#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"

#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <assert.h>
#include <iostream>
#include <cmath>
#include <climits>
#include <sstream>

using namespace std;
using namespace fastjet;
using namespace contrib;

#include <random>
#include <algorithm>
vector<PseudoJet> reshuffle ( const vector<PseudoJet> orig );

/* 
    A helper for geant data
 */
double LookupXsec( TString filename );

/* 
    A helper for Run12 embedding data
 */
double LookupRun12Xsec( TString filename );


/*
   For sorting with a different key
*/
typedef pair<PseudoJet,double> PseudoJetPt;
struct PseudoJetPtGreater {
  bool operator()( PseudoJetPt const & a, PseudoJetPt const & b) { 
    return a.second > b.second ;
  }
};

/*
  To keep original and groomed jets connected
 */
class ResultStruct{
public:
    double pT_lead0;
    double pT_lead3;
    double pT_lead5;
    double pT_lead7;
  PseudoJet orig;
  ResultStruct ( PseudoJet orig, double pT_lead0, double pT_lead3, double pT_lead5, double pT_lead7)  :
    orig(orig),
    pT_lead0(pT_lead0),
    pT_lead3(pT_lead3),
    pT_lead5(pT_lead5),
    pT_lead7(pT_lead7){};
  
  static bool origptgreater( ResultStruct const & a, ResultStruct const & b) { 
    return a.orig.pt()>b.orig.pt();
  };
  
};


/* 
    convenient output
*/
ostream & operator<<(ostream & ostr, const PseudoJet& jet);

/* 
    Helper for chains
 */
void InitializeReader(  std::shared_ptr<TStarJetPicoReader> pReader, const TString InputName, const Long64_t NEvents,
			const int PicoDebugLevel, const double HadronicCorr=0.999999 );

static const Selector NotGhost = !fastjet::SelectorIsPureGhost ();     ///< Helper useful outside the class as well
static const Selector OnlyCharged = NotGhost && ( SelectorChargeRange( -3, -1) || SelectorChargeRange( 1, 3) );     ///< Helper useful outside the class as well
static const Selector OnlyNeutral = NotGhost && SelectorChargeRange( 0, 0);     ///< Helper useful outside the class as well


/*
   The main class
 */
class ppTestAnalysis {

private :

  // These need to be initialized
  // ----------------------------
  ppTestParameters pars;   ///< container to have all analysis parameters in one place
  
  // Internal
  // --------  
  float EtaJetCut;       ///< jet eta
  float EtaGhostCut;     ///< ghost eta

  TDatabasePDG PDGdb;
  
  fastjet::JetDefinition JetDef;       ///< jet definition

  // Relevant jet candidates
  fastjet::Selector select_jet_eta;        ///< jet rapidity selector
  fastjet::Selector select_jet_pt;         ///< jet p<SUB>T</SUB> selector
  fastjet::Selector select_jet;            ///< compound jet selector

  fastjet::GhostedAreaSpec AreaSpec;      ///< ghosted area specification
  fastjet::AreaDefinition AreaDef;        ///< jet area definition

  // Data
  // ----
  Long64_t NEvents=-1;
  TChain* Events = 0;
  TClonesArray* pFullEvent = 0;             ///< Constituents
  TClonesArray* pHardPartons = 0;           ///< For pythia data, the original hard scatter
  TClonesArray* pHardPartonNames = 0;       ///< For pythia data, the original hard scatter names (ubar, g, ...)
  TStarJetVector* pHT = 0;                  ///< the trigger (HT) object, if it exists

  vector<PseudoJet> particles;
  vector<PseudoJet> partons;
  double rho=0;                             ///< background density
  
  std::shared_ptr<TStarJetPicoReader> pReader=0;
  
  Long64_t evi=0;
  
  int PicoDebugLevel=0; /// Control DebugLevel in picoDSTs
  // int PicoDebugLevel=1; /// Control DebugLevel in picoDSTs

  int eventid;
  int runid;
 
  double refmult=0;
 
  double weight=0;    

 
  JetAnalyzer *pJA=0;

  /// Resolution Smearing width for tracks
  /// This function should return sigma (pT - pT_truth) for a given pT
  TF1 * SigmaPt;
  /// Actual smearing factor
  /// Gaussian around 1 with the width from SigmaPt
  TF1 * SmearPt;

  // For matching
  // ------------
  fastjet::Selector SelectClose;

  JetAnalyzer* pJAhi;                      ///< JetAnalyzer object for high pT
  JetAnalyzer* pJAlo;                      ///< JetAnalyzer object for low pT
  
  // --------------------------------------------------------------
  // JetAlgorithm ReclusterJetAlgorithm = fastjet::cambridge_algorithm;
  fastjet::JetAlgorithm ReclusterJetAlgorithm;
  fastjet::JetDefinition ReclusterJetDef;
  //Recluster * recluster=0; 

  vector<ResultStruct> Result; ///< result in a nice structured package
    
public:  

  /* Standard constructor. Parse vectorized argv.
      \param argc: number of arguments
      \param argv: string array of command line options
   */
  ppTestAnalysis ( const int argc, const char** const  );

  /* Destructor. Clean things up
   */
  virtual ~ppTestAnalysis();

  /* Decoupled chain initialization for readability
   */
  bool InitChains ( );

  /* Main routine for one event.
      \return false if at the end of the chain      
   */
  EVENTRESULT RunEvent ( );
 
  // /** This little helper is true if there's at least one 10 GeV jet
  //  */
  // bool Has10Gev;
  // double OverrideJetMin;        ///< quick and dirty way to override the 10 GeV minimum

  // Quick and dirty QA histos - keep them public
  // --------------------------------------------
  TH1D* QA_Tower;
  TH2D* QA_TowerEt;
  TH2D* QA_TowerEta;
  TH2D* QA_TowerPhi;
  
  

  // Getters and Setters
  // -------------------
  inline ppTestParameters& GetPars()         { return pars; };

  /// Get jet radius
  inline double GetR ( )                   { return pars.R; };
  /// Set jet radius
  inline void   SetR ( const double newv ) { pars.R=newv;   };


  /// Handle to pico reader
  inline std::shared_ptr<TStarJetPicoReader> GetpReader()  { return pReader; };

  /// Get the weight of the current event (mainly for PYTHIA)
  inline double GetEventWeight()           { return weight; };

  /// Get the refmult of the current event
  inline double GetRefmult()           { return refmult; };

  /// Get the runid of the current event
  inline double GetRunid()           { return runid; };

  /// Get the eventid of the current event
  inline double GetEventid()           { return eventid; };

  /// Get the Trigger (HT) object if it exists, for matching
  inline TStarJetVector* GetTrigger() const        { return pHT; };
  
  /// The main result of the analysis
  inline const vector<ResultStruct>& GetResult()    { return Result; }
  
  inline double GetRho()                             { return rho; }

  inline const vector<PseudoJet>& GetParticles() const     { return particles; }
  // void SetSigmaPt( TF1* f ){  SigmaPt = f; };

  // /// Get minimum jet p<SUB>T</SUB>
  // inline double GetJet_ptmin ( )                   { return pars.jet_ptmin; };
  // /// Set minimum jet p<SUB>T</SUB>
  // inline void   SetJet_ptmin ( const double newv ) { pars.jet_ptmin=newv;   };

  // /// Get maximum jet p<SUB>T</SUB>
  // inline double GetJet_ptmax ( )                   { return pars.jet_ptmax; };
  // /// Set maximum jet p<SUB>T</SUB>
  // inline void   SetJet_ptmax ( const double newv ) { pars.jet_ptmax=newv;   };

  // /// Get jet rapidity acceptance
  // inline double GetMax_rap ( )                   { return pars.max_rap; };
  // /// Set jet rapidity acceptance
  // inline void   SetMax_rap ( const double newv ) { pars.max_rap=newv;   };

  // /// Get ghosted area rapidity cut, should be >= max_rap + 2*R
  // inline double GetGhost_maxrap ( )                   { return pars.ghost_maxrap; };
  // /// Set ghosted area rapidity cut, should be >= max_rap + 2*R
  // inline void   SetGhost_maxrap ( const double newv ) { pars.ghost_maxrap=newv;   };


  // // Objects will be handed by _reference_! Obviates need for setter
  // /// Handle to jet definition
  // inline fastjet::JetDefinition& GetJetDef () { return JetDef; }
  // /// Handle to selector for low  p<SUB>T</SUB> constituents
  // inline fastjet::Selector& GetLoConsSelector () { return slo; }
  // /// Handle to selector for high  p<SUB>T</SUB> constituents
  // inline fastjet::Selector& GetHiConsSelector () { return shi; }
  
  /// Handle to selector for jet candidates
  inline fastjet::Selector& GetJetSelector () { return select_jet; }

  /// Handle to ghosted area specification
  inline fastjet::GhostedAreaSpec& GetAreaSpec () { return AreaSpec; }
  /// Handle to jet area definition
  inline fastjet::AreaDefinition& GetAreaDef () { return AreaDef; }

  /// Handle to constituents
  inline std::vector<fastjet::PseudoJet> GetConstituents() {return particles; };

   
  //   TStarJetPicoTowerCuts* towerCuts = pReader->GetTowerCuts();
  //   cout << "DumbTest: towerCuts->MaxEt="<< towerCuts->GetMaxEtCut() << endl;
  //   cout << "DumbTest: towerCuts=" << towerCuts << endl;
    
  //   cout << endl;
  //   TStarJetPicoTrackCuts* trackCuts = pReader->GetTrackCuts();
  //   cout << "DumbTest: trackCuts->MaxPt="<< trackCuts->GetMaxPtCut() << endl;
  //   cout << "DumbTest: trackCuts=" << trackCuts << endl;
    
  //   cout << endl;
  //   cout << "DumbTest: pReader->N=" << pReader->GetNOfEvents() << endl;

  //   pReader->ReadEvent( 0 );
    
  //   // for (int i=0; i<10; ++i){
  //   //   pReader->ReadEvent( 200 +i );
  //   //   cout << 200+i << "  " << pReader->GetEvent()->GetHeader()->GetEventId() << endl;
  //   // }
  //   return;
    
  // };

  // /// Handle to JetAnalyzer for high pT constituents
  // inline JetAnalyzer* GetJAhi() {return pJAhi; };
  // /// Handle to JetAnalyzer for low pT constituents
  // inline JetAnalyzer* GetJAlo() {return pJAlo; };

  // /// Handle to unaltered clustering result with high pT constituents
  // inline std::vector<fastjet::PseudoJet> GetJAhiResult() {return JAhiResult; };
  // /// Handle to unaltered clustering result with low pT constituents
  // inline std::vector<fastjet::PseudoJet> GetJAloResult() {return JAloResult; };

  // /// Handle to high pT constituents
  // inline std::vector<fastjet::PseudoJet> GetLoConstituents() {return pLo; };
  // /// Handle to low pT constituents
  // inline std::vector<fastjet::PseudoJet> GetHiConstituents() {return pHi; };

  
  // /// Handle to Dijet result with high pT constituents
  // inline std::vector<fastjet::PseudoJet> GetDiJetsHi() {return DiJetsHi; };
  // /// Handle to Dijet result with low pT constituents
  // inline std::vector<fastjet::PseudoJet> GetDiJetsLo() {return DiJetsLo; };
  
};  

/* Helper to perform the TStarJetPicoReader initialization
 */
TStarJetPicoReader GetReader ( TString ChainPattern="~putschke/Data/Pico_ppHT/*.root", 
			       TString TriggerString="ppHT",
			       TString ChainName="JetTree",
			       const double RefMultCut=0
			       );

/* Slightly different, preferred version of GetReader
 */
shared_ptr<TStarJetPicoReader> SetupReader ( TChain* chain, const ppTestParameters& pars );


/* For use with GeantMc data
 */
void TurnOffCuts ( std::shared_ptr<TStarJetPicoReader> pReader );
  
#endif // __PPTESTANALYSIS_HH
