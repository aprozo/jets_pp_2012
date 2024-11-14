/* @file ppAnalysis.hh
    @author Raghav Kunnawalkam Elayavalli
    @version Revision 1.0
    @brief analysis class
    @details Uses JetAnalyzer objects
    @date March 16, 2022
*/

#ifndef __PPANALYSIS_HH
#define __PPANALYSIS_HH

#include "ppParameters.hh"
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

double LookupRun12Xsec(TString filename);

/*
   For sorting with a different key
*/
typedef pair<PseudoJet, double> PseudoJetPt;
struct PseudoJetPtGreater
{
  bool operator()(PseudoJetPt const &a, PseudoJetPt const &b)
  {
    return a.second > b.second;
  }
};

/*
  To keep original and groomed jets connected
 */
class ResultStruct
{
public:
  PseudoJet orig;
  ResultStruct(PseudoJet orig) : orig(orig) {};

  static bool origptgreater(ResultStruct const &a, ResultStruct const &b)
  {
    return a.orig.pt() > b.orig.pt();
  };
};

/*
    convenient output
*/
ostream &operator<<(ostream &ostr, const PseudoJet &jet);

/*
    Helper for chains
 */
void InitializeReader(std::shared_ptr<TStarJetPicoReader> pReader, const TString InputName, const Long64_t NEvents,
                      const int PicoDebugLevel, const double HadronicCorr = 0.999999);

static const Selector NotGhost = !fastjet::SelectorIsPureGhost();                                           ///< Helper useful outside the class as well
static const Selector OnlyCharged = NotGhost && (SelectorChargeRange(-3, -1) || SelectorChargeRange(1, 3)); ///< Helper useful outside the class as well
static const Selector OnlyNeutral = NotGhost && SelectorChargeRange(0, 0);                                  ///< Helper useful outside the class as well

/*
   The main class
 */
class ppAnalysis
{

private:
  // These need to be initialized
  // ----------------------------
  ppParameters pars; ///< container to have all analysis parameters in one place

  // Internal
  // --------
  float EtaJetCut;   ///< jet eta
  float EtaGhostCut; ///< ghost eta

  TDatabasePDG PDGdb;

  fastjet::JetDefinition JetDef; ///< jet definition

  // Relevant jet candidates
  fastjet::Selector select_jet_eta; ///< jet rapidity selector
  fastjet::Selector select_jet_pt;  ///< jet p<SUB>T</SUB> selector
  fastjet::Selector select_jet;     ///< compound jet selector

  // Data
  // ----
  Long64_t NEvents = -1;
  TChain *Events = 0;
  TClonesArray *pFullEvent = 0; ///< Constituents
  TStarJetVector *pHT = 0;      ///< the trigger (HT) object, if it exists

  vector<PseudoJet> particles;
  vector<PseudoJet> partons;
  double rho = 0; ///< background density

  std::shared_ptr<TStarJetPicoReader> pReader = 0;

  Long64_t evi = 0;

  int PicoDebugLevel = 0; /// Control DebugLevel in picoDSTs
  int eventid;
  int runid;
  int runid1;
  float event_sum_pt;
  double vz;
  int mult;
  double refmult;
  double weight;
  bool is_rejected;
  int njets;

  JetAnalyzer *pJA = 0;

  vector<ResultStruct> Result; ///< result in a nice structured package

public:
  ppAnalysis(const int argc, const char **const);

  /* Destructor. Clean things up
   */
  virtual ~ppAnalysis();

  /* Decoupled chain initialization for readability
   */
  bool InitChains();

  /* Main routine for one event.
      \return false if at the end of the chain
   */
  EVENTRESULT RunEvent();

  // Getters and Setters
  // -------------------
  inline ppParameters &GetPars() { return pars; };

  /// Get jet radius
  inline double GetR() { return pars.R; };

  /// Set jet radius
  inline void SetR(const double newv) { pars.R = newv; };

  /// Handle to pico reader
  inline std::shared_ptr<TStarJetPicoReader> GetpReader() { return pReader; };

  /// Get the weight of the current event (mainly for PYTHIA)
  inline double GetEventWeight() { return weight; };

  /// Get the refmult of the current event
  inline double GetRefmult() { return refmult; };

  /// Get the runid of the current event (this id for geant events can match to bad run ids)
  inline double GetRunid1() { return runid1; };

  /// Get the runid of the current event
  inline double GetRunid() { return runid; };

  /// Get the eventid of the current event
  inline double GetEventid() { return eventid; };

  inline float GetEventSumPt() { return event_sum_pt; };

  /// Get the vz of the current event
  inline double GetVz() { return vz; };

  inline int GetEventMult() { return mult; };

  inline int GetRejectCode() { return is_rejected; };

  /// Get the Trigger (HT) object if it exists, for matching
  inline TStarJetVector *GetTrigger() const { return pHT; };

  /// The main result of the analysis
  inline const vector<ResultStruct> &GetResult() { return Result; }
};

shared_ptr<TStarJetPicoReader> SetupReader(TChain *chain, const ppParameters &pars);

/* For use with GeantMc data
 */
void TurnOffCuts(std::shared_ptr<TStarJetPicoReader> pReader);

#endif // __PPANALYSIS_HH
