/* @file ppAnalysis.hh
    @author Raghav Kunnawalkam Elayavalli
    @brief Stage-1 analysis class: reads TStarJetPico events, finds jets with
           JetAnalyzer and exposes per-event / per-jet results to RunppAna.
*/

#ifndef __PPANALYSIS_HH
#define __PPANALYSIS_HH

#include "JetAnalyzer.hh"
#include "ppParameters.hh"

#include "TChain.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TRandom.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoReader.h"

#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTriggerInfo.h"

#include "TStarJetPicoUtils.h"
#include "TStarJetVector.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetPicoTower.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <assert.h>
#include <climits>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;
using namespace fastjet;
using namespace contrib;

#include <algorithm>
#include <random>

/* For sorting with a different key */
typedef pair<PseudoJet, double> PseudoJetPt;
struct PseudoJetPtGreater {
   bool operator()(PseudoJetPt const &a, PseudoJetPt const &b) { return a.second > b.second; }
};

/* One jet plus the area / UE quantities computed for it */
class ResultStruct {
public:
   PseudoJet orig;
   double area = 0.0;          ///< jet area from ClusterSequenceArea
   double bg_density = 0.0;    ///< off-axis-cones UE density (averaged over the two ±π/2 cones), GeV / unit area
   double pt_corrected = 0.0;  ///< orig.perp() - bg_density * area
   ResultStruct(PseudoJet orig) : orig(orig) {};
   static bool origptgreater(ResultStruct const &a, ResultStruct const &b) { return a.orig.pt() > b.orig.pt(); };
};

ostream &operator<<(ostream &ostr, const PseudoJet &jet);

void InitializeReader(std::shared_ptr<TStarJetPicoReader> pReader, const TString InputName, const Long64_t NEvents,
                      const int PicoDebugLevel, const double HadronicCorr = 0.999999);

// Constituent selectors, also useful outside the class
static const Selector NotGhost = !fastjet::SelectorIsPureGhost();
static const Selector OnlyCharged = NotGhost && (SelectorChargeRange(-3, -1) || SelectorChargeRange(1, 3));
static const Selector OnlyNeutral = NotGhost && SelectorChargeRange(0, 0);

#include "JetQAHistogramManager.hh"

class ppAnalysis {

private:
   ppParameters pars; ///< container to have all analysis parameters in one place

   JetQAHistogramManager QA_hist; ///< QA histograms

   float EtaJetCut;   ///< jet eta
   float EtaGhostCut; ///< ghost eta

   TDatabasePDG PDGdb;

   fastjet::JetDefinition JetDef; ///< jet definition

   fastjet::Selector select_jet_eta; ///< jet rapidity selector
   fastjet::Selector select_jet_pt;  ///< jet p<SUB>T</SUB> selector
   fastjet::Selector select_jet;     ///< compound jet selector

   Long64_t NEvents = -1;
   TChain *Events = 0;
   TClonesArray *pFullEvent = 0; ///< Constituents
   TStarJetVector *pHT = 0;      ///< the trigger (HT) object, if it exists

   vector<PseudoJet> particles;
   vector<PseudoJet> partons;
   double rho = 0; ///< background density

   std::shared_ptr<TStarJetPicoReader> pReader = 0;

   int PicoDebugLevel = 0; /// Control DebugLevel in picoDSTs
   int eventid;
   int runid;
   int runid1;
   float event_sum_pt;
   int mult;
   double refmult;
   double weight;
   int njets;
   bool isTriggerEvent;
   // Per-event JP fire flags = header trigger ids 370601/611/621. In DATA these
   // are the real prescale-accepted hardware bits (the maker copies the nominal
   // trigger-id list verbatim).
   bool firedJP0;
   bool shouldJP0 = false, shouldJP1 = false, shouldJP2 = false;
   bool shouldHwJP0 = false, shouldHwJP1 = false, shouldHwJP2 = false;
   bool shouldHT2 = false;          // full-simulator HT2 decision (bit-8 picos)
   bool haveSimuDecision = false;   // pico carries bit-8 isTrigger() objects
   bool firedJP1;
   bool firedJP2;
   bool firedHT2 = false; // 370531 in the nominal trigger-id list
   bool firedMB = false;  // 370011 or 370001 in the nominal trigger-id list
   double vz;  ///< primary-vertex z (cm); needed for vertex-z reweighting downstream
   double pthat = -1;  ///< event pt-hat, MC picos only (MC header reference-centrality weight); -1 otherwise
   double vx = 0.0;     ///< primary-vertex x
   double vy = 0.0;     ///< primary-vertex y
   double vz_vpd = 0.0; ///< VPD-measured vertex z
   int n_vpd_east = 0;  ///< VPD east hit count (VPDMB-fired proxy: east>=1 && west>=1)
   int n_vpd_west = 0;  ///< VPD west hit count
   // Trigger-simulator ADCs on the DSM scale, for threshold variations downstream:
   // max jet-patch ADC and max high-tower ADC of the event, and the thresholds.
   int jp_adc_max = -1;
   int ht_adc_max = -1;
   int jp_thr[3] = {0, 0, 0};
   int ht_thr[4] = {0, 0, 0, 0};

   JetAnalyzer *pJA = 0;

   vector<ResultStruct> Result; ///< result in a nice structured package

public:
   ppAnalysis(const int argc, const char **const);

   virtual ~ppAnalysis();

   bool InitChains();

   /* Main routine for one event.
       \return false if at the end of the chain
    */
   EVENTRESULT RunEvent();

   inline ppParameters &GetPars() { return pars; };
   inline JetQAHistogramManager &GetHistogramManager() { return QA_hist; };

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

   /// Get the runid of the current event (for geant events this id can collide
   /// with bad run ids)
   inline double GetRunid1() { return runid1; };

   /// Get the runid of the current event
   inline double GetRunid() { return runid; };

   /// Get the eventid of the current event
   inline double GetEventid() { return eventid; };

   inline bool IsTriggerEvent() { return isTriggerEvent; };
   inline bool FiredJP0() { return firedJP0; };
   inline bool FiredJP1() { return firedJP1; };
   inline bool FiredJP2() { return firedJP2; };
   inline bool FiredHT2() { return firedHT2; };
   inline bool FiredMB() { return firedMB; };
   // offline-emulator ("shouldFire") per-event decisions, non-bit-7 objects
   inline bool ShouldJP0() { return shouldJP0; };
   inline bool ShouldJP1() { return shouldJP1; };
   inline bool ShouldJP2() { return shouldJP2; };
   inline bool ShouldHT2() { return shouldHT2; };
   // hardware (kOnline-sim, bit-7) per-event decisions
   inline bool ShouldHwJP0() { return shouldHwJP0; };
   inline bool ShouldHwJP1() { return shouldHwJP1; };
   inline bool ShouldHwJP2() { return shouldHwJP2; };

   inline float GetEventSumPt() { return event_sum_pt; };

   inline int GetEventMult() { return mult; };

   /// Primary-vertex z of the current event (cm)
   inline double GetVz() { return vz; };
   inline double GetPthat() { return pthat; };
   inline double GetVx() { return vx; };
   inline double GetVy() { return vy; };
   inline double GetVpdVz() { return vz_vpd; };
   inline int GetNVpdEast() { return n_vpd_east; };
   inline int GetNVpdWest() { return n_vpd_west; };
   inline int GetJpAdcMax() { return jp_adc_max; };
   inline int GetHtAdcMax() { return ht_adc_max; };
   inline int GetJpThr(int i) { return jp_thr[i]; };
   inline int GetHtThr(int i) { return ht_thr[i]; };

   /// Get the Trigger (HT) object if it exists, for matching
   inline TStarJetVector *GetTrigger() const { return pHT; };

   /// The main result of the analysis
   inline const vector<ResultStruct> &GetResult() { return Result; }
};

shared_ptr<TStarJetPicoReader> SetupReader(TChain *chain, const ppParameters &pars);

/* For use with GeantMc data
 */
void TurnOffCuts(std::shared_ptr<TStarJetPicoReader> pReader);

bool isInsideJetPatch(const int &jetPatch, const float &jetEta, const float &jetPhi);
bool isInsideJetPatchBox(const int &jetPatch, const float &jetEta, const float &jetPhi);

#endif // __PPANALYSIS_HH
