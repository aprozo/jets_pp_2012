/** @file JetAnalyzer.hh
    @author Kauder:Kolja
    @brief Light wrapper around FastJet 3.x, derived from
    fastjet::ClusterSequenceArea. Event selection and track cleanup belong in
    the calling macro; this class is handed PseudoJets and clusters them. It
    also carries the free helpers used with it: the constituent UserInfo, the
    charge and dijet selectors, and the PseudoJet <-> ROOT conversions.
*/

#ifndef JETANALYZER_H
#define JETANALYZER_H

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/Selector.hh"

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

#include <iostream>
#include <sstream>
#include <string>

// For pairs and sorting in the dijet finding
#include <algorithm>
#include <utility>

class JetAnalyzer : public fastjet::ClusterSequenceArea {

private:
   /** Keep a copy of the original constituents. In principle they are reachable
       via this->jets(), but that vector has extra space allocated in it.
    */
   std::vector<fastjet::PseudoJet> &OrigParticles;

   /** Determined by whether we have an area definition */
   bool CanDoBackground;

   /**  For background subtraction.
        Pointer for now, so that it can be 0 if unused.
        Default specs and areas will be supplied.
   */
   fastjet::JetMedianBackgroundEstimator *bkgd_estimator;
   /** Subtractor */
   fastjet::Subtractor *bkgd_subtractor;

   fastjet::BackgroundJetScalarPtDensity *scalarPtDensity;

   /** Background jet cut */
   fastjet::Selector selector_bkgd;
   /**  Background jet definiton */
   fastjet::JetDefinition *jet_def_bkgd;
   /**  Background area definiton */
   fastjet::AreaDefinition *area_def_bkgd;

public:
   /** Standard constructor: passes through to FastJet and sets up the internal
       background estimator with default values. Presence of an AreaDefinition is
       what enables the background capability, hence the second ctor below.
       \param InOrigParticles is the full set of constituent candidates. Passed
      by reference!
       \param JetDef is a fastjet::JetDefinition for the clustering. Passed by
      reference!
       \param AreaDef is a fastjet::AreaDefinition for the clustering
       \param selector_bkgd is a fastjet::Selector for background subtraction
    */
   JetAnalyzer(std::vector<fastjet::PseudoJet> &InOrigParticles, fastjet::JetDefinition &JetDef,
               fastjet::AreaDefinition &AreaDef,
               fastjet::Selector selector_bkgd = fastjet::SelectorAbsRapMax(0.6) * (!fastjet::SelectorNHardest(2)));

   /** Use as ClusterSequence, without area computation and without BG options
       \param InOrigParticles is the full set of constituent candidates. handed
      by reference!
       \param JetDef is a fastjet::JetDefinition for the clustering
    */
   JetAnalyzer(std::vector<fastjet::PseudoJet> &InOrigParticles, fastjet::JetDefinition &JetDef);

   /** Destructor. Take care of all the objects created with new. */
   ~JetAnalyzer()
   {
      if (area_def_bkgd) {
         delete area_def_bkgd;
         area_def_bkgd = 0;
      }
      if (jet_def_bkgd) {
         delete jet_def_bkgd;
         jet_def_bkgd = 0;
      }
      if (bkgd_estimator) {
         delete bkgd_estimator;
         bkgd_estimator = 0;
      }
      if (bkgd_subtractor) {
         delete bkgd_subtractor;
         bkgd_subtractor = 0;
      }
      // scalarPtDensity is deliberately NOT deleted here: ownership appears to
      // pass to the background estimator, and deleting it crashes.
   };

   /** Background functionality.
       Currently, the jet definition is hard-coded to fastjet::kt_algorithm,
      jet_def().R(), and the area definition is computed internally. Expand and
      modify as needed.
    */
   fastjet::Subtractor *GetBackgroundSubtractor();
   /**
      Handle to BackgroundEstimator()
    */
   fastjet::JetMedianBackgroundEstimator *GetBackgroundEstimator()
   {
      if (!bkgd_estimator)
         throw(std::string("No estimator available"));
      return bkgd_estimator;
   };

   /** Set BackgroundEstimator by hand */
   void SetBackgroundEstimator(fastjet::JetMedianBackgroundEstimator *bge) { bkgd_estimator = bge; };

   static const double pi;

   /** Returns an angle between -pi and pi */
   static double phimod2pi(double phi);
};

// The remainder is NOT part of the class.
// =============================================================================
/** Dijet finding as a Selector, fashioned after fastjet::SW_NHardest.
    Searches for and returns dijet pairs within |&phi;1 - &phi;2 - &pi;| <
   &Delta;&phi;. Returns 0 if no pair is found. Only the top two jets are
   compared.
 */
class SelectorDijetWorker : public fastjet::SelectorWorker {
public:
   /** Standard constructor
       \param dPhi: Opening angle, searching for |&phi;1 - &phi;2 - &pi;| <
      &Delta;&phi;
   */
   SelectorDijetWorker(double dPhi) : dPhi(dPhi) {};

   /// the selector's description
   std::string description() const
   {
      std::ostringstream oss;
      oss << "Searches for and returns dijet pairs within |phi1 - phi2 - pi| < " << dPhi;
      return oss.str();
   };

   /// Returns false, we need a jet ensemble.
   bool applies_jet_by_jet() const { return false; };

   /// Never valid on a single jet; throws.
   bool pass(const fastjet::PseudoJet &pj) const
   {
      if (!applies_jet_by_jet())
         throw(std::string("Cannot apply this selector worker to an individual jet"));
      return false;
   };

   /// The relevant method
   void terminator(std::vector<const fastjet::PseudoJet *> &jets) const;

private:
   const double dPhi; ///< Opening angle, searching for |&phi;1 - &phi;2 - &pi;|
                      ///< < &Delta;&phi;
};

/** Helper for sorting pairs by second argument */
struct sort_IntDoubleByDouble {
   /// returns left.second < right.second
   bool operator()(const std::pair<int, double> &left, const std::pair<int, double> &right)
   {
      return left.second < right.second;
   }
};
/** The actual dijet selector.
    \param dPhi: Dijet acceptance angle &Delta;&phi;
 */
fastjet::Selector SelectorDijets(const double dPhi = 0.4);

// =============================================================================
/** Determines whether two vector sets are matched 1 to 1. Enforcing 1-to-1
    avoids pathologies.
 */
bool IsMatched(const std::vector<fastjet::PseudoJet> &jetset1, const std::vector<fastjet::PseudoJet> &jetset2,
               const double Rmax);

/** Check if one of the jets in jetset1 matches the reference. */
bool IsMatched(const std::vector<fastjet::PseudoJet> &jetset1, const fastjet::PseudoJet &reference, const double Rmax);

/** Check if jet1 and jet2 are matched */
bool IsMatched(const fastjet::PseudoJet &jet1, const fastjet::PseudoJet &jet2, const double Rmax);

// =============================================================================
/** vector<PseudoJet> is interfaced with ROOT via TClonesArray<TLorentzVector>,
    so these convert back and forth.
*/
TLorentzVector MakeTLorentzVector(const fastjet::PseudoJet &pj);
fastjet::PseudoJet MakePseudoJet(const TLorentzVector *const lv);

// =============================================================================
/** Constituent- and jet-level payload attached to a PseudoJet, derived from
    PseudoJet::UserInfoBase.
 */
class JetAnalysisUserInfo : public fastjet::PseudoJet::UserInfoBase {
public:
   /// Standard Constructor
   JetAnalysisUserInfo(int quarkcharge = -999, int pid = -9999, std::string tag = "", float number = -1)
      : quarkcharge(quarkcharge), pid(pid), tag(tag), number(number) {};

   /// Charge in units of e
   int GetQuarkCharge() const { return quarkcharge; };

   int GetPID() const { return pid; };

   /// Multi-purpose description
   std::string GetTag() const { return tag; };
   void SetTag(const std::string newtag) { tag = newtag; };

   /// Multi-purpose description
   float GetNumber() const { return number; };
   void SetNumber(const float f) { number = f; };

   bool IsMatchedJP() const { return match_jp; };
   void SetMatchJP(const bool b) { match_jp = b; };

   // Per-threshold JP-patch matches, set independently of each other, for
   // per-trigger response/data combination. IsMatchedJP() is the match for the
   // trigger this job was configured with.
   bool IsMatchedJP0() const { return match_jp0; };
   void SetMatchJP0(const bool b) { match_jp0 = b; };
   bool IsMatchedJP1() const { return match_jp1; };
   void SetMatchJP1(const bool b) { match_jp1 = b; };
   bool IsMatchedJP2() const { return match_jp2; };
   void SetMatchJP2(const bool b) { match_jp2 = b; };

   bool IsMatchedHT() const { return match_ht; };
   void SetMatchHT(const bool b) { match_ht = b; };

   int GetTrackId() const { return trackid; };
   void SetTrackId(const int id) { trackid = id; };

   int GetLeadTowerId() const { return lead_tower_id; };
   void SetLeadTowerId(const int id) { lead_tower_id = id; };

   int GetJpAdc() const { return jp_adc; };
   void SetJpAdc(const int a) { jp_adc = a; };
   int GetHtAdc() const { return ht_adc; };
   void SetHtAdc(const int a) { ht_adc = a; };

   float GetsDCAxy() const { return sdcaxy; };
   void SetsDCAxy(const float s) { sdcaxy = s; };

private:
   const int quarkcharge; ///< Charge in units of e
   const int pid;
   std::string tag; ///< Multi-purpose
   float number;    ///< Multi-purpose
   int trackid;
   bool match_jp;
   bool match_jp0;
   bool match_jp1;
   bool match_jp2;
   bool match_ht;
   int lead_tower_id = -1; ///< BEMC tower id of leading neutral constituent (-1 if none)
   int jp_adc = -1;        ///< near-max JP-patch ADC (>JP0) matched to jet (-1 if none)
   int ht_adc = -1;        ///< max DSM high-tower ADC of the jet's towers (-1 if none above BHT1)
   float sdcaxy = -999;    ///< track signed transverse DCA (charged constituents; for leading-track fake veto)
};

// =============================================================================
/** Selects particles by the charge stored in their JetAnalysisUserInfo.
    Charge is in units of e/3!
*/
class SelectorChargeWorker : public fastjet::SelectorWorker {
public:
   /** Standard constructor
       \param cmin: inclusive lower bound
       \param cmax: inclusive upper bound
   */
   SelectorChargeWorker(const int cmin, const int cmax) : cmin(cmin), cmax(cmax) {};

   /// the selector's description
   std::string description() const
   {
      std::ostringstream oss;
      oss.str("");
      oss << cmin << " <= quark charge <= " << cmax;
      return oss.str();
   };

   /// keeps the ones that have cmin <= quarkcharge <= cmax
   bool pass(const fastjet::PseudoJet &p) const
   {
      // FastJet's `&&` between selectors does NOT short-circuit, so ghosts
      // injected by ClusterSequenceArea reach this predicate even when composed
      // as `NotGhost && SelectorChargeRange(...)` — hence these guards.
      if (p.is_pure_ghost()) return false;
      if (!p.has_user_info<JetAnalysisUserInfo>()) return false;
      const int &quarkcharge = p.user_info<JetAnalysisUserInfo>().GetQuarkCharge();
      return (quarkcharge >= cmin) && (quarkcharge <= cmax);
   };

private:
   const int cmin; ///< inclusive lower bound
   const int cmax; ///< inclusive upper bound
};
/** Builds the selector: Selector sel = SelectorChargeRange( cmin, cmax ); */
fastjet::Selector SelectorChargeRange(const int cmin = -999, const int cmax = 999);

// =============================================================================
/** Helper to get a jet-algorithm enum from a string; some generous spellings
    and abbreviations are accepted.
*/
fastjet::JetAlgorithm AlgoFromString(std::string s);

#endif // JETANALYZER_H