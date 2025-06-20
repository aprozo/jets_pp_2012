/* @file ppAnalysis.hh
    @author Raghav Kunnawalkam Elayavalli
    @version Revision 1.0
    @brief analysis class
    @details Uses JetAnalyzer objects
    @date March 16, 2022
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

// Not needed for analysis per se
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoReader.h"

#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTriggerInfo.h"

#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"
#include "TStarJetVector.h"
#include "TStarJetVectorContainer.h"

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

double LookupRun12Xsec(TString filename);
/*
   For sorting with a different key
*/
typedef pair<PseudoJet, double> PseudoJetPt;
struct PseudoJetPtGreater {
  bool operator()(PseudoJetPt const &a, PseudoJetPt const &b) {
    return a.second > b.second;
  }
};

/*
  To keep original and groomed jets connected
 */
class ResultStruct {
public:
  PseudoJet orig;
  ResultStruct(PseudoJet orig) : orig(orig){};
  static bool origptgreater(ResultStruct const &a, ResultStruct const &b) {
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
void InitializeReader(std::shared_ptr<TStarJetPicoReader> pReader,
                      const TString InputName, const Long64_t NEvents,
                      const int PicoDebugLevel,
                      const double HadronicCorr = 0.999999);

static const Selector NotGhost =
    !fastjet::SelectorIsPureGhost(); ///< Helper useful outside the class as
                                     ///< well
static const Selector OnlyCharged =
    NotGhost &&
    (SelectorChargeRange(-3, -1) ||
     SelectorChargeRange(1, 3)); ///< Helper useful outside the class as well
static const Selector OnlyNeutral =
    NotGhost &&
    SelectorChargeRange(0, 0); ///< Helper useful outside the class as well

// Histograms for  QA histograms

struct HistogramManager {
  void Init() {
    vx = new TH1D("vx", "Primary vertex x; vx, cm; N", 300, -0.5, 0.5);
    vy = new TH1D("vy", "Primary vertex y; vy, cm; N", 300, -0.5, 0.5);
    vz = new TH1D("vz", "Primary vertex z; vz, cm; N", 280, -70, 70);
    vz_vpd = new TH1D("vz_vpd", "Primary vertex z from VPD; vz, cm; N", 280,
                      -70, 70);
    vz_diff = new TH1D("vz_diff", "vertex z VPD-TPC; vz, cm; N", 500, -20, 20);
    pt_sDCAxy_pos = new TH2D("pt_sDCAxy_pos",
                             "sDCAxy positive tracks; sDCA, cm; pt, GeV/c; N",
                             100, -4, 4, 100, 0, 30);
    pt_sDCAxy_neg = new TH2D("pt_sDCAxy_neg",
                             "sDCAxy negative tracks; sDCA, cm;  pt, GeV/c; N",
                             100, -4, 4, 100, 0, 30);
    jetpt_TowerID =
        new TH2D("jetpt_TowerID", "Jet pt vs Tower ID; Jet pt, GeV/c; Tower ID",
                 50, 0, 50, 4801, 0, 4801);
    highjetpt_leadingtower_pt = new TH1D(
        "highjetpt_leadingtower_pt",
        "Leading tower pt in high(>22) pt jets; Leading tower pt, GeV/c; N",
        100, 0, 40);
    highjetpt_leadingtrack_pt = new TH1D(
        "highjetpt_leadingtrack_pt",
        "Leading track pt in high(>22) pt jets; Leading track pt, GeV/c; N",
        100, 0, 40);
    //  problematic jets qa histograms
    sDCAxy = new TH1D("sDCAxy", "sDCAxy; sDCA, cm", 100, -1, 1);
    DCA = new TH1D("DCA", "DCA; DCA, cm", 100, 0, 2);
    Chi2 = new TH1D("Chi2", "Chi2; Chi2, a.u.", 200, 0, 10);
    Chi2PV = new TH1D("Chi2PV", "Chi2PV; Chi2PV, a.u.", 200, 0, 100);
    matched = new TH1D("matched", "Matching with TOF and BEMC", 3, 0, 3);
    matched->GetXaxis()->SetBinLabel(1, "No match");
    matched->GetXaxis()->SetBinLabel(2, "BEMC match");
    matched->GetXaxis()->SetBinLabel(3, "TOF match");
  }
  void FillTrack(TStarJetPicoPrimaryTrack *track) {
    sDCAxy->Fill(track->GetsDCAxy() * track->GetCharge());
    DCA->Fill(track->GetDCA());
    Chi2->Fill(track->GetChi2());
    Chi2PV->Fill(track->GetChi2PV());
    bool matchBEMC = track->GetBemcMatchFlag();
    bool matchTOF = track->GetTofMatchFlag();
    if (matchBEMC)
      matched->Fill(1);
    if (matchTOF)
      matched->Fill(2);
    if (!matchBEMC && !matchTOF)
      matched->Fill(0);
  }

  void Write(TFile *f, TString name = "QA_histograms") {
    f->cd();
    f->mkdir(name);
    f->cd(name);
    vx->Write();
    vy->Write();
    vz->Write();
    vz_vpd->Write();
    vz_diff->Write();
    pt_sDCAxy_pos->Write();
    pt_sDCAxy_neg->Write();
    jetpt_TowerID->Write();
    highjetpt_leadingtower_pt->Write();
    highjetpt_leadingtrack_pt->Write();
    // problematic jets
    sDCAxy->Write();
    DCA->Write();
    Chi2->Write();
    Chi2PV->Write();
    matched->Write();
  }

  TH2D *pt_sDCAxy_pos;             ///< QA for positive DCAxy
  TH2D *pt_sDCAxy_neg;             ///< QA for negative DCAxy
  TH2D *jetpt_TowerID;             ///< QA for jet pt vs Tower ID
  TH1D *highjetpt_leadingtower_pt; ///< QA for leading tower pt in high
                                   ///< pt jets
  TH1D *highjetpt_leadingtrack_pt; ///< QA for leading track pt in high
                                   ///< pt jets
  TH1D *vx;                        ///< QA for primary vertex x
  TH1D *vy;                        ///< QA for primary vertex y
  TH1D *vz;                        ///< QA for primary vertex z
  TH1D *vz_vpd;                    ///< QA for primary vertex z from VPD
  TH1D *vz_diff;                   ///< QA for vertex z VPD-TPC
  // QA for tracks
  TH1D *sDCAxy;
  TH1D *DCA;
  TH1D *Chi2;
  TH1D *Chi2PV;
  TH1D *matched;
};

/*
   The main class
 */
class ppAnalysis {

private:
  // These need to be initialized
  // ----------------------------
  ppParameters pars; ///< container to have all analysis parameters in one place

  HistogramManager QA_hist, QA_hist_problematic; ///< QA histograms

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
  // get HistogramManager
  inline HistogramManager &GetHistogramManager() { return QA_hist; };
  inline HistogramManager &GetHistogramManagerProblematic() {
    return QA_hist_problematic;
  };

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

  /// Get the runid of the current event (this id for geant events can match
  /// to bad run ids)
  inline double GetRunid1() { return runid1; };

  /// Get the runid of the current event
  inline double GetRunid() { return runid; };

  /// Get the eventid of the current event
  inline double GetEventid() { return eventid; };

  inline float GetEventSumPt() { return event_sum_pt; };

  inline int GetEventMult() { return mult; };

  inline int GetRejectCode() { return is_rejected; };

  /// Get the Trigger (HT) object if it exists, for matching
  inline TStarJetVector *GetTrigger() const { return pHT; };

  /// The main result of the analysis
  inline const vector<ResultStruct> &GetResult() { return Result; }
};

shared_ptr<TStarJetPicoReader> SetupReader(TChain *chain,
                                           const ppParameters &pars);

/* For use with GeantMc data
 */
void TurnOffCuts(std::shared_ptr<TStarJetPicoReader> pReader);

#endif // __PPANALYSIS_HH
