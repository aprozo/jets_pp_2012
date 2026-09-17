/* @file ppParameters.hh
   @author Kolja Kauder
   @brief Every Stage-1 cut and analysis parameter in one place, so the same
          values can be included from different macros. Defaults here are the
          production values; some are overridden by ppAnalysis command-line args.
 */

#ifndef PPPARAMETERS_HH
#define PPPARAMETERS_HH

#include "JetAnalyzer.hh"

using fastjet::antikt_algorithm;
using fastjet::cambridge_algorithm;
using fastjet::JetAlgorithm;
using fastjet::kt_algorithm;

enum INTYPE {
   MCTREE,
   INTREE,
   INPICO,
   MCPICO,
   HERWIGTREE
};

// Return values for the main routine
enum class EVENTRESULT {
   PROBLEM,
   ENDOFINPUT,
   NOTACCEPTED,
   NOCONSTS,
   NOJETS,
   JETSFOUND
};

class ppParameters {

public:
   double R = 0.6; ///< Resolution parameter ("radius").

   /// Jet algorithm for the original jets
   JetAlgorithm LargeJetAlgorithm = fastjet::antikt_algorithm;

   /// Alternative reclustering to Cambridge/Aachen
   bool CustomRecluster = false;
   JetAlgorithm CustomReclusterJetAlgorithm;
   bool Recursive = false; ///< Repeat on subjets?

   double PtJetMin = 5.0;    ///< Min jet pT
   double PtJetMax = 1000.0; ///< Max jet pT

   /// Max neutral energy fraction R_T = sum_pT_tower / sum_pT_jet. The value is
   /// stored per jet at Stage-1; the < 0.95 cut is applied at Stage-2
   /// (published analysis, Table 4).
   double MaxJetNEF = 0.95;

   double EtaConsCut = 1.0; ///< Constituent |&eta;| acceptance
   double PtConsMin = 0.2;  ///< Constituent pT minimum (tracks only)
   double PtConsMax = 30;   ///< Constituent pT maximum (tracks only)

   double RefMultCut = 0; ///< Reference multiplicity. Needs to be rethought to
                          ///< accomodate pp and AuAu

   double VzCut = 60;        ///< Vertex z (standard Run-12 window)
   double VzDiffCut = 99999; ///< |Vz(TPC) - Vz(VPD)|

   double DcaCut = 3.0;               ///< flat 3-D DCA cap |dcaGlobal().mag()| < 3 cm, WITH z; applied by reader SetDCACut
   double sDCAxyCut = 99999;          ///< |signed dca_xy| flat cap (off; TdcaPtDep below replaces it)
   double LeadTrackSdcaCut = 99999;   ///< veto a jet whose LEADING charged constituent has |sDCAxy|
                                      ///< above this, to remove single-fake-track jets. Off by default.
   double NMinFit = 12;               ///< minimum number of fit points for tracks
   double FitOverMaxPointsCut = 0.51; ///< NFit / NFitPossible

   /// Minimum track flag: requires flag > 0 (the pico maker keeps flag >= 0).
   int FlagMin = 1;

   /// pT-dependent cut on the FULL 3-D DCA (dcaGlobal().mag() = track->GetDCA(),
   /// WITH z — not the transverse dca_xy), from the published Run-12 pp200 note
   /// (Table 2):
   ///   DCA < 2 cm                 for pt < 0.5 GeV
   ///       < 2.5 cm - pt*(1/GeV)  for 0.5 <= pt < 1.5 GeV   (slope -1)
   ///       < 1 cm                 for pt >= 1.5 GeV
   bool   ApplyTdcaPtDep = true;
   double TdcaPt1     = 0.5;
   double TdcaPt2     = 1.5;
   double TdcaDcaMax1 = 2.0;
   double TdcaDcaMax2 = 1.0;

   double HadronicCorr = 1.0; ///< Fraction of hadronic correction

   double FakeEff = 1.0; ///< fake efficiency for systematics. 0.95 is a reasonable example.

   Int_t IntTowScale = 0;
   /// Tower GAIN: 4.8%
   Float_t fTowUnc = 0.048;
   /// Tower scale for uncertainty;
   float fTowScale = 1.0;

   // High-pT tracks and towers are NOT cut, at the constituent level or at the
   // event level (1000 = no veto). Rejection is per jet instead, via
   // MaxJetTrackPt below; an event-level veto would carve a cliff into the
   // trigger efficiency near reco jet pT ~ 30 GeV.
   double MaxEtCut = 1000;   ///< tower ET cut
   double MaxTrackPt = 1000; ///< track pT cut
   double MaxEventPtCut = 1000; ///< event-level track veto disabled
   double MaxEventEtCut = 1000; ///< the published analysis has no event-level tower cut

   /// Per-jet rejection: drop the jet (but keep the event) if any charged
   /// constituent has pT above this. Applied in the jet loop in ppAnalysis.cxx.
   double MaxJetTrackPt = 30;
   double MinEventEtCut = 0;  ///< min event ET cut for event
   double ManualHtCut = 0.0;  ///< necessary for some embedding picos. Should
                              ///< always equal MinEventEtCut

   /// Geant files have unusable runid / event id; turn this on to synthesize
   /// reproducible ones. On for Geant files, off otherwise.
   bool UseGeantNumbering = false;

   TString InputName = "test.root";
   INTYPE intype = INPICO;        ///< Input type (can be a pico dst, a result tree, an MC tree)
   TString ChainName = "JetTree"; ///< Name of the input chain
   TString TriggerName = "JP2";
   TString OutFileName = "test.root";
};
#endif // PPPARAMETERS_HH
