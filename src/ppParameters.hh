/* @file parameters file
   @author Kolja Kauder
   @brief Common parameters
   @details Used to quickly include the same parameters into different macros.
   @date Mar 23, 2017
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
   // JetAlgorithm LargeJetAlgorithm = fastjet::cambridge_algorithm;

   /// Alternative reclustering to Cambridge/Aachen (at our own risk)
   bool CustomRecluster = false;
   JetAlgorithm CustomReclusterJetAlgorithm;
   // bool CustomRecluster=true;
   // JetAlgorithm CustomReclusterJetAlgorithm = fastjet::antikt_algorithm;
   bool Recursive = false; ///< Repeat on subjets?

   /// Repetitions in the background. Anything other than 1 WILL NOT WORK because
   /// a) we're using explicit ghosts (though we don't have to)
   /// b) more importantly, the background subtractor contains
   /// fastjet::SelectorNHardest(2)
   ///    which doesn't work jet-by-jet and throws an error
   // int GhostRepeat = 1;
   // float GhostArea = 0.005;    ///< ghost area

   double PtJetMin = 5.0;    ///< Min jet pT
   double PtJetMax = 1000.0; ///< Max jet pT

   double MaxJetNEF = 0.95; ///< Jet R_T = sum_pT_tower / sum_pT_jet < 0.95 (analysis note Table 4)

   double EtaConsCut = 1.0; ///< Constituent |&eta;| acceptance || was 1.0
   double PtConsMin = 0.2;  ///< Constituent pT minimum || was 0.2
   double PtConsMax = 30;   ///< Constituent pT maximum || was 30

   double RefMultCut = 0; ///< Reference multiplicity. Needs to be rethought to
                          ///< accomodate pp and AuAu

   double VzCut = 60; ///< Vertex z (matches Dmitry's default min/max_vertex_z = -60..+60)
   // const double VzDiffCut=6;         ///< |Vz(TPC) - Vz(VPD)| <-- NOT WORKING
   // in older data (no VPD)
   double VzDiffCut = 99999; ///< |Vz(TPC) - Vz(VPD)|

   double DcaCut = 3.0;               ///< flat 3-D DCA cap |dcaGlobal().mag()| < 3 cm, WITH z (Dmitry: 3.0); applied by reader SetDCACut
   double sDCAxyCut = 99999;          ///< |signed dca_xy| flat cap (kept off — TdcaPtDep below replaces it)
   double NMinFit = 12;               ///< minimum number of fit points for tracks (Dmitry: 12)
   double FitOverMaxPointsCut = 0.51; ///< NFit / NFitPossible (Dmitry: 0.51)

   /// Minimum track flag. Dmitry's StjTrackCutFlag(0) rejects flag<=0
   /// (i.e. requires flag > 0). Pico maker default is flag>=0.
   int FlagMin = 1;

   /// pT-dependent cut on the FULL 3-D DCA (dcaGlobal().mag() = track->GetDCA()),
   /// matching Dmitry's Run12 StjTrackCutTdcaPtDependent and the PUBLISHED pp200
   /// note (Table 2):
   ///   DCA < 2 cm                 for pt < 0.5 GeV
   ///       < 2.5 cm - pt*(1/GeV)  for 0.5 <= pt < 1.5 GeV   (slope -1)
   ///       < 1 cm                 for pt >= 1.5 GeV
   /// i.e. pt1=0.5, dca1=2.0, pt2=1.5, dca2=1.0. (StjTrackCutTdcaPtDependent's
   /// "Tdca" is the 3-D magnitude, NOT the transverse dca_xy; StjTPCMuDst.cxx:100.
   /// The dca_xy / transverse twin StjTrackCutDcaPtDependent ships pt2=1.0 and is
   /// what Dmitry used for run9 pp200, NOT run12 — this repo follows run12.)
   bool   ApplyTdcaPtDep = true;
   double TdcaPt1     = 0.5;
   double TdcaPt2     = 1.5;
   double TdcaDcaMax1 = 2.0;
   double TdcaDcaMax2 = 1.0;

   double HadronicCorr = 1.0; ///< Fraction of hadronic correction (Dmitry: 1.00)

   double FakeEff = 1.0; ///< fake efficiency for systematics. 0.95 is a reasonable example.

   Int_t IntTowScale = 0;
   /// Tower GAIN: 4.8%
   Float_t fTowUnc = 0.048;
   /// Tower scale for uncertainty;
   float fTowScale = 1.0;

   // ************************************
   // Do NOT cut high tracks and towers!
   // Instead, reject only the affected jet (per-jet, à la Dmitry's
   // make_max_track_pt_cut). Event-level rejection produced a sharp
   // cliff at reco jet pT ≈ 30 GeV in the JP2/all efficiency.
   // ************************************
   double MaxEtCut = 1000;   ///< tower ET cut
   double MaxTrackPt = 1000; ///< track pT cut

   // EVENT rejection cuts — match Dmitry's per-event vetoes
   // (StJetPlots/StJetCut.h:make_max_track_pt_cut, default.nix:24).
   // A charged track with pt > MaxEventPtCut rejects the entire event.
   double MaxEventPtCut = 1000;   ///< Dmitry: 30 GeV
   double MaxEventEtCut = 1000; ///< Dmitry has no event-level tower cut

   // Per-jet rejection: drop the jet (keep the event) if any charged
   // constituent has pT > this. Disabled (set to a huge value) because
   // event-level cut above already removes such events.
   double MaxJetTrackPt = 30;
   double MinEventEtCut = 0;  ///< min event ET cut for event
   double ManualHtCut = 0.0;  ///< necessary for some embedding picos. Should
                              ///< always equal MinEventEtCut

   // Geant files have messed up runid and event id.
   // Switch for fixing. Should be turned on by default for Geant files, off
   // otherwise.
   bool UseGeantNumbering = false;

   TString InputName = "test.root";
   INTYPE intype = INPICO;        ///< Input type (can be a pico dst, a result tree, an MC tree)
   TString ChainName = "JetTree"; ///< Name of the input chain
   TString TriggerName = "JP2";
   TString OutFileName = "test.root";
};
#endif // PPPARAMETERS_HH
