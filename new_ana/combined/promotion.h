#ifndef PROMOTION_CONFIG_H
#define PROMOTION_CONFIG_H
// Promotion combination of the JP triggers (Dmitry's published nominal):
// JP0+JP1+JP2 partitioned into EXCLUSIVE categories by the highest jet-patch
// threshold the event's patches clear, summed with prescale-recovery weights,
// and unfolded ONCE by unregularized matrix inversion on a cell-filtered
// fine-binned response. Everything promotion-specific lives here; everything
// shared with the per-trigger pipeline (bins, C(pt), badRuns, paths) comes
// from ../config.h.
//
// The category of an event is the highest patch threshold it SHOULD fire:
//   cat2 : should_JP2                    jets: trigger_match_JP2 && pt >= 8.4
//   cat1 : should_JP1 && !should_JP2     jets: trigger_match_JP1 && pt >  8.2
//   cat0 : should_JP0 && !should_JP1     jets: trigger_match_JP0 && pt < 22.5
// In data the should bits are the OR of the per-jet hardware trigger_match
// bits; in embedding they are the stored evt_should_JP* (same jet-based
// convention). Data events additionally require the recorded (prescaled)
// hardware accept: cat2 fired_JP0|1|2, cat1 fired_JP0|1, cat0 fired_JP0.
//
// Data enters RAW (no prescale weight) and is normalized by the FULL JP2
// luminosity; the response's measured side carries the sampling probability
// instead:
//   w2 = 1,  w1 = 1/ps0 + 1/ps1 - <1/(ps0*ps1)>  (P(fired JP0 or JP1)),
//   w0 = 1/ps0
// The embedding is run-blind (its runid is a production timestamp, not a data
// anchor run), so the weights are the lumi-weighted expectations over the run
// mix: 1/psN = Leff(JPN)/Leff(JP2) from lumi_zilong_full.root (badRuns-
// consistent, the same source as the normalization) and the union cross term
// lumi-averaged from lists/run_prescales.txt.

#include "../config.h"

#include <TFile.h>
#include <TH1D.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace PromotionConfig {

using CrossSectionConfig::AnalysisConfig;

// ---- category detector-pT windows (GeV, on the UE-subtracted reco pT) ------
// SINGLE definition used by BOTH the data selection (cross_section.cpp) and
// the response reco side (response.cxx) — that identity is the invariant that
// makes the shouldFire turn-on cancel in the unfold.
inline std::string CatWindow(int cat, const std::string &ptCol)
{
   if (cat == 2) return ptCol + " >= 8.4";
   if (cat == 1) return ptCol + " > 8.2";
   return ptCol + " < 22.5";
}

// Per-category jet gate: the jp_match veto + the window, on either side's
// branch naming (prefix "" for data, "reco_" for the matching tree).
inline std::string CatJetGate(int cat, const std::string &prefix, const std::string &ptCol)
{
   return prefix + "trigger_match_JP" + std::to_string(cat) + " && " + CatWindow(cat, ptCol);
}

// =====================================================================
// MEASURED combination-level trigger correction C_JPX(pt), applied to the
// summed data before unfolding (corrections/measure_Cjpx.C, 2026-07-07).
//
// C_JPX = C_sum(pt) = eps_data(pt)/eps_emb(pt), the combination's
// sampling-weighted gate-probability ratio
//   eps = sum_cat s_cat P(cat gate | pt),  s = {w0, w1, 1},
// measured for the EXACT promotion gates (exclusive shouldFire category +
// jp_match veto + windows) on the unbiased fired_JP0 data base vs the
// total_weight-weighted embedding. The response Misses divide by eps_emb;
// dividing the data by C_sum closes the pair by construction. It is ONE
// self-contained in-situ measurement — no inclusive-trigger R/T-hat tables
// (they do not see the exclusive-category shuffle and over-correct the
// turn-on by ~half their size).
//
// Bin-by-bin in the high-statistics turn-on ([8.2,19), the fired_JP0 base
// is deep there); the plateau is frozen to the [19,44) pol0 fit 0.9830 —
// the hardware/simulator ratio is smooth, per-bin base noise (+-2-4%
// above ~30 GeV) must not inject fake structure. NB the plateau level is
// R x T_true with T_true = 0.995 for the combination (per-mil, unlike the
// standalone JP2's ~0.90): dividing it out costs 0.5% normalization and is
// within the trigger-efficiency systematic.
// =====================================================================
inline double JpxTrigEff(double pt)
{
   if (pt <  8.2) return 1.0;    // below every solved bin (never quoted)
   if (pt <  9.7) return 0.8332; // measured C_sum turn-on
   if (pt < 11.5) return 0.8729;
   if (pt < 13.6) return 0.9105;
   if (pt < 16.1) return 0.9490;
   if (pt < 19.0) return 0.9670;
   return 0.9830;                // C_sum plateau fit over [19,44)
}

// ---- fine grid for the response fill ---------------------------------------
// The migration is filled on 0.1 GeV cells (Dmitry's grid extended to cover
// the full analysis range), statistically unreliable cells are removed, and
// only then everything is rebinned to McBins. Every McBins edge is a multiple
// of 0.1 so the rebin is exact. The axes MUST cover the 52-86 feed-down
// buffer: Dmitry's own [0,60] caps were fine for his 52-top grid, but with
// the buffer bin a 60-capped truth axis leaves the buffer column truncated
// (its reco row is real data up to 80) and the last quoted bin (44-52) swings
// by +-8% depending on how the incomplete column is treated.
const int    kNRecoFine = 810;  const double kRecoFineLo = 5.0,  kRecoFineHi = 86.0;
const int    kNMcFine   = 860;  const double kMcFineLo   = 0.0,  kMcFineHi   = 86.0;
const double kFineW     = 0.1;
const int    kBoxRes    = 10; // +-10 fine cells = +-1.0 GeV box
const int    kBoxThresh = 4;  // zero cells whose box entry-sum is <= 4

// ---- floor-restricted inversion block --------------------------------------
// The square block inverted runs from the first McBins bin at or above
// kJpxFloor up to (excluding) kQuoteHi. 9.7 GeV is the validated default; the
// combination is populated down to 6.9 via cat0 but those prescale-deep bins
// can make the unregularized inverse ring — study by editing here.
const double kJpxFloor = 9.7;
const double kQuoteHi  = 86.0; // include the 52-86 feed-down buffer as a real
                               // solved column/row (its reco row is real data,
                               // 52-80); solved, never quoted. With the buffer
                               // excluded (Dmitry's 52 cut) its feed-down is
                               // background-scaled instead and the 44-52 bin
                               // comes out +7% — the standalone square solver
                               // (cross_section_inverse) keeps the buffer too.

// Optional Tikhonov damping of the block solve (0 = plain M^-1 = Dmitry).
// Scale-matched second-difference penalty on the unfolded/embedding-truth
// ratio; 0.014 is the validated light setting if the tail needs damping.
const double kTikhonovLambda = 0.0;

// ---- runtime luminosity (same source + badRuns logic as the per-trigger
// pipeline; the promotion normalizes by the FULL JP2 luminosity) -------------
inline double RuntimeLeff(const AnalysisConfig &cfg, const std::string &trigger)
{
   TFile lf((cfg.workdir + "lumi_zilong_full.root").c_str(), "READ");
   auto *lumi = lf.IsZombie() ? nullptr : (TH1D *)lf.Get(Form("luminosity_%s", trigger.c_str()));
   if (!lumi)
      throw std::runtime_error("luminosity_" + trigger + " missing in lumi_zilong_full.root");
   double sum = 0.0;
   for (int i = 1; i <= lumi->GetNbinsX(); ++i) {
      const char *label = lumi->GetXaxis()->GetBinLabel(i);
      if (!label || !*label) continue;
      if (std::find(cfg.badRuns.begin(), cfg.badRuns.end(), std::stoi(label)) != cfg.badRuns.end())
         continue;
      sum += lumi->GetBinContent(i);
   }
   return sum;
}

// ---- promotion prescale-recovery weights ------------------------------------
// The embedding is run-blind, so the response carries the lumi-weighted
// expected sampling probability of each category over the kept run mix:
//   1/psN        = Leff(JPN)/Leff(JP2)            (exact, lumi-weighted mean)
//   <1/(ps0ps1)> = sum_r L_r/(ps0_r ps1_r)/sum L  (union cross term, from the
//                                                  per-run prescale table)
struct PromotionWeights {
   double w0, w1, w2;
   double LeffFull;
};

inline PromotionWeights Weights(const AnalysisConfig &cfg)
{
   const double L2 = RuntimeLeff(cfg, "JP2"); // JP2 unprescaled == full lumi
   const double L1 = RuntimeLeff(cfg, "JP1");
   const double L0 = RuntimeLeff(cfg, "JP0");
   const double invPs0 = L0 / L2, invPs1 = L1 / L2;

   // Lumi-weighted <1/(ps0*ps1)> over kept runs (per-run ps join per-run L).
   std::map<int, std::pair<double, double>> runPs;
   {
      std::ifstream fp((cfg.workdir + "../lists/run_prescales.txt").c_str());
      std::string line;
      int rid;
      double p0, p1;
      while (std::getline(fp, line)) {
         if (line.empty() || line[0] == '#') continue;
         std::istringstream ss(line);
         if (ss >> rid >> p0 >> p1) runPs[rid] = {p0, p1};
      }
   }
   double cross = invPs0 * invPs1; // fallback: product of the means
   {
      TFile lf((cfg.workdir + "lumi_zilong_full.root").c_str(), "READ");
      auto *lumi = lf.IsZombie() ? nullptr : (TH1D *)lf.Get("luminosity_JP2");
      if (lumi && !runPs.empty()) {
         double num = 0.0, den = 0.0;
         for (int i = 1; i <= lumi->GetNbinsX(); ++i) {
            const char *label = lumi->GetXaxis()->GetBinLabel(i);
            if (!label || !*label) continue;
            const int run = std::stoi(label);
            if (std::find(cfg.badRuns.begin(), cfg.badRuns.end(), run) != cfg.badRuns.end())
               continue;
            auto it = runPs.find(run);
            if (it == runPs.end()) continue;
            const double L = lumi->GetBinContent(i);
            num += L / (it->second.first * it->second.second);
            den += L;
         }
         if (den > 0) cross = num / den;
      }
   }

   PromotionWeights w;
   w.w0 = invPs0;
   w.w1 = invPs0 + invPs1 - cross;
   w.w2 = 1.0;
   w.LeffFull = L2;
   return w;
}

} // namespace PromotionConfig

#endif // PROMOTION_CONFIG_H
