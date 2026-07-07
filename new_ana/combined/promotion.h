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

// Which trigger's measured data-side correction C(pt) applies to a category
// (config.h::TrigEffMeas; JP0 has no measured correction -> 1).
inline const char *CatTrigger(int cat) { return cat == 2 ? "JP2" : (cat == 1 ? "JP1" : "JP0"); }

// ---- fine grid for the response fill (Dmitry's 550x600) --------------------
// The migration is filled on 0.1 GeV cells (detector 550 bins over [5,60],
// particle 600 bins over [0,60]), statistically unreliable cells are removed,
// and only then everything is rebinned to McBins. Every McBins edge <= 52 is
// a multiple of 0.1 so the rebin is exact; the 52-86 feed-down buffer only
// collects fine content up to 60 (it is excluded from the inverted block, its
// feed-down is background-scaled away by the b/matched row factor).
const int    kNRecoFine = 550;  const double kRecoFineLo = 5.0,  kRecoFineHi = 60.0;
const int    kNMcFine   = 600;  const double kMcFineLo   = 0.0,  kMcFineHi   = 60.0;
const double kFineW     = 0.1;
const int    kBoxRes    = 10; // +-10 fine cells = +-1.0 GeV box
const int    kBoxThresh = 4;  // zero cells whose box entry-sum is <= 4

// ---- floor-restricted inversion block --------------------------------------
// The square block inverted runs from the first McBins bin at or above
// kJpxFloor up to (excluding) kQuoteHi. 9.7 GeV is the validated default; the
// combination is populated down to 6.9 via cat0 but those prescale-deep bins
// can make the unregularized inverse ring — study by editing here.
const double kJpxFloor = 9.7;
const double kQuoteHi  = 52.0; // 52-86 buffer never inverted, never quoted

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
