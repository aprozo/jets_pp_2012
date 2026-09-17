// unfold.C — inclusive jet cross section from the eta-block response: matrix inversion (the published
// construction) or Bayesian unfolding of the same ingredients.
//
// Inputs (<results>/): data_blocks_R<R>[_mb][_<variant>].root (build_data.C),
// response_blocks_R<R>[_e23|_mb2023|_mb2021m][_<variant>].root (build_resp.C), the per-run luminosities
// and the chain efficiency (../inputs/) and the bad-run lists.
// Construction (shared by both methods):
//   levels l (exclusive trigger partition 0,1,2; min-bias 3; inclusive single triggers 4 = JP1, 5 = JP2;
//   high tower 6 = HT2), eta blocks i (detector) and j (particle) in {00_05, 05_09};
//   per-level detector-pT windows (common.h): l0 below 22.5 GeV, l1 above 8.2, l2 above 9.7, HT2 above
//   11.5 GeV, the inclusive and min-bias levels unrestricted;
//   data      b_i = sum_l  n_{l,i} / L_l   (L_l = luminosity sampled by level l over the kept runs, pb^-1;
//             min-bias: L_MB x eps_chain, the JP0-sample chain efficiency constant)
//   response  A_ik = sum_l A_{l,ik},  b^emb_i = sum_l b^emb_{l,i},  m_i = sum_k A_ik (pairs in both ranges),
//             x_k = particle spectrum (all generated events);
//   outlier filter (published construction), per pt-hat sample on the level-summed fine arrays: a response
//             cell is dropped when the 20 x 20-cell box around it (0.1 GeV cells) holds <= 4 entries; the
//             dropped pairs leave b^emb and x as well, and so do the unmatched detector (particle) entries
//             whose 20-cell neighbourhood holds <= 4 unmatched entries;
//   fine axes rebinned to the analysis bins.
// method "inv":   M_ik = (b^emb_i / m_i) A_ik / x_k on the 12 published bins (24 x 24 with the blocks);
//                 x_unf = M^-1 b, covariance U B U^T with B_ii = sum_l n_{l,i} / L_l^2 (unit-weight data).
// method "bayes": the same purity (m_i / b^emb_i) applied to the data, then RooUnfoldBayes (nIter
//                 iterations, prior = x) on the fake-free response (A, x), with one detector bin
//                 [5.0, 6.9) and two particle bins [3.5, 5.2), [5.2, 6.9) below the first published bin,
//                 so that the migrations across 6.9 GeV sit inside the matrix.
// mode: <trigger>[_e23][_cal][_mbins][_<variant>], trigger = jp (exclusive partition), jp1 / jp2
//       (inclusive single trigger), ht2 (high tower), mb2021m / mb2023 (min-bias level of the 2021 / 2023
//       min-bias pass); e23: the 2023 jet-patch-pass response; cal: that response with its detector pT
//       calibrated to the 2021 jet energy scale per particle bin; mbins: the 5-60 GeV bins of the R_AA
//       reference; variant: a name of variants.h, which selects the data and response files and, for the
//       Bayesian method, the prior tilt. The variant "embstat" adds the replica loop of the embedding
//       statistics (1000 Poisson draws of the entry counts of the fine response, fakes and misses, each
//       times the cell's average weight, re-unfolded: canonical_embstat carries the nominal values with
//       the per-bin RMS as error).
// Result: the |eta| < 0.5 particle block divided by the bin width, in pb/GeV, compared with the published
// table.
// Output <results>/xsec_inversion_R<R>[_<mode>][_bayes<n>].root: canonical (12 bins), canonical_05_09,
// reference (stat), reference_syst.
// root -l -b -q -e 'gSystem->Load("libRooUnfold");' 'unfold.C+("0.5", "jp", "inv", 0)'
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "common.h"
#include "variants.h"

using namespace CrossSectionConfig;
using namespace Stage2;

// ---- trigger efficiency of the high-tower level ------------------------------------------------------

// HT2 trigger efficiency of the embedding relative to the data, per detector-pT bin and eta block: the
// fraction of JP2-patch-tagged jets that also carry an HT2 tower, data / 2023 embedding (the simulated
// raw single-tower ADC in charged-hadron jets is too high). The response HT2 level is scaled by it; the
// +-3 % on the factor is a component of the HT2 band (systematics.C).
static double Ht2Eff(double detPt, int block)
{
   static const double edge[11] = {9.7, 11.5, 13.6, 16.1, 19.0, 22.5, 26.6, 31.4, 37.2, 44.0, 52.0};
   static const double eps[2][10] = {{0.987, 0.957, 0.941, 0.925, 0.937, 0.926, 0.954, 0.969, 0.933, 1.054},
                                     {1.014, 0.766, 0.870, 0.881, 0.873, 0.894, 0.894, 1.013, 0.937, 0.944}};
   for (int k = 0; k < 10; ++k)
      if (detPt >= edge[k] && detPt < edge[k + 1]) return eps[block][k];
   return 1.0;
}

// ---- luminosity --------------------------------------------------------------------------------------

// luminosity sampled by one jet-patch / high-tower level over the kept runs, in pb^-1
static double LevelLuminosity(const AnalysisConfig &cfg, const char *table, const std::vector<int> &badRuns)
{
   TFile lumiFile((cfg.workdir + "inputs/lumi_zilong_full.root").c_str(), "READ");
   auto *lumi = (TH1D *)lumiFile.Get(Form("luminosity_%s", table));
   if (!lumi) throw std::runtime_error("luminosity histogram missing");
   double sum = 0;
   for (int bin = 1; bin <= lumi->GetNbinsX(); ++bin) {
      const char *label = lumi->GetXaxis()->GetBinLabel(bin);
      if (!label || !*label) continue;
      if (IsBadRun(badRuns, std::stoi(label))) continue;
      sum += lumi->GetBinContent(bin);
   }
   return sum;
}

// The chain efficiency eps_chain = P(VPDMB fired and |vz_vpd - vz| < 6 cm | jet event) is an EVENT-level
// quantity, so one number serves every radius: the error-weighted JP0-sample constant over the R = 0.3,
// 0.4, 0.5 measurements, which agree within 0.5 sigma. R = 0.2 (the poorest statistics) sits 2.6 sigma
// below them and is left out; the 0.050-0.062 spread is the systematic on eps_chain.
static double ChainEfficiency(const AnalysisConfig &cfg)
{
   double sumWeights = 0;
   double sumWeightedValues = 0;
   for (const char *radius : {"0.3", "0.4", "0.5"}) {
      TFile epsFile((cfg.workdir + "inputs/eps_chain3_R" + radius + ".root").c_str(), "READ");
      auto *denominator = (TH1D *)epsFile.Get("den0");
      auto *chain = (TH1D *)epsFile.Get("chain0");
      if (!denominator || !chain) throw std::runtime_error("eps_chain3 histograms missing");
      for (int bin = 1; bin <= denominator->GetNbinsX(); ++bin) {
         if (denominator->GetBinContent(bin) <= 0 || chain->GetBinContent(bin) <= 0) continue;
         const double value = chain->GetBinContent(bin) / denominator->GetBinContent(bin);
         const double error = chain->GetBinError(bin) / denominator->GetBinContent(bin);
         sumWeights += 1 / (error * error);
         sumWeightedValues += value / (error * error);
      }
   }
   return sumWeights > 0 ? sumWeightedValues / sumWeights : 0;
}

// min-bias luminosity over the kept runs times the chain efficiency
static double MinBiasLuminosity(const AnalysisConfig &cfg, const std::vector<int> &badRuns)
{
   double lumiSum = 0;
   {
      TFile lumiFile((cfg.workdir + "inputs/lumi_VPDMB_true.root").c_str(), "READ");
      auto *lumi = (TH1D *)lumiFile.Get("luminosity_MBtrue");
      if (!lumi) throw std::runtime_error("min-bias luminosity histogram missing");
      for (int bin = 1; bin <= lumi->GetNbinsX(); ++bin) {
         const char *label = lumi->GetXaxis()->GetBinLabel(bin);
         if (!label || !*label || lumi->GetBinContent(bin) <= 0) continue;
         if (IsBadRun(badRuns, std::stoi(label))) continue;
         lumiSum += lumi->GetBinContent(bin);
      }
   }
   const double eps = ChainEfficiency(cfg);
   printf("[unfold] min-bias: L_MB %.5f pb^-1 x eps_chain %.4f (radius-independent, R=0.3-0.5 JP0 samples) "
          "= %.6f pb^-1\n",
          lumiSum, eps, lumiSum * eps);
   return lumiSum * eps;
}

// ---- fine arrays of one pt-hat sample ----------------------------------------------------------------

// The level-summed response of one sample on the 0.1 GeV cells, with the detector window of every level
// already applied. Each quantity comes weighted (the cross section) and as a plain entry count, which is
// what the outlier filter and the replica loop need.
struct Fine {
   static const int ND = kNDetFine;
   static const int NP = kNParFine;
   std::vector<double> A[2][2];  // response (detector cell, particle cell), per eta block pair
   std::vector<double> Ae[2][2]; // its entry count
   std::vector<double> b[2];     // detector spectrum of the level(s), per detector block
   std::vector<double> be[2];    // its entry count
   std::vector<double> x[2];     // particle spectrum, per particle block
   std::vector<double> xe[2];    // its entry count

   Fine()
   {
      for (int i = 0; i < 2; ++i) {
         b[i].assign(ND, 0);
         be[i].assign(ND, 0);
         x[i].assign(NP, 0);
         xe[i].assign(NP, 0);
         for (int j = 0; j < 2; ++j) {
            A[i][j].assign(ND * NP, 0);
            Ae[i][j].assign(ND * NP, 0);
         }
      }
   }

   void Add(const Fine &other)
   {
      for (int i = 0; i < 2; ++i) {
         for (int f = 0; f < ND; ++f) {
            b[i][f] += other.b[i][f];
            be[i][f] += other.be[i][f];
         }
         for (int g = 0; g < NP; ++g) {
            x[i][g] += other.x[i][g];
            xe[i][g] += other.xe[i][g];
         }
         for (int j = 0; j < 2; ++j)
            for (int c = 0; c < ND * NP; ++c) {
               A[i][j][c] += other.A[i][j][c];
               Ae[i][j][c] += other.Ae[i][j][c];
            }
      }
   }
};

// ---- the published outlier filter --------------------------------------------------------------------

// the filter window: cells [i-9, i+10] on each axis, i.e. 20 cells = 2 GeV, zero outside the array
static const int kBoxBelow = 9;
static const int kBoxAbove = 11;
static const double kBoxMinEntries = 4.0; // a cell survives when its box holds MORE than this

// running sum of v over the filter window, one entry per cell
static std::vector<double> BoxSum1D(const std::vector<double> &v)
{
   const int n = (int)v.size();
   std::vector<double> cumulative(n + 1, 0);
   std::vector<double> box(n, 0);
   for (int i = 0; i < n; ++i) cumulative[i + 1] = cumulative[i] + v[i];
   for (int i = 0; i < n; ++i) {
      const int lo = std::max(0, i - kBoxBelow);
      const int hi = std::min(n, i + kBoxAbove);
      box[i] = cumulative[hi] - cumulative[lo];
   }
   return box;
}

// the same over the two axes of the response, by a summed-area table
static std::vector<double> BoxSum2D(const std::vector<double> &v, int nDet, int nPar)
{
   std::vector<double> cumulative((nDet + 1) * (nPar + 1), 0);
   std::vector<double> box(nDet * nPar, 0);
   auto at = [&](int f, int g) -> double & { return cumulative[f * (nPar + 1) + g]; };
   for (int f = 0; f < nDet; ++f)
      for (int g = 0; g < nPar; ++g) at(f + 1, g + 1) = v[f * nPar + g] + at(f, g + 1) + at(f + 1, g) - at(f, g);
   for (int f = 0; f < nDet; ++f)
      for (int g = 0; g < nPar; ++g) {
         const int f0 = std::max(0, f - kBoxBelow);
         const int f1 = std::min(nDet, f + kBoxAbove);
         const int g0 = std::max(0, g - kBoxBelow);
         const int g1 = std::min(nPar, g + kBoxAbove);
         box[f * nPar + g] = at(f1, g1) - at(f0, g1) - at(f1, g0) + at(f0, g0);
      }
   return box;
}

// which response cells the filter throws away, one flag per cell and eta-block pair
struct DroppedCells {
   std::vector<char> flag[2][2];
};

// mark the response cells whose 20 x 20-cell box holds too few entries, and count what that costs
static DroppedCells MarkThinResponseCells(const Fine &fine, double &nCellsRemoved, double &weightRemoved)
{
   const int nCells = Fine::ND * Fine::NP;
   DroppedCells dropped;
   for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j) {
         const std::vector<double> box = BoxSum2D(fine.Ae[i][j], Fine::ND, Fine::NP);
         dropped.flag[i][j].assign(nCells, 0);
         for (int c = 0; c < nCells; ++c) {
            if (box[c] > kBoxMinEntries) continue;
            dropped.flag[i][j][c] = 1;
            if (fine.Ae[i][j][c] <= 0) continue;
            ++nCellsRemoved;
            weightRemoved += fine.A[i][j][c];
         }
      }
   return dropped;
}

// The detector spectrum loses the pairs of the dropped cells, and the fakes (detector jets with no
// particle partner in range) of every thinly populated detector cell. [parLo, parHi) is the particle
// range of the analysis bins: a pair outside it is a fake, not a migration.
static void FilterDetectorSpectrum(Fine &fine, const DroppedCells &dropped, int parLo, int parHi)
{
   const int ND = Fine::ND;
   const int NP = Fine::NP;
   for (int i = 0; i < 2; ++i) {
      std::vector<double> matchedEntries(ND, 0);
      std::vector<double> matchedWeight(ND, 0);
      std::vector<double> subtractWeight(ND, 0);
      std::vector<double> subtractEntries(ND, 0);
      for (int j = 0; j < 2; ++j)
         for (int f = 0; f < ND; ++f)
            for (int g = parLo; g < parHi; ++g) {
               const int c = f * NP + g;
               matchedEntries[f] += fine.Ae[i][j][c];
               matchedWeight[f] += fine.A[i][j][c];
               if (!dropped.flag[i][j][c]) continue;
               subtractWeight[f] += fine.A[i][j][c];
               subtractEntries[f] += fine.Ae[i][j][c];
            }
      std::vector<double> fakes(ND);
      for (int f = 0; f < ND; ++f) fakes[f] = fine.be[i][f] - matchedEntries[f];
      const std::vector<double> box = BoxSum1D(fakes);
      for (int f = 0; f < ND; ++f) {
         if (box[f] <= kBoxMinEntries) {
            subtractWeight[f] += fine.b[i][f] - matchedWeight[f];
            subtractEntries[f] += fakes[f];
         }
         fine.b[i][f] -= subtractWeight[f];
         fine.be[i][f] -= subtractEntries[f];
      }
   }
}

// The same on the particle side: the pairs of the dropped cells, and the misses (particle jets with no
// detector partner in range) of every thinly populated particle cell.
static void FilterParticleSpectrum(Fine &fine, const DroppedCells &dropped, int detLo, int detHi)
{
   const int NP = Fine::NP;
   for (int j = 0; j < 2; ++j) {
      std::vector<double> matchedEntries(NP, 0);
      std::vector<double> matchedWeight(NP, 0);
      std::vector<double> subtractWeight(NP, 0);
      std::vector<double> subtractEntries(NP, 0);
      for (int i = 0; i < 2; ++i)
         for (int f = detLo; f < detHi; ++f)
            for (int g = 0; g < NP; ++g) {
               const int c = f * NP + g;
               matchedEntries[g] += fine.Ae[i][j][c];
               matchedWeight[g] += fine.A[i][j][c];
               if (!dropped.flag[i][j][c]) continue;
               subtractWeight[g] += fine.A[i][j][c];
               subtractEntries[g] += fine.Ae[i][j][c];
            }
      std::vector<double> misses(NP);
      for (int g = 0; g < NP; ++g) misses[g] = fine.xe[j][g] - matchedEntries[g];
      const std::vector<double> box = BoxSum1D(misses);
      for (int g = 0; g < NP; ++g) {
         if (box[g] <= kBoxMinEntries) {
            subtractWeight[g] += fine.x[j][g] - matchedWeight[g];
            subtractEntries[g] += misses[g];
         }
         fine.x[j][g] -= subtractWeight[g];
         fine.xe[j][g] -= subtractEntries[g];
      }
   }
}

// empty the dropped cells, now that the spectra have been corrected for them
static void ClearDroppedCells(Fine &fine, const DroppedCells &dropped)
{
   const int nCells = Fine::ND * Fine::NP;
   for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
         for (int c = 0; c < nCells; ++c) {
            if (!dropped.flag[i][j][c]) continue;
            fine.A[i][j][c] = 0;
            fine.Ae[i][j][c] = 0;
         }
}

// Drop the response cells that sit in a thinly populated neighbourhood, and take the pairs they hold out
// of the detector and particle spectra as well, so that the three stay consistent. Unmatched entries
// (fakes on the detector side, misses on the particle side) are dropped by the same criterion.
// [detLo, detHi) and [parLo, parHi) are the fine-index ranges of the analysis bins.
static void ApplyOutlierFilter(Fine &fine, int detLo, int detHi, int parLo, int parHi, double &nCellsRemoved,
                               double &weightRemoved)
{
   const DroppedCells dropped = MarkThinResponseCells(fine, nCellsRemoved, weightRemoved);
   FilterDetectorSpectrum(fine, dropped, parLo, parHi);
   FilterParticleSpectrum(fine, dropped, detLo, detHi);
   ClearDroppedCells(fine, dropped);
}

// ---- reading the response file -----------------------------------------------------------------------

// one pt-hat sample of the response file into the fine arrays, summed over the requested levels with the
// detector window of each level applied
static void ReadSample(TFile &file, const std::string &suffix, Fine &fine, const std::vector<int> &levels)
{
   const int ND = Fine::ND;
   const int NP = Fine::NP;
   for (int j = 0; j < 2; ++j) {
      auto *parSpectrum = (TH1D *)file.Get(Form("x_%s%s", kEtaBlockName[j], suffix.c_str()));
      auto *parEntries = (TH1D *)file.Get(Form("x_%s_entries%s", kEtaBlockName[j], suffix.c_str()));
      if (!parSpectrum || !parEntries) throw std::runtime_error("particle histograms missing" + suffix);
      for (int g = 0; g < NP; ++g) {
         fine.x[j][g] = parSpectrum->GetBinContent(g + 1);
         fine.xe[j][g] = parEntries->GetBinContent(g + 1);
      }
   }
   for (int i = 0; i < 2; ++i)
      for (int level : levels) {
         const char *name = kLevelName[level];
         auto *detSpectrum = (TH1D *)file.Get(Form("b_%s_%s%s", name, kEtaBlockName[i], suffix.c_str()));
         auto *detEntries = (TH1D *)file.Get(Form("b_%s_%s_entries%s", name, kEtaBlockName[i], suffix.c_str()));
         if (!detSpectrum || !detEntries)
            throw std::runtime_error(std::string("response histograms missing for level ") + name +
                                     " (rebuild the response)");
         for (int f = 0; f < ND; ++f) {
            const double detPt = detSpectrum->GetXaxis()->GetBinCenter(f + 1);
            if (!KeepDet(level, detPt)) continue;
            const double efficiency = level == kLevelHighTower ? Ht2Eff(detPt, i) : 1.0;
            fine.b[i][f] += detSpectrum->GetBinContent(f + 1) * efficiency;
            fine.be[i][f] += detEntries->GetBinContent(f + 1);
         }
         for (int j = 0; j < 2; ++j) {
            auto *response = (TH2D *)file.Get(
               Form("A_%s_%s_%s%s", name, kEtaBlockName[i], kEtaBlockName[j], suffix.c_str()));
            auto *responseEntries = (TH2D *)file.Get(
               Form("A_%s_%s_%s_entries%s", name, kEtaBlockName[i], kEtaBlockName[j], suffix.c_str()));
            if (!response || !responseEntries)
               throw std::runtime_error(std::string("response matrices missing for level ") + name);
            for (int f = 0; f < ND; ++f) {
               const double detPt = response->GetXaxis()->GetBinCenter(f + 1);
               if (!KeepDet(level, detPt)) continue;
               const double efficiency = level == kLevelHighTower ? Ht2Eff(detPt, i) : 1.0;
               for (int g = 0; g < NP; ++g) {
                  fine.A[i][j][f * NP + g] += response->GetBinContent(f + 1, g + 1) * efficiency;
                  fine.Ae[i][j][f * NP + g] += responseEntries->GetBinContent(f + 1, g + 1);
               }
            }
         }
      }
}

// the pt-hat sample tags stored in a response file, from the names of its particle-entry histograms
static std::vector<std::string> SampleTags(TFile &file)
{
   std::vector<std::string> samples;
   const std::string prefix = "x_00_05_entries_";
   TIter keys(file.GetListOfKeys());
   while (auto *key = (TKey *)keys()) {
      const std::string name = key->GetName();
      if (name.rfind(prefix, 0) == 0) samples.push_back(name.substr(prefix.size()));
   }
   return samples;
}

// ---- analysis binning ---------------------------------------------------------------------------------

// The two axes of the unfolding. Everything is accumulated on the 0.1 GeV cells and rebinned here; the
// detector and particle axes carry the eta blocks unrolled after each other, so an index is
// block * nDetBins + bin.
struct Binning {
   std::vector<double> quotedEdges; // the bins the cross section is quoted in
   std::vector<double> detEdges;    // detector-axis bin edges (may start below the first quoted bin)
   std::vector<double> parEdges;    // particle-axis bin edges (may start below the first quoted bin)
   std::vector<int> detFine;        // fine-cell index of every detector edge
   std::vector<int> parFine;        // fine-cell index of every particle edge
   int nQuoted = 0;                 // quoted bins
   int nDetBins = 0;                // detector bins per eta block
   int nParBins = 0;                // particle bins per eta block
   int nBlocks = 0;                 // eta blocks
   int nDetTotal = 0;               // detector bins over all blocks (the rows of A)
   int nParTotal = 0;               // particle bins over all blocks (the columns of A)
   int detOffset = 0;               // index of the first quoted bin on the detector axis
   int parOffset = 0;               // index of the first quoted bin on the particle axis
};

// the published 12 bins (6.9-52 GeV), or the 5-60 GeV bins of the R_AA reference with "mbins". The
// Bayesian method adds one detector bin and two particle bins below the first quoted bin, so that the
// migrations across 6.9 GeV sit inside the matrix; the inversion stays on the published square binning.
static Binning MakeBinning(bool bayes, bool mbins)
{
   if (mbins && !bayes) throw std::runtime_error("mbins needs the Bayesian method");
   std::vector<double> published(McBins().begin(), McBins().end());
   while (!published.empty() && published.back() > 52.5) published.pop_back();

   Binning bins;
   bins.quotedEdges = mbins ? std::vector<double>{5, 10, 15, 20, 25, 30, 35, 40, 50, 60} : published;
   bins.detEdges = mbins ? std::vector<double>{5.0,  6.9,  8.2,  9.7,  11.5, 13.6, 16.1, 19.0,
                                               22.5, 26.6, 31.4, 37.2, 44.0, 52.0, 60.0}
                         : published;
   bins.parEdges = bins.quotedEdges;
   if (bayes) {
      if (!mbins) bins.detEdges.insert(bins.detEdges.begin(), 5.0);
      const std::vector<double> lowPar = mbins ? std::vector<double>{3.5} : std::vector<double>{3.5, 5.2};
      bins.parEdges.insert(bins.parEdges.begin(), lowPar.begin(), lowPar.end());
   }
   bins.nQuoted = (int)bins.quotedEdges.size() - 1;
   bins.nDetBins = (int)bins.detEdges.size() - 1;
   bins.nParBins = (int)bins.parEdges.size() - 1;
   bins.nBlocks = kNEtaBlocks;
   bins.nDetTotal = bins.nBlocks * bins.nDetBins;
   bins.nParTotal = bins.nBlocks * bins.nParBins;
   bins.detOffset = bins.nDetBins - bins.nQuoted;
   bins.parOffset = bins.nParBins - bins.nQuoted;
   bins.detFine.resize(bins.nDetBins + 1);
   bins.parFine.resize(bins.nParBins + 1);
   for (int k = 0; k <= bins.nDetBins; ++k) bins.detFine[k] = FineIndexOfEdge(bins.detEdges[k], kDetFineMin);
   for (int k = 0; k <= bins.nParBins; ++k) bins.parFine[k] = FineIndexOfEdge(bins.parEdges[k], kParFineMin);
   return bins;
}

// ---- the mode string ------------------------------------------------------------------------------------

// what a mode string asks for
struct ModeSpec {
   std::string trigger;         // jp | jp1 | jp2 | ht2 | mb2021m | mb2023
   std::vector<int> levels;     // the response / data levels the trigger is built from
   bool minBias = false;        // the min-bias data and luminosity
   bool use2023 = false;        // the 2023 jet-patch-pass response
   bool calibrateScale = false; // calibrate its detector pT to the 2021 energy scale
   bool wideBins = false;       // the 5-60 GeV bins of the R_AA reference
   std::string variant;         // the systematic variant name ("" nominal)
};

static ModeSpec ParseMode(const std::string &mode)
{
   const std::vector<std::string> tokens = SplitWords(mode, '_');
   if (tokens.empty()) throw std::runtime_error("empty mode");
   ModeSpec spec;
   spec.trigger = tokens[0];
   spec.use2023 = spec.trigger == "mb2023";
   for (size_t k = 1; k < tokens.size(); ++k) {
      if (tokens[k] == "e23") {
         spec.use2023 = true;
      } else if (tokens[k] == "e23cal") {
         spec.use2023 = true;
         spec.calibrateScale = true;
      } else if (tokens[k] == "cal") {
         spec.calibrateScale = true;
      } else if (tokens[k] == "mbins") {
         spec.wideBins = true;
      } else {
         spec.variant += (spec.variant.empty() ? "" : "_") + tokens[k];
      }
   }
   spec.minBias = spec.trigger.rfind("mb", 0) == 0;
   if (spec.trigger == "jp1") spec.levels = {4};
   else if (spec.trigger == "jp2") spec.levels = {5};
   else if (spec.trigger == "ht2") spec.levels = {kLevelHighTower};
   else if (spec.minBias) spec.levels = {kLevelMinBias};
   else if (spec.trigger == "jp") spec.levels = {0, 1, 2};
   else throw std::runtime_error("unknown trigger in mode " + mode);
   return spec;
}

// the response file suffix the mode asks for ("" = the nominal 2021 jet-patch-pass response)
static std::string ResponseSuffix(const ModeSpec &spec)
{
   if (spec.trigger == "mb2023") return "_mb2023";
   if (spec.trigger == "mb2021m") return "_mb2021m";
   return spec.use2023 ? "_e23" : "";
}

// ---- data spectra ---------------------------------------------------------------------------------------

// The level-weighted data vector b_i = sum_l n_{l,i} / L_l, its variance sum_l n_{l,i} / L_l^2, and the
// plain entry count over the levels (which the publication uses as its error estimator).
static void ReadDataSpectra(TFile &file, const std::vector<int> &levels, const Binning &bins,
                            const double luminosity[kNLevels], TVectorD &b, TVectorD &bVar, TVectorD &nEntries)
{
   for (int i = 0; i < bins.nBlocks; ++i)
      for (int level : levels) {
         auto *spectrum = (TH1D *)file.Get(Form("d550_%s_%s", kLevelName[level], kEtaBlockName[i]));
         if (!spectrum) throw std::runtime_error("data histogram missing (run build_data for this mode)");
         for (int k = 0; k < bins.nDetBins; ++k)
            for (int f = bins.detFine[k]; f < bins.detFine[k + 1]; ++f) {
               const double detPt = spectrum->GetXaxis()->GetBinCenter(f + 1);
               if (!KeepDet(level, detPt)) continue;
               const double n = spectrum->GetBinContent(f + 1);
               b[i * bins.nDetBins + k] += n / luminosity[level];
               bVar[i * bins.nDetBins + k] += n / (luminosity[level] * luminosity[level]);
               nEntries[i * bins.nDetBins + k] += n;
            }
      }
}

// ---- response assembly -----------------------------------------------------------------------------------

// every pt-hat sample of one response file, filtered and summed. keepSamples, when given, collects the
// filtered samples for the replica loop.
static Fine ReadFilteredResponse(TFile &file, const std::vector<int> &levels, const Binning &bins, bool filter,
                                 double &nCellsRemoved, double &weightRemoved, std::vector<Fine> *keepSamples,
                                 size_t *nSamples)
{
   Fine total;
   const std::vector<std::string> samples = SampleTags(file);
   for (const std::string &sample : samples) {
      Fine fine;
      ReadSample(file, "_" + sample, fine, levels);
      if (filter)
         ApplyOutlierFilter(fine, bins.detFine[0], bins.detFine[bins.nDetBins], bins.parFine[0],
                            bins.parFine[bins.nParBins], nCellsRemoved, weightRemoved);
      total.Add(fine);
      if (keepSamples) keepSamples->push_back(fine);
   }
   if (nSamples) *nSamples = samples.size();
   return total;
}

// rebin the fine arrays to the analysis bins, with the eta blocks unrolled along both axes
static void RebinToAnalysisBins(const Fine &fine, const Binning &bins, TMatrixD &A, TVectorD &bEmb, TVectorD &x)
{
   A.ResizeTo(bins.nDetTotal, bins.nParTotal);
   bEmb.ResizeTo(bins.nDetTotal);
   x.ResizeTo(bins.nParTotal);
   A.Zero();
   bEmb.Zero();
   x.Zero();
   for (int i = 0; i < bins.nBlocks; ++i) {
      for (int k = 0; k < bins.nDetBins; ++k)
         for (int f = bins.detFine[k]; f < bins.detFine[k + 1]; ++f) bEmb[i * bins.nDetBins + k] += fine.b[i][f];
      for (int j = 0; j < bins.nBlocks; ++j)
         for (int k = 0; k < bins.nDetBins; ++k)
            for (int f = bins.detFine[k]; f < bins.detFine[k + 1]; ++f)
               for (int q = 0; q < bins.nParBins; ++q)
                  for (int g = bins.parFine[q]; g < bins.parFine[q + 1]; ++g)
                     A(i * bins.nDetBins + k, j * bins.nParBins + q) += fine.A[i][j][f * Fine::NP + g];
   }
   for (int j = 0; j < bins.nBlocks; ++j)
      for (int q = 0; q < bins.nParBins; ++q)
         for (int g = bins.parFine[q]; g < bins.parFine[q + 1]; ++g) x[j * bins.nParBins + q] += fine.x[j][g];
}

// The prior tilt of the Bayesian regularisation component: truth events reweighted by (pT/20 GeV)^tilt,
// i.e. the columns of A and x. Fakes (b^emb - m) are not truth events and keep their weight, so b^emb
// follows the change of the matched part only.
static const double kTiltPivot = 20.0; // GeV

static void ApplyPriorTilt(const Binning &bins, double tilt, TMatrixD &A, TVectorD &bEmb, TVectorD &x)
{
   if (tilt == 0) return;
   TVectorD matchedBefore(bins.nDetTotal);
   for (int r = 0; r < bins.nDetTotal; ++r) {
      matchedBefore[r] = 0;
      for (int c = 0; c < bins.nParTotal; ++c) matchedBefore[r] += A(r, c);
   }
   for (int c = 0; c < bins.nParTotal; ++c) {
      const int q = c % bins.nParBins;
      const double weight = std::pow(0.5 * (bins.parEdges[q] + bins.parEdges[q + 1]) / kTiltPivot, tilt);
      x[c] *= weight;
      for (int r = 0; r < bins.nDetTotal; ++r) A(r, c) *= weight;
   }
   for (int r = 0; r < bins.nDetTotal; ++r) {
      double matchedAfter = 0;
      for (int c = 0; c < bins.nParTotal; ++c) matchedAfter += A(r, c);
      bEmb[r] += matchedAfter - matchedBefore[r];
   }
}

// m_r = the matched (in-range) part of detector bin r
static TVectorD MatchedPerDetBin(const TMatrixD &A, const Binning &bins)
{
   TVectorD matched(bins.nDetTotal);
   for (int r = 0; r < bins.nDetTotal; ++r) {
      matched[r] = 0;
      for (int c = 0; c < bins.nParTotal; ++c) matched[r] += A(r, c);
   }
   return matched;
}

// ---- the two unfolding methods ----------------------------------------------------------------------------

// x_unf = M^-1 b with the covariance U B U^T. B is the exact variance sum_l n_l / L_l^2 of the
// level-weighted data; Bp is the publication's estimator b^2 / N_entries, which ignores the level weights
// (the JP0 level carries 48 x the JP1 weight) and is smaller wherever a high-weight level adds few entries.
static void SolveByInversion(const TMatrixD &M, const TVectorD &b, const TVectorD &bVar, const TVectorD &nEntries,
                             const Binning &bins, TVectorD &xUnfolded, TMatrixD &cov, TMatrixD &covPub)
{
   if (bins.nDetTotal != bins.nParTotal) throw std::runtime_error("inversion needs a square matrix");
   TMatrixD unfoldingMatrix(M);
   unfoldingMatrix.Invert();
   xUnfolded = unfoldingMatrix * b;
   TMatrixD dataCov(bins.nDetTotal, bins.nDetTotal);
   for (int r = 0; r < bins.nDetTotal; ++r) dataCov(r, r) = bVar[r];
   cov = unfoldingMatrix * dataCov * TMatrixD(TMatrixD::kTransposed, unfoldingMatrix);
   TMatrixD dataCovPub(bins.nDetTotal, bins.nDetTotal);
   for (int r = 0; r < bins.nDetTotal; ++r)
      dataCovPub(r, r) = nEntries[r] > 0 ? b[r] * b[r] / nEntries[r] : 0.0;
   covPub = unfoldingMatrix * dataCovPub * TMatrixD(TMatrixD::kTransposed, unfoldingMatrix);
}

// the per-bin purity m/b^emb applied to the data, then RooUnfoldBayes on the fake-free response (A, x)
static void SolveByBayes(const TMatrixD &A, const TVectorD &bEmb, const TVectorD &x, const TVectorD &matched,
                         const TVectorD &b, const TVectorD &bVar, const Binning &bins, int nIter,
                         TVectorD &xUnfolded, TMatrixD &cov)
{
   const int nDet = bins.nDetTotal;
   const int nPar = bins.nParTotal;
   TH1D measured("hmeas", "", nDet, 0, nDet);
   TH1D measuredMc("hmeasMC", "", nDet, 0, nDet);
   TH1D truth("htru", "", nPar, 0, nPar);
   TH2D response("hA", "", nDet, 0, nDet, nPar, 0, nPar);
   for (int r = 0; r < nDet; ++r) {
      const double purity = bEmb[r] > 0 ? matched[r] / bEmb[r] : 0.0;
      measured.SetBinContent(r + 1, b[r] * purity);
      measured.SetBinError(r + 1, std::sqrt(std::max(0.0, bVar[r])) * purity);
      measuredMc.SetBinContent(r + 1, matched[r]);
      for (int c = 0; c < nPar; ++c) response.SetBinContent(r + 1, c + 1, A(r, c));
   }
   for (int c = 0; c < nPar; ++c) truth.SetBinContent(c + 1, x[c]);
   RooUnfoldResponse rooResponse(&measuredMc, &truth, &response);
   rooResponse.UseOverflow(false);
   RooUnfoldBayes bayes(&rooResponse, &measured, nIter);
   bayes.SetVerbose(0);
   TH1 *unfolded = bayes.Hreco(RooUnfold::kCovariance);
   TMatrixD bayesCov = bayes.Ereco(RooUnfold::kCovariance);
   for (int c = 0; c < nPar; ++c) {
      xUnfolded[c] = unfolded->GetBinContent(c + 1);
      for (int d = 0; d < nPar; ++d) cov(c, d) = bayesCov(c, d);
   }
}

// Unfold one response (A, b^emb, x) with the data b. M_ik = (b^emb_i / m_i) A_ik / x_k is the migration
// matrix of the published construction; the Bayesian method uses it only for the fold check.
static void Solve(const TMatrixD &A, const TVectorD &bEmb, const TVectorD &x, const TVectorD &b,
                  const TVectorD &bVar, const TVectorD &nEntries, const Binning &bins, bool bayes, int nIter,
                  TVectorD &xUnfolded, TMatrixD &cov, TMatrixD &covPub, TMatrixD &M)
{
   const TVectorD matched = MatchedPerDetBin(A, bins);
   M.ResizeTo(bins.nDetTotal, bins.nParTotal);
   xUnfolded.ResizeTo(bins.nParTotal);
   cov.ResizeTo(bins.nParTotal, bins.nParTotal);
   covPub.ResizeTo(bins.nParTotal, bins.nParTotal);
   cov.Zero();
   covPub.Zero();
   for (int r = 0; r < bins.nDetTotal; ++r)
      for (int c = 0; c < bins.nParTotal; ++c)
         M(r, c) = (matched[r] > 0 && x[c] > 0) ? (bEmb[r] / matched[r]) * A(r, c) / x[c] : 0.0;
   if (bayes)
      SolveByBayes(A, bEmb, x, matched, b, bVar, bins, nIter, xUnfolded, cov);
   else
      SolveByInversion(M, b, bVar, nEntries, bins, xUnfolded, cov, covPub);
}

// ---- the jet-energy-scale calibration of the 2023 response ---------------------------------------------------

// <pT_det / pT_part> of one particle bin, over the |eta| < 0.5 blocks of the response
static double MeanDetOverPar(const Fine &fine, const Binning &bins, int parBin)
{
   double sumWeight = 0;
   double sumRatio = 0;
   for (int g = bins.parFine[parBin]; g < bins.parFine[parBin + 1]; ++g) {
      const double parPt = (g + 0.5) * kFineCell;
      for (int f = 0; f < Fine::ND; ++f) {
         const double weight = fine.A[0][0][f * Fine::NP + g];
         if (weight <= 0) continue;
         sumWeight += weight;
         sumRatio += weight * (kDetFineMin + (f + 0.5) * kFineCell) / parPt;
      }
   }
   return sumWeight > 0 ? sumRatio / sumWeight : 1.0;
}

// add a weight at a fractional detector-cell position, split linearly between the two nearest cells
static void AddSplitBetweenCells(std::vector<double> &target, int stride, int parCell, double cellPosition,
                                 double weight)
{
   const int cell = (int)std::floor(cellPosition);
   const double fraction = cellPosition - cell;
   if (cell >= 0 && cell < Fine::ND) target[cell * stride + parCell] += weight * (1 - fraction);
   if (cell + 1 >= 0 && cell + 1 < Fine::ND) target[(cell + 1) * stride + parCell] += weight * fraction;
}

// Move every detector-side cell of the 2023 response to d' = d x r, r being the ratio of the 2021 and the
// 2023 <pT_det / pT_part> in the particle bin the cell belongs to. This puts the 2023 sample on the 2021
// jet energy scale without touching its resolution or efficiency.
static void CalibrateDetectorScale(Fine &response, const Fine &reference, const Binning &bins)
{
   std::vector<double> scale(bins.nParBins, 1.0);
   printf("[unfold] JES calibration to the 2021 sample, per particle bin:");
   for (int q = 0; q < bins.nParBins; ++q) {
      const double referenceRatio = MeanDetOverPar(reference, bins, q);
      const double responseRatio = MeanDetOverPar(response, bins, q);
      scale[q] = (referenceRatio > 0 && responseRatio > 0) ? referenceRatio / responseRatio : 1.0;
      printf(" %.4f", scale[q]);
   }
   printf("\n");
   // the scale of a particle cell, and of a detector cell taken at the same pT
   auto scaleOfParCell = [&](int g) {
      for (int q = 0; q < bins.nParBins; ++q)
         if (g >= bins.parFine[q] && g < bins.parFine[q + 1]) return scale[q];
      return g < bins.parFine[0] ? scale[0] : scale[bins.nParBins - 1];
   };
   auto scaleOfDetCell = [&](int f) {
      const double detPt = kDetFineMin + (f + 0.5) * kFineCell;
      for (int q = 0; q < bins.nParBins; ++q)
         if (detPt >= bins.parEdges[q] && detPt < bins.parEdges[q + 1]) return scale[q];
      return detPt < bins.parEdges[0] ? scale[0] : scale[bins.nParBins - 1];
   };
   Fine calibrated = response;
   for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) std::fill(calibrated.A[i][j].begin(), calibrated.A[i][j].end(), 0.0);
      std::fill(calibrated.b[i].begin(), calibrated.b[i].end(), 0.0);
      for (int j = 0; j < 2; ++j)
         for (int f = 0; f < Fine::ND; ++f)
            for (int g = 0; g < Fine::NP; ++g) {
               const double weight = response.A[i][j][f * Fine::NP + g];
               if (weight <= 0) continue;
               const double detPt = kDetFineMin + (f + 0.5) * kFineCell;
               const double moved = detPt * scaleOfParCell(g);
               AddSplitBetweenCells(calibrated.A[i][j], Fine::NP, g, (moved - kDetFineMin) / kFineCell - 0.5, weight);
            }
      for (int f = 0; f < Fine::ND; ++f) {
         const double weight = response.b[i][f];
         if (weight <= 0) continue;
         const double detPt = kDetFineMin + (f + 0.5) * kFineCell;
         const double moved = detPt * scaleOfDetCell(f);
         AddSplitBetweenCells(calibrated.b[i], 1, 0, (moved - kDetFineMin) / kFineCell - 0.5, weight);
      }
   }
   response = calibrated;
}

// ---- the replica loop of the embedding statistics -------------------------------------------------------------

static const int kNReplicas = 1000;
static const UInt_t kReplicaSeed = 31337;

// One Poisson replica of one filtered sample: every response cell, fake and miss entry count is redrawn
// and scaled by its average weight. [detLo, detHi) x [parLo, parHi) is the analysis range; pairs outside
// it count as fakes / misses.
static void DrawPoissonReplica(const Fine &fine, Fine &replica, TRandom3 &rng, int detLo, int detHi, int parLo,
                               int parHi)
{
   const int ND = Fine::ND;
   const int NP = Fine::NP;
   for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
         for (int c = 0; c < ND * NP; ++c) {
            const double entries = fine.Ae[i][j][c];
            if (entries <= 0) {
               replica.A[i][j][c] = 0;
               replica.Ae[i][j][c] = 0;
               continue;
            }
            const double drawn = rng.Poisson(entries);
            replica.Ae[i][j][c] = drawn;
            replica.A[i][j][c] = drawn * fine.A[i][j][c] / entries;
         }
   for (int i = 0; i < 2; ++i)
      for (int f = 0; f < ND; ++f) {
         double matchedEntries = 0;
         double matchedWeight = 0;
         double replicaWeight = 0;
         for (int j = 0; j < 2; ++j)
            for (int g = parLo; g < parHi; ++g) {
               const int c = f * NP + g;
               matchedEntries += fine.Ae[i][j][c];
               matchedWeight += fine.A[i][j][c];
               replicaWeight += replica.A[i][j][c];
            }
         const double fakeEntries = fine.be[i][f] - matchedEntries;
         const double fakeWeight = fine.b[i][f] - matchedWeight;
         const double drawn = fakeEntries > 0 ? rng.Poisson(fakeEntries) : 0;
         replica.b[i][f] = replicaWeight + (fakeEntries > 0 ? drawn * fakeWeight / fakeEntries : 0);
         replica.be[i][f] = 0;
      }
   for (int j = 0; j < 2; ++j)
      for (int g = 0; g < NP; ++g) {
         double matchedEntries = 0;
         double matchedWeight = 0;
         double replicaWeight = 0;
         for (int i = 0; i < 2; ++i)
            for (int f = detLo; f < detHi; ++f) {
               const int c = f * NP + g;
               matchedEntries += fine.Ae[i][j][c];
               matchedWeight += fine.A[i][j][c];
               replicaWeight += replica.A[i][j][c];
            }
         const double missEntries = fine.xe[j][g] - matchedEntries;
         const double missWeight = fine.x[j][g] - matchedWeight;
         const double drawn = missEntries > 0 ? rng.Poisson(missEntries) : 0;
         replica.x[j][g] = replicaWeight + (missEntries > 0 ? drawn * missWeight / missEntries : 0);
         replica.xe[j][g] = 0;
      }
}

// per-bin RMS of the unfolded result over the replicas
static TVectorD ReplicaRms(const std::vector<Fine> &samples, const Binning &bins, const TVectorD &b,
                           const TVectorD &bVar, const TVectorD &nEntries, bool bayes, int nIter, double tilt)
{
   TRandom3 rng(kReplicaSeed);
   TVectorD sum(bins.nParTotal);
   TVectorD sumSquares(bins.nParTotal);
   sum.Zero();
   sumSquares.Zero();
   Fine replica;
   for (int r = 0; r < kNReplicas; ++r) {
      Fine total;
      for (const Fine &sample : samples) {
         DrawPoissonReplica(sample, replica, rng, bins.detFine[0], bins.detFine[bins.nDetBins], bins.parFine[0],
                            bins.parFine[bins.nParBins]);
         total.Add(replica);
      }
      TMatrixD A;
      TVectorD bEmb;
      TVectorD x;
      RebinToAnalysisBins(total, bins, A, bEmb, x);
      ApplyPriorTilt(bins, tilt, A, bEmb, x);
      TVectorD xUnfolded;
      TMatrixD cov;
      TMatrixD covPub;
      TMatrixD M;
      Solve(A, bEmb, x, b, bVar, nEntries, bins, bayes, nIter, xUnfolded, cov, covPub, M);
      for (int c = 0; c < bins.nParTotal; ++c) {
         sum[c] += xUnfolded[c];
         sumSquares[c] += xUnfolded[c] * xUnfolded[c];
      }
      if ((r + 1) % 100 == 0) {
         printf("[unfold] replica %d / %d\n", r + 1, kNReplicas);
         fflush(stdout);
      }
   }
   TVectorD rms(bins.nParTotal);
   for (int c = 0; c < bins.nParTotal; ++c) {
      const double mean = sum[c] / kNReplicas;
      rms[c] = std::sqrt(std::max(0.0, sumSquares[c] / kNReplicas - mean * mean));
   }
   return rms;
}

// ---- printouts -------------------------------------------------------------------------------------------

// how much of the detector spectrum is unmatched, and how much of the particle spectrum is reconstructed
static void PrintFakeAndEfficiency(const TMatrixD &A, const TVectorD &bEmb, const TVectorD &x, const Binning &bins)
{
   const TVectorD matched = MatchedPerDetBin(A, bins);
   printf("[unfold] fake fraction 1 - m/b_emb (detector) and efficiency m_col/x (particle), |eta|<0.5 block:\n");
   for (int k = 0; k < std::max(bins.nDetBins, bins.nParBins); ++k) {
      char detText[64] = "";
      char parText[64] = "";
      if (k < bins.nDetBins)
         snprintf(detText, 64, "[%5.1f,%5.1f)  fake %.3f", bins.detEdges[k], bins.detEdges[k + 1],
                  bEmb[k] > 0 ? 1 - matched[k] / bEmb[k] : 0);
      if (k < bins.nParBins) {
         double column = 0;
         for (int r = 0; r < bins.nDetTotal; ++r) column += A(r, k);
         snprintf(parText, 64, "[%5.1f,%5.1f)  eff %.3f", bins.parEdges[k], bins.parEdges[k + 1],
                  x[k] > 0 ? column / x[k] : 0);
      }
      printf("   %-32s %s\n", detText, parText);
   }
}

// fold the published cross section back through M and compare with the measured detector spectrum
static void PrintFoldCheck(const TMatrixD &M, const TVectorD &b, const TVectorD &xUnfolded, TH1D *published,
                           const Binning &bins)
{
   printf("[unfold] fold check  b_data / (M x published), detector |eta|<0.5 bins:\n");
   TVectorD truth(bins.nParTotal);
   for (int c = 0; c < bins.nParTotal; ++c) truth[c] = xUnfolded[c];
   for (int k = 0; k < bins.nQuoted; ++k) {
      const double centre = 0.5 * (bins.quotedEdges[k] + bins.quotedEdges[k + 1]);
      truth[bins.parOffset + k] =
         published->GetBinContent(published->FindBin(centre)) * (bins.quotedEdges[k + 1] - bins.quotedEdges[k]);
   }
   const TVectorD predicted = M * truth;
   for (int k = 0; k < bins.nQuoted; ++k) {
      const int r = bins.detOffset + k;
      printf("   [%5.1f,%5.1f)  data %10.4g  pred %10.4g  ratio %.3f\n", bins.quotedEdges[k],
             bins.quotedEdges[k + 1], b[r], predicted[r], predicted[r] > 0 ? b[r] / predicted[r] : 0);
   }
}

// bin-by-bin comparison with the published cross section
static void PrintRatioToPublished(const TH1D *result, TH1D *published, const Binning &bins, bool bayes)
{
   const char *method = bayes ? "bayes" : "inversion";
   printf("[ratio %s/published] |eta|<0.5:\n", method);
   double chi2 = 0;
   int ndf = 0;
   for (int k = 0; k < bins.nQuoted; ++k) {
      const double centre = 0.5 * (bins.quotedEdges[k] + bins.quotedEdges[k + 1]);
      const double value = result->GetBinContent(k + 1);
      const double error = result->GetBinError(k + 1);
      const double reference = published ? published->GetBinContent(published->FindBin(centre)) : 0;
      const double referenceError = published ? published->GetBinError(published->FindBin(centre)) : 0;
      printf("   pT=%5.2f  ours=%9.4g +- %.2g%%  published=%9.4g  ratio=%.3f\n", centre, value,
             value > 0 ? 100 * error / value : 0, reference, reference > 0 ? value / reference : 0);
      if (reference <= 0) continue;
      const double variance = error * error + referenceError * referenceError;
      if (variance <= 0) continue;
      chi2 += (value - reference) * (value - reference) / variance;
      ++ndf;
   }
   printf("[ratio %s/published] chi2/ndf = %.2f / %d\n", method, chi2, ndf);
}

// ---- the cross section ------------------------------------------------------------------------------------

// one eta block of the result: the unfolded particle spectrum divided by the bin width, in pb/GeV
static TH1D *MakeResultHistogram(const char *name, const TVectorD &xUnfolded, const TMatrixD &cov,
                                 const Binning &bins, int block)
{
   TH1D *result = new TH1D(name, ";p_{T} [GeV/c];d^{2}#sigma/dp_{T}d#eta [pb/(GeV/c)]", bins.nQuoted,
                           bins.quotedEdges.data());
   result->SetDirectory(0);
   for (int k = 0; k < bins.nQuoted; ++k) {
      const int c = block * bins.nParBins + bins.parOffset + k;
      const double width = bins.quotedEdges[k + 1] - bins.quotedEdges[k];
      result->SetBinContent(k + 1, xUnfolded[c] / width);
      result->SetBinError(k + 1, std::sqrt(std::max(0.0, cov(c, c))) / width);
   }
   return result;
}

// mode / method / nIter: see the header
void unfold(const char *jetR = "0.5", const char *mode = "jp", const char *method = "inv", int nIter = 4)
{
   const bool filter = true;
   AnalysisConfig cfg;
   const std::vector<int> badRuns = AllBadRuns(cfg);
   const std::string modeText = mode;
   const bool bayes = std::string(method) == "bayes";

   // 1. what to build, from which files
   const ModeSpec spec = ParseMode(modeText);
   const Syst::Variant variant = Syst::Parse(spec.variant);
   const Binning bins = MakeBinning(bayes, spec.wideBins);
   const std::string radiusDir = RadiusDir(jetR);
   const std::string outSuffix =
      (modeText == "jp" ? "" : "_" + modeText) + (bayes ? Form("_bayes%d", nIter) : "");
   const std::string dataFile =
      radiusDir + "data_blocks_R" + jetR + (spec.minBias ? "_mb" : "") + variant.dataTag() + ".root";
   const std::string responseFile =
      radiusDir + "response_blocks_R" + jetR + ResponseSuffix(spec) + variant.respTag() + ".root";
   TFile dataInput(dataFile.c_str(), "READ");
   TFile responseInput(responseFile.c_str(), "READ");
   if (dataInput.IsZombie() || responseInput.IsZombie())
      throw std::runtime_error("inputs missing (" + dataFile + ", " + responseFile +
                               "): run build_data / build_resp first");
   printf("[unfold] mode %s  method %s%s  levels:", mode, method, bayes ? Form(" (%d iterations)", nIter) : "");
   for (int level : spec.levels) printf(" %s", kLevelName[level]);
   printf("\n[unfold] data %s\n[unfold] response %s\n", dataFile.c_str(), responseFile.c_str());
   if (!spec.variant.empty())
      printf("[unfold] variant %s%s%s\n", spec.variant.c_str(),
             variant.priorTilt != 0 ? Form(" (prior tilt %+.2f)", variant.priorTilt) : "",
             variant.embStat ? " (replica loop)" : "");

   // 2. the luminosity each level sampled over the kept runs
   double luminosity[kNLevels] = {0, 0, 0, 0, 0, 0, 0};
   for (int level : spec.levels)
      luminosity[level] =
         level == kLevelMinBias ? MinBiasLuminosity(cfg, badRuns) : LevelLuminosity(cfg, kLumiName[level], badRuns);
   printf("[unfold] luminosity per level (kept runs):");
   for (int level : spec.levels) printf("  %s %.5f pb^-1", kLevelName[level], luminosity[level]);
   printf("\n");

   // 3. the data vector and its variance
   TVectorD b(bins.nDetTotal);
   TVectorD bVar(bins.nDetTotal);
   TVectorD nEntries(bins.nDetTotal); // plain entry count: the publication's error estimator
   ReadDataSpectra(dataInput, spec.levels, bins, luminosity, b, bVar, nEntries);

   // 4. the response, filtered per pt-hat sample and summed
   std::vector<Fine> samples; // kept for the replica loop
   double nCellsRemoved = 0;
   double weightRemoved = 0;
   size_t nSamples = 0;
   Fine total = ReadFilteredResponse(responseInput, spec.levels, bins, filter, nCellsRemoved, weightRemoved,
                                     variant.embStat ? &samples : nullptr, &nSamples);
   printf("[unfold] %zu samples; outlier filter %s: %.0f filled cells removed (%.3g pb)\n", nSamples,
          filter ? "on" : "off", nCellsRemoved, weightRemoved);

   // 5. optionally put the 2023 response on the 2021 jet energy scale
   if (spec.calibrateScale) {
      TFile referenceInput((radiusDir + "response_blocks_R" + jetR + ".root").c_str(), "READ");
      double refCells = 0;
      double refWeight = 0;
      const Fine reference =
         ReadFilteredResponse(referenceInput, spec.levels, bins, filter, refCells, refWeight, nullptr, nullptr);
      CalibrateDetectorScale(total, reference, bins);
   }

   // 6. rebin, tilt the prior if the variant asks for it, and unfold
   TMatrixD A;
   TVectorD bEmb;
   TVectorD x;
   RebinToAnalysisBins(total, bins, A, bEmb, x);
   const double tilt = bayes ? variant.priorTilt : 0.0;
   ApplyPriorTilt(bins, tilt, A, bEmb, x);
   PrintFakeAndEfficiency(A, bEmb, x, bins);
   TVectorD xUnfolded;
   TMatrixD cov;
   TMatrixD covPub; // covariance with the publication's data variance b^2/N_entries (inversion only)
   TMatrixD M;
   Solve(A, bEmb, x, b, bVar, nEntries, bins, bayes, nIter, xUnfolded, cov, covPub, M);
   if (!bayes) {
      printf("[unfold] data variance, exact / publication estimator, detector |eta|<0.5 bins:");
      for (int k = 0; k < bins.nDetBins; ++k)
         printf(" %.2f", covPub(k, k) > 0 ? std::sqrt(cov(k, k) / covPub(k, k)) : 0.0);
      printf("\n");
   }

   // 7. the embedding-statistics component, when this variant asks for it
   TVectorD embeddingRms(bins.nParTotal);
   embeddingRms.Zero();
   if (variant.embStat) {
      embeddingRms = ReplicaRms(samples, bins, b, bVar, nEntries, bayes, nIter, tilt);
      printf("[unfold] embedding statistics, RMS / value per published bin (|eta|<0.5):");
      for (int k = 0; k < bins.nQuoted; ++k) {
         const int c = bins.parOffset + k;
         printf(" %.3f", xUnfolded[c] > 0 ? embeddingRms[c] / xUnfolded[c] : 0);
      }
      printf("\n");
   }

   // 8. the result histograms
   TH1D *result[2] = {nullptr, nullptr};
   TH1D *resultEmbStat[2] = {nullptr, nullptr};
   for (int j = 0; j < bins.nBlocks; ++j) {
      result[j] = MakeResultHistogram(j == 0 ? "canonical" : "canonical_05_09", xUnfolded, cov, bins, j);
      if (!variant.embStat) continue;
      resultEmbStat[j] = (TH1D *)result[j]->Clone(j == 0 ? "canonical_embstat" : "canonical_05_09_embstat");
      resultEmbStat[j]->SetDirectory(0);
      for (int k = 0; k < bins.nQuoted; ++k)
         resultEmbStat[j]->SetBinError(k + 1, embeddingRms[j * bins.nParBins + bins.parOffset + k] /
                                                 (bins.quotedEdges[k + 1] - bins.quotedEdges[k]));
   }
   TH1D *resultPubErr = nullptr;
   if (!bayes) {
      resultPubErr = (TH1D *)result[0]->Clone("canonical_pubErr");
      resultPubErr->SetDirectory(0);
      for (int k = 0; k < bins.nQuoted; ++k) {
         const int c = bins.parOffset + k;
         resultPubErr->SetBinError(k + 1, std::sqrt(std::max(0.0, covPub(c, c))) /
                                             (bins.quotedEdges[k + 1] - bins.quotedEdges[k]));
      }
   }

   // 9. the published table, the fold check and the comparison
   TH1D *publishedStat = nullptr;
   TH1D *publishedSyst = nullptr;
   {
      const std::string path = cfg.workdir + "inputs/jet_cross_section_publishedR" + jetR + ".root";
      if (!gSystem->AccessPathName(path.c_str())) {
         TFile publishedFile(path.c_str(), "READ");
         if (auto *h = (TH1D *)publishedFile.Get("crossSection_statistic")) {
            publishedStat = (TH1D *)h->Clone("reference");
            publishedStat->SetDirectory(0);
         }
         if (auto *h = (TH1D *)publishedFile.Get("crossSection_systematic")) {
            publishedSyst = (TH1D *)h->Clone("reference_syst");
            publishedSyst->SetDirectory(0);
         }
      }
   }
   const bool comparable = publishedStat != nullptr && !spec.wideBins;
   if (comparable) PrintFoldCheck(M, b, xUnfolded, publishedStat, bins);
   PrintRatioToPublished(result[0], comparable ? publishedStat : nullptr, bins, bayes);

   // 10. write everything out
   TFile out((radiusDir + "xsec_inversion_R" + jetR + outSuffix + ".root").c_str(), "RECREATE");
   result[0]->Write();
   if (result[1]) result[1]->Write();
   if (resultPubErr) resultPubErr->Write();
   if (resultEmbStat[0]) resultEmbStat[0]->Write();
   if (resultEmbStat[1]) resultEmbStat[1]->Write();
   if (publishedStat) publishedStat->Write("reference");
   if (publishedSyst) publishedSyst->Write("reference_syst");
   M.Write("M");
   b.Write("b_data");
   x.Write("x_emb");
   bEmb.Write("b_emb");
   cov.Write("cov"); // statistical covariance of the unfolded spectrum (both blocks)
   out.Close();
   printf("[unfold] wrote %sxsec_inversion_R%s%s.root\n", radiusDir.c_str(), jetR, outSuffix.c_str());
}
