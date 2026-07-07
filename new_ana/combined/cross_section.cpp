// Promotion-combined (JPX) cross section — Dmitry's published nominal method.
//
// DATA: one pass over the Stage-1 data merge partitions events into the
// exclusive promotion categories (hardware shouldFire proxy = OR of the
// per-jet hardware trigger_match bits), applies the per-category detector
// windows + recorded-accept gates, divides each category by its measured
// trigger correction C(pt) (config.h::TrigEffMeas, the hardware->simulator
// ruler bridge), and sums — RAW counts, full JP2 luminosity.
//
// RESPONSE: reads the fine ingredients written by response.cxx, removes
// statistically unreliable cells (raw box-entry sum <= 4 in a +-1 GeV box,
// Dmitry's filtered()), routes the removed content consistently out of b and
// x, re-coarsens to McBins, and builds
//   M_ij = (b_i / matched_i) * A_ij / x_j
// (background/fake boost x migration x matching+trigger efficiency).
//
// SOLVE: floor-restricted square block [kJpxFloor, 52), unregularized
// inversion x = M^-1 b (optional scale-matched second-difference Tikhonov,
// promotion.h::kTikhonovLambda). Covariance X = R diag(err_b^2) R^T.
//
// Outputs (config.h kWorkDir):
//   xsec_JPX_R<R>.root                    hist "canonical" + "reference"
//   comparison_with_dmitriy_R<R>_JPX.pdf

#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TPad.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorD.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "promotion.h"

using namespace CrossSectionConfig;
using namespace PromotionConfig;
AnalysisConfig cfg;

static void AppendBadRunsFromFile(const std::string &path)
{
   std::ifstream fp(path);
   if (!fp) {
      std::cerr << "[badruns] file missing: " << path << " — skipped\n";
      return;
   }
   int run, n = 0;
   while (fp >> run) {
      cfg.badRuns.push_back(run);
      ++n;
   }
   std::cout << "[badruns] appended " << n << " runs from " << path << std::endl;
}

// ---------------------------------------------------------------------------
// Dmitry's box sums (np.pad + cumsum): out[i] = sum of v over the asymmetric
// window [i-res, i+res-1] (width 2*res).
static std::vector<double> boxconv1d(const std::vector<double> &v, int res)
{
   const int n = (int)v.size();
   const int off = 2 * res;
   std::vector<double> padded(n + 2 * res, 0.0);
   for (int i = 0; i < n; ++i)
      padded[i + res] = v[i];
   std::vector<double> S(padded.size() + 1, 0.0);
   for (size_t k = 0; k < padded.size(); ++k)
      S[k + 1] = S[k] + padded[k];
   std::vector<double> out(n, 0.0);
   for (int i = 0; i < n; ++i)
      out[i] = S[i + off] - S[i];
   return out;
}

static std::vector<std::vector<double>> boxconv2d(const std::vector<std::vector<double>> &M, int res)
{
   const int nr = (int)M.size();
   if (nr == 0) return {};
   const int nc = (int)M[0].size();
   std::vector<std::vector<double>> tmp(nr, std::vector<double>(nc, 0.0));
   for (int j = 0; j < nc; ++j) {
      std::vector<double> col(nr);
      for (int i = 0; i < nr; ++i) col[i] = M[i][j];
      std::vector<double> bc = boxconv1d(col, res);
      for (int i = 0; i < nr; ++i) tmp[i][j] = bc[i];
   }
   std::vector<std::vector<double>> out(nr, std::vector<double>(nc, 0.0));
   for (int i = 0; i < nr; ++i)
      out[i] = boxconv1d(tmp[i], res);
   return out;
}

// ---------------------------------------------------------------------------
// DATA: category-partitioned, C(pt)-corrected, summed raw spectrum on McBins.
// The three PRE-correction category histograms are cached on disk
// (data_JPX_R<R>_categories.root): the 16 GB tree read happens once, and
// floor/lambda/trigEff studies reuse the cache. Delete the file (or run after
// a selection change — the driver does not do this for you) to force a
// re-read; the C(pt) division and any trigEffScale are applied AFTER loading,
// so every variation is served correctly from the same cache.
static std::unique_ptr<TH1D> RawCombined(const std::string &jetR, const Systematic &syst)
{
   const std::vector<double> &bins = McBins();
   const std::string cachePath =
      Form("%sdata_JPX_R%s_categories.root", cfg.workdir.c_str(), jetR.c_str());

   // The UE-fraction variation changes every jet's pT (windows + spectrum), so
   // it cannot be served from the nominal cache: bypass (and do not overwrite).
   const bool ueVar = std::abs(syst.ueFraction - 1.0) > 1e-9;

   TH1D *cat[3] = {nullptr, nullptr, nullptr};
   if (!ueVar && !gSystem->AccessPathName(cachePath.c_str())) {
      TFile fc(cachePath.c_str(), "READ");
      for (int c = 0; c < 3; ++c) {
         auto *h = (TH1D *)fc.Get(Form("jpx_cat%d", c));
         if (h) {
            cat[c] = (TH1D *)h->Clone(Form("jpx_cat%d_c", c));
            cat[c]->SetDirectory(0);
         }
      }
      if (cat[0] && cat[1] && cat[2])
         std::cout << "[jpx][data] using cached category histograms " << cachePath
                   << " (delete after any selection change!)" << std::endl;
   }

   if (!cat[0] || !cat[1] || !cat[2]) {
      const std::string dataFile = Form("%smerged_data_R%s.root", cfg.datapath.c_str(), jetR.c_str());
      ROOT::RDataFrame df("ResultTree", dataFile.c_str());

      auto anyOf = [](const ROOT::VecOps::RVec<bool> &v) {
         for (auto b : v)
            if (b) return true;
         return false;
      };

      auto dfn = df.Filter(
                      [](int runIndex) {
                         return runIndex >= 0 && std::find(cfg.badRuns.begin(), cfg.badRuns.end(),
                                                           runIndex) == cfg.badRuns.end();
                      },
                      {"runid1"})
                    .Define("evt_sf0", anyOf, {"trigger_match_JP0"})
                    .Define("evt_sf1", anyOf, {"trigger_match_JP1"})
                    .Define("evt_sf2", anyOf, {"trigger_match_JP2"})
                    .Define("rcat_evt", "evt_sf2 ? 2 : (evt_sf1 ? 1 : (evt_sf0 ? 0 : -1))")
                    // UE-fraction variation: pt_corrected = pt_raw - area*rho, so
                    // the varied pT is pt_corrected + (1-f)*area*rho.
                    .Define("ptv", ueVar ? std::string(Form("pt_corrected + (%.6f)*jet_area*bg_density",
                                                            1.0 - syst.ueFraction))
                                         : std::string("pt_corrected"));
      if (ueVar)
         std::cout << "[jpx][data] UE-fraction variation: detector fraction " << syst.ueFraction
                   << " (cache bypassed)" << std::endl;

      // Per-category jet masks: the SAME CatJetGate as the response reco side
      // (promotion.h — one definition); data additionally requires the
      // recorded hardware accept.
      const std::string base = "abs(det_eta) < 0.5 && neutral_fraction <= 0.95 && ";
      auto d2 = dfn.Define("sel2", base + "rcat_evt == 2 && " + CatJetGate(2, "", "ptv") +
                                      " && (fired_JP0 || fired_JP1 || fired_JP2)")
                   .Define("pt2", "ptv[sel2]");
      auto d1 = d2.Define("sel1", base + "rcat_evt == 1 && " + CatJetGate(1, "", "ptv") +
                                     " && (fired_JP0 || fired_JP1)")
                   .Define("pt1", "ptv[sel1]");
      auto d0 = d1.Define("sel0", base + "rcat_evt == 0 && " + CatJetGate(0, "", "ptv") +
                                     " && fired_JP0")
                   .Define("pt0", "ptv[sel0]");

      auto h2 = d0.Histo1D({"jpx_cat2", "", (int)bins.size() - 1, bins.data()}, "pt2");
      auto h1 = d0.Histo1D({"jpx_cat1", "", (int)bins.size() - 1, bins.data()}, "pt1");
      auto h0 = d0.Histo1D({"jpx_cat0", "", (int)bins.size() - 1, bins.data()}, "pt0");

      if (!ueVar) {
         TFile fc(cachePath.c_str(), "RECREATE");
         h0->Write("jpx_cat0");
         h1->Write("jpx_cat1");
         h2->Write("jpx_cat2");
         fc.Close();
         std::cout << "[jpx][data] category cache written: " << cachePath << std::endl;
      }

      cat[0] = (TH1D *)h0->Clone("jpx_cat0_c");
      cat[1] = (TH1D *)h1->Clone("jpx_cat1_c");
      cat[2] = (TH1D *)h2->Clone("jpx_cat2_c");
      for (int c = 0; c < 3; ++c) cat[c]->SetDirectory(0);
   }

   std::cout << "[jpx][data] category jets (raw): cat0=" << cat[0]->Integral()
             << "  cat1=" << cat[1]->Integral() << "  cat2=" << cat[2]->Integral() << std::endl;

   TH1D *sum = (TH1D *)cat[2]->Clone(Form("JPX_raw_R%s", jetR.c_str()));
   sum->SetDirectory(0);
   sum->Add(cat[1]);
   sum->Add(cat[0]);
   for (int c = 0; c < 3; ++c) delete cat[c];

   // Measured combination-level trigger correction: divide the summed
   // hardware-gated data by C_JPX(pt) = That_sum (turn-on) then R (plateau
   // ruler) — promotion.h::JpxTrigEff, measured for the exact promotion gates.
   std::cout << "[jpx][data] dividing by C_JPX(pt) (promotion.h JpxTrigEff), trigEffScale "
             << syst.trigEffScale << std::endl;
   for (int i = 1; i <= sum->GetNbinsX(); ++i) {
      const double p = JpxTrigEff(sum->GetXaxis()->GetBinCenter(i)) * syst.trigEffScale;
      if (p > 0) {
         sum->SetBinContent(i, sum->GetBinContent(i) / p);
         sum->SetBinError(i, sum->GetBinError(i) / p);
      }
   }
   return std::unique_ptr<TH1D>(sum);
}

// ---------------------------------------------------------------------------
// RESPONSE: fine file -> cell filter -> b/x subtractions -> coarse Ac, bc, xc
// (+ their propagated variances, for the embedding-statistics toys).
struct CoarseResponse {
   std::vector<std::vector<double>> A, Aerr2; // [reco][mc]
   std::vector<double> b, x, berr2, xerr2;
};

static CoarseResponse FilterAndCoarsen(const std::string &jetR, const Systematic &syst)
{
   // A shape systematic (jesShift/jerSmear) reads its own rebuilt fine file;
   // everything else reuses the nominal ingredients.
   const std::string respTag = syst.needsResponse() ? SystTag(syst) : std::string("");
   const TString fineName =
      Form("%sresponse_JPX_R%s_fine%s.root", cfg.workdir.c_str(), jetR.c_str(), respTag.c_str());
   TFile fin(fineName, "READ");
   if (fin.IsZombie()) {
      std::cerr << "[FATAL] " << fineName << " missing — run response.cxx first" << std::endl;
      gSystem->Exit(2);
   }
   auto *A_w = (TH2D *)fin.Get("A_fine");
   auto *A_e = (TH2D *)fin.Get("A_entries_fine");
   auto *A_x = (TH2D *)fin.Get("A_xfine");
   auto *b_w = (TH1D *)fin.Get("b_fine");
   auto *b_e = (TH1D *)fin.Get("b_entries_fine");
   auto *x_w = (TH1D *)fin.Get("x_fine");
   auto *x_e = (TH1D *)fin.Get("x_entries_fine");
   if (!A_w || !A_e || !A_x || !b_w || !b_e || !x_w || !x_e) {
      std::cerr << "[FATAL] fine ingredients incomplete in " << fineName << std::endl;
      gSystem->Exit(2);
   }

   const std::vector<double> &grid = McBins();
   const int nb = (int)grid.size() - 1;
   const double winLo = grid.front(), winHi = grid.back();
   auto mcInWindow = [&](int j) {
      const double c = kMcFineLo + (j + 0.5) * kFineW;
      return c >= winLo && c < winHi;
   };
   auto recoInWindow = [&](int i) {
      const double c = kRecoFineLo + (i + 0.5) * kFineW;
      return c >= winLo && c < winHi;
   };

   // Local copies (Aw is the filtered matrix used downstream; the ORIGINAL
   // TH2s are kept for the b/x subtraction projections).
   std::vector<std::vector<double>> Aw(kNRecoFine, std::vector<double>(kNMcFine, 0.0));
   std::vector<std::vector<double>> Ae(kNRecoFine, std::vector<double>(kNMcFine, 0.0));
   for (int i = 0; i < kNRecoFine; ++i)
      for (int j = 0; j < kNMcFine; ++j) {
         Aw[i][j] = A_w->GetBinContent(i + 1, j + 1);
         Ae[i][j] = A_e->GetBinContent(i + 1, j + 1);
      }

   // Pre-filter matched entry projections (for the term-2 "unmatched" split).
   std::vector<double> matched_b_e(kNRecoFine, 0.0), matched_x_e(kNMcFine, 0.0);
   for (int i = 0; i < kNRecoFine; ++i)
      for (int j = 0; j < kNMcFine; ++j) {
         if (mcInWindow(j)) matched_b_e[i] += Ae[i][j];
         if (recoInWindow(i)) matched_x_e[j] += Ae[i][j];
      }

   // --- Box filter: zero cells whose raw-entry box sum is <= kBoxThresh ----
   std::vector<std::vector<double>> boxed = boxconv2d(Ae, kBoxRes);
   long long nOut = 0;
   double wRemoved = 0.0;
   std::vector<std::vector<bool>> isOut(kNRecoFine, std::vector<bool>(kNMcFine, false));
   for (int i = 0; i < kNRecoFine; ++i)
      for (int j = 0; j < kNMcFine; ++j)
         if (boxed[i][j] <= (double)kBoxThresh) {
            isOut[i][j] = true;
            if (Ae[i][j] > 0) {
               ++nOut;
               wRemoved += Aw[i][j];
            }
            Aw[i][j] = 0.0;
         }
   std::cout << "[jpx][filter] outlier cells zeroed: " << nOut << " (weight removed " << wRemoved << ")"
             << std::endl;

   // --- b / x subtractions (both terms of Dmitry's filtered()) -------------
   std::vector<double> bw(kNRecoFine, 0.0), xw(kNMcFine, 0.0);
   for (int i = 0; i < kNRecoFine; ++i) bw[i] = b_w->GetBinContent(i + 1);
   for (int j = 0; j < kNMcFine; ++j) xw[j] = x_w->GetBinContent(j + 1);

   // term 1: removed-A projections, ORIGINAL contents. The b (reco) side uses
   // the meas_w-weighted A; the x (truth) side uses the total_weight-weighted
   // A_xfine so the truth normalization stays exact.
   for (int i = 0; i < kNRecoFine; ++i)
      for (int j = 0; j < kNMcFine; ++j) {
         if (!isOut[i][j]) continue;
         if (mcInWindow(j)) bw[i] -= A_w->GetBinContent(i + 1, j + 1);
         if (recoInWindow(i)) xw[j] -= A_x->GetBinContent(i + 1, j + 1);
      }

   // term 2: sparse UNMATCHED fine bins (fake/miss tails), 1D box on entries.
   std::vector<double> unmatched_b_e(kNRecoFine, 0.0), unmatched_x_e(kNMcFine, 0.0);
   for (int i = 0; i < kNRecoFine; ++i) unmatched_b_e[i] = b_e->GetBinContent(i + 1) - matched_b_e[i];
   for (int j = 0; j < kNMcFine; ++j) unmatched_x_e[j] = x_e->GetBinContent(j + 1) - matched_x_e[j];
   std::vector<double> b_box = boxconv1d(unmatched_b_e, kBoxRes);
   std::vector<double> x_box = boxconv1d(unmatched_x_e, kBoxRes);
   for (int i = 0; i < kNRecoFine; ++i)
      if (b_box[i] <= (double)kBoxThresh) {
         double matched_b_w = 0.0;
         for (int j = 0; j < kNMcFine; ++j)
            if (mcInWindow(j)) matched_b_w += A_w->GetBinContent(i + 1, j + 1);
         bw[i] -= (b_w->GetBinContent(i + 1) - matched_b_w);
      }
   for (int j = 0; j < kNMcFine; ++j)
      if (x_box[j] <= (double)kBoxThresh) {
         double matched_x_w = 0.0;
         for (int i = 0; i < kNRecoFine; ++i)
            if (recoInWindow(i)) matched_x_w += A_x->GetBinContent(i + 1, j + 1);
         xw[j] -= (x_w->GetBinContent(j + 1) - matched_x_w);
      }
   for (auto &v : bw)
      if (v < 0) v = 0.0;
   for (auto &v : xw)
      if (v < 0) v = 0.0;

   // --- Re-coarsen to the square analysis grid ------------------------------
   auto coarseIdx = [&](double center) -> int {
      for (int k = 0; k < nb; ++k)
         if (center >= grid[k] && center < grid[k + 1]) return k;
      return -1;
   };
   CoarseResponse R;
   R.A.assign(nb, std::vector<double>(nb, 0.0));
   R.Aerr2.assign(nb, std::vector<double>(nb, 0.0));
   R.b.assign(nb, 0.0);
   R.x.assign(nb, 0.0);
   R.berr2.assign(nb, 0.0);
   R.xerr2.assign(nb, 0.0);
   for (int i = 0; i < kNRecoFine; ++i) {
      const double rc = kRecoFineLo + (i + 0.5) * kFineW;
      const int ri = coarseIdx(rc);
      if (ri < 0) continue;
      R.b[ri] += bw[i];
      R.berr2[ri] += b_w->GetBinError(i + 1) * b_w->GetBinError(i + 1);
      for (int j = 0; j < kNMcFine; ++j) {
         if (Aw[i][j] == 0.0) continue;
         const double mc = kMcFineLo + (j + 0.5) * kFineW;
         const int mi = coarseIdx(mc);
         if (mi < 0) continue;
         R.A[ri][mi] += Aw[i][j];
         R.Aerr2[ri][mi] += A_w->GetBinError(i + 1, j + 1) * A_w->GetBinError(i + 1, j + 1);
      }
   }
   for (int j = 0; j < kNMcFine; ++j) {
      const double mc = kMcFineLo + (j + 0.5) * kFineW;
      const int mi = coarseIdx(mc);
      if (mi < 0) continue;
      R.x[mi] += xw[j];
      R.xerr2[mi] += x_w->GetBinError(j + 1) * x_w->GetBinError(j + 1);
   }
   return R;
}

// ---------------------------------------------------------------------------
// systName selects a preset from config.h::Systematics(); default "nominal" is
// the physics result. A variation writes xsec_JPX_R<R>_<name>.root; the Dmitry
// comparison PDF is drawn only for the nominal.
//
// lambdaOverride / floorOverride (explicit STUDY arguments, not flags): pass
// >= 0 to try a different Tikhonov damping or block floor without editing
// promotion.h. Study results overwrite the same output files — rerun the
// defaults afterwards. The defaults (<0) use promotion.h.
void cross_section(const char *systName = "nominal", double lambdaOverride = -1.0,
                   double floorOverride = -1.0)
{
   ROOT::EnableImplicitMT(kImtThreads);
   DefineCustomColors();
   gStyle->SetOptStat(0);
   TH1::SetDefaultSumw2();

   const Systematic syst = FindSystematic(systName);
   std::cout << "[jpx] systematic = " << syst.name << std::endl;

   // Damping precedence: explicit study argument > the variation's jpxLambda
   // (the "jpxDamp" unfolding systematic) > the promotion.h default.
   const double tikLambda =
      (lambdaOverride >= 0.0) ? lambdaOverride : (syst.jpxLambda >= 0.0 ? syst.jpxLambda : kTikhonovLambda);
   const double jpxFloor = (floorOverride >= 0.0) ? floorOverride : kJpxFloor;

   AppendBadRunsFromFile(cfg.workdir + "../lists/dmitry_extras.list");

   for (const auto &jetR : cfg.jetRs) {
      const std::string dataPath = Form("%smerged_data_R%s.root", cfg.datapath.c_str(), jetR.c_str());
      if (gSystem->AccessPathName(dataPath.c_str())) {
         std::cerr << "[warn] " << dataPath << " missing; skipping JPX R=" << jetR << std::endl;
         continue;
      }

      // ---- inputs ----------------------------------------------------------
      std::unique_ptr<TH1D> hData = RawCombined(jetR, syst);
      CoarseResponse cr = FilterAndCoarsen(jetR, syst);

      const std::vector<double> &grid = McBins();
      const int nb = (int)grid.size() - 1;

      // ---- floor-restricted square block ------------------------------------
      int i0 = 0;
      for (int k = 0; k < nb; ++k)
         if (grid[k] >= jpxFloor - 1e-6) {
            i0 = k;
            break;
         }
      int i1 = nb - 1;
      for (int k = nb - 1; k >= 0; --k)
         if (grid[k] < kQuoteHi - 1e-6) {
            i1 = k;
            break;
         }
      const int nfloor = i1 - i0 + 1;
      std::cout << "[jpx][solve] block bins [" << grid[i0] << "," << grid[i1 + 1] << ") — " << nfloor
                << "x" << nfloor << ", lambda=" << tikLambda << std::endl;

      // M_ij = (b_i / matched_i) * A_ij / x_j, matched_i = in-block row sum
      // (matched content with truth OUTSIDE the block — below-floor feed-up —
      // is treated as background by the b/matched boost).
      auto buildM = [&](const CoarseResponse &c) -> TMatrixD {
         TMatrixD Mm(nfloor, nfloor);
         for (int a = 0; a < nfloor; ++a) {
            const int i = i0 + a;
            const double bi = c.b[i];
            double mi = 0.0;
            for (int k = 0; k < nfloor; ++k) mi += c.A[i][i0 + k];
            const double rowScale = (mi > 0) ? bi / mi : 0.0;
            for (int k = 0; k < nfloor; ++k) {
               const int j = i0 + k;
               const double xj = c.x[j];
               Mm(a, k) = (xj > 0) ? rowScale * c.A[i][j] / xj : 0.0;
            }
         }
         return Mm;
      };
      TMatrixD M = buildM(cr);

      TVectorD bdata(nfloor), bdataErr2(nfloor);
      for (int a = 0; a < nfloor; ++a) {
         bdata[a] = hData->GetBinContent(i0 + a + 1);
         bdataErr2[a] = hData->GetBinError(i0 + a + 1) * hData->GetBinError(i0 + a + 1);
      }

      // QA: fold the REFERENCE truth through M and compare with the data
      // vector, row by row — localizes any residual as data-side (shows up
      // here) vs solve-side (does not). Approximate in the buffer row (the
      // reference stops at 52).
      {
         const std::string refPath =
            Form((cfg.workdir + "jet_cross_section_dmitriyR%s.root").c_str(), jetR.c_str());
         if (!gSystem->AccessPathName(refPath.c_str())) {
            TFile rf(refPath.c_str(), "READ");
            auto *r = (TH1D *)rf.Get("crossSection_systematic");
            if (r) {
               const double Lfold = RuntimeLeff(cfg, "JP2");
               TVectorD xref(nfloor);
               for (int a = 0; a < nfloor; ++a) {
                  const int j = i0 + a;
                  const double lo = grid[j], hi = grid[j + 1];
                  const int rb = r->FindBin(0.5 * (lo + hi));
                  const double v = r->GetBinContent(rb); // dsigma/dpt/deta
                  // back to raw counts: x_j = v * dpt * (2*deta) * L
                  xref[a] = v * (hi - lo) * 2.0 * (1.0 - std::stod(jetR)) * Lfold;
               }
               TVectorD bpred = M * xref;
               std::cout << "[jpx][fold-QA] b_data / (M x Dmitry-truth)  per reco bin:\n";
               for (int a = 0; a < nfloor; ++a)
                  printf("    [%5.1f,%5.1f)  data=%10.4g  pred=%10.4g  ratio=%6.3f\n", grid[i0 + a],
                         grid[i0 + a + 1], bdata[a], bpred[a], bpred[a] > 0 ? bdata[a] / bpred[a] : 0.0);
            }
         }
      }

      // ---- solve (shared by the nominal and the embedding-stat toys) ---------
      auto solveR = [&](const TMatrixD &Mm, const CoarseResponse &c) -> TMatrixD {
         TMatrixD Rr(nfloor, nfloor);
         if (tikLambda <= 0.0) {
            Double_t det = 0.0;
            TMatrixD Minv(Mm);
            Minv.Invert(&det);
            if (det == 0.0)
               std::cerr << "[jpx][WARN] filtered M is singular — inversion unreliable" << std::endl;
            Rr = Minv;
            return Rr;
         }
         // Scale-matched scale-invariant 2nd-difference (curvature) damping:
         // R = (M^T W M + l^2 D^-1 L^T L D^-1)^-1 M^T W, W = diag(1/b_data^2),
         // D = diag(xref) with xref = embedding truth forward-fold-matched to
         // the data scale.
         TVectorD xref(nfloor);
         for (int a = 0; a < nfloor; ++a) xref[a] = std::max(0.0, c.x[i0 + a]);
         TVectorD Mxr = Mm * xref;
         double num = 0.0, den = 0.0;
         for (int a = 0; a < nfloor; ++a) {
            num += bdata[a] * Mxr[a];
            den += Mxr[a] * Mxr[a];
         }
         const double C = (den > 0) ? num / den : 1.0;
         for (int a = 0; a < nfloor; ++a) xref[a] *= C;
         TMatrixD Dinv(nfloor, nfloor);
         Dinv.Zero();
         for (int a = 0; a < nfloor; ++a) Dinv(a, a) = (xref[a] > 0) ? 1.0 / xref[a] : 0.0;
         TMatrixD W(nfloor, nfloor);
         W.Zero();
         for (int a = 0; a < nfloor; ++a) W(a, a) = (bdata[a] > 0) ? 1.0 / (bdata[a] * bdata[a]) : 0.0;
         const int nL = (nfloor >= 3) ? nfloor - 2 : 0;
         TMatrixD L(nL > 0 ? nL : 1, nfloor);
         L.Zero();
         for (int r = 0; r < nL; ++r) {
            L(r, r) = 1.0;
            L(r, r + 1) = -2.0;
            L(r, r + 2) = 1.0;
         }
         TMatrixD Lt(TMatrixD::kTransposed, L);
         TMatrixD Mt(TMatrixD::kTransposed, Mm);
         TMatrixD MtW = Mt * W;
         TMatrixD Areg = MtW * Mm;
         TMatrixD pen = Dinv * (Lt * L) * Dinv;
         Areg += (tikLambda * tikLambda) * pen;
         Double_t det = 0.0;
         TMatrixD AregInv(Areg);
         AregInv.Invert(&det);
         if (det == 0.0)
            std::cerr << "[jpx][WARN] regularized normal matrix singular" << std::endl;
         Rr = AregInv * MtW;
         return Rr;
      };

      TMatrixD R = solveR(M, cr);
      TVectorD xunf = R * bdata;
      TMatrixD B(nfloor, nfloor);
      B.Zero();
      for (int a = 0; a < nfloor; ++a) B(a, a) = bdataErr2[a];
      TMatrixD RT(TMatrixD::kTransposed, R);
      TMatrixD X = R * B * RT;

      // ---- embedding-statistics toys (the response-statistics systematic) ----
      // Dmitry's simu_stat term: resample the response ingredients within their
      // statistical errors, re-solve, take the per-bin spread. The variation
      // writes canonical = nominal +- 1 sigma(toys).
      if (syst.embStatSign != 0) {
         const int kToys = 200;
         gRandom->SetSeed(20260707); // fixed: embStatUp/Down see identical toys
         std::vector<double> sum(nfloor, 0.0), sum2(nfloor, 0.0);
         int used = 0;
         for (int t = 0; t < kToys; ++t) {
            CoarseResponse ct = cr;
            for (int i = 0; i < nb; ++i) {
               ct.b[i] = std::max(0.0, cr.b[i] + gRandom->Gaus(0.0, std::sqrt(cr.berr2[i])));
               ct.x[i] = std::max(0.0, cr.x[i] + gRandom->Gaus(0.0, std::sqrt(cr.xerr2[i])));
               for (int j = 0; j < nb; ++j)
                  if (cr.A[i][j] > 0)
                     ct.A[i][j] =
                        std::max(0.0, cr.A[i][j] + gRandom->Gaus(0.0, std::sqrt(cr.Aerr2[i][j])));
            }
            TMatrixD Mt = buildM(ct);
            TMatrixD Rt = solveR(Mt, ct);
            TVectorD xt = Rt * bdata;
            ++used;
            for (int a = 0; a < nfloor; ++a) {
               sum[a] += xt[a];
               sum2[a] += xt[a] * xt[a];
            }
         }
         std::cout << "[jpx][embStat] " << used << " toys, per-bin sigma/nominal:";
         for (int a = 0; a < nfloor; ++a) {
            const double mean = sum[a] / used;
            const double sig = std::sqrt(std::max(0.0, sum2[a] / used - mean * mean));
            printf(" %.1f%%", xunf[a] != 0 ? 100.0 * sig / std::abs(xunf[a]) : 0.0);
            xunf[a] += syst.embStatSign * sig;
         }
         std::cout << std::endl;
      }

      TH1D *h = new TH1D(Form("JPX_unfolded_R%s", jetR.c_str()), ";p_{T} [GeV/c];", nb, grid.data());
      h->SetDirectory(0);
      for (int a = 0; a < nfloor; ++a) {
         h->SetBinContent(i0 + a + 1, xunf[a]);
         h->SetBinError(i0 + a + 1, std::sqrt(std::max(0.0, X(a, a))));
      }

      // ---- normalize ---------------------------------------------------------
      h->Scale(1.0 / (2.0 * (1.0 - std::stod(jetR)))); // eta acceptance
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
         const double w = h->GetBinWidth(i);
         h->SetBinContent(i, h->GetBinContent(i) / w);
         h->SetBinError(i, h->GetBinError(i) / w);
      }
      const double Leff = RuntimeLeff(cfg, "JP2") * syst.lumiScale; // full lumi (JP2 unprescaled)
      std::cout << "[jpx][lumi] Leff(full, JP2) = " << Leff << " pb^-1 (lumiScale " << syst.lumiScale
                << ")" << std::endl;
      h->Scale(1.0 / Leff);

      // ---- compare + write ---------------------------------------------------
      const std::string refPath =
         Form((cfg.workdir + "jet_cross_section_dmitriyR%s.root").c_str(), jetR.c_str());
      TH1D *ref = nullptr;
      if (!gSystem->AccessPathName(refPath.c_str())) {
         TFile rf(refPath.c_str(), "READ");
         auto *r = (TH1D *)rf.Get("crossSection_systematic");
         if (r) {
            ref = (TH1D *)r->Clone(Form("ref_R%s", jetR.c_str()));
            ref->SetDirectory(0);
         }
      }
      if (ref && syst.name == "nominal") {
         std::cout << "[ratio JPX/Dmitry] (solved bins)\n";
         for (int b = 1; b <= ref->GetNbinsX(); ++b) {
            const double x = ref->GetXaxis()->GetBinCenter(b);
            const double rv = ref->GetBinContent(b);
            const double mine = h->GetBinContent(h->FindBin(x));
            if (rv > 0 && mine != 0)
               printf("    pT=%5.2f  JPX=%9.3g  Dmitry=%9.3g  ratio=%5.3f\n", x, mine, rv, mine / rv);
         }

         TCanvas *c = new TCanvas("c_jpx", "", 900, 900);
         c->Divide(1, 2);
         auto *p1 = (TPad *)c->cd(1);
         p1->SetPad(0, 0.5, 1, 1);
         p1->SetBottomMargin(0);
         p1->SetLogy();
         ref->SetTitle("Promotion combination JP0+JP1+JP2 (matrix inversion);p_{T} [GeV/c];"
                       "d^{2}#sigma/dp_{T}d#eta [pb/(GeV/c)]");
         ref->SetLineColor(kViolet);
         ref->SetFillColorAlpha(kViolet, 0.20);
         ref->SetMarkerSize(0);
         ref->GetYaxis()->SetRangeUser(1.2, 1e7);
         ref->Draw("E2");
         h->SetMarkerStyle(20);
         h->SetMarkerColor(kAzure - 1);
         h->SetLineColor(kAzure - 1);
         h->Draw("E1 SAME");
         TLegend *leg = new TLegend(0.55, 0.65, 0.88, 0.88);
         leg->SetBorderSize(0);
         leg->AddEntry(ref, "Dmitry Table III", "f");
         leg->AddEntry(h, "JPX (promotion, M^{-1})", "lep");
         leg->Draw();

         auto *p2 = (TPad *)c->cd(2);
         p2->SetPad(0, 0, 1, 0.5);
         p2->SetTopMargin(0);
         p2->SetBottomMargin(0.30);
         TH1D *frame = (TH1D *)ref->Clone("jpx_ratio_frame");
         frame->Reset();
         frame->SetDirectory(0);
         frame->SetTitle(";p_{T} [GeV/c];JPX / Dmitry");
         frame->GetYaxis()->SetRangeUser(0.5, 1.5);
         frame->Draw();
         TH1D *band = (TH1D *)ref->Clone("jpx_ref_band");
         band->Reset();
         band->SetDirectory(0);
         for (int b = 1; b <= band->GetNbinsX(); ++b) {
            const double cv = ref->GetBinContent(b), ev = ref->GetBinError(b);
            band->SetBinContent(b, 1.0);
            band->SetBinError(b, cv > 0 ? ev / cv : 0.0);
         }
         band->SetFillColorAlpha(kViolet, 0.20);
         band->SetFillStyle(1001);
         band->SetLineColor(kViolet);
         band->SetMarkerSize(0);
         band->Draw("E2 SAME");
         TLine *line = new TLine(ref->GetXaxis()->GetXmin(), 1.0, ref->GetXaxis()->GetXmax(), 1.0);
         line->SetLineColor(kViolet);
         line->SetLineWidth(2);
         line->Draw("SAME");
         TH1D *r = new TH1D("jpx_ratio", "", ref->GetNbinsX(), ref->GetXaxis()->GetXbins()->GetArray());
         r->SetDirectory(0);
         for (int b = 1; b <= r->GetNbinsX(); ++b) {
            const double x = r->GetXaxis()->GetBinCenter(b);
            const double rv = ref->GetBinContent(b);
            const int bh = h->FindBin(x);
            r->SetBinContent(b, rv > 0 ? h->GetBinContent(bh) / rv : 0.0);
            r->SetBinError(b, rv > 0 ? h->GetBinError(bh) / rv : 0.0);
         }
         r->SetMarkerStyle(20);
         r->SetMarkerColor(kAzure - 1);
         r->SetLineColor(kAzure - 1);
         r->Draw("E1 SAME");
         c->SaveAs(Form("%scomparison_with_dmitriy_R%s_JPX.pdf", cfg.workdir.c_str(), jetR.c_str()));
      }

      TFile fout(Form("%sxsec_JPX_R%s%s.root", cfg.workdir.c_str(), jetR.c_str(), SystTag(syst).c_str()),
                 "RECREATE");
      h->Write("canonical");
      if (ref) ref->Write("reference");
      StampProvenance(syst);
      fout.Close();
      std::cout << "[jpx] wrote " << cfg.workdir << "xsec_JPX_R" << jetR << SystTag(syst) << ".root"
                << std::endl;
   }
}
