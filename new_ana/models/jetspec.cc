// jetspec.cc — the particle-level inclusive jet spectra of one generator, from its pT-hat-binned
// particle trees, with the Stage-1 particle-level jet definition (src/ppAnalysis.cxx at intype MCPICO):
// constituents |eta| < 2.5 without a pT cut; anti-kT with an active area from explicit ghosts; jets
// |eta| < 0.5 above the truth-pass floor of 1.5 GeV; the underlying event subtracted with the off-axis
// cones (rho = mean pT density of the two cones of radius R at +-pi/2 in phi, pt_corr = pt - rho x area).
//
// Every pt<lo>_<hi>.root of <tree-dir> is one pT-hat bin, weighted sigma_pb / N_generated; for pythia6*
// the three softest bins carry the cross-section corrections of the published analysis. Each sample then
// passes the published outlier filter in its one-dimensional form: a 0.1 GeV cell is dropped when the
// 20-cell window around it holds no more than 4 unweighted entries of that sample (a single high-weight
// soft event that made a hard jet).
//
// Input : <tree-dir>/pt<lo>_<hi>.root, the particle trees of tree.h.
// Output: <out.root>, weighted counts in 0.1 GeV bins per radius (compare.C rebins and divides by the bin
//         width; |eta| < 0.5 gives d eta = 1):
//           ue_R<R>, noue_R<R>            pt_corr / raw pt, weight sigma / N   [+ _entries: unit weights]
//           ue_R<R>_soft, noue_R<R>_soft  the same times the soft reweight of the published analysis
//                                         (soft_reweight.h; the Pythia 6 spectra are quoted with it)
//           jetpt_vs_pthat_R0.5           raw jet pT against the hard-process pT, before the filter
// Usage : jetspec <generator> <tree-dir> <out.root> [radii, default "0.2 0.3 0.4 0.5"]
#include "tree.h"
#include "../soft_reweight.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TVector2.h>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/Selector.hh>
#include <cstdio>
#include <glob.h>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Cross-section corrections of the three softest pT-hat bins, as in the published analysis
// (stage2/build_resp.C); they apply to the Pythia 6 samples only.
static const std::map<std::string, double> kSigmaCorr = {
   {"pt2_3", 1. / 1.228}, {"pt3_4", 1. / 1.051}, {"pt4_5", 1. / 1.014}};

// The spectra are filled in fine 0.1 GeV cells and rebinned downstream.
static const int kNFine = 1000;
static const double kFineMax = 100.0;

// The particle-level jet definition of the Stage-1 truth pass.
static const double kEtaConstituent = 2.5; // constituents kept, no pT cut
static const double kEtaJet = 0.5;         // |eta| acceptance of the published measurement
static const double kPtJetMin = 1.5;       // raw jet floor of the truth pass (container.sh PJMIN=1.5)
static const double kGhostArea = 0.04;     // active-area ghosts
static const int kGhostRepeat = 1;

// The published outlier filter: the window [i-9, i+10] of 20 cells around cell i must hold more than
// four unweighted entries, otherwise the cell is a single reweighted soft event and is dropped.
static const int kFilterCellsBelow = 9;
static const int kFilterCellsAbove = 11;
static const double kFilterMinEntries = 4.0;

// The two jet definitions written side by side, and the three copies of each.
enum JetDefinitionIndex { kUeSubtracted = 0, kRawJet = 1, kNJetDefs = 2 };
enum SpectrumVariant { kWeighted = 0, kWeightedSoft = 1, kEntries = 2, kNVariants = 3 };

// Position of one histogram in the flat [radius][jet definition][variant] vectors.
static size_t HistIndex(size_t radius, int jetDef, int variant)
{
   return (radius * kNJetDefs + jetDef) * kNVariants + variant;
}

// Mean pT density of the two cones of radius R placed at +-pi/2 in phi from the jet: the underlying-event
// estimate of the published analysis.
static double OffAxisRho(const fastjet::PseudoJet &jet, const std::vector<fastjet::PseudoJet> &particles, double R)
{
   const double phiPlus = jet.phi() + M_PI / 2;
   const double phiMinus = jet.phi() - M_PI / 2;
   const double R2 = R * R;
   double sumPlus = 0, sumMinus = 0;
   for (const auto &p : particles) {
      const double dEta2 = (p.eta() - jet.eta()) * (p.eta() - jet.eta());
      const double dPhiPlus = TVector2::Phi_mpi_pi(p.phi() - phiPlus);
      const double dPhiMinus = TVector2::Phi_mpi_pi(p.phi() - phiMinus);
      if (dEta2 + dPhiPlus * dPhiPlus < R2) sumPlus += p.perp();
      if (dEta2 + dPhiMinus * dPhiMinus < R2) sumMinus += p.perp();
   }
   return 0.5 * (sumPlus + sumMinus) / (M_PI * R2);
}

// Zero every cell whose 20-cell neighbourhood holds too few unweighted entries, in the entry histogram
// and in the weighted ones that go with it; accumulate the weight thrown away.
static void ApplyOutlierFilter(TH1D *entries, std::vector<TH1D *> weighted, double &removed)
{
   std::vector<double> cumulative(kNFine + 1, 0);
   for (int i = 0; i < kNFine; ++i) cumulative[i + 1] = cumulative[i] + entries->GetBinContent(i + 1);
   for (int i = 0; i < kNFine; ++i) {
      const int lo = std::max(0, i - kFilterCellsBelow);
      const int hi = std::min(kNFine, i + kFilterCellsAbove);
      if (cumulative[hi] - cumulative[lo] > kFilterMinEntries) continue;
      for (TH1D *h : weighted) {
         removed += h->GetBinContent(i + 1);
         h->SetBinContent(i + 1, 0);
         h->SetBinError(i + 1, 0);
      }
      entries->SetBinContent(i + 1, 0);
      entries->SetBinError(i + 1, 0);
   }
}

// "0.2 0.3 0.4 0.5" -> the list of jet radii.
static std::vector<double> ParseRadii(const char *text)
{
   std::vector<double> radii;
   std::istringstream stream(text);
   double R;
   while (stream >> R) radii.push_back(R);
   return radii;
}

// Book the running totals and the per-sample working histograms, one set per radius and jet definition.
static void BookSpectra(const std::vector<double> &radii, std::vector<TH1D *> &total, std::vector<TH1D *> &sample)
{
   const char *defName[kNJetDefs] = {"ue", "noue"};
   const char *variantSuffix[kNVariants] = {"", "_soft", "_entries"};
   for (size_t k = 0; k < radii.size(); ++k) {
      for (int d = 0; d < kNJetDefs; ++d) {
         const std::string base = Form("%s_R%.1f", defName[d], radii[k]);
         for (int v = 0; v < kNVariants; ++v) {
            const size_t i = HistIndex(k, d, v);
            total[i] = new TH1D((base + variantSuffix[v]).c_str(), "", kNFine, 0, kFineMax);
            sample[i] = new TH1D((base + variantSuffix[v] + "_sample").c_str(), "", kNFine, 0, kFineMax);
            total[i]->Sumw2();
            sample[i]->Sumw2();
         }
      }
   }
}

// The final state of one event, as fastjet input.
static void ReadParticles(const Models::Event &event, std::vector<fastjet::PseudoJet> &particles)
{
   particles.clear();
   for (int i = 0; i < event.n; ++i) {
      fastjet::PseudoJet p(event.px[i], event.py[i], event.pz[i], event.e[i]);
      if (std::fabs(p.eta()) <= kEtaConstituent) particles.push_back(p);
   }
}

// Cluster one event at every radius and fill the per-sample spectra, with and without the UE subtraction.
static void FillEvent(const std::vector<fastjet::PseudoJet> &particles, const std::vector<double> &radii,
                      double weight, double weightSoft, std::vector<TH1D *> &sample, TH2D &jetPtVsPtHat,
                      double ptHat)
{
   for (size_t k = 0; k < radii.size(); ++k) {
      const double R = radii[k];
      fastjet::JetDefinition definition(fastjet::antikt_algorithm, R);
      fastjet::AreaDefinition area(fastjet::active_area_explicit_ghosts,
                                   fastjet::GhostedAreaSpec(1.0 + 2.0 * R, kGhostRepeat, kGhostArea));
      fastjet::ClusterSequenceArea cluster(particles, definition, area);
      for (const auto &jet : fastjet::SelectorAbsEtaMax(kEtaJet)(cluster.inclusive_jets(kPtJetMin))) {
         const double pt = jet.perp();
         const double ptCorrected = pt - OffAxisRho(jet, particles, R) * jet.area();
         const double jetPt[kNJetDefs] = {ptCorrected, pt};
         for (int d = 0; d < kNJetDefs; ++d) {
            sample[HistIndex(k, d, kWeighted)]->Fill(jetPt[d], weight);
            sample[HistIndex(k, d, kWeightedSoft)]->Fill(jetPt[d], weightSoft);
            sample[HistIndex(k, d, kEntries)]->Fill(jetPt[d]);
         }
         if (R == 0.5) jetPtVsPtHat.Fill(pt, ptHat, weight);
      }
   }
}

int main(int argc, char **argv)
{
   if (argc < 4 || argc > 5) {
      fprintf(stderr, "usage: jetspec <generator> <tree-dir> <out.root> [radii]\n");
      return 1;
   }
   const std::string generator = argv[1];
   const std::string treeDir = argv[2];
   const std::vector<double> radii = ParseRadii(argc > 4 ? argv[4] : "0.2 0.3 0.4 0.5");

   // 1. book the spectra: running totals and the working set of the sample being read
   const size_t nHist = radii.size() * kNJetDefs * kNVariants;
   std::vector<TH1D *> total(nHist), sample(nHist);
   BookSpectra(radii, total, sample);
   TH2D jetPtVsPtHat("jetpt_vs_pthat_R0.5", ";jet p_{T} [GeV];p-hat [GeV]", 200, 0, 100, 140, 0, 70);
   jetPtVsPtHat.Sumw2();

   glob_t files;
   if (glob((treeDir + "/pt*_*.root").c_str(), 0, nullptr, &files) != 0 || files.gl_pathc == 0) {
      fprintf(stderr, "jetspec: no pt*_*.root in %s\n", treeDir.c_str());
      return 2;
   }

   long long nEvents = 0;
   double removedAll = 0, filledAll = 0;
   // 2. one pT-hat bin at a time: cluster its events, filter it, add it to the totals
   for (size_t f = 0; f < files.gl_pathc; ++f) {
      const std::string file = files.gl_pathv[f];
      const std::string base = file.substr(file.rfind('/') + 1);
      const std::string sampleName = base.substr(0, base.size() - 5); // strip ".root"
      Models::Reader reader(file);
      const double correction =
         (generator.rfind("pythia6", 0) == 0 && kSigmaCorr.count(sampleName)) ? kSigmaCorr.at(sampleName) : 1.0;
      const double weight = reader.nev > 0 ? reader.sigma_pb * correction / (double)reader.nev : 0.0;
      for (TH1D *h : sample) h->Reset();

      std::vector<fastjet::PseudoJet> particles;
      while (reader.Next()) {
         ++nEvents;
         ReadParticles(reader.ev, particles);
         const double weightSoft = weight * SoftReweight::weight(reader.ev.pthat);
         FillEvent(particles, radii, weight, weightSoft, sample, jetPtVsPtHat, reader.ev.pthat);
      }

      double removed = 0, filled = 0;
      for (size_t k = 0; k < radii.size(); ++k) {
         for (int d = 0; d < kNJetDefs; ++d) {
            filled += sample[HistIndex(k, d, kWeighted)]->Integral();
            ApplyOutlierFilter(sample[HistIndex(k, d, kEntries)],
                               {sample[HistIndex(k, d, kWeighted)], sample[HistIndex(k, d, kWeightedSoft)]}, removed);
            for (int v = 0; v < kNVariants; ++v) total[HistIndex(k, d, v)]->Add(sample[HistIndex(k, d, v)]);
         }
      }
      removedAll += removed;
      filledAll += filled;
      printf("jetspec: %-10s %9lld events  sigma %.4g pb  corr %.3f  weight %.4g pb/event  "
             "filter removed %.2g%% of the weight\n",
             sampleName.c_str(), reader.nev, reader.sigma_pb, correction, weight,
             filled > 0 ? 100 * removed / filled : 0);
   }

   // 3. write the spectra of the generator
   TFile out(argv[3], "RECREATE");
   for (TH1D *h : total) h->Write();
   jetPtVsPtHat.Write();
   TNamed("summary", Form("%s: %zu pT-hat bins, %lld events, outlier filter removed %.3g%% of the weighted jets",
                          generator.c_str(), files.gl_pathc, nEvents, filledAll > 0 ? 100 * removedAll / filledAll : 0))
      .Write();
   out.Close();
   globfree(&files);
   printf("jetspec: %s -> %s (%lld events, filter removed %.3g%% of the weight)\n", generator.c_str(), argv[3],
          nEvents, filledAll > 0 ? 100 * removedAll / filledAll : 0);
   return 0;
}
