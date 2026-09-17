// compare.C — the unfolded cross section against the particle-level generator spectra of jetspec.cc,
// per radius and per jet definition:
//   ue    the nominal definition (underlying event subtracted at both levels) in the published bins:
//         JP1 (Bayes 4 it., 2023 embedding) with its systematic band where one exists (R = 0.5), the
//         min-bias level, the combined jet-patch levels and, at R = 0.5, the published Table III;
//   noue  raw jets at both levels in the R_AA bins 5-60 GeV: the reference set of results/Michal
//         (JP1 with its band, the combined levels, the min-bias level).
// Generators: pythia6 (the embedding generator, quoted with the soft reweight and the sample
// corrections of the published analysis), pythia8 (Monash), pythia8detroit (the RHIC underlying-event
// tune), herwig.
// Inputs:  results/R<R>/ (analysis/run.sh), results/Michal/, results/models/spectra_<generator>.root.
// Outputs: results/models/models_R<R>_<def>.txt (the table with the ratios model / JP1) and .root
//          (every histogram). The figures are drawn by analysis/plot.C.
//   root -l -b -q 'compare.C+("0.5")'
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "../config.h"

using namespace CrossSectionConfig;

namespace {

// one generator sample of new_ana/models
struct Gen {
   const char *name; // sample name: the spectra_<name>.root file and the column of the table
   bool soft;        // true when the sample is quoted with the soft reweight of the published analysis
};
const std::vector<Gen> kGens = {{"pythia6", true}, {"pythia8", false}, {"pythia8detroit", false}, {"herwig", false}};

// a clone detached from its file, nullptr when the file or the histogram is absent
TH1D *Get(const std::string &file, const char *name, const char *as)
{
   if (gSystem->AccessPathName(file.c_str())) return nullptr;
   TFile f(file.c_str());
   auto *h = (TH1D *)f.Get(name);
   if (!h) return nullptr;
   h = (TH1D *)h->Clone(as);
   h->SetDirectory(0);
   return h;
}

// the fine (0.1 GeV) generator spectrum summed into the target bins and divided by their width
TH1D *Rebin(const TH1D *fine, const std::vector<double> &edges, const char *name)
{
   auto *h = new TH1D(name, "", edges.size() - 1, edges.data());
   h->SetDirectory(0);
   h->Sumw2();
   for (int k = 1; k <= h->GetNbinsX(); ++k) {
      double sum = 0;
      double errorSquared = 0;
      for (int i = 1; i <= fine->GetNbinsX(); ++i) {
         const double centre = fine->GetBinCenter(i);
         if (centre < edges[k - 1] || centre >= edges[k]) continue;
         sum += fine->GetBinContent(i);
         errorSquared += fine->GetBinError(i) * fine->GetBinError(i);
      }
      const double width = edges[k] - edges[k - 1];
      h->SetBinContent(k, sum / width);
      h->SetBinError(k, std::sqrt(errorSquared) / width);
   }
   return h;
}

// num / den bin by bin, with both statistical errors
TH1D *Ratio(const TH1D *num, const TH1D *den, const char *name)
{
   auto *h = (TH1D *)num->Clone(name);
   h->SetDirectory(0);
   for (int k = 1; k <= h->GetNbinsX(); ++k) {
      const double n = num->GetBinContent(k);
      const double d = den->GetBinContent(k);
      if (n <= 0 || d <= 0) {
         h->SetBinContent(k, 0);
         h->SetBinError(k, 0);
         continue;
      }
      const double relNum = num->GetBinError(k) / n;
      const double relDen = den->GetBinError(k) / d;
      h->SetBinContent(k, n / d);
      h->SetBinError(k, n / d * std::sqrt(relNum * relNum + relDen * relDen));
   }
   return h;
}

// the relative statistical error of one bin, in percent
double Pct(const TH1D *h, int k)
{
   const double v = h->GetBinContent(k);
   return v > 0 ? 100 * h->GetBinError(k) / v : 0;
}

} // namespace

void compare_def(const char *jetR, const char *def)
{
   const bool ue = std::string(def) == "ue";
   const std::string rdir = RadiusDir(jetR);
   const std::string mdir = kWorkDir + "results/models/";
   gSystem->mkdir(mdir.c_str(), true);
   // the published truth bins without the feed-down buffer, or the R_AA bins of the no-UE set
   const std::vector<double> edges =
      ue ? std::vector<double>(McBins().begin(), McBins().end() - 1) : std::vector<double>{5, 10, 15, 20, 25, 30, 35, 40, 50, 60};

   // 1. the data: JP1 with its band, the combined jet-patch levels, the min-bias level
   TH1D *jp1 = nullptr;
   TH1D *up = nullptr;
   TH1D *dn = nullptr;
   TH1D *mb = nullptr;
   TH1D *comb = nullptr;
   TH1D *pub = nullptr;
   TH1D *pubs = nullptr;
   if (ue) {
      jp1 = Get(rdir + "xsec_inversion_R" + jetR + "_jp1_e23_bayes4.root", "canonical", "jp1");
      up = Get(rdir + "xsec_jp1_e23_R" + jetR + "_syst_bayes4.root", "syst_up", "jp1_syst_up");
      dn = Get(rdir + "xsec_jp1_e23_R" + jetR + "_syst_bayes4.root", "syst_dn", "jp1_syst_dn");
      mb = Get(rdir + "xsec_inversion_R" + jetR + "_mb2023_bayes4.root", "canonical", "mb");
      comb = Get(rdir + "xsec_inversion_R" + jetR + "_jp_e23_bayes4.root", "canonical", "combined");
      if (std::string(jetR) == "0.5") {
         pub = Get(kWorkDir + "inputs/jet_cross_section_publishedR0.5.root", "crossSection_statistic", "published");
         pubs = Get(kWorkDir + "inputs/jet_cross_section_publishedR0.5.root", "crossSection_systematic", "published_syst");
      }
   } else {
      const std::string f = kWorkDir + "results/Michal/xsec_R" + jetR + "_noUE_bins5-60.root";
      jp1 = Get(f, "xsec_jp1", "jp1");
      up = Get(f, "xsec_jp1_syst_up", "jp1_syst_up");
      dn = Get(f, "xsec_jp1_syst_dn", "jp1_syst_dn");
      mb = Get(f, "xsec_mb", "mb");
      comb = Get(f, "xsec_combined", "combined");
   }
   if (!jp1 || !mb || !comb) {
      printf("[compare] R = %s %s: data missing\n", jetR, def);
      return;
   }
   if (jp1->GetNbinsX() != (int)edges.size() - 1) {
      printf("[compare] R = %s %s: data binning differs from the target bins\n", jetR, def);
      return;
   }
   if (!up) printf("[compare] R = %s %s: no JP1 systematic band (statistical only)\n", jetR, def);

   // 2. the generators, rebinned to the data bins, and their ratios to JP1
   std::vector<TH1D *> models(kGens.size(), nullptr);
   std::vector<TH1D *> ratios(kGens.size(), nullptr);
   for (size_t g = 0; g < kGens.size(); ++g) {
      const std::string hn = std::string(def) + "_R" + jetR + (kGens[g].soft ? "_soft" : "");
      TH1D *fine = Get(mdir + "spectra_" + kGens[g].name + ".root", hn.c_str(), (std::string(kGens[g].name) + "_fine").c_str());
      if (!fine) {
         printf("[compare] no %s spectrum (%s)\n", kGens[g].name, hn.c_str());
         continue;
      }
      models[g] = Rebin(fine, edges, kGens[g].name);
      ratios[g] = Ratio(models[g], jp1, (std::string(kGens[g].name) + "_over_jp1").c_str());
      // below 9.7 GeV the model bins are carried by a few jets of the softest pT-hat samples: not quoted
      for (int k = 1; k <= models[g]->GetNbinsX(); ++k) {
         if (models[g]->GetXaxis()->GetBinUpEdge(k) > (ue ? 9.7 : 10.0) + 1e-6) continue;
         models[g]->SetBinContent(k, 0);
         models[g]->SetBinError(k, 0);
         ratios[g]->SetBinContent(k, 0);
         ratios[g]->SetBinError(k, 0);
      }
      delete fine;
   }
   TH1D *pubRatio = pub ? Ratio(pub, jp1, "published_over_jp1") : nullptr;

   // 3. the table
   const std::string stem = mdir + "models_R" + jetR + "_" + def;
   FILE *t = fopen((stem + ".txt").c_str(), "w");
   fprintf(t, "# STAR Run-12 pp 200 GeV inclusive jets, anti-kT R = %s, |eta| < 0.5, d2sigma/dpT deta [pb/GeV], particle level, %s.\n", jetR,
           ue ? "underlying event subtracted at both levels (the nominal definition)" : "WITHOUT underlying-event subtraction (raw jets at both levels)");
   fprintf(t, "# Data: JP1 (Bayes 4 it., 2023 embedding; %s) with stat%% and its systematic band syst+%%/syst-%% (luminosity 5.6 %% not included%s);\n",
           ue ? "quoted from 8.2 GeV, the first bin is an extrapolation" : "quoted from 10 GeV, the 5-10 GeV bin is an extrapolation", up ? "" : "; NO band at this radius");
   fprintf(t, "#       combined = jet-patch levels JP0+JP1+JP2, mb = min-bias level (to ~20 GeV), statistical only%s.\n", pub ? "; published = Table III (stat%, syst%)" : "");
   fprintf(t, "# Models are quoted from %s GeV: below, the bins are carried by a few jets of the softest pT-hat samples (0 = not quoted).\n", ue ? "9.7" : "10");
   fprintf(t, "# Models (statistical errors, %%): pythia6 = PYTHIA 6.4.29 Perugia 2012 with PARP(90)=0.213 (the embedding generator: QCD 2->2 from\n");
   fprintf(t, "#   pT-hat 2 GeV, the soft reweight and the softest-bin corrections of the published analysis), pythia8 = PYTHIA 8.313 Monash 2013,\n");
   fprintf(t, "#   pythia8detroit = PYTHIA 8.313 with the RHIC Detroit tune (PRD 105, 016011), herwig = Herwig 7.3.0 default tune; all with the 13\n");
   fprintf(t, "#   particles of the embedding request undecayed. Ratios = model / JP1.\n");
   fprintf(t, "# %-11s %12s %7s %7s %7s %12s %7s %12s %7s", "pT [GeV]", "jp1", "stat%", "syst+%", "syst-%", "combined", "stat%", "mb", "stat%");
   if (pub) fprintf(t, " %12s %7s %7s", "published", "stat%", "syst%");
   for (const Gen &g : kGens) fprintf(t, " %12s %7s", g.name, "stat%");
   for (const Gen &g : kGens) fprintf(t, " %14s", (std::string(g.name) + "/jp1").c_str());
   if (pub) fprintf(t, " %9s", "pub/jp1");
   fprintf(t, "\n");
   for (int k = 1; k <= jp1->GetNbinsX(); ++k) {
      const double v = jp1->GetBinContent(k);
      fprintf(t, "  %4.1f-%-6.1f %12.4e %7.2f %7.2f %7.2f", edges[k - 1], edges[k], v, Pct(jp1, k), (up && v > 0) ? 100 * up->GetBinContent(k) / v : 0, (dn && v > 0) ? 100 * dn->GetBinContent(k) / v : 0);
      fprintf(t, " %12.4e %7.2f %12.4e %7.2f", comb->GetBinContent(k), Pct(comb, k), mb->GetBinContent(k), Pct(mb, k));
      if (pub) fprintf(t, " %12.4e %7.2f %7.2f", pub->GetBinContent(k), Pct(pub, k), pubs ? Pct(pubs, k) : 0);
      for (size_t g = 0; g < kGens.size(); ++g) fprintf(t, " %12.4e %7.2f", models[g] ? models[g]->GetBinContent(k) : 0, models[g] ? Pct(models[g], k) : 0);
      for (size_t g = 0; g < kGens.size(); ++g) fprintf(t, " %14.3f", ratios[g] ? ratios[g]->GetBinContent(k) : 0);
      if (pubRatio) fprintf(t, " %9.3f", pubRatio->GetBinContent(k));
      fprintf(t, "\n");
   }
   fclose(t);

   // 4. the same numbers as histograms, for the figures of analysis/plot.C
   TFile fo((stem + ".root").c_str(), "RECREATE");
   for (TH1D *h : {jp1, up, dn, mb, comb, pub, pubs, pubRatio}) {
      if (h) h->Write();
   }
   for (size_t g = 0; g < kGens.size(); ++g) {
      if (!models[g]) continue;
      models[g]->Write();
      ratios[g]->Write();
   }
   fo.Close();
   printf("[compare] wrote %s.{txt,root}\n", stem.c_str());
}

void compare(const char *jetR = "0.5")
{
   compare_def(jetR, "ue");
   compare_def(jetR, "noue");
}
