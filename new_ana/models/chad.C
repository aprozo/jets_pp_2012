// chad.C — the hadronisation correction of the fixed-order comparison, as in the STAR paper
// (arXiv:2603.28695, Sec. V): C_had(pT) = sigma_particle / sigma_parton per bin, from Pythia 6
// Perugia 2012 (PARP(90) = 0.213), both levels with the off-axis-cone underlying-event subtraction
// and the parton level taken from the record before hadronisation (the samples of
// "run_models.sh chad": pythia6_t<tune> and pythia6_t<tune>_partons, the same events).
// The uncertainty is the quadrature sum of three variations of the Perugia 2012 tune: radiation
// (371 radHi / 372 radLo, the larger deviation), fragmentation (376 FL / 377 FT, the larger
// deviation) and the Innsbruck hadronisation tune (373).
// Inputs:  results/models/spectra_pythia6_t<tune>[_partons].root.
// Outputs: results/models/chad_R0.5.txt and .root, in the published bins (the nominal definition)
//          and in the 5-60 GeV bins. The figure is drawn by analysis/plot.C.
//   root -l -b -q 'chad.C+'
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "../config.h"

using namespace CrossSectionConfig;

namespace {

// the Perugia 2012 tunes: the nominal one and the five variations that make the uncertainty
const int kTuneNominal = 370;
const int kTuneRadHi = 371;
const int kTuneRadLo = 372;
const int kTuneInnsbruck = 373;
const int kTuneFragL = 376;
const int kTuneFragT = 377;

// the jet spectrum of one tune at one level (particle or parton), R = 0.5
TH1D *Load(const std::string &mdir, int tune, bool partons, const char *def)
{
   const std::string f = mdir + "spectra_pythia6_t" + std::to_string(tune) + (partons ? "_partons" : "") + ".root";
   if (gSystem->AccessPathName(f.c_str())) return nullptr;
   TFile in(f.c_str());
   auto *h = (TH1D *)in.Get(Form("%s_R0.5", def));
   if (!h) return nullptr;
   h = (TH1D *)h->Clone(Form("%s_t%d%s", def, tune, partons ? "_partons" : ""));
   h->SetDirectory(0);
   return h;
}

// the fine spectrum summed into the target bins; no division by the width, the ratio does not need it
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
      h->SetBinContent(k, sum);
      h->SetBinError(k, std::sqrt(errorSquared));
   }
   return h;
}

// C_had of one tune: the particle-level spectrum over the parton-level one of the same events
TH1D *CHad(const std::string &mdir, int tune, const char *def, const std::vector<double> &edges, const char *name)
{
   TH1D *particleFine = Load(mdir, tune, false, def);
   TH1D *partonFine = Load(mdir, tune, true, def);
   if (!particleFine || !partonFine) {
      printf("[chad] tune %d: missing spectra\n", tune);
      return nullptr;
   }
   TH1D *particle = Rebin(particleFine, edges, "a");
   TH1D *parton = Rebin(partonFine, edges, "b");
   auto *chad = (TH1D *)particle->Clone(name);
   chad->SetDirectory(0);
   for (int k = 1; k <= chad->GetNbinsX(); ++k) {
      const double x = particle->GetBinContent(k);
      const double y = parton->GetBinContent(k);
      // the two levels are the same events: the statistical error comes from the particle level alone
      chad->SetBinContent(k, y > 0 ? x / y : 0);
      chad->SetBinError(k, (x > 0 && y > 0) ? x / y * particle->GetBinError(k) / x : 0);
   }
   delete particle;
   delete parton;
   delete particleFine;
   delete partonFine;
   return chad;
}

// the deviation of one variation from the nominal C_had in one bin
double Deviation(const TH1D *variation, double nominal, int k)
{
   return variation ? std::fabs(variation->GetBinContent(k) - nominal) : 0.0;
}

} // namespace

// C_had and its uncertainty in one set of bins: the table rows and the histograms of one block.
// tag names the histograms inside the output file ("pub" for the published bins, "raa" for 5-60 GeV).
void chad_def(const char *def, const std::vector<double> &edges, const char *tag, FILE *t, TFile &fo)
{
   const std::string mdir = kWorkDir + "results/models/";
   TH1D *nominal = CHad(mdir, kTuneNominal, def, edges, Form("chad_%s", tag));
   if (!nominal) return;
   TH1D *radiation[2] = {CHad(mdir, kTuneRadHi, def, edges, Form("chad_%s_radHi", tag)),
                         CHad(mdir, kTuneRadLo, def, edges, Form("chad_%s_radLo", tag))};
   TH1D *fragmentation[2] = {CHad(mdir, kTuneFragL, def, edges, Form("chad_%s_FL", tag)),
                             CHad(mdir, kTuneFragT, def, edges, Form("chad_%s_FT", tag))};
   TH1D *innsbruck = CHad(mdir, kTuneInnsbruck, def, edges, Form("chad_%s_IBK", tag));

   auto *uncertainty = (TH1D *)nominal->Clone(Form("chad_%s_syst", tag));
   uncertainty->SetDirectory(0);
   fprintf(t, "# %s: C_had = sigma_particle / sigma_parton, Pythia 6 Perugia 2012 (PARP(90)=0.213), R = 0.5, |eta| < 0.5, UE subtracted at both levels\n", def);
   fprintf(t, "# %-11s %8s %7s %8s %8s %8s %8s\n", "pT [GeV]", "C_had", "stat", "syst", "rad", "frag", "IBK");
   for (int k = 1; k <= nominal->GetNbinsX(); ++k) {
      const double value = nominal->GetBinContent(k);
      const double dRad = std::max(Deviation(radiation[0], value, k), Deviation(radiation[1], value, k));
      const double dFrag = std::max(Deviation(fragmentation[0], value, k), Deviation(fragmentation[1], value, k));
      const double dInnsbruck = Deviation(innsbruck, value, k);
      const double total = std::sqrt(dRad * dRad + dFrag * dFrag + dInnsbruck * dInnsbruck);
      uncertainty->SetBinContent(k, total);
      uncertainty->SetBinError(k, 0);
      fprintf(t, "  %4.1f-%-6.1f %8.4f %7.4f %8.4f %8.4f %8.4f %8.4f\n", edges[k - 1], edges[k], value, nominal->GetBinError(k), total, dRad, dFrag, dInnsbruck);
   }

   fo.cd();
   nominal->Write();
   uncertainty->Write();
   for (TH1D *h : {radiation[0], radiation[1], fragmentation[0], fragmentation[1], innsbruck}) {
      if (h) h->Write();
   }
}

void chad()
{
   const std::string mdir = kWorkDir + "results/models/";
   gSystem->mkdir(mdir.c_str(), true);
   // the published truth bins without the feed-down buffer, and the R_AA bins of the no-UE set
   const std::vector<double> publishedBins(McBins().begin(), McBins().end() - 1);
   const std::vector<double> raaBins = {5, 10, 15, 20, 25, 30, 35, 40, 50, 60};

   FILE *t = fopen((mdir + "chad_R0.5.txt").c_str(), "w");
   fprintf(t, "# Hadronisation correction for the fixed-order (parton-level) predictions, as in arXiv:2603.28695 Sec. V:\n");
   fprintf(t, "# C_had = sigma_particle / sigma_parton per bin from Pythia 6 Perugia 2012 with PARP(90) = 0.213 (the embedding generator),\n");
   fprintf(t, "# the parton level = the record before hadronisation, the off-axis-cone UE subtraction applied at both levels.\n");
   fprintf(t, "# syst = quadrature of rad (max of radHi 371 / radLo 372), frag (max of FL 376 / FT 377) and IBK (373); stat = MC statistics.\n");
   TFile fo((mdir + "chad_R0.5.root").c_str(), "RECREATE");
   chad_def("ue", publishedBins, "pub", t, fo);
   fprintf(t, "#\n");
   chad_def("ue", raaBins, "raa", t, fo);
   fclose(t);
   fo.Close();
   printf("[chad] wrote %schad_R0.5.{txt,root}\n", mdir.c_str());
}
