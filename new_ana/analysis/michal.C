// michal.C — the pp reference for the R_AA comparison: particle-level inclusive jet cross sections
// WITHOUT the underlying-event subtraction (raw jets at detector and particle level), |eta| < 0.5, in the
// coarse bins 5 10 15 20 25 30 35 40 50 60 GeV.
//
// Reference column: the JP1 trigger with its systematic band (systematics.C, mode jp1_e23_noue_mbins:
// energy scale, tracking, embedding statistics, underlying event, regularisation; up and down added in
// quadrature; the 5.6 % luminosity uncertainty is left outside). It is quoted from 10 GeV; below that the
// JP1 matching efficiency is 2 %, so the 5-10 GeV bin is read off the combined column.
// Comparison columns, statistical errors only: "combined" = the exclusive jet-patch levels JP0+JP1+JP2
// over the whole range, "mb" = the min-bias level up to 20 GeV.
//
// Input : <results>/R<R>/xsec_inversion_R<R>_<mode>_noue_mbins[_bayes<n>].root (canonical) and the band
//         file xsec_jp1_e23_noue_mbins_R<R>_syst[_bayes<n>].root (syst_up, syst_dn).
// Output: <results>/Michal/xsec_R<R>_noUE_bins5-60.root (the three spectra + the band) and .txt (the table).
// Run   : root -l -b -q 'michal.C+("0.5")'
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <cstdio>
#include <string>
#include <vector>
#include "../config.h"
using namespace CrossSectionConfig;

// Lowest pT at which the JP1 reference is quoted; below it the trigger is still in its turn-on.
static const double kJp1QuoteLoGeV = 10.0;
// Highest pT the min-bias level reaches with useful statistics.
static const double kMbQuoteHiGeV = 20.0;

// One column of the table: the unfolding mode that produced it and the name it is written under.
struct Column {
   const char *label; // short name in the table header and in the ROOT file
   const char *mode;  // unfolding mode of stage2/unfold.C
};

// The unfolded no-UE spectrum of one mode, detached from its file and titled for the output.
static TH1D *ReadSpectrum(const std::string &radiusDir, const char *jetR, const Column &col, const std::string &bayesTag)
{
   const std::string path =
      radiusDir + "xsec_inversion_R" + jetR + "_" + col.mode + "_noue_mbins" + bayesTag + ".root";
   if (gSystem->AccessPathName(path.c_str())) {
      printf("[michal] missing %s\n", path.c_str());
      return nullptr;
   }
   TFile in(path.c_str());
   TH1D *spectrum = (TH1D *)((TH1D *)in.Get("canonical"))->Clone(Form("xsec_%s", col.label));
   spectrum->SetDirectory(nullptr);
   spectrum->SetTitle(Form("d^{2}#sigma/dp_{T}d#eta, |#eta|<0.5, anti-k_{T} R=%s, "
                           "no UE subtraction, %s;p_{T} [GeV/c];pb/(GeV/c)",
                           jetR, col.label));
   return spectrum;
}

// The systematic band of the JP1 reference; absent until analysis/run.sh michal has run for this radius.
static void ReadJp1Band(const std::string &radiusDir, const char *jetR, const std::string &bayesTag, TH1D *&up,
                        TH1D *&down)
{
   up = nullptr;
   down = nullptr;
   const std::string path = radiusDir + "xsec_jp1_e23_noue_mbins_R" + jetR + "_syst" + bayesTag + ".root";
   if (gSystem->AccessPathName(path.c_str())) {
      printf("[michal] no systematic band yet for R = %s\n", jetR);
      return;
   }
   TFile in(path.c_str());
   up = (TH1D *)((TH1D *)in.Get("syst_up"))->Clone("syst_up");
   up->SetDirectory(nullptr);
   down = (TH1D *)((TH1D *)in.Get("syst_dn"))->Clone("syst_dn");
   down->SetDirectory(nullptr);
}

// One quantity as a percentage of the cross section in the same bin, guarded against empty bins.
static double PercentOfValue(double value, double reference)
{
   return reference > 0 ? 100 * value / reference : 0;
}

// The plain-text table: one row per pT bin, the JP1 reference with its band first, then the other columns.
static void WriteTable(const std::string &path, const char *jetR, int nIter, const std::vector<Column> &columns,
                       const std::vector<TH1D *> &spectra, const TH1D *up, const TH1D *down)
{
   FILE *table = fopen(path.c_str(), "w");
   fprintf(table,
           "# STAR Run-12 pp 200 GeV inclusive jets, anti-kT R = %s, |eta| < 0.5, d2sigma/dpT deta in pb/GeV, "
           "particle level,\n",
           jetR);
   fprintf(table,
           "# WITHOUT underlying-event subtraction (raw jets at detector and particle level). "
           "Bayesian unfolding (%d it.), 2023 embedding.\n",
           nIter);
   fprintf(table, "# Statistical (data) and systematic uncertainties in percent of the value; "
                  "the luminosity uncertainty (5.6 %%) is not included.\n");
   fprintf(table, "# jp1 = the reference (JP1 trigger, from %.0f GeV) with its systematic band (syst+%% / syst-%%);\n",
           kJp1QuoteLoGeV);
   fprintf(table,
           "# combined = jet-patch levels JP0+JP1+JP2 (the 5-%.0f GeV bin), mb = min-bias level (to %.0f GeV), "
           "both statistical only.\n",
           kJp1QuoteLoGeV, kMbQuoteHiGeV);

   fprintf(table, "# %-11s %12s %8s %8s %8s", "pT [GeV]", "jp1", "stat%", "syst+%", "syst-%");
   for (const Column &col : columns) {
      if (std::string(col.label) != "jp1") fprintf(table, " %12s %8s", col.label, "stat%");
   }
   fprintf(table, "\n");

   const TH1D *jp1 = spectra[1];
   for (int bin = 1; bin <= spectra[0]->GetNbinsX(); ++bin) {
      fprintf(table, "  %4.0f-%-6.0f", spectra[0]->GetBinLowEdge(bin), spectra[0]->GetBinLowEdge(bin + 1));
      const double reference = jp1->GetBinContent(bin);
      fprintf(table, " %12.4e %8.2f %8.2f %8.2f", reference, PercentOfValue(jp1->GetBinError(bin), reference),
              up ? PercentOfValue(up->GetBinContent(bin), reference) : 0,
              down ? PercentOfValue(down->GetBinContent(bin), reference) : 0);
      for (size_t i = 0; i < spectra.size(); ++i) {
         if (i == 1) continue; // the JP1 reference is already printed
         const double value = spectra[i]->GetBinContent(bin);
         fprintf(table, " %12.4e %8.2f", value, PercentOfValue(spectra[i]->GetBinError(bin), value));
      }
      fprintf(table, "\n");
   }
   fclose(table);
}

void michal(const char *jetR = "0.5", const char *method = "bayes", int nIter = 4)
{
   const bool bayes = std::string(method) == "bayes";
   const std::string bayesTag = bayes ? Form("_bayes%d", nIter) : "";
   const std::string radiusDir = RadiusDir(jetR);
   const std::string outDir = AnalysisConfig().workdir + "results/Michal/";
   gSystem->mkdir(outDir.c_str(), true);

   // 1. the three spectra, all unfolded from the no-UE variant on the coarse 5-60 GeV bins
   const std::vector<Column> columns = {{"combined", "jp_e23"}, {"jp1", "jp1_e23"}, {"mb", "mb2023"}};
   std::vector<TH1D *> spectra;
   for (const Column &col : columns) {
      TH1D *spectrum = ReadSpectrum(radiusDir, jetR, col, bayesTag);
      if (!spectrum) return;
      spectra.push_back(spectrum);
   }

   // 2. the systematic band of the JP1 reference
   TH1D *up = nullptr;
   TH1D *down = nullptr;
   ReadJp1Band(radiusDir, jetR, bayesTag, up, down);

   // 3. write the spectra and the band
   const std::string stem = outDir + "xsec_R" + jetR + "_noUE_bins5-60";
   TFile out((stem + ".root").c_str(), "RECREATE");
   for (TH1D *spectrum : spectra) spectrum->Write();
   if (up) {
      up->Write("xsec_jp1_syst_up");
      down->Write("xsec_jp1_syst_dn");
   }
   out.Close();

   // 4. the same numbers as a table
   WriteTable(stem + ".txt", jetR, nIter, columns, spectra, up, down);
   printf("[michal] wrote %s.{txt,root}\n", stem.c_str());
}
