// ue_onoff.C — the inclusive jet cross section with and without the underlying-event subtraction, per
// trigger, so that the size of the UE correction can be read off bin by bin.
//
// "with UE sub." is the nominal jet definition (the off-axis-cone UE subtracted at detector and particle
// level); "without" is the variant "noue" of variants.h (raw jet pT at both levels, i.e. jets including
// their underlying event). The two are a different observable, not a systematic, and are unfolded the
// same way (Bayes with nIter iterations, or matrix inversion).
//
// Input : <results>/R<R>/xsec_inversion_R<R>_<trigger>[_noue][_bayes<n>].root (canonical).
// Output: <results>/R<R>/ue_onoff_R<R>[_bayes<n>].txt — per trigger and pT bin the two cross sections
//         (pb/GeV), their statistical errors in percent, and the ratio without / with.
// Run   : root -l -b -q 'ue_onoff.C+("0.5")'
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <cstdio>
#include <string>
#include <vector>
#include "../config.h"
using namespace CrossSectionConfig;

// Split a space-separated argument ("jp1_e23 jp2_e23 ...") into its words.
static std::vector<std::string> SplitWords(const std::string &text)
{
   std::vector<std::string> words;
   for (size_t begin = 0; begin < text.size();) {
      const size_t space = text.find(' ', begin);
      const size_t count = (space == std::string::npos) ? std::string::npos : space - begin;
      words.push_back(text.substr(begin, count));
      if (space == std::string::npos) break;
      begin = space + 1;
   }
   return words;
}

// Statistical error in percent of the cross section, guarded against empty bins.
static double StatErrorPercent(double value, double reference)
{
   return reference > 0 ? 100 * value / reference : 0;
}

void ue_onoff(const char *jetR = "0.5", const char *triggers = "jp1_e23 jp2_e23 ht2_e23 mb2023",
              const char *method = "bayes", int nIter = 4)
{
   const bool bayes = std::string(method) == "bayes";
   const std::string bayesTag = bayes ? Form("_bayes%d", nIter) : "";
   const std::string radiusDir = RadiusDir(jetR);
   const std::vector<std::string> triggerList = SplitWords(triggers);

   // 1. the table header
   const std::string outPath = radiusDir + "ue_onoff_R" + jetR + bayesTag + ".txt";
   FILE *table = fopen(outPath.c_str(), "w");
   fprintf(table,
           "# d2sigma/dpT deta [pb/GeV], |eta| < 0.5, R = %s, %s: jets with the UE subtracted (nominal) "
           "and without (raw jets, both levels)\n",
           jetR, bayes ? Form("Bayes %d iterations", nIter) : "inversion");
   fprintf(table, "# %-9s %-12s %12s %7s %12s %7s %8s\n", "trigger", "pT [GeV]", "with UE sub.", "stat%", "without",
           "stat%", "ratio");

   // 2. one block of rows per trigger
   for (const std::string &trigger : triggerList) {
      const std::string pathWith = radiusDir + "xsec_inversion_R" + jetR + "_" + trigger + bayesTag + ".root";
      const std::string pathWithout = radiusDir + "xsec_inversion_R" + jetR + "_" + trigger + "_noue" + bayesTag + ".root";
      if (gSystem->AccessPathName(pathWith.c_str()) || gSystem->AccessPathName(pathWithout.c_str())) {
         printf("[ue] %s: missing input\n", trigger.c_str());
         continue;
      }
      TFile fileWith(pathWith.c_str());
      TFile fileWithout(pathWithout.c_str());
      TH1D *withUe = (TH1D *)fileWith.Get("canonical");
      TH1D *withoutUe = (TH1D *)fileWithout.Get("canonical");
      for (int bin = 1; bin <= withUe->GetNbinsX(); ++bin) {
         const double valueWith = withUe->GetBinContent(bin);
         const double valueWithout = withoutUe->GetBinContent(bin);
         fprintf(table, "  %-9s %5.1f-%-6.1f %12.4e %7.2f %12.4e %7.2f %8.3f\n", trigger.c_str(),
                 withUe->GetBinLowEdge(bin), withUe->GetBinLowEdge(bin + 1), valueWith,
                 StatErrorPercent(withUe->GetBinError(bin), valueWith), valueWithout,
                 StatErrorPercent(withoutUe->GetBinError(bin), valueWithout),
                 valueWith > 0 ? valueWithout / valueWith : 0);
      }
   }
   fclose(table);
   printf("[ue] wrote %s\n", outPath.c_str());
}
