// systematics.C — systematic uncertainties of an unfolded cross section from the variant results of
// unfold.C, and the trigger-to-trigger and radius-to-radius ratios with the same variants applied to both.
//
// combine(R, mode, method, nIter): reads <results>/xsec_inversion_R<R>_<mode>[_<member>][_bayes<n>].root
//   for the nominal and every member of every component (variants.h), builds per component the envelope
//   around the nominal (up = max(members) - nominal, down = nominal - min(members), each floored at zero;
//   the embedding-statistics component takes the per-bin RMS of the replica loop, symmetric), adds the up
//   and the down deviations of the components in quadrature separately, and writes
//   <results>/xsec_<mode>_R<R>_syst[_bayes<n>].root: canonical (nominal, statistical errors), syst_up /
//   syst_dn (absolute), syst_<component>_up / _dn, canonical_syst (TGraphAsymmErrors), and the text table
//   <results>/syst_<mode>_R<R>[_bayes<n>].txt with max(up, down) / value per bin rounded up to 0.1 %.
//   Luminosity (5.6 %) is not in the band.
// ratios(R, base, others, method, nIter): the same for r = sigma_T / sigma_base, T in others: every member
//   recomputed as the ratio of the two variant results, so that the shifts common to both cancel; output
//   <results>/ratio_<T>_over_<base>_R<R>[_bayes<n>].root and .txt.
// ratios_radii(mode, radii, baseR, method, nIter): the same for sigma(R) / sigma(baseR) of one mode.
// list_variants(what): the variant names the shell drivers (lib.sh) loop over.
// root -l -b -q 'systematics.C+("0.5", "jp", "inv", 0)'
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "common.h"
#include "variants.h"

using namespace CrossSectionConfig;
using namespace Stage2;

namespace {

// the luminosity uncertainty is quoted separately and is not part of the band
const double kLuminosityPercent = 5.6;
// the +-3 % on the HT2 trigger-efficiency factor of unfold.C::Ht2Eff
const double kHt2EfficiencyUncertainty = 0.03;

// ---- one result file ---------------------------------------------------------------------------------

// the file of one unfolding result: the member appended to the mode, except the iteration members, which
// are the same mode at their own iteration count
std::string XsecFile(const std::string &radiusDir, const char *jetR, const std::string &mode,
                     const std::string &member, bool bayes, int nIter)
{
   std::string modeWithMember = mode;
   int iterations = nIter;
   if (member == "iter2") iterations = 2;
   else if (member == "iter6") iterations = 6;
   else if (!member.empty()) modeWithMember += "_" + member;
   // unfold.C writes the bare mode "jp" without a mode suffix
   const std::string modeSuffix = modeWithMember == "jp" ? "" : "_" + modeWithMember;
   const std::string iterSuffix = bayes ? Form("_bayes%d", iterations) : "";
   return radiusDir + "xsec_inversion_R" + jetR + modeSuffix + iterSuffix + ".root";
}

// one histogram out of a result file, detached from it
TH1D *Get(const std::string &file, const char *name, bool required = true)
{
   if (gSystem->AccessPathName(file.c_str())) {
      if (required) throw std::runtime_error("missing " + file);
      return nullptr;
   }
   TFile input(file.c_str(), "READ");
   auto *histogram = (TH1D *)input.Get(name);
   if (!histogram) {
      if (required) throw std::runtime_error(std::string("missing ") + name + " in " + file);
      return nullptr;
   }
   histogram = (TH1D *)histogram->Clone();
   histogram->SetDirectory(0);
   return histogram;
}

// ---- one component of the band -----------------------------------------------------------------------

// the absolute up and down deviation of one component, per bin
struct Band {
   std::vector<double> up;
   std::vector<double> dn;
};

typedef std::pair<std::string, Band> NamedBand;

// Envelope of one component around the nominal: the highest and the lowest member, each floored at the
// nominal so that a one-sided component stays one-sided. A component given as an RMS (the embedding
// statistics) carries its deviation as the error of its own histogram and is symmetric.
Band Envelope(const TH1D *nominal, const std::vector<TH1D *> &members, const TH1D *rmsHist)
{
   const int nBins = nominal->GetNbinsX();
   Band band;
   band.up.assign(nBins, 0);
   band.dn.assign(nBins, 0);
   if (rmsHist) {
      for (int k = 1; k <= nBins; ++k) {
         band.up[k - 1] = rmsHist->GetBinError(k);
         band.dn[k - 1] = rmsHist->GetBinError(k);
      }
      return band;
   }
   for (int k = 1; k <= nBins; ++k) {
      const double nominalValue = nominal->GetBinContent(k);
      double highest = nominalValue;
      double lowest = nominalValue;
      for (const TH1D *member : members) {
         highest = std::max(highest, member->GetBinContent(k));
         lowest = std::min(lowest, member->GetBinContent(k));
      }
      band.up[k - 1] = highest - nominalValue;
      band.dn[k - 1] = nominalValue - lowest;
   }
   return band;
}

// the quoted percentage: max(up, down) / value, rounded up to 0.1 %
double Pct(double up, double dn, double value)
{
   return value > 0 ? std::ceil(1000.0 * std::max(up, dn) / value) / 10.0 : 0;
}

// The pT range a mode is quoted over: JP1 from 8.2 GeV (its detector window), JP2 from 9.7, HT2 from
// 11.5, min-bias up to 22.5 GeV; on the 5-60 GeV bins JP1/JP2/HT2 from 10 GeV and min-bias up to 20 GeV.
// Bins outside are kept in the .root output and left out of the table.
void QuoteRange(const std::string &mode, double &lo, double &hi)
{
   const std::string trigger = mode.substr(0, mode.find('_'));
   const bool wideBins = mode.find("mbins") != std::string::npos;
   lo = wideBins ? 5.0 : 6.9;
   hi = wideBins ? 60.0 : 52.0;
   if (trigger == "jp1") lo = wideBins ? 10.0 : 8.2;
   else if (trigger == "jp2") lo = wideBins ? 10.0 : 9.7;
   else if (trigger == "ht2") lo = wideBins ? 10.0 : 11.5;
   else if (trigger.rfind("mb", 0) == 0) hi = wideBins ? 20.0 : 22.5;
}

// component names become histogram names, so spaces turn into underscores
std::string HistName(const std::string &componentName)
{
   std::string name = componentName;
   std::replace(name.begin(), name.end(), ' ', '_');
   return name;
}

// ---- output ---------------------------------------------------------------------------------------------

// the up and down deviations of the components added in quadrature, separately
Band TotalBand(const std::vector<NamedBand> &components, int nBins)
{
   Band total;
   total.up.assign(nBins, 0);
   total.dn.assign(nBins, 0);
   for (const NamedBand &component : components) {
      for (int k = 0; k < nBins; ++k) {
         total.up[k] += component.second.up[k] * component.second.up[k];
         total.dn[k] += component.second.dn[k] * component.second.dn[k];
      }
   }
   for (int k = 0; k < nBins; ++k) {
      total.up[k] = std::sqrt(total.up[k]);
      total.dn[k] = std::sqrt(total.dn[k]);
   }
   return total;
}

// one deviation vector as a histogram on the binning of the nominal
void WriteDeviation(const TH1D *nominal, const char *name, const std::vector<double> &values)
{
   TH1D *histogram = (TH1D *)nominal->Clone(name);
   for (int k = 1; k <= nominal->GetNbinsX(); ++k) {
      histogram->SetBinContent(k, values[k - 1]);
      histogram->SetBinError(k, 0);
   }
   histogram->Write(name);
}

// the .root output: the nominal, the total and per-component deviations, and the band as a graph
void WriteRootBand(const std::string &path, const TH1D *nominal, const std::vector<NamedBand> &components,
                   const Band &total, double quoteLo, double quoteHi)
{
   const int nBins = nominal->GetNbinsX();
   TFile out(path.c_str(), "RECREATE");
   nominal->Write("canonical");
   WriteDeviation(nominal, "syst_up", total.up);
   WriteDeviation(nominal, "syst_dn", total.dn);
   for (const NamedBand &component : components) {
      const std::string name = HistName(component.first);
      WriteDeviation(nominal, ("syst_" + name + "_up").c_str(), component.second.up);
      WriteDeviation(nominal, ("syst_" + name + "_dn").c_str(), component.second.dn);
   }
   TGraphAsymmErrors graph(nBins);
   for (int k = 1; k <= nBins; ++k) {
      graph.SetPoint(k - 1, nominal->GetBinCenter(k), nominal->GetBinContent(k));
      graph.SetPointError(k - 1, nominal->GetBinWidth(k) / 2, nominal->GetBinWidth(k) / 2, total.dn[k - 1],
                          total.up[k - 1]);
   }
   graph.Write("canonical_syst");
   TParameter<double>("quote_lo", quoteLo).Write();
   TParameter<double>("quote_hi", quoteHi).Write();
   out.Close();
}

// the .txt table: per bin the value, the statistical error, the total and every component, in percent
void WriteTextBand(const std::string &path, const std::string &title, const TH1D *nominal,
                   const std::vector<NamedBand> &components, const Band &total, double quoteLo, double quoteHi)
{
   const int nBins = nominal->GetNbinsX();
   FILE *table = fopen(path.c_str(), "w");
   fprintf(table, "# %s: systematic uncertainties in percent of the value, max(up, down) rounded up to 0.1 %%\n",
           title.c_str());
   fprintf(table, "# quoted range %.1f-%.1f GeV; the luminosity uncertainty (%.1f %%) is not included\n", quoteLo,
           quoteHi, kLuminosityPercent);
   fprintf(table, "# %-12s %10s %8s %8s", "pT [GeV]", "value", "stat%", "total%");
   for (const NamedBand &component : components) fprintf(table, " %14s", HistName(component.first).c_str());
   fprintf(table, " %8s %8s\n", "up%", "down%");
   for (int k = 1; k <= nBins; ++k) {
      if (nominal->GetBinLowEdge(k) < quoteLo - 1e-6 || nominal->GetBinLowEdge(k + 1) > quoteHi + 1e-6) continue;
      const double value = nominal->GetBinContent(k);
      fprintf(table, "  %5.1f-%-6.1f %10.4g %8.2f %8.1f", nominal->GetBinLowEdge(k), nominal->GetBinLowEdge(k + 1),
              value, value > 0 ? 100 * nominal->GetBinError(k) / value : 0,
              Pct(total.up[k - 1], total.dn[k - 1], value));
      for (const NamedBand &component : components)
         fprintf(table, " %14.1f", Pct(component.second.up[k - 1], component.second.dn[k - 1], value));
      fprintf(table, " %8.2f %8.2f\n", value > 0 ? 100 * total.up[k - 1] / value : 0,
              value > 0 ? 100 * total.dn[k - 1] / value : 0);
   }
   fclose(table);
}

// the same table on the terminal, so that a driver run shows the band it just built
void PrintBand(const TH1D *nominal, const std::vector<NamedBand> &components, const Band &total)
{
   for (int k = 1; k <= nominal->GetNbinsX(); ++k) {
      const double value = nominal->GetBinContent(k);
      printf("   %5.1f-%-5.1f total %5.1f%%", nominal->GetBinLowEdge(k), nominal->GetBinLowEdge(k + 1),
             Pct(total.up[k - 1], total.dn[k - 1], value));
      for (const NamedBand &component : components)
         printf("  %s %5.1f%%", component.first.c_str(),
                Pct(component.second.up[k - 1], component.second.dn[k - 1], value));
      printf("\n");
   }
}

void WriteOut(const std::string &radiusDir, const std::string &stem, TH1D *nominal,
              const std::vector<NamedBand> &components, const std::string &title, double quoteLo = 0,
              double quoteHi = 1e9)
{
   const Band total = TotalBand(components, nominal->GetNbinsX());
   WriteRootBand(radiusDir + stem + ".root", nominal, components, total, quoteLo, quoteHi);
   WriteTextBand(radiusDir + stem + ".txt", title, nominal, components, total, quoteLo, quoteHi);
   printf("[syst] %s\n", (radiusDir + stem + ".txt").c_str());
   PrintBand(nominal, components, total);
}

// the method tag of the output names
std::string MethodTag(bool bayes, int nIter)
{
   return bayes ? Form("_bayes%d", nIter) : "";
}

} // namespace

// the variant lists for the drivers (stage2/lib.sh): what = response | data | databuilds | members_inv |
// members_bayes
void list_variants(const char *what = "response")
{
   const std::string which = what;
   std::vector<std::string> names;
   if (which == "response") names = Syst::ResponseVariants();
   else if (which == "data") names = Syst::DataVariants();
   else if (which == "databuilds") names = Syst::DataBuilds();
   else if (which.rfind("members_", 0) == 0)
      for (const Syst::Component &component : Syst::Components(which == "members_bayes"))
         for (const std::string &member : component.members) names.push_back(member);
   printf("VARIANTS");
   for (const std::string &name : names) printf(" %s", name.c_str());
   printf("\n");
}

// the band of one unfolded cross section
void combine(const char *jetR = "0.5", const char *mode = "jp", const char *method = "inv", int nIter = 0)
{
   const bool bayes = std::string(method) == "bayes";
   const std::string radiusDir = RadiusDir(jetR);
   TH1D *nominal = Get(XsecFile(radiusDir, jetR, mode, "", bayes, nIter), "canonical");

   std::vector<NamedBand> components;
   for (const Syst::Component &component : Syst::Components(bayes)) {
      std::vector<TH1D *> members;
      TH1D *rms = nullptr;
      for (const std::string &member : component.members) {
         const std::string file = XsecFile(radiusDir, jetR, mode, member, bayes, nIter);
         TH1D *histogram = Get(file, component.rms ? "canonical_embstat" : "canonical", false);
         if (!histogram) {
            printf("[syst] WARNING member %s missing (%s) - component %s incomplete\n", member.c_str(),
                   file.c_str(), component.name.c_str());
            continue;
         }
         if (component.rms) rms = histogram;
         else members.push_back(histogram);
      }
      components.push_back({component.name, Envelope(nominal, members, rms)});
   }

   // the high-tower modes carry the uncertainty of their trigger-efficiency factor on top
   if (std::string(mode).rfind("ht2", 0) == 0) {
      Band band;
      for (int k = 1; k <= nominal->GetNbinsX(); ++k) {
         band.up.push_back(kHt2EfficiencyUncertainty * nominal->GetBinContent(k));
         band.dn.push_back(kHt2EfficiencyUncertainty * nominal->GetBinContent(k));
      }
      components.push_back({"trigger efficiency", band});
   }

   double quoteLo = 0;
   double quoteHi = 0;
   QuoteRange(mode, quoteLo, quoteHi);
   const std::string stem = std::string("xsec_") + mode + "_R" + jetR + "_syst" + MethodTag(bayes, nIter);
   const std::string title =
      std::string("mode ") + mode + ", " + method + (bayes ? Form(" %d iterations", nIter) : "");
   WriteOut(radiusDir, stem, nominal, components, title, quoteLo, quoteHi);
}

// ratio of every trigger in others to the base trigger, with each member recomputed on both, so that the
// shifts common to the two triggers cancel
void ratios(const char *jetR = "0.5", const char *base = "jp1_e23", const char *others = "mb2023 jp2_e23 ht2_e23",
            const char *method = "bayes", int nIter = 4)
{
   const bool bayes = std::string(method) == "bayes";
   const std::string radiusDir = RadiusDir(jetR);
   const std::vector<std::string> triggers = SplitWords(others);
   TH1D *baseNominal = Get(XsecFile(radiusDir, jetR, base, "", bayes, nIter), "canonical");

   for (const std::string &trigger : triggers) {
      TH1D *triggerNominal = Get(XsecFile(radiusDir, jetR, trigger, "", bayes, nIter), "canonical");
      TH1D *nominal = RatioWithErrors(triggerNominal, baseNominal, "canonical");
      std::vector<NamedBand> components;
      for (const Syst::Component &component : Syst::Components(bayes)) {
         std::vector<TH1D *> members;
         TH1D *rms = nullptr;
         for (const std::string &member : component.members) {
            const std::string triggerFile = XsecFile(radiusDir, jetR, trigger, member, bayes, nIter);
            const std::string baseFile = XsecFile(radiusDir, jetR, base, member, bayes, nIter);
            if (component.rms) {
               // the replicas of the two triggers are drawn with the same seed, so their RMS values are
               // propagated as if uncorrelated rather than differenced
               TH1D *triggerRms = Get(triggerFile, "canonical_embstat");
               TH1D *baseRms = Get(baseFile, "canonical_embstat");
               rms = RatioWithErrors(triggerRms, baseRms, "rms");
            } else {
               members.push_back(RatioWithErrors(Get(triggerFile, "canonical"), Get(baseFile, "canonical"),
                                                 ("m_" + member).c_str()));
            }
         }
         components.push_back({component.name, Envelope(nominal, members, rms)});
      }
      double triggerLo = 0;
      double triggerHi = 0;
      double baseLo = 0;
      double baseHi = 0;
      QuoteRange(trigger, triggerLo, triggerHi);
      QuoteRange(base, baseLo, baseHi);
      const std::string stem = "ratio_" + trigger + "_over_" + base + "_R" + jetR + MethodTag(bayes, nIter);
      WriteOut(radiusDir, stem, nominal, components, "ratio " + trigger + " / " + base + ", " + method,
               std::max(triggerLo, baseLo), std::min(triggerHi, baseHi));
   }
}

// ratio of one mode at radius R to the same mode at the base radius, every member recomputed on both radii
// (the energy-scale members largely cancel)
// output results/R<R>/rratio_<mode>_R<R>_over_R<base>[_bayes<n>].{root,txt}
void ratios_radii(const char *mode = "jp1_e23", const char *radii = "0.2 0.3 0.4", const char *baseR = "0.5",
                  const char *method = "bayes", int nIter = 4)
{
   const bool bayes = std::string(method) == "bayes";
   const std::vector<std::string> radiusList = SplitWords(radii);
   const std::string baseDir = RadiusDir(baseR);
   TH1D *baseNominal = Get(XsecFile(baseDir, baseR, mode, "", bayes, nIter), "canonical");

   for (const std::string &radius : radiusList) {
      const std::string radiusDir = RadiusDir(radius.c_str());
      TH1D *nominal = RatioWithErrors(Get(XsecFile(radiusDir, radius.c_str(), mode, "", bayes, nIter), "canonical"),
                                      baseNominal, "canonical");
      std::vector<NamedBand> components;
      for (const Syst::Component &component : Syst::Components(bayes)) {
         std::vector<TH1D *> members;
         TH1D *rms = nullptr;
         const char *name = component.rms ? "canonical_embstat" : "canonical";
         for (const std::string &member : component.members) {
            TH1D *atRadius = Get(XsecFile(radiusDir, radius.c_str(), mode, member, bayes, nIter), name, false);
            TH1D *atBase = Get(XsecFile(baseDir, baseR, mode, member, bayes, nIter), name, false);
            if (!atRadius || !atBase) {
               printf("[syst] WARNING member %s missing at R=%s or R=%s - component %s incomplete\n",
                      member.c_str(), radius.c_str(), baseR, component.name.c_str());
               continue;
            }
            if (component.rms) rms = RatioWithErrors(atRadius, atBase, "rms");
            else members.push_back(RatioWithErrors(atRadius, atBase, ("m_" + member).c_str()));
         }
         components.push_back({component.name, Envelope(nominal, members, rms)});
      }
      double quoteLo = 0;
      double quoteHi = 0;
      QuoteRange(mode, quoteLo, quoteHi);
      const std::string stem =
         std::string("rratio_") + mode + "_R" + radius + "_over_R" + baseR + MethodTag(bayes, nIter);
      const std::string title = std::string("ratio R=") + radius + " / R=" + baseR + ", " + mode + ", " + method;
      WriteOut(radiusDir, stem, nominal, components, title, quoteLo, quoteHi);
   }
}
