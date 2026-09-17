// plot.C — the figures of the published-analysis reproduction: the cross section against Table III,
// the ladder from the reproduction to the analysis result, and the systematic band against Table I.
// Inputs: results/R<R>/xsec_inversion_R<R>*.root and results/R<R>/xsec_jp_R<R>_syst.root, both
// written by published/run.sh, and the published table in inputs/.
// A figure whose input files are absent is skipped with a one-line message.
//   root -l -b -q 'plot.C+("0.5")'      -> results/figures/published_*.pdf
#include "../plot_style.h"

#include <string>
#include <vector>

using namespace PlotStyle;

// the number of quoted pT bins, i.e. the length of the published Table I
const int kTableIBins = 12;

// Table I: the total systematic uncertainty and its energy-scale part, in percent, bin by bin
const double kPublishedTotalPct[kTableIBins] = {17.1, 25.4, 12.5, 17.3, 10.7, 8.8, 13.0, 9.9, 12.3, 12.9, 16.2, 26.4};
const double kPublishedScalePct[kTableIBins] = {15.9, 24.6, 10.3, 16.6, 9.8, 7.9, 12.5, 9.2, 12.0, 12.6, 14.4, 25.3};

// one component of the systematic band, as systematics.C names it inside the syst file
struct Component {
   const char *name;  // histogram prefix: syst_<name>_up and syst_<name>_dn
   const char *label; // legend entry
   int color;         // line colour
};
const std::vector<Component> kComponents = {{"energy_scale", "energy scale", kAzure},
                                            {"tracking", "tracking", kRed},
                                            {"embedding_statistics", "embedding statistics", kGreen - 2},
                                            {"underlying_event", "underlying event", kOrange + 7}};

// the larger of the two deviations, in percent of the nominal cross section, bin by bin
static TH1D *PercentOfNominal(const TH1D *nominal, const TH1D *up, const TH1D *dn, const char *newName)
{
   TH1D *pct = (TH1D *)nominal->Clone(newName);
   pct->SetDirectory(nullptr);
   pct->Reset();
   pct->GetXaxis()->SetRangeUser(kPtLo, kPtHi); // only the quoted bins are shown
   for (int k = 1; k <= kTableIBins; ++k) {
      const double v = nominal->GetBinContent(k);
      if (v <= 0) continue;
      pct->SetBinContent(k, 100 * std::max(up->GetBinContent(k), dn->GetBinContent(k)) / v);
   }
   return pct;
}

// the twelve Table I percentages in the binning of the figure, as points
static TH1D *PublishedPercent(const TH1D *shape, const double values[kTableIBins], const char *newName,
                              int color, int marker)
{
   TH1D *h = (TH1D *)shape->Clone(newName);
   h->SetDirectory(nullptr);
   h->Reset();
   for (int k = 1; k <= kTableIBins; ++k) h->SetBinContent(k, values[k - 1]);
   h->SetMarkerStyle(marker);
   h->SetMarkerColor(color);
   h->SetMarkerSize(kBigMarkerSize);
   h->SetLineColor(color);
   return h;
}

// Figure published_tableIII_R<R>: the reproduced cross section — exclusive jet-patch levels, matrix
// inversion, 2021 embedding — on the published Table III, with the ratio to it below. The violet
// band is the published systematic uncertainty, the points are this analysis.
// Reads results/R<R>/xsec_inversion_R<R>.root ("canonical")
// and inputs/jet_cross_section_publishedR0.5.root (Table III, statistical and systematic).
static void TableIII(const char *jetR)
{
   TH1D *ours = Take(RadiusDir(jetR) + "xsec_inversion_R" + jetR + ".root", "canonical", "inv");
   TH1D *reference = Take(Published(), "crossSection_statistic", "ref");
   TH1D *band = RefBand(false, "refband");
   TH1D *ratioBand = RefBand(true, "refrel");
   if (!ours || !reference || !band || !ratioBand) return;

   TPad *top = nullptr;
   TPad *bottom = nullptr;
   TCanvas *c = TwoPads("tableIII", top, bottom);

   top->cd();
   top->SetLogy();
   band->SetTitle("Inclusive jet cross section, R = 0.5, |#eta| < 0.5;"
                  " p_{T} [GeV/c]; d^{2}#sigma / dp_{T} d#eta [pb/(GeV/c)]");
   Axes(band, 0.045, 0.05, true);
   band->GetXaxis()->SetRangeUser(kPtLo, kPtHi);
   band->Draw("E2");
   Points(ours, "JPX");
   ours->Draw("P SAME");
   TLegend *leg = Leg(0.40, 0.70, 0.94, 0.90);
   leg->AddEntry(band, "published Table III (stat., syst.)", "f");
   leg->AddEntry(ours, "this analysis: inversion, 2021 embedding", "p");
   leg->Draw();

   bottom->cd();
   ratioBand->SetTitle("; p_{T} [GeV/c]; this / published");
   Axes(ratioBand, 0.055, 0.06);
   ratioBand->GetXaxis()->SetRangeUser(kPtLo, kPtHi);
   ratioBand->GetYaxis()->SetRangeUser(0.8, 1.2);
   ratioBand->GetYaxis()->SetNdivisions(505);
   ratioBand->Draw("E2");
   Guides(kPtLo, kPtHi);
   TH1D *ratio = Ratio(ours, reference, "r", true);
   Points(ratio, "JPX");
   ratio->Draw("P SAME");

   Save(c, std::string("published_tableIII_R") + jetR);
}

// Figure published_ladder_R<R>: the three steps from the reproduction to the analysis result, each
// as a ratio to Table III — matrix inversion on the 2021 embedding, Bayes with four iterations on
// the same embedding, and Bayes on the 2023 embedding. The violet band is the published systematic.
// Reads results/R<R>/xsec_inversion_R<R>.root, _bayes4.root and _jp_e23_bayes4.root ("canonical")
// and inputs/jet_cross_section_publishedR0.5.root.
static void Ladder(const char *jetR)
{
   const std::string dir = RadiusDir(jetR);
   TH1D *reference = Take(Published(), "crossSection_statistic", "ref2");
   TH1D *ratioBand = RefBand(true, "refrel2");
   TH1D *inversion = Take(dir + "xsec_inversion_R" + jetR + ".root", "canonical", "l1");
   TH1D *bayes2021 = Take(dir + "xsec_inversion_R" + jetR + "_bayes4.root", "canonical", "l2");
   TH1D *bayes2023 = Take(dir + "xsec_inversion_R" + jetR + "_jp_e23_bayes4.root", "canonical", "l4");
   if (!reference || !ratioBand || !inversion || !bayes2021 || !bayes2023) return;

   TCanvas *c = OnePad("ladder");
   ratioBand->SetTitle("The ladder: exclusive JP levels, ratio to Table III; p_{T} [GeV/c]; this / published");
   Axes(ratioBand);
   ratioBand->GetXaxis()->SetRangeUser(kPtLo, kPtHi);
   ratioBand->GetYaxis()->SetRangeUser(0.8, 1.2);
   ratioBand->Draw("E2");
   Guides(kPtLo, kPtHi);

   TH1D *rInversion = Ratio(inversion, reference, "r1", true);
   TH1D *rBayes2021 = Ratio(bayes2021, reference, "r2", true);
   TH1D *rBayes2023 = Ratio(bayes2023, reference, "r4", true);
   Points(rInversion, "JPX");
   Marker(rBayes2021, kAzure, 20, kMarkerSize);
   Marker(rBayes2023, kGreen - 2, 33, kStarMarkerSize);
   rInversion->Draw("P SAME");
   rBayes2021->Draw("P SAME");
   rBayes2023->Draw("P SAME");

   TLegend *leg = Leg(0.14, 0.66, 0.70, 0.90);
   leg->AddEntry(ratioBand, "published systematic band", "f");
   leg->AddEntry(rInversion, "inversion, 2021 embedding (the reproduction)", "p");
   leg->AddEntry(rBayes2021, "Bayes 4 it., 2021 embedding", "p");
   leg->AddEntry(rBayes2023, "Bayes 4 it., 2023 embedding", "p");
   leg->Draw();

   Save(c, std::string("published_ladder_R") + jetR);
}

// Figure published_syst_R<R>: the systematic band of the reproduction against the published one.
// Our total and its four components are lines in percent of the cross section; the published total
// and its energy-scale part (Table I) are points, over the twelve quoted bins.
// Reads results/R<R>/xsec_jp_R<R>_syst.root ("canonical", "syst_up", "syst_dn" and the
// "syst_<component>_up" / "_dn" pairs); the published percentages are the tables above.
static void Syst(const char *jetR)
{
   const std::string file = RadiusDir(jetR) + "xsec_jp_R" + jetR + "_syst.root";
   TH1D *nominal = Take(file, "canonical", "snom");
   TH1D *up = Take(file, "syst_up", "sup");
   TH1D *dn = Take(file, "syst_dn", "sdn");
   if (!nominal || !up || !dn) return;

   TCanvas *c = OnePad("syst");
   TH1D *total = PercentOfNominal(nominal, up, dn, "tot");
   total->SetTitle("Systematic band of the inversion (2021 embedding) against the published Table I;"
                   " p_{T} [GeV/c]; systematic uncertainty [%]");
   Axes(total);
   total->GetXaxis()->SetRangeUser(kPtLo, kPtHi);
   total->GetYaxis()->SetRangeUser(0, 35);
   Line(total, kBlack, 3);
   total->Draw("HIST");

   TLegend *leg = Leg(0.14, 0.58, 0.60, 0.90, 0.034);
   leg->AddEntry(total, "total, this analysis", "l");

   TH1D *publishedTotal = PublishedPercent(total, kPublishedTotalPct, "pubtot", kRefColor, 20);
   publishedTotal->Draw("P SAME");
   leg->AddEntry(publishedTotal, "total, published Table I", "p");

   for (const Component &comp : kComponents) {
      const std::string prefix = std::string("syst_") + comp.name;
      TH1D *componentUp = Take(file, (prefix + "_up").c_str(), (std::string("cu") + comp.name).c_str());
      TH1D *componentDn = Take(file, (prefix + "_dn").c_str(), (std::string("cd") + comp.name).c_str());
      if (!componentUp || !componentDn) continue;
      TH1D *pct = PercentOfNominal(nominal, componentUp, componentDn, comp.name);
      Line(pct, comp.color);
      pct->Draw("HIST SAME");
      leg->AddEntry(pct, comp.label, "l");
   }

   TH1D *publishedScale = PublishedPercent(total, kPublishedScalePct, "pubes", kAzure, 24);
   publishedScale->Draw("P SAME");
   leg->AddEntry(publishedScale, "energy scale, published Table I", "p");
   leg->Draw();

   Save(c, std::string("published_syst_R") + jetR);
}

void plot(const char *jetR = "0.5")
{
   Init();
   TableIII(jetR);
   Ladder(jetR);
   Syst(jetR);
}
