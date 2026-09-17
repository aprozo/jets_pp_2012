// plot_ladder.C — the ladder that leads from the reproduction of the published cross section to the
// analysis result: every rung is a ratio to the published Table III, the reproduction is also shown
// as a spectrum, and two figures compare the statistical errors and the bin-to-bin correlations of
// matrix inversion and Bayesian unfolding.
// Inputs: results/R<R>/xsec_inversion_R<R>*.root, all written by "published/run.sh ladder"; each of
// them carries the published reference ("reference", "reference_syst") next to our "canonical".
// Output: results/figures/published_ladder_<rung>_R<R>.pdf, copied to docs/talk/ladder_<rung>.pdf,
// which is the name docs/talk/ladder.tex includes. Every rung also prints its ratios bin by bin.
// A rung whose result files are absent is skipped with a one-line message.
//   root -l -b -q 'plot_ladder.C+("0.5")'
#include "../plot_style.h"

#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TMatrixD.h>

#include <string>
#include <vector>

using namespace PlotStyle;

// the pT window of the ladder figures: the published bins reach 52 GeV
const double kLadderPtLo = 5.0;
const double kLadderPtHi = 52.0;

// the deck palette in the roles of this figure set
const int kColReproduction = SeriesColor(0); // matrix inversion on the 2021 embedding
const int kColVariant1 = SeriesColor(1);     // the first alternative of a rung
const int kColVariant2 = SeriesColor(2);     // the second
const int kColVariant3 = SeriesColor(3);     // the third
const int kColVariant4 = SeriesColor(4);     // the fourth (only the error figure needs so many)
const int kColVariant5 = SeriesColor(5);     // the fifth
const int kColPublished = kBlack;            // the published points and the unity line

// marker styles: the filled and open variants that keep overlapping series apart
const int kMkCircle = 20;
const int kMkSquare = 21;
const int kMkDiamond = 23;
const int kMkCircleOpen = 24;
const int kMkSquareOpen = 25;
const int kMkTriangleOpen = 26;

// a legend entry up to this many characters still fits two to a row
const size_t kShortLabel = 40;

// the two embedding samples, as the talk names them
const char *kEmb2021 = "2021 embedding (Cold QCD request)";
const char *kEmb2023 = "2023 embedding (HP request)";

// one curve of a rung: which result file it comes from and how it is drawn
struct Series {
   std::string file;  // result file name, inside results/R<R>/
   std::string label; // legend entry
   int marker;        // marker style
   int color;         // marker and line colour
   double shift;      // horizontal offset of the points [GeV], so that series do not overlap
};

// one curve of the error figure: which result file and which of its histograms
struct ErrorSeries {
   std::string file;  // result file name, inside results/R<R>/
   const char *hist;  // histogram to read: the result, or the publication's own error estimator
   std::string label; // legend entry
   int marker;        // marker style
   int color;         // marker and line colour
};

// a ratio of two of our own results, without any published reference
struct DirectRatio {
   std::string numerator;   // result file of the numerator
   std::string denominator; // result file of the denominator
   std::string label;       // legend entry
   int marker;              // marker style
   int color;               // marker and line colour
   double shift;            // horizontal offset of the points [GeV]
};

// a canvas whose legend lives in a strip above the plot, where it can never cover a point
struct StripCanvas {
   TCanvas *canvas;  // the canvas itself
   TPad *plot;       // the pad the figure is drawn in (current after MakeStripCanvas)
   TPad *strip;      // the narrow pad above it that holds the legend
   TLegend *legend;  // the legend, drawn by DrawLegendStrip
};

// the result file of one unfolding variant: "" is the reproduction, "_bayes4" a Bayes run, ...
static std::string LadderFile(const std::string &jetR, const std::string &variant)
{
   return "xsec_inversion_R" + jetR + variant + ".root";
}

// one histogram of a ladder result file
static TH1D *Load(const std::string &jetR, const std::string &file, const char *hist)
{
   return Take(RadiusDir(jetR) + file, hist, hist);
}

// label and title sizes for the talk-shaped pads of this macro
static void PadAxes(TH1 *h, double titleSize = 0.062, double labelSize = 0.052)
{
   h->GetXaxis()->SetTitleSize(titleSize);
   h->GetYaxis()->SetTitleSize(titleSize);
   h->GetXaxis()->SetLabelSize(labelSize);
   h->GetYaxis()->SetLabelSize(labelSize);
   h->GetXaxis()->SetTitleOffset(1.0);
   h->GetYaxis()->SetTitleOffset(0.85);
   h->GetXaxis()->SetNdivisions(510);
   h->GetYaxis()->SetNdivisions(505);
}

// Write one figure to results/figures/ and copy it to docs/talk/ under the short name that
// docs/talk/ladder.tex includes.
static void SaveRung(TCanvas *c, const std::string &name, const std::string &jetR)
{
   Save(c, "published_ladder_" + name + "_R" + jetR);
   const std::string figure = FigDir() + "published_ladder_" + name + "_R" + jetR + ".pdf";
   const std::string talk = kWorkDir + "../docs/talk/ladder_" + name + ".pdf";
   gSystem->Exec(("cp -f " + figure + " " + talk).c_str());
}

// A canvas with a legend strip on top and the plot below. The strip grows with the number of
// legend rows, so that the plot never has to make room for the legend inside its frame.
static StripCanvas MakeStripCanvas(const std::string &name, int nEntries, int nColumns = 2)
{
   gStyle->SetPaperSize(24.0, 15.0);
   StripCanvas sc;
   sc.canvas = new TCanvas(("c_" + name).c_str(), "", 900, 560);
   const int rows = (nEntries + nColumns - 1) / nColumns;
   const double stripHeight = 0.068 * rows + 0.03;
   sc.strip = new TPad(("pl_" + name).c_str(), "", 0, 1 - stripHeight, 1, 1);
   sc.strip->Draw();
   sc.plot = new TPad(("pp_" + name).c_str(), "", 0, 0, 1, 1 - stripHeight);
   sc.plot->Draw();
   sc.plot->SetLeftMargin(0.10);
   sc.plot->SetRightMargin(0.03);
   sc.plot->SetBottomMargin(0.15);
   sc.plot->SetTopMargin(0.03);
   sc.legend = new TLegend(0.08, 0.02, 0.98, 0.98);
   sc.legend->SetNColumns(nColumns);
   sc.legend->SetBorderSize(0);
   sc.legend->SetFillStyle(0);
   sc.legend->SetTextSize(std::min(0.42, 0.80 / rows));
   sc.legend->SetMargin(0.10);
   sc.plot->cd();
   return sc;
}

// draw the legend into its strip and return to the canvas
static void DrawLegendStrip(const StripCanvas &sc)
{
   sc.strip->cd();
   sc.legend->Draw();
   sc.canvas->cd();
}

// Ratio of h to the published reference. A bin wider than the reference bins — the deck merges the
// bins above 22.5 GeV — is compared with the width-weighted mean of the reference over its range.
static TGraphErrors *RatioToPublished(const TH1D *h, const TH1D *reference, double xshift, double ptMin, double ptMax)
{
   auto *g = new TGraphErrors();
   if (!h || !reference) return g;
   for (int b = 1; b <= h->GetNbinsX(); ++b) {
      const double lo = h->GetXaxis()->GetBinLowEdge(b);
      const double hi = h->GetXaxis()->GetBinUpEdge(b);
      const double x = 0.5 * (lo + hi);
      if (x < ptMin || x > ptMax || h->GetBinContent(b) == 0) continue;
      double weighted = 0;
      double width = 0;
      for (int r = 1; r <= reference->GetNbinsX(); ++r) {
         const double refLo = reference->GetXaxis()->GetBinLowEdge(r);
         const double refHi = reference->GetXaxis()->GetBinUpEdge(r);
         const double overlap = std::min(hi, refHi) - std::max(lo, refLo);
         if (overlap <= 1e-6) continue;
         weighted += reference->GetBinContent(r) * overlap;
         width += overlap;
      }
      if (width <= 0 || weighted <= 0) continue;
      const double refValue = weighted / width;
      g->SetPoint(g->GetN(), x + xshift, h->GetBinContent(b) / refValue);
      g->SetPointError(g->GetN() - 1, 0, h->GetBinError(b) / refValue);
   }
   return g;
}

// the relative statistical error of a result, in percent, as a graph over pT
static TGraphErrors *RelativeError(const TH1D *h, int color, int marker, int width)
{
   auto *g = new TGraphErrors();
   for (int b = 1; b <= h->GetNbinsX(); ++b) {
      if (h->GetBinContent(b) <= 0) continue;
      g->SetPoint(g->GetN(), h->GetBinCenter(b), 100 * h->GetBinError(b) / h->GetBinContent(b));
   }
   g->SetMarkerStyle(marker);
   g->SetMarkerColor(color);
   g->SetLineColor(color);
   g->SetMarkerSize(kBigMarkerSize);
   g->SetLineWidth(width);
   return g;
}

// the deck marker attributes of one ladder series
static void StyleSeries(TGraphErrors *g, int marker, int color, double size = kStarMarkerSize)
{
   g->SetMarkerStyle(marker);
   g->SetMarkerColor(color);
   g->SetLineColor(color);
   g->SetMarkerSize(size);
   g->SetLineWidth(kLineWidth);
}

// the ratios of one series, bin by bin, for the talk's tables
static void PrintSeries(const std::string &rung, const std::string &label, const TGraphErrors *g)
{
   printf("%-8s %-52s", rung.c_str(), label.c_str());
   for (int k = 0; k < g->GetN(); ++k) printf(" %.3f", g->GetY()[k]);
   printf("\n");
}

// The uncertainty of the published reference around unity: the systematic as the violet band of the
// deck, the statistical as a narrower grey one. Both are drawn behind the points of the rung.
static void ReferenceBands(const TH1D *reference, const TH1D *referenceSyst, const std::string &name, TLegend *leg)
{
   if (referenceSyst) {
      TH1D *syst = (TH1D *)referenceSyst->Clone(("bsys_" + name).c_str());
      syst->Reset();
      syst->SetDirectory(nullptr);
      for (int b = 1; b <= referenceSyst->GetNbinsX(); ++b) {
         const double value = referenceSyst->GetBinContent(b);
         syst->SetBinContent(b, 1.0);
         syst->SetBinError(b, value > 0 ? referenceSyst->GetBinError(b) / value : 0);
      }
      Band(syst, kRefColor, kRefAlpha);
      syst->SetLineColor(0);
      syst->Draw("E2 SAME");
      leg->AddEntry(syst, "published, syst. unc.", "f");
   }
   TH1D *stat = (TH1D *)reference->Clone(("band_" + name).c_str());
   stat->Reset();
   stat->SetDirectory(nullptr);
   for (int b = 1; b <= reference->GetNbinsX(); ++b) {
      const double value = reference->GetBinContent(b);
      stat->SetBinContent(b, 1.0);
      stat->SetBinError(b, value > 0 ? reference->GetBinError(b) / value : 0);
   }
   Band(stat, kGray + 2, 0.45);
   stat->SetLineColor(0);
   stat->Draw("E2 SAME");
   leg->AddEntry(stat, "published, stat. unc.", "f");
}

// The number of legend columns that still lets every entry be read: two short labels fit side by
// side, one long one needs the whole width.
static int LegendColumns(const std::vector<Series> &series)
{
   size_t longest = 0;
   for (const Series &s : series) longest = std::max(longest, s.label.size());
   return longest > kShortLabel ? 1 : 2;
}

// One rung: every series of the list as a ratio to the published Table III, over the published
// uncertainty bands, with the legend in the strip above the plot. The ratios are also printed.
static void Rung(const std::string &jetR, const std::string &name, const std::vector<Series> &series, double ylo,
                 double yhi, double ptMin = kLadderPtLo, double ptMax = kLadderPtHi,
                 const char *referenceFile = nullptr)
{
   const std::string file = referenceFile ? referenceFile : LadderFile(jetR, "");
   TH1D *reference = Load(jetR, file, "reference");
   TH1D *referenceSyst = Load(jetR, file, "reference_syst");
   if (!reference) return;

   StripCanvas sc = MakeStripCanvas(name, (int)series.size() + 2, LegendColumns(series));
   TH1D *frame = new TH1D(("frame_" + name).c_str(), ";jet p_{T} [GeV];this analysis / published", 100, ptMin, ptMax);
   frame->SetDirectory(nullptr);
   frame->GetYaxis()->SetRangeUser(ylo, yhi);
   PadAxes(frame);
   frame->Draw();
   ReferenceBands(reference, referenceSyst, name, sc.legend);
   Baseline(ptMin, ptMax);

   for (const Series &s : series) {
      TH1D *h = Load(jetR, s.file, "canonical");
      TGraphErrors *g = RatioToPublished(h, reference, s.shift, ptMin, ptMax);
      StyleSeries(g, s.marker, s.color);
      g->Draw("P SAME");
      sc.legend->AddEntry(g, s.label.c_str(), "lep");
      PrintSeries(name, s.label, g);
   }
   DrawLegendStrip(sc);
   SaveRung(sc.canvas, name, jetR);
}

// The reproduction as a spectrum: the published Table III and our inversion on top of each other,
// with their ratio below. This is the first rung, the one that shows the absolute normalisation.
static void Spectrum(const std::string &jetR, const std::string &name, const std::string &file, const char *label)
{
   TH1D *reference = Load(jetR, file, "reference");
   TH1D *referenceSyst = Load(jetR, file, "reference_syst");
   TH1D *ours = Load(jetR, file, "canonical");
   if (!reference || !ours) return;

   gStyle->SetPaperSize(24.0, 15.0);
   auto *c = new TCanvas(("c_" + name).c_str(), "", 900, 560);
   auto *top = new TPad("pt", "", 0, 0.42, 1, 1);
   top->Draw();
   top->SetBottomMargin(0.02);
   top->SetTopMargin(0.04);
   top->SetLeftMargin(0.10);
   top->SetRightMargin(0.03);
   top->SetLogy();
   auto *bottom = new TPad("pb", "", 0, 0, 1, 0.42);
   bottom->Draw();
   bottom->SetTopMargin(0.03);
   bottom->SetBottomMargin(0.30);
   bottom->SetLeftMargin(0.10);
   bottom->SetRightMargin(0.03);

   top->cd();
   auto *frameTop = new TH1D("ft", ";;d^{2}#sigma / dp_{T} d#eta  [pb / GeV]", 100, kLadderPtLo, kLadderPtHi);
   frameTop->SetDirectory(nullptr);
   frameTop->GetYaxis()->SetRangeUser(1.5, 3e7);
   PadAxes(frameTop, 0.075, 0.065);
   frameTop->GetXaxis()->SetLabelSize(0);
   frameTop->GetYaxis()->SetTitleOffset(0.62);
   frameTop->Draw();

   TH1D *published = (TH1D *)reference->Clone("refS");
   Marker(published, kColPublished, kMkCircle, kMarkerSize);
   published->Draw("E1 SAME");
   TH1D *reproduction = (TH1D *)ours->Clone("ours");
   Marker(reproduction, kColReproduction, kMkCircleOpen, kStarMarkerSize);
   reproduction->Draw("E1 SAME");

   TLegend *leg = Leg(0.45, 0.68, 0.96, 0.92, 0.062);
   leg->AddEntry(published, "STAR Run-12, published Table III", "lep");
   leg->AddEntry(reproduction, label, "lep");
   leg->Draw();
   TLatex caption;
   caption.SetNDC();
   caption.SetTextSize(0.062);
   caption.DrawLatex(0.14, 0.10, "pp #sqrt{s} = 200 GeV, anti-k_{T} R = 0.5, |#eta| < 0.5");

   bottom->cd();
   auto *frameBottom = new TH1D("fb", ";jet p_{T} [GeV];this / published", 100, kLadderPtLo, kLadderPtHi);
   frameBottom->SetDirectory(nullptr);
   frameBottom->GetYaxis()->SetRangeUser(0.90, 1.10);
   PadAxes(frameBottom, 0.105, 0.09);
   frameBottom->GetYaxis()->SetTitleOffset(0.45);
   frameBottom->GetXaxis()->SetTitleOffset(1.15);
   frameBottom->GetYaxis()->SetNdivisions(404);
   frameBottom->Draw();
   TLegend *unused = Leg(0.0, 0.0, 0.01, 0.01); // the ratio pad shows the bands but repeats no legend
   ReferenceBands(reference, referenceSyst, name, unused);
   Baseline(kLadderPtLo, kLadderPtHi);
   TGraphErrors *ratio = RatioToPublished(ours, reference, 0.0, kLadderPtLo, kLadderPtHi);
   StyleSeries(ratio, kMkCircleOpen, kColReproduction);
   ratio->Draw("P SAME");

   SaveRung(c, name, jetR);
}

// The relative statistical error of every unfolding setting against the published one: matrix
// inversion carries the exact variance of the level-weighted data, the Bayesian errors grow with
// the number of iterations. The published error estimator of the inversion is shown as well.
static void Errors(const std::string &jetR, const std::string &name, const std::vector<ErrorSeries> &series)
{
   TH1D *reference = Load(jetR, LadderFile(jetR, ""), "reference");
   if (!reference) return;

   StripCanvas sc = MakeStripCanvas(name, (int)series.size() + 1, 3);
   sc.plot->SetLogy();
   TH1D *frame = new TH1D(("frame_" + name).c_str(), ";jet p_{T} [GeV];statistical error [%]", 100, kLadderPtLo, kLadderPtHi);
   frame->SetDirectory(nullptr);
   frame->GetYaxis()->SetRangeUser(0.08, 40);
   PadAxes(frame);
   frame->Draw();

   TGraphErrors *publishedError = RelativeError(reference, kColPublished, kMkCircle, 3);
   publishedError->Draw("PL SAME");
   sc.legend->AddEntry(publishedError, "published", "lp");
   for (const ErrorSeries &s : series) {
      TH1D *h = Load(jetR, s.file, s.hist);
      if (!h) continue;
      TGraphErrors *g = RelativeError(h, s.color, s.marker, kLineWidth);
      g->Draw("PL SAME");
      sc.legend->AddEntry(g, s.label.c_str(), "lp");
   }
   DrawLegendStrip(sc);
   SaveRung(sc.canvas, name, jetR);
}

// The statistical correlation matrix of one result, in the published bins of |eta| < 0.5. The
// covariance of the eta blocks is stored as one matrix; block 0 is the one that is quoted, and the
// Bayes axis carries extra low bins in front of it.
static TH2D *CorrelationMatrix(const std::string &jetR, const std::string &file, const char *newName)
{
   TFile f((RadiusDir(jetR) + file).c_str());
   if (f.IsZombie()) return nullptr;
   auto *covariance = (TMatrixD *)f.Get("cov");
   auto *result = (TH1D *)f.Get("canonical");
   if (!covariance || !result) return nullptr;
   const int nBins = result->GetNbinsX();
   const int nBlockBins = covariance->GetNrows() / 2;
   const int offset = nBlockBins - nBins;

   std::vector<double> edges(nBins + 1);
   for (int b = 0; b <= nBins; ++b) edges[b] = result->GetXaxis()->GetBinLowEdge(b + 1);
   auto *matrix = new TH2D(newName, ";jet p_{T} [GeV];jet p_{T} [GeV]", nBins, edges.data(), nBins, edges.data());
   matrix->SetDirectory(nullptr);
   for (int i = 0; i < nBins; ++i) {
      for (int j = 0; j < nBins; ++j) {
         const double diagonal = (*covariance)(offset + i, offset + i) * (*covariance)(offset + j, offset + j);
         matrix->SetBinContent(i + 1, j + 1, diagonal > 0 ? (*covariance)(offset + i, offset + j) / std::sqrt(diagonal) : 0);
      }
   }
   return matrix;
}

// The correlation matrices of two results side by side: matrix inversion leaves neighbouring bins
// anti-correlated, Bayesian unfolding smooths them into a positively correlated block.
static void Correlations(const std::string &jetR, const std::string &name, const std::string &fileA, const char *labelA,
                         const std::string &fileB, const char *labelB)
{
   TH2D *matrices[2] = {CorrelationMatrix(jetR, fileA, "corrA"), CorrelationMatrix(jetR, fileB, "corrB")};
   const char *labels[2] = {labelA, labelB};
   if (!matrices[0] || !matrices[1]) {
      printf("[plot] no covariance in %s or %s\n", fileA.c_str(), fileB.c_str());
      return;
   }
   gStyle->SetPaperSize(24.0, 12.0);
   gStyle->SetPalette(kTemperatureMap);
   gStyle->SetNumberContours(40);
   auto *c = new TCanvas(("c_" + name).c_str(), "", 1000, 500);
   c->Divide(2, 1, 0.005, 0.005);
   for (int i = 0; i < 2; ++i) {
      c->cd(i + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetRightMargin(0.16);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.10);
      matrices[i]->GetZaxis()->SetRangeUser(-1, 1);
      PadAxes(matrices[i], 0.055, 0.045);
      matrices[i]->GetYaxis()->SetTitleOffset(1.15);
      matrices[i]->GetZaxis()->SetLabelSize(0.045);
      matrices[i]->Draw("COLZ");
      TLatex caption;
      caption.SetNDC();
      caption.SetTextSize(0.055);
      caption.DrawLatex(0.13, 0.93, labels[i]);
   }
   SaveRung(c, name, jetR);
}

// The ratio of two of our own results, without the published reference: the two share the data, so
// the ratio isolates one ingredient (the method at fixed embedding, or the embedding at fixed
// method) and carries no meaningful error bar.
static void Direct(const std::string &jetR, const std::string &name, const std::vector<DirectRatio> &ratios, double ylo,
                   double yhi, const char *ytitle)
{
   StripCanvas sc = MakeStripCanvas(name, (int)ratios.size(), 1);
   TH1D *frame = new TH1D(("frame_" + name).c_str(), (std::string(";jet p_{T} [GeV];") + ytitle).c_str(), 100,
                          kLadderPtLo, kLadderPtHi);
   frame->SetDirectory(nullptr);
   frame->GetYaxis()->SetRangeUser(ylo, yhi);
   PadAxes(frame);
   frame->Draw();
   Baseline(kLadderPtLo, kLadderPtHi);

   for (const DirectRatio &r : ratios) {
      TH1D *numerator = Load(jetR, r.numerator, "canonical");
      TH1D *denominator = Load(jetR, r.denominator, "canonical");
      if (!numerator || !denominator) continue;
      auto *g = new TGraphErrors();
      for (int k = 1; k <= numerator->GetNbinsX(); ++k) {
         if (denominator->GetBinContent(k) <= 0) continue;
         g->SetPoint(g->GetN(), numerator->GetBinCenter(k) + r.shift,
                     numerator->GetBinContent(k) / denominator->GetBinContent(k));
      }
      StyleSeries(g, r.marker, r.color);
      g->Draw("PL SAME");
      sc.legend->AddEntry(g, r.label.c_str(), "lp");
      PrintSeries(name, r.label, g);
   }
   DrawLegendStrip(sc);
   SaveRung(sc.canvas, name, jetR);
}

void plot_ladder(const char *jetR = "0.5")
{
   Init();
   const std::string R = jetR;

   // rung 1: the reproduction against the published spectrum
   Spectrum(R, "rung1", LadderFile(R, ""), "this analysis: inversion, 2021 embedding");

   // rung 2: the unfolding method at fixed embedding
   Rung(R, "rung2",
        {{LadderFile(R, ""), "inversion", kMkCircle, kColReproduction, -0.45},
         {LadderFile(R, "_bayes2"), "Bayes, 2 iterations", kMkSquare, kColVariant1, -0.15},
         {LadderFile(R, "_bayes4"), "Bayes, 4 iterations", kMkSquareOpen, kColVariant2, 0.15},
         {LadderFile(R, "_bayes10"), "Bayes, 10 iterations", kMkTriangleOpen, kColVariant3, 0.45}},
        0.90, 1.10);

   // the statistical error of every method, including the publication's own estimator
   Errors(R, "errors",
          {{LadderFile(R, ""), "canonical", "inversion", kMkCircle, kColReproduction},
           {LadderFile(R, ""), "canonical_pubErr", "inversion, published estimator", kMkCircleOpen, kColReproduction},
           {LadderFile(R, "_bayes2"), "canonical", "Bayes 2 it.", kMkSquare, kColVariant1},
           {LadderFile(R, "_bayes4"), "canonical", "Bayes 4 it.", kMkSquareOpen, kColVariant2},
           {LadderFile(R, "_bayes10"), "canonical", "Bayes 10 it.", kMkTriangleOpen, kColVariant3},
           {LadderFile(R, "_bayes30"), "canonical", "Bayes 30 it.", kMkDiamond, kColVariant4},
           {LadderFile(R, "_bayes100"), "canonical", "Bayes 100 it.", kMkCircleOpen, kColVariant5}});

   // the bin-to-bin correlations the two methods leave behind
   Correlations(R, "corr", LadderFile(R, ""), "inversion", LadderFile(R, "_bayes4"), "Bayes, 4 iterations");

   // the two ingredients on their own: the method at fixed embedding, the embedding at fixed method
   Direct(R, "method",
          {{LadderFile(R, "_bayes4"), LadderFile(R, ""), "Bayes 4 it. / inversion, 2021 embedding", kMkSquare,
            kColReproduction, -0.25},
           {LadderFile(R, "_jp_e23_bayes4"), LadderFile(R, "_jp_e23"), "Bayes 4 it. / inversion, 2023 embedding",
            kMkSquareOpen, kColVariant1, 0.25}},
          0.90, 1.10, "Bayes / inversion");
   Direct(R, "embedding",
          {{LadderFile(R, "_jp_e23"), LadderFile(R, ""), "2023 / 2021, inversion", kMkCircle, kColReproduction, -0.25},
           {LadderFile(R, "_jp_e23_bayes4"), LadderFile(R, "_bayes4"), "2023 / 2021, Bayes 4 it.", kMkSquareOpen,
            kColVariant1, 0.25}},
          0.90, 1.15, "2023 / 2021");

   // rung 3: the embedding sample, at fixed method, with and without the energy-scale calibration
   Rung(R, "rung3",
        {{LadderFile(R, ""), std::string("inversion, ") + kEmb2021, kMkCircle, kColReproduction, -0.3},
         {LadderFile(R, "_jp_e23"), std::string("inversion, ") + kEmb2023, kMkSquare, kColVariant1, 0.0},
         {LadderFile(R, "_jp_e23cal"), "inversion, 2023, energy scale calibrated to 2021", kMkSquareOpen, kColVariant1,
          0.3}},
        0.80, 1.20);
   Rung(R, "rung3b",
        {{LadderFile(R, "_bayes4"), std::string("Bayes 4 it., ") + kEmb2021, kMkSquareOpen, kColVariant2, -0.3},
         {LadderFile(R, "_jp_e23_bayes4"), std::string("Bayes 4 it., ") + kEmb2023, kMkSquare, kColVariant1, 0.0},
         {LadderFile(R, "_jp_e23cal_bayes4"), "Bayes 4 it., 2023, energy scale calibrated to 2021", kMkTriangleOpen,
          kColVariant1, 0.3}},
        0.80, 1.20);

   // rung 4: the reproduction against the analysis choice (Bayes on the 2023 embedding)
   Rung(R, "rung4",
        {{LadderFile(R, ""), std::string("inversion, ") + kEmb2021, kMkCircle, kColReproduction, -0.3},
         {LadderFile(R, "_jp_e23_bayes4"), std::string("Bayes 4 it., ") + kEmb2023, kMkSquare, kColVariant1, 0.0},
         {LadderFile(R, "_jp_e23cal_bayes4"), "Bayes 4 it., 2023, energy scale calibrated", kMkSquareOpen, kColVariant1,
          0.3}},
        0.85, 1.15);

   // rung 5: one jet-patch trigger on its own instead of the exclusive level construction
   Rung(R, "rung5a",
        {{LadderFile(R, ""), std::string("all levels, inversion, ") + kEmb2021, kMkCircle, kColReproduction, -0.3},
         {LadderFile(R, "_jp1_e23"), "JP1 alone, inversion, 2023", kMkSquare, kColVariant1, 0.0},
         {LadderFile(R, "_jp1_e23_bayes4"), "JP1 alone, Bayes 4 it., 2023", kMkSquareOpen, kColVariant4, 0.3}},
        0.75, 1.25, kLadderPtLo, 46.);
   Rung(R, "rung5b",
        {{LadderFile(R, ""), std::string("all levels, inversion, ") + kEmb2021, kMkCircle, kColReproduction, -0.3},
         {LadderFile(R, "_jp2_e23"), "JP2 alone, inversion, 2023", kMkSquare, kColVariant1, 0.0},
         {LadderFile(R, "_jp2_e23_bayes4"), "JP2 alone, Bayes 4 it., 2023", kMkSquareOpen, kColVariant4, 0.3}},
        0.75, 1.25, kLadderPtLo, 46.);

   // rung 6: the min-bias level, the one stream that shares nothing with the jet-patch triggers
   Rung(R, "rung6",
        {{LadderFile(R, ""), std::string("all levels, inversion, ") + kEmb2021, kMkCircle, kColReproduction, -0.3},
         {LadderFile(R, "_mb2021m"), "min-bias, inversion, 2021", kMkSquare, kColVariant1, -0.1},
         {LadderFile(R, "_mb2021m_bayes4"), "min-bias, Bayes 4 it., 2021", kMkSquareOpen, kColVariant1, 0.1},
         {LadderFile(R, "_mb2023_bayes4"), "min-bias, Bayes 4 it., 2023", kMkTriangleOpen, kColVariant4, 0.3}},
        0.70, 1.40, kLadderPtLo, 27.);
}
