// plot_style.h — one visual language for every figure of the deck.
//
// The rules: a LINEAR jet pT axis from kPtLo to kPtHi, the published Table III as a violet band,
// data as points in the colour and marker of their trigger, simulation and systematic uncertainties
// as bands, no grid, ticks on all four sides, ratio pads with guides at 1 and at 0.9 / 1.1.
// Every figure goes to FigDir() as <deck>_<figure>[_R<R>].pdf.
//
// Used by: published/plot.C, published/plot_ladder.C, published/embedding_diff_plot.C, analysis/plot.C.
// Those macros read the Stage-2 result files under Results(); every path comes from config.h.
//
//   #include "../plot_style.h"
//   using namespace PlotStyle;
#ifndef PLOT_STYLE_H
#define PLOT_STYLE_H

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPolyLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "config.h"

namespace PlotStyle {
using namespace CrossSectionConfig;

// ---------------------------------------------------------------------------------------------
// Constants of the deck
// ---------------------------------------------------------------------------------------------

const double kPtLo = 5.0;    // left edge of the jet pT axis [GeV/c]
const double kPtHi = 55.0;   // right edge of the jet pT axis [GeV/c]
const double kGuideLo = 0.9; // lower dashed guide of a ratio pad
const double kGuideHi = 1.1; // upper dashed guide of a ratio pad

const int kRefColor = kViolet; // the published analysis: violet everywhere

const double kRefAlpha = 0.20;   // published systematic band
const double kSystAlpha = 0.25;  // systematic band of our own data points
const double kSoftAlpha = 0.30;  // legend swatch of a point series, and the JP1 reference band
const double kBandAlpha = 0.35;  // simulation bands
const double kFaintAlpha = 0.15; // a simulation band that only guides the eye

const double kMarkerSize = 1.3;     // round and square markers
const double kBigMarkerSize = 1.4;  // open and triangular markers
const double kStarMarkerSize = 1.6; // diamonds and stars, which read smaller at equal size

const int kLineWidth = 2;

// ---------------------------------------------------------------------------------------------
// Paths (everything through config.h: no absolute path in a macro)
// ---------------------------------------------------------------------------------------------

// where every figure of the deck is written; created on demand
inline std::string FigDir()
{
   const std::string dir = kWorkDir + "results/figures/";
   gSystem->mkdir(dir.c_str(), true);
   return dir;
}

// the root of the Stage-2 result files (results/R<R>/, results/models/, results/Michal/)
inline std::string Results()
{
   return kWorkDir + "results/";
}

// the published Table III: the cross section and its systematic uncertainty, R = 0.5
inline std::string Published()
{
   return kWorkDir + "inputs/jet_cross_section_publishedR0.5.root";
}

// the generator comparison of models/compare.C for one radius and jet definition ("ue" or "noue")
inline std::string ModelFile(const std::string &jetR, const std::string &definition)
{
   return Results() + "models/models_R" + jetR + "_" + definition + ".root";
}

// ---------------------------------------------------------------------------------------------
// Palettes: one colour and one marker per trigger, per radius, per generator, per extra series
// ---------------------------------------------------------------------------------------------

// the colour of a trigger; "JPX" (the combined jet-patch levels) and anything unknown come out magenta
inline int TrigColor(const std::string &trigger)
{
   if (trigger == "JP2") return kAzure;
   if (trigger == "HT2") return kRed;
   if (trigger == "JP1") return kOrange + 7;
   if (trigger == "JP0") return kGreen - 2;
   if (trigger == "MB") return kBlack;
   return kMagenta + 1;
}

// the marker of a trigger, in the same order as TrigColor
inline int TrigMarker(const std::string &trigger)
{
   if (trigger == "JP2") return 20;
   if (trigger == "HT2") return 21;
   if (trigger == "JP1") return 22;
   if (trigger == "JP0") return 33;
   return 34;
}

// one jet radius as it appears in a figure
struct Radius {
   const char *name; // "0.2" ... "0.5": the R<name> of the result directory
   int color;        // line and marker colour
   int marker;       // marker style
};
const std::vector<Radius> kRadii = {
   {"0.2", kGreen - 2, 33}, {"0.3", kAzure, 20}, {"0.4", kOrange + 7, 22}, {"0.5", kMagenta + 1, 34}};

// one generator of new_ana/models as it appears in a figure
struct Model {
   const char *name;  // histogram name inside the models_R<R>_<def>.root files
   const char *label; // legend entry
   int color;         // band and line colour
};
const std::vector<Model> kModels = {{"pythia6", "PYTHIA 6 Perugia 2012", kGreen + 2},
                                    {"pythia8", "PYTHIA 8 Monash", kAzure},
                                    {"pythia8detroit", "PYTHIA 8 Detroit", kRed},
                                    {"herwig", "Herwig 7.3", kGray + 2}};

// colour of the i-th curve of a figure that compares methods rather than triggers (unfolding
// iterations, embedding samples, jet observables); it wraps around
inline int SeriesColor(int i)
{
   const int palette[6] = {kAzure, kOrange + 7, kGreen - 2, kMagenta + 1, kRed, kGray + 2};
   return palette[i % 6];
}

// marker of the i-th curve of such a figure, in the same order as SeriesColor
inline int SeriesMarker(int i)
{
   const int palette[6] = {20, 21, 22, 33, 34, 24};
   return palette[i % 6];
}

// ---------------------------------------------------------------------------------------------
// Global ROOT style
// ---------------------------------------------------------------------------------------------

// the deck style: no stats box, a centred title, no grid, ticks on all four sides
inline void Init()
{
   TH1::AddDirectory(kFALSE);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleFontSize(0.045);
   gStyle->SetTitleX(0.5);
   gStyle->SetTitleAlign(23);
   gStyle->SetTitleBorderSize(0);
   gStyle->SetTitleFillColor(0);
   gStyle->SetPadGridX(0);
   gStyle->SetPadGridY(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetEndErrorSize(4);
}

// ---------------------------------------------------------------------------------------------
// Reading result files
// ---------------------------------------------------------------------------------------------

// One histogram of a result file, cloned and detached from the file; nullptr with a one-line
// message when the file or the histogram is absent, so a macro can skip a figure whose inputs
// were never produced.
inline TH1D *Take(const std::string &file, const char *name, const char *newName)
{
   TFile f(file.c_str(), "READ");
   if (f.IsZombie()) {
      printf("[plot] cannot open %s\n", file.c_str());
      return nullptr;
   }
   auto *h = (TH1D *)f.Get(name);
   if (!h) {
      printf("[plot] %s missing in %s\n", name, file.c_str());
      return nullptr;
   }
   // a histogram that was drawn before being written carries a stats box owned by the file
   h->GetListOfFunctions()->Delete();
   h = (TH1D *)h->Clone(newName);
   h->SetDirectory(nullptr);
   return h;
}

// Empty the bins outside the quoted pT range of a result. systematics.C writes that range into the
// header of the companion table ("# quoted range lo-hi GeV"); the variations follow the value.
inline void Quoted(const std::string &file, TH1D *h, TH1D *up = nullptr, TH1D *dn = nullptr)
{
   double lo = 0.0;
   double hi = 1e9;
   std::string table = file;
   table.replace(table.rfind(".root"), 5, ".txt");
   std::ifstream in(table.c_str());
   std::string line;
   while (std::getline(in, line)) {
      if (sscanf(line.c_str(), "# quoted range %lf-%lf GeV", &lo, &hi) == 2) break;
   }
   for (int k = 1; k <= h->GetNbinsX(); ++k) {
      if (h->GetBinLowEdge(k) >= lo - 1e-6 && h->GetBinLowEdge(k + 1) <= hi + 1e-6) continue;
      h->SetBinContent(k, 0);
      h->SetBinError(k, 0);
      if (up) up->SetBinContent(k, 0);
      if (dn) dn->SetBinContent(k, 0);
   }
}

// empty every bin that starts below a pT floor (a trigger turn-on, or a soft region not quoted)
inline void Below(TH1D *h, double ptFloor)
{
   for (int k = 1; k <= h->GetNbinsX(); ++k) {
      if (h->GetBinLowEdge(k) >= ptFloor - 1e-6) continue;
      h->SetBinContent(k, 0);
      h->SetBinError(k, 0);
   }
}

// ---------------------------------------------------------------------------------------------
// Attributes: points, lines, bands
// ---------------------------------------------------------------------------------------------

// a histogram as points of one colour; its fill is the legend swatch of its systematic band
inline void Marker(TH1D *h, int color, int marker, double size = kMarkerSize, double alpha = kSoftAlpha)
{
   h->SetLineColor(color);
   h->SetMarkerColor(color);
   h->SetMarkerStyle(marker);
   h->SetMarkerSize(size);
   h->SetLineWidth(kLineWidth);
   h->SetFillColorAlpha(color, alpha);
}

// a histogram as the points of one trigger: its colour and its marker
inline void Points(TH1D *h, const std::string &trigger, double size = kMarkerSize)
{
   Marker(h, TrigColor(trigger), TrigMarker(trigger), size);
}

// a histogram as a plain line without markers (percentages, efficiencies, generator curves)
inline void Line(TH1D *h, int color, int width = kLineWidth, int style = 1)
{
   h->SetLineColor(color);
   h->SetLineWidth(width);
   h->SetLineStyle(style);
   h->SetMarkerSize(0);
}

// a histogram as a filled band of half-width equal to its bin errors, to be drawn with option "E2"
inline void Band(TH1D *h, int color, double alpha = kBandAlpha)
{
   h->SetLineColor(color);
   h->SetLineWidth(kLineWidth);
   h->SetFillColorAlpha(color, alpha);
   h->SetFillStyle(1001);
   h->SetMarkerSize(0);
}

// A band drawn as pad primitives over the FILLED bins only, so that an empty bin leaves a gap:
// opt "2" is one box per bin (the systematic band of data points), opt "3" a filled polygon through
// the bin centres (a smooth simulation band). With line = true the central value is drawn on top.
// The returned box carries the attributes and is meant to be the legend entry of the band.
inline TBox *DrawBand(TH1D *h, int color, double alpha, bool line, const char *opt = "3", int lineStyle = 1)
{
   auto *key = new TBox(0, 0, 0, 0);
   key->SetFillColorAlpha(color, alpha);
   key->SetFillStyle(1001);
   key->SetLineColor(color);
   key->SetLineWidth(kLineWidth);
   key->SetLineStyle(lineStyle);

   std::vector<double> x;
   std::vector<double> upper;
   std::vector<double> lower;
   for (int k = 1; k <= h->GetNbinsX(); ++k) {
      if (h->GetBinContent(k) == 0) continue;
      const double y = h->GetBinContent(k);
      const double e = h->GetBinError(k);
      if (std::string(opt) == "2") {
         auto *box = new TBox(h->GetBinLowEdge(k), y - e, h->GetBinLowEdge(k + 1), y + e);
         box->SetFillColorAlpha(color, alpha);
         box->SetFillStyle(1001);
         box->SetLineColor(color);
         box->SetLineWidth(0);
         box->Draw("SAME");
      }
      x.push_back(h->GetBinCenter(k));
      upper.push_back(y + e);
      lower.push_back(y - e);
   }
   const int n = (int)x.size();
   if (n == 0) return key;

   if (std::string(opt) == "3") {
      auto *poly = new TPolyLine(2 * n);
      for (int i = 0; i < n; ++i) {
         poly->SetPoint(i, x[i], upper[i]);
         poly->SetPoint(2 * n - 1 - i, x[i], lower[i]);
      }
      poly->SetFillColorAlpha(color, alpha);
      poly->SetLineWidth(0);
      poly->Draw("f SAME");
   }
   if (line) {
      auto *central = new TPolyLine(n);
      for (int i = 0; i < n; ++i) central->SetPoint(i, x[i], (upper[i] + lower[i]) / 2);
      central->SetLineColor(color);
      central->SetLineWidth(kLineWidth);
      central->SetLineStyle(lineStyle);
      central->Draw("SAME");
   }
   return key;
}

// a dashed grey horizontal reference line, for a panel that carries no ratio guides
inline void Baseline(double xlo, double xhi, double y = 1.0)
{
   auto *line = new TLine(xlo, y, xhi, y);
   line->SetLineColor(kGray + 1);
   line->SetLineStyle(2);
   line->Draw("SAME");
}

// the guides of a ratio pad: unity in the colour of the published reference, +-10% dashed and grey
inline void Guides(double xlo, double xhi, double lo = kGuideLo, double hi = kGuideHi)
{
   auto *unity = new TLine(xlo, 1.0, xhi, 1.0);
   unity->SetLineColor(kRefColor);
   unity->SetLineWidth(kLineWidth);
   unity->Draw("SAME");
   for (double y : {lo, hi}) {
      auto *guide = new TLine(xlo, y, xhi, y);
      guide->SetLineColor(kGray + 1);
      guide->SetLineStyle(2);
      guide->Draw("SAME");
   }
}

// ---------------------------------------------------------------------------------------------
// Derived histograms
// ---------------------------------------------------------------------------------------------

// The ratio a / b bin by bin. With numeratorErrorOnly the denominator counts as exact, which is the
// case when our result is compared with the published one (their errors are not ours to combine).
inline TH1D *Ratio(const TH1D *a, const TH1D *b, const char *newName, bool numeratorErrorOnly = false)
{
   TH1D *r = (TH1D *)a->Clone(newName);
   r->SetDirectory(nullptr);
   r->Reset();
   for (int k = 1; k <= a->GetNbinsX(); ++k) {
      const double va = a->GetBinContent(k);
      const double vb = b->GetBinContent(k);
      if (va <= 0 || vb <= 0) continue;
      const double ea = a->GetBinError(k) / va;
      const double eb = numeratorErrorOnly ? 0.0 : b->GetBinError(k) / vb;
      r->SetBinContent(k, va / vb);
      r->SetBinError(k, va / vb * std::sqrt(ea * ea + eb * eb));
   }
   return r;
}

// The systematic band of a result as a histogram whose error is the half-width of the band: around
// the central value (unity = false) or around 1 for a ratio pad (unity = true). The up and down
// histograms hold the ABSOLUTE deviations written by systematics.C.
inline TH1D *RelBand(const TH1D *centre, const TH1D *up, const TH1D *dn, const char *newName, bool unity)
{
   TH1D *band = (TH1D *)centre->Clone(newName);
   band->SetDirectory(nullptr);
   for (int k = 1; k <= band->GetNbinsX(); ++k) {
      const double v = centre->GetBinContent(k);
      if (v <= 0) {
         band->SetBinContent(k, 0);
         band->SetBinError(k, 0);
         continue;
      }
      const double relUp = up->GetBinContent(k) / v;
      const double relDn = dn->GetBinContent(k) / v;
      const double scale = unity ? 1.0 : v;
      band->SetBinContent(k, scale * (1 + (relUp - relDn) / 2));
      band->SetBinError(k, scale * (relUp + relDn) / 2);
   }
   return band;
}

// the published Table III as a violet band: around 1 (its relative systematic) or around its value
inline TH1D *RefBand(bool unity, const char *newName)
{
   TH1D *value = Take(Published(), "crossSection_statistic", (std::string(newName) + "_c").c_str());
   TH1D *syst = Take(Published(), "crossSection_systematic", (std::string(newName) + "_s").c_str());
   if (!value || !syst) return nullptr;
   TH1D *band = (TH1D *)value->Clone(newName);
   band->SetDirectory(nullptr);
   for (int k = 1; k <= band->GetNbinsX(); ++k) {
      const double v = value->GetBinContent(k);
      band->SetBinContent(k, unity ? 1.0 : v);
      band->SetBinError(k, v > 0 ? (unity ? syst->GetBinError(k) / v : syst->GetBinError(k)) : 0.0);
   }
   Band(band, kRefColor, kRefAlpha);
   return band;
}

// ---------------------------------------------------------------------------------------------
// A measured spectrum or ratio with its systematic band: the standard content of a result file
// ---------------------------------------------------------------------------------------------

struct Measured {
   TH1D *value = nullptr; // "canonical": the cross section, or a ratio of two of them
   TH1D *up = nullptr;    // "syst_up":   absolute upward systematic deviation
   TH1D *dn = nullptr;    // "syst_dn":   absolute downward systematic deviation
};

// true when all three histograms were found, i.e. the figure can be drawn
inline bool Complete(const Measured &m)
{
   return m.value && m.up && m.dn;
}

// Read canonical + syst_up + syst_dn from one result file and blank the bins outside its quoted
// range. nameTag makes the three ROOT names unique within one figure.
inline Measured TakeMeasured(const std::string &file, const std::string &nameTag)
{
   Measured m;
   m.value = Take(file, "canonical", (nameTag + "_val").c_str());
   m.up = Take(file, "syst_up", (nameTag + "_up").c_str());
   m.dn = Take(file, "syst_dn", (nameTag + "_dn").c_str());
   if (Complete(m)) Quoted(file, m.value, m.up, m.dn);
   return m;
}

// Draw a measurement the way the deck always draws one: the systematic band as one box per filled
// bin, then the statistical error bars and the points on top. Returns the value histogram, which
// carries the colour and the fill of the band and is the legend entry to use with option "pf".
inline TH1D *DrawMeasured(const Measured &m, const std::string &nameTag, int color, int marker,
                          double size = kMarkerSize, double alpha = kSystAlpha)
{
   TH1D *band = RelBand(m.value, m.up, m.dn, (nameTag + "_band").c_str(), false);
   DrawBand(band, color, alpha, false, "2");
   Marker(m.value, color, marker, size, alpha);
   m.value->Draw("P SAME");
   return m.value;
}

// ---------------------------------------------------------------------------------------------
// Canvases, frames, axes, legends
// ---------------------------------------------------------------------------------------------

// a canvas with a spectrum pad on top and a ratio pad below, sharing the pT axis
inline TCanvas *TwoPads(const char *name, TPad *&top, TPad *&bottom)
{
   auto *c = new TCanvas(name, name, 900, 900);
   c->Divide(1, 2);
   top = (TPad *)c->cd(1);
   top->SetPad(0, 0.42, 1, 1);
   top->SetBottomMargin(0.02);
   top->SetLeftMargin(0.13);
   top->SetRightMargin(0.04);
   top->SetTopMargin(0.06);
   bottom = (TPad *)c->cd(2);
   bottom->SetPad(0, 0, 1, 0.42);
   bottom->SetTopMargin(0.03);
   bottom->SetBottomMargin(0.28);
   bottom->SetLeftMargin(0.13);
   bottom->SetRightMargin(0.04);
   return c;
}

// a single wide pad, the shape of every ratio-only figure of the deck
inline TCanvas *OnePad(const char *name)
{
   auto *c = new TCanvas(name, name, 900, 560);
   c->SetLeftMargin(0.12);
   c->SetRightMargin(0.04);
   c->SetBottomMargin(0.14);
   c->SetTopMargin(0.06);
   return c;
}

// label and title sizes of one pad; hideX blanks the pT axis of the upper pad of a two-pad canvas
inline void Axes(TH1D *h, double label = 0.045, double title = 0.05, bool hideX = false)
{
   h->GetXaxis()->SetLabelSize(hideX ? 0 : label);
   h->GetXaxis()->SetTitleSize(hideX ? 0 : title);
   h->GetXaxis()->SetTitleOffset(1.1);
   h->GetYaxis()->SetLabelSize(label);
   h->GetYaxis()->SetTitleSize(title);
   h->GetYaxis()->SetTitleOffset(title >= 0.055 ? 0.9 : 1.2);
}

// An empty histogram carrying only the axes of a pad, titled "<title>;<x title>;<y title>", over
// the pT range xlo..xhi and the y range ylo..yhi. The curves go on top of it with "SAME".
inline TH1D *Frame(const char *name, const char *titles, double xlo, double xhi, double ylo, double yhi)
{
   auto *frame = new TH1D(name, titles, 1, xlo, xhi);
   frame->SetDirectory(nullptr);
   frame->SetStats(0);
   frame->SetMinimum(ylo);
   frame->SetMaximum(yhi);
   Axes(frame);
   return frame;
}

// a legend without border or fill, placed by its pad coordinates
inline TLegend *Leg(double x1 = 0.55, double y1 = 0.55, double x2 = 0.94, double y2 = 0.90, double textSize = 0.038)
{
   auto *l = new TLegend(x1, y1, x2, y2);
   l->SetBorderSize(0);
   l->SetFillStyle(0);
   l->SetTextSize(textSize);
   return l;
}

// ---------------------------------------------------------------------------------------------
// Normalised distributions: a shape comparison (vertex z, multiplicity, background density)
// ---------------------------------------------------------------------------------------------

struct Curve {
   TH1D *hist;        // the distribution; the panel scales it to unit integral
   std::string label; // legend entry
   int color;         // line colour
   int style;         // line style: 1 solid, 2 dashed
};

// One pad of normalised distributions: every curve is scaled to unit integral and drawn as a line,
// the first one carrying the axes, with room above the highest curve for the legend. The caller
// owns the pad (margins, log scale) and the legend, which is filled and drawn here.
inline void NormalisedPanel(const std::vector<Curve> &curves, const char *axisTitles, TLegend *leg,
                            double headroom = 1.3)
{
   double maximum = 0.0;
   for (const Curve &c : curves) {
      if (!c.hist || c.hist->Integral() <= 0) return;
      c.hist->Scale(1.0 / c.hist->Integral());
      Line(c.hist, c.color, kLineWidth, c.style);
      maximum = std::max(maximum, c.hist->GetMaximum());
   }
   TH1D *first = curves.front().hist;
   first->SetTitle(axisTitles);
   Axes(first, 0.05, 0.06);
   first->GetYaxis()->SetTitleOffset(1.1);
   first->SetMaximum(headroom * maximum);
   first->Draw("HIST");
   for (size_t i = 1; i < curves.size(); ++i) curves[i].hist->Draw("HIST SAME");
   for (const Curve &c : curves) leg->AddEntry(c.hist, c.label.c_str(), "l");
   leg->Draw();
}

// ---------------------------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------------------------

// write one canvas to FigDir() as <stem>.pdf and say where it went
inline void Save(TCanvas *c, const std::string &stem)
{
   const std::string path = FigDir() + stem + ".pdf";
   c->SaveAs(path.c_str());
   printf("[plot] %s\n", path.c_str());
}

} // namespace PlotStyle

#endif // PLOT_STYLE_H
