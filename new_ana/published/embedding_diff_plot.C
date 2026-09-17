// embedding_diff_plot.C — what separates the 2023 embedding sample from the 2021 one. Everything is
// shown as the ratio 2023 / 2021 at the same PARTICLE jet pT, so that the two samples are compared
// at equal truth and any difference belongs to the detector simulation or to the generated event:
// the jet energy scale and resolution, the matching efficiency, the fragmentation of the detector
// and of the generated jets, the background density, the vertex.
// Inputs: results/R0.5/embedding_diff_e21.root and _e23.root, written by embedding_diff.C; the
// histograms of a sample carry its tag as a name suffix.
// Outputs: results/figures/published_embedding_diff_{detector,generator,event}.pdf, plus a printed
// table of every ratio bin by bin and the vertex-gate fraction.
//   root -l -b -q 'embedding_diff_plot.C+'
#include "../plot_style.h"

#include <TGraphErrors.h>
#include <TLatex.h>
#include <TProfile.h>

#include <string>
#include <vector>

using namespace PlotStyle;

// the particle-jet pT window of every panel
const double kPartPtLo = 5.0;
const double kPartPtHi = 52.0;

// the two embedding samples, opened once by the entry point
TFile *gEmb2021 = nullptr; // Cold QCD request, pi0/eta/Sigma0 decayed in the stored generator record
TFile *gEmb2023 = nullptr; // HP request, the same particles left undecayed

// where a panel's legend sits inside its frame, in pad coordinates
struct LegendBox {
   double x1; // left
   double y1; // bottom
   double x2; // right
   double y2; // top
};
const LegendBox kLegendTop = {0.16, 0.72, 0.95, 0.89};    // the default: a strip across the top
const LegendBox kLegendMiddle = {0.50, 0.40, 0.95, 0.60}; // when the curves themselves fill the top

// one curve of a ratio panel
struct RatioCurve {
   TGraphErrors *graph; // 2023 / 2021 versus particle jet pT
   const char *label;   // legend entry
};

// the suffix that names the histograms of one sample inside its file
static const char *Tag(TFile *f)
{
   return f == gEmb2021 ? "e21" : "e23";
}

// one profile of one sample
static TProfile *Profile(TFile *f, const char *name)
{
   return (TProfile *)f->Get(Form("%s_%s", name, Tag(f)));
}

// one histogram of one sample
static TH1D *Hist(TFile *f, const char *name)
{
   return (TH1D *)f->Get(Form("%s_%s", name, Tag(f)));
}

// the ratio 2023 / 2021 of the means of one profile, with the statistical errors of both samples
static TGraphErrors *ProfileRatio(const char *name)
{
   auto *g = new TGraphErrors();
   TProfile *newer = Profile(gEmb2023, name);
   TProfile *older = Profile(gEmb2021, name);
   if (!newer || !older) return g;
   for (int k = 1; k <= newer->GetNbinsX(); ++k) {
      const double a = newer->GetBinContent(k);
      const double b = older->GetBinContent(k);
      if (b <= 0) continue;
      const double ea = newer->GetBinError(k);
      const double eb = older->GetBinError(k);
      g->SetPoint(g->GetN(), newer->GetBinCenter(k), a / b);
      g->SetPointError(g->GetN() - 1, 0, a / b * std::sqrt(ea * ea / (a * a) + eb * eb / (b * b)));
   }
   return g;
}

// The ratio 2023 / 2021 of the SPREADS of one profile: sqrt(<x^2> - <x>^2) from the profile of the
// quantity and the profile of its square. That spread is the jet energy resolution.
static TGraphErrors *SpreadRatio(const char *meanName, const char *meanSquareName)
{
   auto *g = new TGraphErrors();
   TProfile *newer = Profile(gEmb2023, meanName);
   TProfile *older = Profile(gEmb2021, meanName);
   TProfile *newerSquare = Profile(gEmb2023, meanSquareName);
   TProfile *olderSquare = Profile(gEmb2021, meanSquareName);
   if (!newer || !older || !newerSquare || !olderSquare) return g;
   for (int k = 1; k <= newer->GetNbinsX(); ++k) {
      const double meanNew = newer->GetBinContent(k);
      const double meanOld = older->GetBinContent(k);
      const double spreadNew = std::sqrt(std::max(0.0, newerSquare->GetBinContent(k) - meanNew * meanNew));
      const double spreadOld = std::sqrt(std::max(0.0, olderSquare->GetBinContent(k) - meanOld * meanOld));
      if (spreadOld <= 0) continue;
      g->SetPoint(g->GetN(), newer->GetBinCenter(k), spreadNew / spreadOld);
   }
   return g;
}

// the ratio 2023 / 2021 of one histogram, bin by bin
static TGraphErrors *HistRatio(const char *name)
{
   auto *g = new TGraphErrors();
   TH1D *newer = Hist(gEmb2023, name);
   TH1D *older = Hist(gEmb2021, name);
   if (!newer || !older) return g;
   for (int k = 1; k <= newer->GetNbinsX(); ++k) {
      if (older->GetBinContent(k) <= 0) continue;
      g->SetPoint(g->GetN(), newer->GetBinCenter(k), newer->GetBinContent(k) / older->GetBinContent(k));
   }
   return g;
}

// matched jets over generated jets of one sample, per particle pT bin, for one jet selection
static TH1D *MatchingEfficiency(TFile *f, const char *selection)
{
   TH1D *matched = Hist(f, Form("nM_%s", selection));
   TH1D *generated = Hist(f, "nP");
   if (!matched || !generated) return nullptr;
   auto *efficiency = (TH1D *)matched->Clone(Form("eff_%s_%s", selection, Tag(f)));
   efficiency->SetDirectory(nullptr);
   efficiency->Divide(generated);
   return efficiency;
}

// the ratio 2023 / 2021 of two matching efficiencies
static TGraphErrors *EfficiencyRatio(const char *selection)
{
   auto *g = new TGraphErrors();
   TH1D *newer = MatchingEfficiency(gEmb2023, selection);
   TH1D *older = MatchingEfficiency(gEmb2021, selection);
   if (!newer || !older) return g;
   for (int k = 1; k <= newer->GetNbinsX(); ++k) {
      if (older->GetBinContent(k) <= 0) continue;
      g->SetPoint(g->GetN(), newer->GetBinCenter(k), newer->GetBinContent(k) / older->GetBinContent(k));
   }
   return g;
}

// the deck margins every panel of this macro uses
static void PreparePad(TVirtualPad *pad)
{
   pad->cd();
   gPad->SetLeftMargin(0.14);
   gPad->SetRightMargin(0.03);
   gPad->SetBottomMargin(0.15);
   gPad->SetTopMargin(0.10);
}

// One panel of 2023 / 2021 ratios over the particle jet pT, with a dashed line at unity, the
// legend inside the frame above the curves and the panel title over the frame.
static void RatioPanel(TCanvas *c, int pad, const char *title, const std::vector<RatioCurve> &curves, double ylo,
                       double yhi, const LegendBox &box = kLegendTop)
{
   PreparePad(c->cd(pad));
   TH1D *frame = Frame(Form("fr%d_%s", pad, title), ";particle jet p_{T} [GeV];2023 / 2021", kPartPtLo, kPartPtHi, ylo, yhi);
   Axes(frame, 0.05, 0.06);
   frame->GetYaxis()->SetTitleOffset(1.1);
   frame->Draw("AXIS");
   Baseline(kPartPtLo, kPartPtHi);

   TLegend *leg = Leg(box.x1, box.y1, box.x2, box.y2, 0.05);
   for (size_t i = 0; i < curves.size(); ++i) {
      TGraphErrors *g = curves[i].graph;
      g->SetMarkerStyle(SeriesMarker((int)i));
      g->SetMarkerColor(SeriesColor((int)i));
      g->SetLineColor(SeriesColor((int)i));
      g->SetMarkerSize(kMarkerSize);
      g->SetLineWidth(kLineWidth);
      g->Draw("P SAME");
      leg->AddEntry(g, curves[i].label, "lp");
   }
   leg->Draw();

   TLatex caption;
   caption.SetNDC();
   caption.SetTextSize(0.06);
   caption.DrawLatex(0.14, 0.93, title);
}

// Figure published_embedding_diff_detector: what the detector simulation does differently in the
// two samples, at the same particle pT — the jet energy scale <pT_det / pT_part>, its spread (the
// resolution), the fraction of generated jets that are reconstructed and matched, and the
// fragmentation of the detector jet (neutral energy fraction, constituents, leading fraction).
// Reads the r_*, r2_*, nM_*, nP, nefD_*, nconD_* and leadD_* histograms of both samples.
static void DetectorFigure()
{
   gStyle->SetPaperSize(24, 16); // a canvas of 900 x 600 fills this page
   auto *c = new TCanvas("c1", "", 900, 600);
   c->Divide(2, 2, 0.004, 0.004);
   RatioPanel(c, 1, "jet energy scale  #LTp_{T}^{det}/p_{T}^{part}#GT",
              {{ProfileRatio("r_jp"), "jets of the JP levels"}, {ProfileRatio("r_all"), "all jets, no trigger"}}, 0.97,
              1.03);
   RatioPanel(c, 2, "jet energy resolution  #sigma(p_{T}^{det}/p_{T}^{part})",
              {{SpreadRatio("r_jp", "r2_jp"), "jets of the JP levels"},
               {SpreadRatio("r_all", "r2_all"), "all jets, no trigger"}},
              0.90, 1.10);
   RatioPanel(c, 3, "matching efficiency  (matched, 5 #leq p_{T}^{det} < 60) / generated",
              {{EfficiencyRatio("jp"), "jets of the JP levels"}, {EfficiencyRatio("all"), "all jets, no trigger"}},
              0.90, 1.10);
   RatioPanel(c, 4, "detector jet: neutral fraction, constituents, leading fraction",
              {{ProfileRatio("nefD_all"), "neutral energy fraction"},
               {ProfileRatio("nconD_all"), "number of constituents"},
               {ProfileRatio("leadD_all"), "leading constituent / jet p_{T}"}},
              0.90, 1.10);
   Save(c, "published_embedding_diff_detector");
}

// Figure published_embedding_diff_generator: what the GENERATED events already differ in, before
// any detector — the fragmentation of the particle jets at the same pT, the background density at
// both levels, the weighted particle-jet spectrum, and the vertex z of the thrown and the
// reconstructed events. Reads nefP, nconP, leadP, bgP, bgD_all, nP, vzThrown and vzReco.
static void GeneratorFigure()
{
   gStyle->SetPaperSize(24, 16);
   auto *c = new TCanvas("c2", "", 900, 600);
   c->Divide(2, 2, 0.004, 0.004);
   RatioPanel(c, 1, "generated jets: fragmentation at the same p_{T}",
              {{ProfileRatio("nefP"), "neutral energy fraction"},
               {ProfileRatio("nconP"), "number of constituents"},
               {ProfileRatio("leadP"), "leading particle / jet p_{T}"}},
              0.90, 1.10, kLegendMiddle);
   RatioPanel(c, 2, "background density #rho (jet-area method)",
              {{ProfileRatio("bgP"), "particle level"}, {ProfileRatio("bgD_all"), "detector level"}}, 0.5, 1.5);
   RatioPanel(c, 3, "generated jet spectrum (weighted), 2023 / 2021",
              {{HistRatio("nP"), "particle jets, |#eta| < 0.5"}}, 0.90, 1.10);

   PreparePad(c->cd(4));
   TLegend *leg = Leg(0.16, 0.68, 0.6, 0.89, 0.05);
   NormalisedPanel({{Hist(gEmb2021, "vzReco"), "2021 reconstructed (weighted)", SeriesColor(0), 1},
                    {Hist(gEmb2023, "vzReco"), "2023 reconstructed (weighted)", SeriesColor(1), 1},
                    {Hist(gEmb2021, "vzThrown"), "2021 thrown", SeriesColor(0), 2},
                    {Hist(gEmb2023, "vzThrown"), "2023 thrown", SeriesColor(1), 2}},
                   ";vertex z [cm];fraction of events", leg);
   TLatex caption;
   caption.SetNDC();
   caption.SetTextSize(0.06);
   caption.DrawLatex(0.14, 0.93, "vertex z");

   Save(c, "published_embedding_diff_generator");
}

// Figure published_embedding_diff_event: the reconstructed background density per event, the
// event-level measure of the underlying event of a sample, normalised to unit area, mean in the legend.
// Reads the bgEv histogram of both samples.
static void EventFigure()
{
   auto *c = OnePad("c3");
   PreparePad(c);
   gPad->SetLogy();
   TH1D *older = Hist(gEmb2021, "bgEv");
   TH1D *newer = Hist(gEmb2023, "bgEv");
   if (!older || !newer) return;
   TLegend *leg = Leg(0.55, 0.72, 0.95, 0.89, 0.05);
   NormalisedPanel({{older, Form("2021, mean %.3f GeV", older->GetMean()), SeriesColor(0), 1},
                    {newer, Form("2023, mean %.3f GeV", newer->GetMean()), SeriesColor(1), 1}},
                   ";detector #rho per event [GeV];fraction of events", leg);
   Save(c, "published_embedding_diff_event");
}

// the ratios of the detector and generator figures, bin by bin, and the event-level means
static void PrintTable()
{
   printf("bin      JES23/21   JER23/21   eff23/21   nefD   nconD   leadD  |  nefP   nconP   leadP  | rhoP  rhoD\n");
   TGraphErrors *columns[11] = {ProfileRatio("r_jp"),        SpreadRatio("r_jp", "r2_jp"),
                                EfficiencyRatio("jp"),       ProfileRatio("nefD_all"),
                                ProfileRatio("nconD_all"),   ProfileRatio("leadD_all"),
                                ProfileRatio("nefP"),        ProfileRatio("nconP"),
                                ProfileRatio("leadP"),       ProfileRatio("bgP"),
                                ProfileRatio("bgD_all")};
   for (int k = 0; k < columns[0]->GetN(); ++k) {
      printf("%5.1f  ", columns[0]->GetX()[k]);
      for (int i = 0; i < 11; ++i) printf(" %6.3f ", k < columns[i]->GetN() ? columns[i]->GetY()[k] : 0.0);
      printf("\n");
   }

   // meta bin 1 = generated events, bin 2 = events with a reconstructed vertex inside the gate
   auto *meta2021 = (TH1D *)gEmb2021->Get("meta");
   auto *meta2023 = (TH1D *)gEmb2023->Get("meta");
   if (meta2021 && meta2023) {
      printf("vertex-gate fraction: 2021 %.4f  2023 %.4f\n", meta2021->GetBinContent(2) / meta2021->GetBinContent(1),
             meta2023->GetBinContent(2) / meta2023->GetBinContent(1));
   }
}

void embedding_diff_plot()
{
   Init();
   const std::string dir = RadiusDir("0.5");
   gEmb2021 = TFile::Open((dir + "embedding_diff_e21.root").c_str());
   gEmb2023 = TFile::Open((dir + "embedding_diff_e23.root").c_str());
   if (!gEmb2021 || gEmb2021->IsZombie() || !gEmb2023 || gEmb2023->IsZombie()) {
      printf("[plot] no embedding_diff_e21/e23.root in %s: run published/embedding_diff.C first\n", dir.c_str());
      return;
   }
   // the printed table first: the figures normalise the event distributions in place
   PrintTable();
   DetectorFigure();
   GeneratorFigure();
   EventFigure();
}
