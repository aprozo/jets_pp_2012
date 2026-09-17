// plot.C — the figures of the analysis deck: every trigger alone against Table III, the ratios to
// JP1, the radius ratios with and without the generators, the underlying-event on/off ratio, the
// no-UE reference set, the generators against JP1 at every radius, and the hadronisation correction.
// Inputs: results/R<R>/ (analysis/run.sh), results/Michal/ (analysis/run.sh michal) and
// results/models/ (models/run_models.sh compare). A figure whose inputs are absent is skipped.
//   root -l -b -q 'plot.C+'      -> results/figures/analysis_*.pdf
#include "../plot_style.h"

#include <string>
#include <vector>

using namespace PlotStyle;

// one trigger as it appears in the two trigger figures
struct TriggerSeries {
   const char *mode;       // the mode in the result file names: xsec_<mode>_R<R>_syst_bayes4.root
   const char *trigger;    // key of the trigger palette (colour and marker)
   const char *label;      // legend entry of the cross-section figure
   const char *ratioLabel; // legend entry of the ratio-to-JP1 figure; nullptr for JP1 itself
};
const std::vector<TriggerSeries> kTriggers = {
   {"mb2023", "MB", "min-bias", "min-bias / JP1"},
   {"jp1_e23", "JP1", "JP1", nullptr},
   {"jp2_e23", "JP2", "JP2", "JP2 / JP1"},
   {"ht2_e23", "HT2", "HT2, efficiency corrected", "HT2 / JP1, efficiency corrected"}};

// the radius the other radii are divided by in every radius ratio
const char *kReferenceRadius = "0.5";

// the pT below which a generator bin is carried by a handful of jets of the softest pT-hat samples
const double kModelFloor = 9.7;

// the reconstruction floor of the jet-patch levels: below it no unfolded point is shown
const double kRecoFloor = 8.2;

// Figure analysis_triggers_R<R>: each trigger unfolded on its own (Bayes, four iterations, 2023
// embedding). At R = 0.5 as a ratio to the published Table III inside the published systematic band;
// at the other radii as the cross sections with their own bands. Only the quoted range is drawn.
// Reads results/R<R>/xsec_<mode>_R<R>_syst_bayes4.root ("canonical" + its .txt quoted range)
// and inputs/jet_cross_section_publishedR0.5.root.
static void Triggers(const char *jetR)
{
   const std::string dir = RadiusDir(jetR);
   const bool published = std::string(jetR) == kReferenceRadius; // Table III exists at R = 0.5 only
   TCanvas *c = OnePad((std::string("triggers_R") + jetR).c_str());
   TLegend *leg = Leg(0.52, 0.13, 0.94, 0.40, 0.034);

   // R = 0.5: the ratio to the published Table III inside its band
   TH1D *reference = nullptr;
   if (published) {
      reference = Take(Published(), "crossSection_statistic", "ref3");
      TH1D *ratioBand = RefBand(true, "refrel3");
      if (!reference || !ratioBand) return;
      ratioBand->SetTitle("Every trigger alone, Bayes 4 it., 2023 embedding; p_{T} [GeV/c]; this / published");
      Axes(ratioBand);
      ratioBand->GetXaxis()->SetRangeUser(kPtLo, kPtHi);
      ratioBand->GetYaxis()->SetRangeUser(0.7, 1.3);
      ratioBand->Draw("E2");
      Guides(kPtLo, kPtHi);
      leg->AddEntry(ratioBand, "published systematic band", "f");
   } else {
      // other radii: the cross sections themselves with their bands
      TH1D *frame = Frame((std::string("frame_triggers_") + jetR).c_str(),
                          Form("Every trigger alone, R = %s, Bayes 4 it., 2023 embedding; p_{T} [GeV/c];"
                               " d^{2}#sigma / dp_{T} d#eta [pb/(GeV/c)]", jetR),
                          kPtLo, kPtHi, 0.5, 5e7);
      c->SetLogy();
      frame->Draw();
      leg->SetY1(0.62);
      leg->SetY2(0.90);
   }

   for (const TriggerSeries &t : kTriggers) {
      const std::string file = dir + "xsec_" + t.mode + "_R" + jetR + "_syst_bayes4.root";
      if (published) {
         TH1D *xsec = Take(file, "canonical", t.mode);
         if (!xsec) continue;
         Quoted(file, xsec);
         TH1D *ratio = Ratio(xsec, reference, (std::string("rt_") + t.mode).c_str(), true);
         Points(ratio, t.trigger);
         ratio->Draw("P SAME");
         leg->AddEntry(ratio, t.label, "p");
      } else {
         Measured m = TakeMeasured(file, std::string("tr_") + t.mode);
         if (!Complete(m)) continue;
         TH1D *drawn = DrawMeasured(m, std::string("tr_") + t.mode, TrigColor(t.trigger), TrigMarker(t.trigger));
         leg->AddEntry(drawn, Form("%s (stat., syst.)", t.label), "pf");
      }
   }
   leg->Draw();
   Save(c, std::string("analysis_triggers_R") + jetR);
}

// Figure analysis_ratios_to_jp1_R<R>: min-bias, JP2 and HT2 divided by JP1, each with its own
// systematic band (the correlated parts cancel in the ratio), Bayes 4 it. on the 2023 embedding.
// This is the internal consistency check of the trigger stitching.
// Reads results/R<R>/ratio_<mode>_over_jp1_e23_R<R>_bayes4.root
// ("canonical", "syst_up", "syst_dn" + its .txt quoted range).
static void RatiosToJp1(const char *jetR)
{
   const std::string dir = RadiusDir(jetR);
   TCanvas *c = OnePad((std::string("ratios_jp1_R") + jetR).c_str());
   TH1D *frame = Frame((std::string("frame_rjp1_") + jetR).c_str(),
                       "Ratios to JP1, Bayes 4 it., 2023 embedding; p_{T} [GeV/c]; trigger / JP1",
                       kPtLo, kPtHi, 0.6, 1.5);
   frame->Draw("AXIS");
   Guides(kPtLo, kPtHi);

   TLegend *leg = Leg(0.14, 0.70, 0.70, 0.90);
   for (const TriggerSeries &t : kTriggers) {
      if (!t.ratioLabel) continue;
      const std::string file = dir + "ratio_" + t.mode + "_over_jp1_e23_R" + jetR + "_bayes4.root";
      const Measured ratio = TakeMeasured(file, std::string("rr_") + t.mode);
      if (!Complete(ratio)) continue;
      TH1D *drawn = DrawMeasured(ratio, std::string("rr_") + t.mode, TrigColor(t.trigger), TrigMarker(t.trigger));
      leg->AddEntry(drawn, t.ratioLabel, "pf");
   }
   leg->Draw();
   Save(c, std::string("analysis_ratios_to_jp1_R") + jetR);
}

// Figure analysis_radii: sigma(R) / sigma(R = 0.5) for R = 0.2, 0.3 and 0.4 (JP1, Bayes 4 it., 2023
// embedding) with their systematic bands, compared with the same ratio in PYTHIA 6 Perugia 2012.
// The narrower the cone, the more of the jet is lost, so the ratio falls with decreasing R.
// Reads results/R<R>/rratio_jp1_e23_R<R>_over_R0.5_bayes4.root ("canonical", "syst_up", "syst_dn")
// and the PYTHIA 6 particle-level spectra results/models/models_R<R>_ue.root.
static void Radii()
{
   TCanvas *c = OnePad("radii");
   TH1D *frame = Frame("frame_radii",
                       "Radius ratios, JP1, Bayes 4 it., 2023 embedding, |#eta| < 0.5;"
                       " p_{T} [GeV/c]; #sigma(R) / #sigma(R = 0.5)",
                       kPtLo, kPtHi, 0, 1.05);
   frame->Draw("AXIS");

   TLegend *leg = Leg(0.50, 0.14, 0.94, 0.42, 0.034);
   TH1D *modelReference = Take(ModelFile(kReferenceRadius, "ue"), "pythia6", "p6ref");
   for (const Radius &r : kRadii) {
      if (std::string(r.name) == kReferenceRadius) continue;
      const std::string file =
         RadiusDir(r.name) + "rratio_jp1_e23_R" + r.name + "_over_R" + kReferenceRadius + "_bayes4.root";
      const Measured ratio = TakeMeasured(file, std::string("rr") + r.name);
      if (!Complete(ratio)) continue;
      TH1D *drawn = DrawMeasured(ratio, std::string("rr") + r.name, r.color, r.marker, kBigMarkerSize);
      leg->AddEntry(drawn, Form("R = %s / R = %s (stat., syst.)", r.name, kReferenceRadius), "pf");

      TH1D *model = Take(ModelFile(r.name, "ue"), "pythia6", (std::string("p6") + r.name).c_str());
      if (!model || !modelReference) continue;
      TH1D *modelRatio = Ratio(model, modelReference, (std::string("p6r") + r.name).c_str());
      Below(modelRatio, kModelFloor);
      TBox *key = DrawBand(modelRatio, r.color, kFaintAlpha, true, "3", 2);
      if (std::string(r.name) == kRadii.front().name) leg->AddEntry(key, "PYTHIA 6 Perugia 2012 (band: MC stat.)", "lf");
   }
   leg->Draw();
   Save(c, "analysis_radii");
}

// Figure analysis_radii_models_R<R>: one radius ratio R / 0.5 at a time (JP1, Bayes 4 it., 2023
// embedding) against all four generators, which differ here mainly in their parton shower and
// their underlying event. Both levels are UE subtracted, so the ratio tests the jet shape.
// Reads results/R<R>/rratio_jp1_e23_R<R>_over_R0.5_bayes4.root ("canonical", "syst_up", "syst_dn")
// and results/models/models_R<R>_ue.root / models_R0.5_ue.root for every generator.
static void RadiiModels(const char *jetR)
{
   const std::string file =
      RadiusDir(jetR) + "rratio_jp1_e23_R" + jetR + "_over_R" + kReferenceRadius + "_bayes4.root";
   const Measured ratio = TakeMeasured(file, std::string("drr") + jetR);
   if (!Complete(ratio)) return;

   TCanvas *c = OnePad((std::string("rrm") + jetR).c_str());
   const double ymax = std::string(jetR) == "0.2" ? 0.6 : std::string(jetR) == "0.3" ? 0.85 : 1.05;
   TH1D *frame = Frame((std::string("frame_rrm") + jetR).c_str(),
                       Form("Radius ratio R = %s / R = %s, JP1 and generators, |#eta| < 0.5, UE subtracted;"
                            " p_{T} [GeV/c]; #sigma(R = %s) / #sigma(R = %s)",
                            jetR, kReferenceRadius, jetR, kReferenceRadius),
                       kPtLo, kPtHi, 0, ymax);
   frame->Draw("AXIS");

   TLegend *leg = Leg(0.50, 0.14, 0.94, 0.44, 0.034);
   TH1D *drawn = DrawMeasured(ratio, std::string("drr") + jetR, TrigColor("JP1"), TrigMarker("JP1"),
                              kMarkerSize, kSoftAlpha);
   leg->AddEntry(drawn, "JP1, Bayes 4 it., 2023 embedding (stat., syst.)", "pf");

   for (const Model &m : kModels) {
      TH1D *here = Take(ModelFile(jetR, "ue"), m.name, (std::string("ma") + m.name + jetR).c_str());
      TH1D *reference = Take(ModelFile(kReferenceRadius, "ue"), m.name, (std::string("mb") + m.name + jetR).c_str());
      if (!here || !reference) continue;
      TH1D *modelRatio = Ratio(here, reference, (std::string("mr") + m.name + jetR).c_str());
      Below(modelRatio, kModelFloor);
      leg->AddEntry(DrawBand(modelRatio, m.color, kSoftAlpha, true, "3"), m.label, "lf");
   }
   leg->Draw();
   Save(c, std::string("analysis_radii_models_R") + jetR);
}

// Figure analysis_ue_onoff: how much the off-axis-cone underlying-event subtraction removes, as
// the ratio of the unsubtracted to the subtracted cross section for the four radii (JP1, Bayes
// 4 it.). The effect grows with the cone area and falls with pT.
// Reads results/R<R>/xsec_inversion_R<R>_jp1_e23_bayes4.root and
// results/R<R>/xsec_inversion_R<R>_jp1_e23_noue_bayes4.root ("canonical").
static void UE()
{
   TCanvas *c = OnePad("ue");
   TH1D *frame = Frame("frame_ue",
                       "Underlying event: JP1 without / with the off-axis subtraction, Bayes 4 it.;"
                       " p_{T} [GeV/c]; #sigma(no UE subtraction) / #sigma(UE subtracted)",
                       kPtLo, kPtHi, 0.95, 1.5);
   frame->Draw("AXIS");
   Baseline(kPtLo, kPtHi);

   TLegend *leg = Leg(0.55, 0.60, 0.94, 0.90);
   for (const Radius &r : kRadii) {
      const std::string dir = RadiusDir(r.name);
      TH1D *subtracted = Take(dir + "xsec_inversion_R" + r.name + "_jp1_e23_bayes4.root", "canonical",
                              (std::string("w") + r.name).c_str());
      TH1D *raw = Take(dir + "xsec_inversion_R" + r.name + "_jp1_e23_noue_bayes4.root", "canonical",
                       (std::string("o") + r.name).c_str());
      if (!subtracted || !raw) continue;
      TH1D *ratio = Ratio(raw, subtracted, (std::string("ue") + r.name).c_str());
      Below(ratio, kRecoFloor);
      Marker(ratio, r.color, r.marker, kBigMarkerSize);
      ratio->Draw("P SAME");
      leg->AddEntry(ratio, Form("R = %s", r.name), "p");
   }
   leg->Draw();
   Save(c, "analysis_ue_onoff");
}

// Figure analysis_noue_R<R>: the reference set without any underlying-event subtraction, in the
// 5-60 GeV binning — JP1 with its systematic band from 10 GeV, the combined jet-patch levels and
// the min-bias level up to 20 GeV — as a spectrum and as a ratio to JP1.
// Reads results/Michal/xsec_R<R>_noUE_bins5-60.root ("xsec_jp1", "xsec_jp1_syst_up",
// "xsec_jp1_syst_dn", "xsec_combined", "xsec_mb").
static void NoUeSet(const char *jetR)
{
   const double xlo = 4.0;  // the 5-60 GeV binning needs a little room on both sides
   const double xhi = 61.0;
   const double jp1Floor = 10.0; // JP1 is quoted from 10 GeV in this set
   const double mbCeiling = 20.0;

   const std::string file = Results() + "Michal/xsec_R" + jetR + "_noUE_bins5-60.root";
   TH1D *jp1 = Take(file, "xsec_jp1", "mjp1");
   TH1D *up = Take(file, "xsec_jp1_syst_up", "mu");
   TH1D *dn = Take(file, "xsec_jp1_syst_dn", "md");
   TH1D *combined = Take(file, "xsec_combined", "mcomb");
   TH1D *minbias = Take(file, "xsec_mb", "mmb");
   if (!jp1 || !up || !dn || !combined || !minbias) return;

   Below(jp1, jp1Floor);
   for (int k = 1; k <= minbias->GetNbinsX(); ++k) {
      if (minbias->GetBinLowEdge(k + 1) <= mbCeiling + 1e-6) continue;
      minbias->SetBinContent(k, 0);
      minbias->SetBinError(k, 0);
   }

   TPad *top = nullptr;
   TPad *bottom = nullptr;
   TCanvas *c = TwoPads((std::string("michal") + jetR).c_str(), top, bottom);

   top->cd();
   top->SetLogy();
   TH1D *band = RelBand(jp1, up, dn, "mband", false);
   band->SetTitle(Form("No UE subtraction, R = %s, |#eta| < 0.5, Bayes 4 it., 2023 embedding;"
                       " p_{T} [GeV/c]; d^{2}#sigma / dp_{T} d#eta [pb/(GeV/c)]",
                       jetR));
   Band(band, TrigColor("JP1"), kSoftAlpha);
   Axes(band, 0.045, 0.05, true);
   band->GetXaxis()->SetRangeUser(xlo, xhi);
   band->SetMinimum(0.1);
   band->Draw("E2");
   Points(jp1, "JP1");
   jp1->Draw("P SAME");
   Points(combined, "JPX");
   combined->Draw("P SAME");
   Points(minbias, "MB");
   minbias->Draw("P SAME");
   TLegend *leg = Leg(0.50, 0.62, 0.94, 0.90);
   leg->AddEntry(jp1, "JP1, from 10 GeV/c (stat., syst.)", "pf");
   leg->AddEntry(combined, "JP0+JP1+JP2 combined (stat.)", "p");
   leg->AddEntry(minbias, "min-bias, to 20 GeV/c (stat.)", "p");
   leg->Draw();

   bottom->cd();
   TH1D *ratioBand = RelBand(jp1, up, dn, "mrb", true);
   ratioBand->SetTitle("; p_{T} [GeV/c]; ratio to JP1");
   Band(ratioBand, TrigColor("JP1"), kSoftAlpha);
   Axes(ratioBand, 0.055, 0.06);
   ratioBand->GetXaxis()->SetRangeUser(xlo, xhi);
   ratioBand->GetYaxis()->SetRangeUser(0.7, 1.4);
   ratioBand->GetYaxis()->SetNdivisions(505);
   ratioBand->Draw("E2");
   Guides(xlo, xhi);
   TH1D *rCombined = Ratio(combined, jp1, "rc");
   TH1D *rMinbias = Ratio(minbias, jp1, "rm");
   Points(rCombined, "JPX");
   Points(rMinbias, "MB");
   rCombined->Draw("P SAME");
   rMinbias->Draw("P SAME");

   Save(c, std::string("analysis_noue_R") + jetR);
}

// Figure analysis_models_<def>_R<R>: the four generators against the JP1 cross section, as a
// spectrum and as a ratio to JP1; at R = 0.5 with the UE subtraction the published Table III is
// shown as well. "ue" is the nominal definition in the published bins, "noue" the 5-60 GeV set.
// Reads results/models/models_R<R>_<def>.root ("jp1", "jp1_syst_up", "jp1_syst_dn", one histogram
// and one "<generator>_over_jp1" per generator, and at R = 0.5 "published" / "published_syst").
static void Models(const char *definition, const char *jetR)
{
   const bool ue = std::string(definition) == "ue";
   const double xlo = ue ? kPtLo : 4.0;
   const double xhi = ue ? kPtHi : 61.0;
   const double ptFloor = ue ? kRecoFloor : 10.0;

   const std::string file = ModelFile(jetR, definition);
   TH1D *jp1 = Take(file, "jp1", "gjp1");
   TH1D *up = Take(file, "jp1_syst_up", "gu");
   TH1D *dn = Take(file, "jp1_syst_dn", "gd");
   if (!jp1 || !up || !dn) return;
   TH1D *published = (ue && std::string(jetR) == "0.5") ? Take(file, "published", "gpub") : nullptr;
   TH1D *publishedSyst = published ? Take(file, "published_syst", "gpubS") : nullptr;
   Below(jp1, ptFloor);

   const std::string stem = std::string("analysis_models_") + definition + "_R" + jetR;
   TPad *top = nullptr;
   TPad *bottom = nullptr;
   TCanvas *c = TwoPads(stem.c_str(), top, bottom);

   top->cd();
   top->SetLogy();
   TH1D *band = RelBand(jp1, up, dn, "gband", false);
   band->SetTitle(Form("Generators against the data, R = %s, |#eta| < 0.5, %s;"
                       " p_{T} [GeV/c]; d^{2}#sigma / dp_{T} d#eta [pb/(GeV/c)]",
                       jetR, ue ? "UE subtracted" : "no UE subtraction"));
   Band(band, TrigColor("JP1"), kSoftAlpha);
   Axes(band, 0.045, 0.05, true);
   band->GetXaxis()->SetRangeUser(xlo, xhi);
   band->SetMinimum(ue ? 0.5 : 0.1);
   band->Draw("E2");

   TLegend *leg = Leg(0.42, 0.55, 0.94, 0.90, 0.034);
   leg->AddEntry(band, ue ? "JP1, Bayes 4 it., 2023 embedding (stat., syst.)" : "JP1, from 10 GeV/c (stat., syst.)", "pf");
   if (published && publishedSyst) {
      TH1D *publishedBand = (TH1D *)published->Clone("pb");
      for (int k = 1; k <= publishedBand->GetNbinsX(); ++k) publishedBand->SetBinError(k, publishedSyst->GetBinError(k));
      Band(publishedBand, kRefColor, kRefAlpha);
      publishedBand->Draw("E2 SAME");
      leg->AddEntry(publishedBand, "published Table III (stat., syst.)", "f");
   }

   std::vector<TH1D *> modelRatios;
   for (const Model &m : kModels) {
      TH1D *spectrum = Take(file, m.name, m.name);
      TH1D *ratio = Take(file, (std::string(m.name) + "_over_jp1").c_str(), (std::string("r") + m.name).c_str());
      if (!spectrum || !ratio) continue;
      leg->AddEntry(DrawBand(spectrum, m.color, kBandAlpha, true, "3"), m.label, "lf");
      ratio->SetLineColor(m.color);
      modelRatios.push_back(ratio);
   }
   Points(jp1, "JP1");
   jp1->Draw("P SAME");
   leg->Draw();

   bottom->cd();
   TH1D *ratioBand = RelBand(jp1, up, dn, "grb", true);
   ratioBand->SetTitle("; p_{T} [GeV/c]; model / JP1");
   Band(ratioBand, TrigColor("JP1"), kSoftAlpha);
   Axes(ratioBand, 0.055, 0.06);
   ratioBand->GetXaxis()->SetRangeUser(xlo, xhi);
   ratioBand->GetYaxis()->SetRangeUser(0.4, 1.7);
   ratioBand->GetYaxis()->SetNdivisions(505);
   ratioBand->Draw("E2");
   Guides(xlo, xhi);
   for (TH1D *ratio : modelRatios) DrawBand(ratio, ratio->GetLineColor(), kBandAlpha, true, "3");
   if (published) {
      TH1D *publishedRatio = Take(file, "published_over_jp1", "pr");
      if (publishedRatio) {
         Marker(publishedRatio, kRefColor, 20, 1.2);
         publishedRatio->SetLineWidth(1); // thin bars: the published points are a reference, not a result
         publishedRatio->Draw("P SAME");
      }
   }

   Save(c, stem);
}

// Figure analysis_chad_R0.5: the hadronisation correction C_had = sigma_particle / sigma_parton of
// PYTHIA 6 Perugia 2012, the factor a fixed-order (parton-level) prediction has to be multiplied by
// before it can be compared with the measured particle-level cross section.
// Reads results/models/chad_R0.5.root ("chad_pub" and its uncertainty "chad_pub_syst", the
// quadrature sum of the radiation, fragmentation and Innsbruck tune variations).
static void Chad()
{
   const std::string file = Results() + "models/chad_R0.5.root";
   TH1D *chad = Take(file, "chad_pub", "chad");
   TH1D *uncertainty = Take(file, "chad_pub_syst", "chadS");
   if (!chad || !uncertainty) return;

   TCanvas *c = OnePad("chad");
   TH1D *band = (TH1D *)chad->Clone("chadband");
   for (int k = 1; k <= band->GetNbinsX(); ++k) {
      const double u = uncertainty->GetBinContent(k);
      band->SetBinError(k, u > 0 ? u : uncertainty->GetBinError(k));
   }
   band->SetTitle("Hadronisation correction, PYTHIA 6 Perugia 2012, R = 0.5, |#eta| < 0.5;"
                  " p_{T} [GeV/c]; C_{had} = #sigma_{particle} / #sigma_{parton}");
   Band(band, kModels.front().color, kSoftAlpha);
   Axes(band);
   band->GetXaxis()->SetRangeUser(kPtLo, kPtHi);
   band->GetYaxis()->SetRangeUser(0.3, 1.1);
   band->Draw("E2");
   Baseline(kPtLo, kPtHi);

   Marker(chad, kModels.front().color, 20, kMarkerSize);
   chad->Draw("P SAME");

   TLegend *leg = Leg(0.45, 0.16, 0.94, 0.34);
   leg->AddEntry(chad, "C_{had} (stat.)", "p");
   leg->AddEntry(band, "tune variations in quadrature", "f");
   leg->Draw();

   Save(c, "analysis_chad_R0.5");
}

void plot()
{
   Init();
   for (const char *jetR : {"0.5", "0.2", "0.3", "0.4"}) {
      Triggers(jetR);
      RatiosToJp1(jetR);
      NoUeSet(jetR);
      Models("ue", jetR);
      Models("noue", jetR);
   }
   Radii();
   UE();
   for (const char *jetR : {"0.2", "0.3", "0.4"}) RadiiModels(jetR);
   Chad();
}
