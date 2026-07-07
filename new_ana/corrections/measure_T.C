// measure_T.C ----------------------------------------------------------------
// Measure the trigger-turn-on DATA/EMBEDDING correction T(pt) for one STAR
// jet-patch trigger (JP1 or JP2), pp200 Run-12 inclusive jets, anti-kT R=0.5.
//
// WHAT IS T?  For a jet of given pT ask: "does the SIMULATED trigger match this
// jet?"  That probability (the turn-on curve) is measured twice:
//   P_data(pt) = P(sim trigger match | jet pt) in real DATA, evaluated on the
//                JP0-fired base sample.  JP0's threshold is far below JP1/JP2,
//                so at low pt the base is unbiased w.r.t. the JP1/JP2 decision.
//   P_emb(pt)  = the same probability in EMBEDDING (matched MC tree, native
//                per-jet reco_trigger_match_<trig> gate — the gate the
//                response build itself uses).
//   T(pt)      = P_data(pt) / P_emb(pt).
// The unfolding response already models the turn-on itself: reco jets failing
// the trigger gate become Miss entries.  What the Misses CANNOT model is any
// INFIDELITY of the simulated calorimeter: the simulated BEMC fires slightly
// "hotter" than the real one (root cause: the pi0/eta generator decay setting
// differs between embedding requests), so the embedding turn-on rises faster
// than the data one.  T(pt) is exactly that residual data/sim ratio, in situ.
//
// That(pt) = the turn-on SHAPE only: normalize T at the [16.1,19) bin and clamp
// to 1.0 at and above it. NOTE (2026-07-07): the plateau of T is NOT ~1 — it
// sits ~0.90 (JP2) / ~0.95 (JP1), i.e. T = R x T_true with T_true =
// P(sw|data)/P(sw|emb) a REAL data/embedding emulator-efficiency difference.
// That difference is deliberately LEFT UNCORRECTED (Dmitry-faithful: he unfolds
// software-data by software-embedding and carries it as an embedding systematic);
// the clamp discards it here, and the plateau RULER bridge is supplied instead by
// R (corrections/measure_R.C) in config.h::TrigEffMeas. Dividing by the un-clamped
// T would over-lift JP2 ~+10%.
//
// Usage (ROOT>=6):  root -l -b -q 'measure_T.C("<merged_data_R0.5.root>",
//                     "<merged_matching_R0.5.root>", 29, "JP1")'
// adcRegisterPlus1 = DSM threshold register + 1 = the native sim patch-ADC
// gate on the STORED data jp_patch_adc (JP1: 28+1=29, JP2: 36+1=37).
// -----------------------------------------------------------------------------
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <cstdio>
#include <vector>

void measure_T(const char *jp0DataTree, const char *matchedTree, int adcRegisterPlus1, const char *trig)
{
   const std::vector<double> QB = {6.9, 8.2, 9.7, 11.5, 13.6, 16.1, 19, 22.5, 26.6, 31.4, 37.2, 44, 52};
   const int nb = (int)QB.size() - 1; // the analysis' quoted truth bins (GeV)

   // ---- DATA: ResultTree = one row per EVENT; per-jet branches are RVecs, so
   // "pt_corrected[sel]" keeps only the jets passing the per-jet mask `sel`
   // (|detector eta|<0.5, neutral fraction<0.95 = the analysis jet selection).
   // `hit` additionally asks the stored jet-patch ADC to clear the native sim
   // threshold: "would the simulated trigger have matched this jet?".
   ROOT::RDataFrame d("ResultTree", jp0DataTree);
   auto dd = d.Filter("fired_JP0") // unbiased base: the event fired JP0
               .Define("sel", "abs(det_eta)<0.5 && neutral_fraction<0.95")
               .Define("hit", Form("sel && jp_patch_adc>=%d", adcRegisterPlus1))
               .Define("ptAll", "pt_corrected[sel]")
               .Define("ptHit", "pt_corrected[hit]");
   auto hA = dd.Histo1D({"hA", "", nb, QB.data()}, "ptAll"); // denominator
   auto hH = dd.Histo1D({"hH", "", nb, QB.data()}, "ptHit"); // numerator

   // ---- EMBEDDING: MatchedTree = one row per truth-reco pair. Same reco-side
   // jet selection (reco_pt>-500 means "a reco partner exists"); the numerator
   // is the native per-jet trigger gate used by the response build itself.
   ROOT::RDataFrame m("MatchedTree", matchedTree);
   auto mb = m.Filter("reco_pt>-500 && std::abs(reco_det_eta)<0.5 && reco_neutral_fraction<0.95");
   auto hE = mb.Histo1D({"hE", "", nb, QB.data()}, "reco_pt");
   auto hF = mb.Filter(Form("reco_trigger_match_%s", trig)).Histo1D({"hF", "", nb, QB.data()}, "reco_pt");

   // ---- turn the four counts into efficiencies (binomial errors), then T.
   TH1D *pD = (TH1D *)hH->Clone("pD"); pD->Divide(hH.GetPtr(), hA.GetPtr(), 1, 1, "B");
   TH1D *pE = (TH1D *)hF->Clone("pE"); pE->Divide(hF.GetPtr(), hE.GetPtr(), 1, 1, "B");
   TH1D *T  = (TH1D *)pD->Clone("T");  T->Divide(pE);
   TH1D *Th = (TH1D *)T->Clone("Th");  // That = shape-only version of T:
   const double Tnorm = T->GetBinContent(T->FindFixBin(17.0)); // T at [16.1,19)
   for (int i = 1; i <= nb; ++i) {
      if (QB[i - 1] >= 16.1 - 1e-6) { Th->SetBinContent(i, 1.0); Th->SetBinError(i, 0); } // clamp plateau
      else { Th->SetBinContent(i, T->GetBinContent(i) / Tnorm); Th->SetBinError(i, T->GetBinError(i) / Tnorm); }
   }
   printf("\n%s: T(pt) = P(sim match | pt, data JP0 base) / P(sim match | pt, embedding)\n", trig);
   for (int i = 1; i <= nb; ++i)
      printf("  [%5.1f,%5.1f)  P_data=%.4f  P_emb=%.4f  T=%.4f  That=%.4f\n", QB[i - 1], QB[i],
             pD->GetBinContent(i), pE->GetBinContent(i), T->GetBinContent(i), Th->GetBinContent(i));

   // ---- one clear figure: turn-ons on top, their ratio T (+ applied That) below.
   gStyle->SetOptStat(0);
   TCanvas c("c", "", 700, 800); c.Divide(1, 2, 0, 0);
   c.cd(1); gPad->SetGrid();
   pD->SetTitle(Form("%s turn-on;;P(sim trigger match)", trig));
   pD->GetYaxis()->SetRangeUser(0.0, 1.09); pD->SetMarkerStyle(20); pD->SetLineColor(kBlue + 1);
   pD->SetMarkerColor(kBlue + 1); pD->Draw("E1"); pE->SetMarkerStyle(24); pE->SetMarkerColor(kRed + 1);
   pE->SetLineColor(kRed + 1); pE->Draw("E1 SAME");
   TLegend lg(0.55, 0.15, 0.88, 0.35); lg.SetBorderSize(0);
   lg.AddEntry(pD, "DATA (JP0-fired base)", "lp"); lg.AddEntry(pE, "EMBEDDING (matched)", "lp"); lg.Draw();
   c.cd(2); gPad->SetGrid();
   T->SetTitle(";jet p_{T} [GeV];T = P_{data}/P_{emb}");
   T->GetYaxis()->SetRangeUser(0.35, 1.14); T->SetMarkerStyle(21); T->SetMarkerColor(kBlack); T->Draw("E1");
   Th->SetLineColor(kGreen + 2); Th->SetLineWidth(2); Th->Draw("HIST SAME");
   TLine one(QB.front(), 1, QB.back(), 1); one.SetLineStyle(2); one.Draw();
   TLegend lg2(0.55, 0.15, 0.88, 0.32); lg2.SetBorderSize(0);
   lg2.AddEntry(T, "T (measured)", "lp"); lg2.AddEntry(Th, "#hat{T} (normalized+clamped)", "l"); lg2.Draw();
   c.SaveAs(Form("measure_T_%s.pdf", trig));
}
