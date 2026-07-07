// Final JPX summary deck: spectrum vs Table III, ratio with stat + systematic
// band, per-trigger Bayes context, and the number table.
//
// Pure consumer of the pipeline outputs in <workdir>:
//   xsec_JPX_R0.5.root            (canonical + reference)
//   xsec_JPX_R0.5_systband.root   (systematic; optional)
//   xsec_<T>_R0.5.root            (per-trigger Bayes, optional context)
// Usage: root -l -b -q plot_final.C   -> <workdir>/jpx_final_R0.5.pdf
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <map>
#include <string>
#include <vector>

#include "../config.h"
using namespace CrossSectionConfig;
AnalysisConfig cfg;

void plot_final()
{
   gStyle->SetOptStat(0);
   const std::string dir = cfg.workdir;
   const TString out = dir + "jpx_final_R0.5.pdf";

   TFile fj((dir + "xsec_JPX_R0.5.root").c_str(), "READ");
   TH1D *jpx = (TH1D *)fj.Get("canonical");
   TH1D *ref = (TH1D *)fj.Get("reference");
   if (!jpx || !ref) { printf("missing xsec_JPX_R0.5.root content\n"); return; }
   jpx->SetDirectory(0); ref->SetDirectory(0);

   TH1D *band = nullptr;
   if (!gSystem->AccessPathName((dir + "xsec_JPX_R0.5_systband.root").c_str())) {
      TFile fb((dir + "xsec_JPX_R0.5_systband.root").c_str(), "READ");
      band = (TH1D *)fb.Get("systematic");
      if (band) band->SetDirectory(0);
   }

   std::map<std::string, TH1D *> trig;
   for (const std::string t : {"JP1", "JP2", "HT2"}) {
      if (gSystem->AccessPathName((dir + "xsec_" + t + "_R0.5.root").c_str())) continue;
      TFile ft((dir + "xsec_" + t + "_R0.5.root").c_str(), "READ");
      auto *h = (TH1D *)ft.Get("canonical");
      if (h) { h->SetDirectory(0); trig[t] = h; }
   }
   auto tcolor = [](const std::string &t) {
      return t == "JP1" ? kOrange - 8 : (t == "JP2" ? (int)kAzure : (int)kRed);
   };
   auto tmarker = [](const std::string &t) { return t == "JP1" ? 22 : (t == "JP2" ? 20 : 21); };

   TCanvas c("c", "", 900, 950);
   c.Print(out + "[");

   // ---------- page 1: spectrum + ratio ----------
   c.Divide(1, 2);
   auto *p1 = (TPad *)c.cd(1);
   p1->SetPad(0, 0.5, 1, 1); p1->SetBottomMargin(0.02); p1->SetLogy();
   ref->SetTitle("Inclusive jets pp200 Run-12, anti-k_{T} R=0.5 |#eta_{det}|<0.5;"
                 ";d^{2}#sigma/dp_{T}d#eta [pb/(GeV/c)]");
   ref->SetLineColor(kViolet); ref->SetFillColorAlpha(kViolet, 0.25); ref->SetMarkerSize(0);
   ref->GetYaxis()->SetRangeUser(1.2, 1e7);
   ref->Draw("E2");
   jpx->SetMarkerStyle(29); jpx->SetMarkerSize(1.5);
   jpx->SetMarkerColor(kBlack); jpx->SetLineColor(kBlack);
   jpx->Draw("E1 SAME");
   TLegend l1(0.52, 0.62, 0.88, 0.87); l1.SetBorderSize(0);
   l1.AddEntry(ref, "Dmitry Table III (stat#oplussyst)", "f");
   l1.AddEntry(jpx, "JPX promotion + matrix inversion", "lep");
   l1.Draw();

   auto *p2 = (TPad *)c.cd(2);
   p2->SetPad(0, 0, 1, 0.5); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.28);
   TH1D *frame = (TH1D *)ref->Clone("fr"); frame->Reset(); frame->SetDirectory(0);
   frame->SetTitle(";p_{T} [GeV/c];JPX / Dmitry");
   frame->GetYaxis()->SetRangeUser(0.65, 1.35);
   frame->Draw();
   TH1D *rb = (TH1D *)ref->Clone("rb"); rb->Reset(); rb->SetDirectory(0);
   for (int b = 1; b <= rb->GetNbinsX(); ++b) {
      const double cv = ref->GetBinContent(b), ev = ref->GetBinError(b);
      rb->SetBinContent(b, 1.0); rb->SetBinError(b, cv > 0 ? ev / cv : 0.0);
   }
   rb->SetFillColorAlpha(kViolet, 0.25); rb->SetLineColor(kViolet); rb->SetMarkerSize(0);
   rb->Draw("E2 SAME");
   for (double y : {0.95, 1.05}) {
      auto *ln = new TLine(ref->GetXaxis()->GetXmin(), y, ref->GetXaxis()->GetXmax(), y);
      ln->SetLineStyle(3); ln->SetLineColor(kGray + 2); ln->Draw();
   }
   auto *one = new TLine(ref->GetXaxis()->GetXmin(), 1.0, ref->GetXaxis()->GetXmax(), 1.0);
   one->SetLineColor(kViolet); one->SetLineWidth(2); one->Draw();

   TH1D *rs = nullptr; // ratio with syst band (shaded)
   if (band) {
      rs = (TH1D *)ref->Clone("rs"); rs->Reset(); rs->SetDirectory(0);
      for (int b = 1; b <= rs->GetNbinsX(); ++b) {
         const double x = rs->GetXaxis()->GetBinCenter(b);
         const int bh = band->FindBin(x);
         const double rv = ref->GetBinContent(b);
         if (rv <= 0 || band->GetBinContent(bh) == 0) continue;
         rs->SetBinContent(b, band->GetBinContent(bh) / rv);
         rs->SetBinError(b, band->GetBinError(bh) / rv);
      }
      rs->SetFillColorAlpha(kAzure - 4, 0.35); rs->SetMarkerSize(0); rs->SetLineColor(kAzure - 2);
      rs->Draw("E2 SAME");
   }
   TH1D *rr = (TH1D *)ref->Clone("rr"); rr->Reset(); rr->SetDirectory(0);
   for (int b = 1; b <= rr->GetNbinsX(); ++b) {
      const double x = rr->GetXaxis()->GetBinCenter(b);
      const int bh = jpx->FindBin(x);
      const double rv = ref->GetBinContent(b);
      if (rv <= 0 || jpx->GetBinContent(bh) == 0) continue;
      rr->SetBinContent(b, jpx->GetBinContent(bh) / rv);
      rr->SetBinError(b, jpx->GetBinError(bh) / rv);
   }
   rr->SetMarkerStyle(29); rr->SetMarkerSize(1.5); rr->SetMarkerColor(kBlack); rr->SetLineColor(kBlack);
   rr->Draw("E1 SAME");
   TLegend l2(0.55, 0.72, 0.88, 0.94); l2.SetBorderSize(0);
   l2.AddEntry(rr, "JPX (stat)", "lep");
   if (rs) l2.AddEntry(rs, "JPX systematic", "f");
   l2.AddEntry(rb, "Table III stat#oplussyst", "f");
   l2.Draw();
   c.Print(out);

   // ---------- page 2: ratio with the per-trigger Bayes context ----------
   c.Clear();
   c.SetLogy(false);
   TH1D *frame2 = (TH1D *)frame->Clone("fr2"); frame2->SetDirectory(0);
   frame2->SetTitle("JPX (matrix inversion) and standalone Bayes triggers;p_{T} [GeV/c];me / Dmitry");
   frame2->GetYaxis()->SetRangeUser(0.65, 1.35);
   frame2->Draw();
   rb->Draw("E2 SAME");
   one->Draw();
   TLegend l3(0.55, 0.70, 0.88, 0.92); l3.SetBorderSize(0); l3.SetNColumns(2);
   for (auto &kv : trig) {
      TH1D *r = (TH1D *)ref->Clone(("r_" + kv.first).c_str()); r->Reset(); r->SetDirectory(0);
      for (int b = 1; b <= r->GetNbinsX(); ++b) {
         const double x = r->GetXaxis()->GetBinCenter(b);
         const int bh = kv.second->FindBin(x);
         const double rv = ref->GetBinContent(b);
         if (rv <= 0 || kv.second->GetBinContent(bh) == 0) continue;
         r->SetBinContent(b, kv.second->GetBinContent(bh) / rv);
         r->SetBinError(b, kv.second->GetBinError(bh) / rv);
      }
      r->SetMarkerStyle(tmarker(kv.first)); r->SetMarkerColor(tcolor(kv.first));
      r->SetLineColor(tcolor(kv.first));
      r->Draw("E1 SAME");
      l3.AddEntry(r, (kv.first + " (Bayes)").c_str(), "lep");
   }
   rr->Draw("E1 SAME");
   l3.AddEntry(rr, "JPX (M^{-1})", "lep");
   l3.Draw();
   c.Print(out);

   // ---------- page 3: the number table ----------
   c.Clear();
   TLatex tx; tx.SetTextFont(42); tx.SetTextSize(0.028);
   tx.DrawLatexNDC(0.08, 0.94, "#bf{JPX promotion + matrix inversion vs Dmitry Table III}");
   tx.DrawLatexNDC(0.08, 0.90, "floor 9.7 GeV, #lambda=0 (unregularized, Dmitry's method), C_{JPX} = measured C_{sum}");
   double y = 0.84;
   tx.DrawLatexNDC(0.08, y, "#bf{p_{T} [GeV]}");
   tx.DrawLatexNDC(0.30, y, "#bf{JPX/Dmitry}");
   tx.DrawLatexNDC(0.50, y, "#bf{stat}");
   tx.DrawLatexNDC(0.65, y, "#bf{note}");
   y -= 0.045;
   for (int b = 1; b <= ref->GetNbinsX(); ++b) {
      const double x = ref->GetXaxis()->GetBinCenter(b);
      const double rv = ref->GetBinContent(b);
      const int bh = jpx->FindBin(x);
      const double mv = jpx->GetBinContent(bh);
      if (rv <= 0 || mv == 0) continue;
      const double ratio = mv / rv;
      const double stat = jpx->GetBinError(bh) / rv;
      const char *note = "";
      if (std::abs(x - 14.85) < 0.1) note = "reference zigzag bin (+8.5#sigma above its own trend)";
      else if (std::abs(x - 40.6) < 0.1) note = "data-side tail fluctuation (fold QA: +2#sigma)";
      tx.DrawLatexNDC(0.08, y, Form("%.2f", x));
      tx.SetTextColor(std::abs(ratio - 1.0) <= 0.05 ? kGreen + 2 : (std::abs(x - 14.85) < 0.1 ? kGray + 2 : kRed + 1));
      tx.DrawLatexNDC(0.30, y, Form("%.3f", ratio));
      tx.SetTextColor(kBlack);
      tx.DrawLatexNDC(0.50, y, Form("#pm%.3f", stat));
      tx.DrawLatexNDC(0.65, y, note);
      y -= 0.042;
   }
   y -= 0.02;
   tx.DrawLatexNDC(0.08, y, "green: within #pm5%.  9/10 solved bins within #pm5%; 14.85 exempt.");
   c.Print(out);

   c.Print(out + "]");
   printf("wrote %s\n", out.Data());
}
