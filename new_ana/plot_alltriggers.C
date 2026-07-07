// Overlay every trigger's unfolded cross section against Dmitry's Table III.
// 2-panel figure: log-y spectra vs the paper band (top) + me/paper ratio (bottom).
//
// Pure consumer: reads the per-trigger summary files xsec_<T>_R0.5.root that
// cross_section.cpp wrote, each holding TH1D "canonical" (my spectrum) + TH1D
// "reference" (Dmitry's paper band). Independent per-trigger files can be
// produced with their own knobs and overlaid here without re-unfolding.
//
// Usage: root -l -b -q 'plot_alltriggers.C("<dir>","comparison_alltriggers_R0.5.pdf")'
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <map>
#include <string>
#include <vector>

void plot_alltriggers(const char *dir = ".", const char *outname = "comparison_alltriggers_R0.5.pdf")
{
   const std::vector<std::string> trigs = {"JPX", "JP2", "HT2", "JP1", "JP0"};
   auto color_for = [](const std::string &t) -> int {
      if (t == "JPX") return kBlack;
      if (t == "JP2") return kAzure;    if (t == "HT2") return kRed;
      if (t == "JP1") return kOrange-8; if (t == "JP0") return kGreen-2;
      return kMagenta + 1;
   };
   auto marker_for = [](const std::string &t) -> int {
      if (t == "JPX") return 29;
      if (t == "JP2") return 20; if (t == "HT2") return 21;
      if (t == "JP1") return 22; if (t == "JP0") return 33;
      return 34;
   };
   gStyle->SetOptStat(0);
   TH1D *ref = nullptr;
   std::map<std::string, TH1D *> spectra;
   for (const auto &t : trigs) {
      TFile *f = TFile::Open(Form("%s/xsec_%s_R0.5.root", dir, t.c_str()));
      if (!f || f->IsZombie()) { printf("skip %s (no file)\n", t.c_str()); continue; }
      auto *h = (TH1D *)f->Get("canonical");
      if (!h) { printf("skip %s (no canonical)\n", t.c_str()); f->Close(); continue; }
      h->SetDirectory(0); spectra[t] = h;
      if (!ref) { ref = (TH1D *)f->Get("reference"); ref->SetDirectory(0); }
   }
   if (!ref) { printf("no inputs found\n"); return; }

   TCanvas *c = new TCanvas("c_alltrig", "", 900, 1000);
   c->Divide(1, 2);
   auto *p1 = (TPad *)c->cd(1); p1->SetPad(0, 0.5, 1, 1); p1->SetBottomMargin(0); p1->SetLogy();
   ref->SetTitle("Inclusive jet cross section R=0.5 |#eta|<0.5; p_{T} [GeV/c]; "
                 "d^{2}#sigma/dp_{T}d#eta [pb/(GeV/c)]");
   ref->SetLineColor(kViolet); ref->SetFillColorAlpha(kViolet, 0.20); ref->SetMarkerSize(0);
   ref->Draw("E2");
   TLegend *leg = new TLegend(0.60, 0.50, 0.88, 0.88); leg->SetBorderSize(0); leg->SetNColumns(2);
   leg->AddEntry(ref, "Dmitry (paper)", "f");
   for (const auto &t : trigs) {
      if (!spectra.count(t)) continue;
      TH1D *h = spectra[t];
      h->SetMarkerStyle(marker_for(t)); h->SetMarkerColor(color_for(t)); h->SetLineColor(color_for(t));
      h->Draw("E1 SAME"); leg->AddEntry(h, t.c_str(), "lep");
   }
   leg->Draw();

   auto *p2 = (TPad *)c->cd(2); p2->SetPad(0, 0, 1, 0.5); p2->SetTopMargin(0); p2->SetBottomMargin(0.30);
   TH1D *frame = (TH1D *)ref->Clone("ratio_frame"); frame->Reset(); frame->SetDirectory(0);
   frame->SetTitle("; p_{T} [GeV/c]; me / paper"); frame->GetYaxis()->SetRangeUser(0.8, 1.4); frame->Draw();
   TH1D *band = (TH1D *)ref->Clone("ref_band"); band->Reset(); band->SetDirectory(0);
   for (int b = 1; b <= band->GetNbinsX(); ++b) {
      const double cv = ref->GetBinContent(b), ev = ref->GetBinError(b);
      band->SetBinContent(b, 1.0); band->SetBinError(b, cv > 0 ? ev / cv : 0.0);
   }
   band->SetFillColorAlpha(kViolet, 0.20); band->SetFillStyle(1001); band->SetLineColor(kViolet);
   band->SetMarkerSize(0); band->Draw("E2 SAME");
   TLine *line = new TLine(ref->GetXaxis()->GetXmin(), 1.0, ref->GetXaxis()->GetXmax(), 1.0);
   line->SetLineColor(kViolet); line->SetLineWidth(2); line->Draw("SAME");
   for (const auto &t : trigs) {
      if (!spectra.count(t)) continue;
      TH1D *h = spectra[t];
      TH1D *r = new TH1D(Form("ratio_%s", t.c_str()), "", ref->GetNbinsX(),
                         ref->GetXaxis()->GetXbins()->GetArray());
      r->SetDirectory(0);
      for (int b = 1; b <= r->GetNbinsX(); ++b) {
         const double x = r->GetXaxis()->GetBinCenter(b); const int bh = h->FindBin(x);
         const double rv = ref->GetBinContent(b);
         r->SetBinContent(b, rv > 0 ? h->GetBinContent(bh) / rv : 0.0);
         r->SetBinError(b, rv > 0 ? h->GetBinError(bh) / rv : 0.0);
      }
      r->SetMarkerStyle(marker_for(t)); r->SetMarkerColor(color_for(t)); r->SetMarkerSize(1.4);
      r->SetLineColor(color_for(t));
      r->Draw("E1 SAME");
   }
   c->SaveAs(Form("%s/%s", dir, outname));
   printf("wrote %s/%s (%zu triggers)\n", dir, outname, spectra.size());
}
