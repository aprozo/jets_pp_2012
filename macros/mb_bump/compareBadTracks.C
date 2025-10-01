#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>

void draw(const TString &trigger = "MB", const TString &trigger_title = "MB") {

  TCanvas *c1 = new TCanvas("c1", "QA comparison", 1200, 800);
  c1->cd();
  gPad->SetGridx();
  gPad->SetLogy();
  c1->SaveAs(trigger + "comparison.pdf[");
  TLegend *legend = new TLegend(0.6, 0.65, 0.99, 0.92);

  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.06);
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.07);

  TFile *inputFile = TFile::Open("~/jets_" + trigger + ".root");
  //   read all histograms

  const int nSize = 5;
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening file!" << std::endl;
    return;
  }

  TString dir[] = {"QA_histograms/", "QA_histograms_problematic/"};
  TString suffix[nSize] = {"sDCAxy", "DCA", "Chi2", "Chi2PV", "matched"};

  TTree *tree = (TTree *)inputFile->Get("ResultTree");
  Double_t pt;
  if (!tree) {
    std::cerr << "Error: Tree ResultTree not found in file!" << std::endl;
    return;
  }
  tree->SetBranchAddress("pt", &pt);

  vector<float> pt_reco_bins = {5.0,  6.0,  7.0,  8.0,  9.0,  10.0, 11.0,
                                12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
                                19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0,
                                26.0, 27.0, 28.0, 29.0, 30.0, 31.0};

  TH1D *h_pt = new TH1D("h_pt", ";p_{t} (GeV/c);dN/dp", pt_reco_bins.size() - 1,
                        &pt_reco_bins[0]);
  h_pt->SetLineColor(kBlue);
  h_pt->SetLineWidth(2);
  h_pt->SetMarkerColor(kBlue);
  h_pt->SetMarkerStyle(21);
  h_pt->GetXaxis()->SetTitleSize(0.07);
  h_pt->GetYaxis()->SetTitleSize(0.05);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    h_pt->Fill(pt);
  }

  // normalize the histogram by pt bin width
  for (int i = 1; i <= h_pt->GetNbinsX(); i++) {
    double binWidth = h_pt->GetBinWidth(i);
    if (binWidth > 0) {
      h_pt->SetBinContent(i, h_pt->GetBinContent(i) / binWidth);
      h_pt->SetBinError(i, h_pt->GetBinError(i) / binWidth);
    } else {
      std::cerr << "Warning: Bin width is zero for bin " << i
                << ". Skipping normalization for this bin." << std::endl;
    }
  }

  // $$ \frac{d\sigma}{dp_T} = A \cdot \left(1 + \frac{p_T}{p_0} \right)^{-n} $$
  TF1 *fitFunc = new TF1("fitFunc", "[0] * (1 + x/[1])^(-[2])", 5, 20);
  TF1 *fitFuncExt = new TF1("fitFunc", "[0] * (1 + x/[1])^(-[2])", 5, 31);
  fitFunc->SetParameters(6e7, 1.7, 8);
  fitFunc->SetParNames("A", "p0", "n");
  h_pt->Fit(fitFunc, "QRN");
  fitFuncExt->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1),
                            fitFunc->GetParameter(2));
  fitFuncExt->SetLineColor(kRed);
  c1->cd();
  h_pt->Draw("E1");
  fitFuncExt->Draw("same");

  latex->DrawLatexNDC(0.3, 0.2, trigger_title);
  c1->SaveAs(trigger + "comparison.pdf");

  // TFile *hist = TFile::Open("~/jets_MB_leadingtrack_highpt.root");
  TH1D *QA_highjetpt_leadingtower_pt =
      (TH1D *)inputFile->Get(dir[0] + "highjetpt_leadingtower_pt");
  TH1D *QA_highjetpt_leadingtrack_pt =
      (TH1D *)inputFile->Get(dir[0] + "highjetpt_leadingtrack_pt");

  if (!QA_highjetpt_leadingtower_pt || !QA_highjetpt_leadingtrack_pt) {
    std::cerr << "Error: Histograms highjetpt_leadingtower_pt or "
                 "highjetpt_leadingtrack_pt not found!"
              << std::endl;
    return;
  }
  // c1->Divide(2, 1);
  // c1->cd(1);
  // gPad->SetLogy();
  // QA_highjetpt_leadingtower_pt->Draw();

  // c1->cd(2);
  c1->cd();
  gPad->SetLogy();
  QA_highjetpt_leadingtrack_pt->GetXaxis()->SetRangeUser(0, 35);
  QA_highjetpt_leadingtrack_pt->GetYaxis()->SetTitle("Counts");
  QA_highjetpt_leadingtrack_pt->GetXaxis()->SetTitleSize(0.07);
  QA_highjetpt_leadingtrack_pt->GetYaxis()->SetTitleSize(0.05);
  QA_highjetpt_leadingtrack_pt->SetTitle("");

  QA_highjetpt_leadingtrack_pt->Draw();

  TLine *line =
      new TLine(22, 0, 22, 1.05 * QA_highjetpt_leadingtrack_pt->GetMaximum());
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw("same");
  latex->SetTextColor(kRed);
  latex->DrawLatexNDC(0.6, 0.75, "problematic");
  latex->SetTextColor(kBlack);
  c1->cd();
  latex->DrawLatexNDC(0.3, 0.92, trigger_title);

  c1->SaveAs(trigger + "comparison.pdf");
  c1->Clear();

  TH1D *QA_problematic[nSize];
  TH1D *QA_good[nSize];
  for (int i = 0; i < nSize; i++) {
    QA_good[i] = (TH1D *)inputFile->Get(dir[0] + suffix[i]);
    if (!QA_good[i]) {
      std::cerr << "Error: Histogram " << dir[0] + suffix[i] << " not found!"
                << std::endl;
      continue;
    }
    QA_good[i]->SetLineColor(kBlue);
    QA_good[i]->SetLineWidth(2);
    QA_good[i]->SetMarkerColor(kBlue);
    QA_good[i]->SetMarkerStyle(21);

    QA_problematic[i] = (TH1D *)inputFile->Get(dir[1] + suffix[i]);
    if (!QA_problematic[i]) {
      std::cerr << "Error: Histogram " << dir[1] + suffix[i] << " not found!"
                << std::endl;
      continue;
    }
    QA_problematic[i]->SetLineColor(kRed);
    QA_problematic[i]->SetLineWidth(2);
    QA_problematic[i]->SetMarkerColor(kRed);
    QA_problematic[i]->SetMarkerStyle(20);
    // increas title labels
    QA_problematic[i]->GetYaxis()->SetTitleSize(0.05);
    QA_good[i]->GetYaxis()->SetTitleSize(0.05);
    QA_problematic[i]->GetXaxis()->SetTitleSize(0.07);
    QA_good[i]->GetXaxis()->SetTitleSize(0.07);
  }

  TH2D *QA_pt_sDCAxy_pos = (TH2D *)inputFile->Get(dir[0] + "pt_sDCAxy_pos");
  TH2D *QA_pt_sDCAxy_neg = (TH2D *)inputFile->Get(dir[0] + "pt_sDCAxy_neg");

  if (!QA_pt_sDCAxy_pos || !QA_pt_sDCAxy_neg) {
    std::cerr << "Error: Histograms pt_sDCAxy_pos or pt_sDCAxy_neg not found!"
              << std::endl;
    return;
  }

  vector<float> pt_bins = {10, 15, 20, 30};

  for (size_t i = 0; i < pt_bins.size() - 1; i++) {
    float pt_min = pt_bins[i];
    float pt_max = pt_bins[i + 1];
    int bin_min = QA_pt_sDCAxy_pos->GetYaxis()->FindBin(pt_min);
    int bin_max = QA_pt_sDCAxy_pos->GetYaxis()->FindBin(pt_max);

    TH1D *project_pos =
        QA_pt_sDCAxy_pos->ProjectionX("pt_sDCAxy_pos_proj", bin_min, bin_max);
    TH1D *project_neg =
        QA_pt_sDCAxy_neg->ProjectionX("pt_sDCAxy_neg_proj", bin_min, bin_max);

    // increase title size
    project_pos->GetYaxis()->SetTitleSize(0.05);
    project_neg->GetYaxis()->SetTitleSize(0.05);
    project_pos->GetXaxis()->SetTitleSize(0.07);
    project_neg->GetXaxis()->SetTitleSize(0.07);

    project_pos->SetLineColor(kRed);
    project_neg->SetLineColor(kBlue);
    project_pos->SetLineWidth(2);
    project_neg->SetLineWidth(2);

    project_pos->SetTitle(";sDCAxy (cm);Normalized Counts");
    project_pos->GetYaxis()->SetTitleOffset(0.9);
    gPad->SetLogy();
    project_pos->GetXaxis()->SetRangeUser(-1.5, 1.5);
    project_pos->DrawNormalized();
    project_neg->DrawNormalized("same");
    if (i == 0) {
      legend->AddEntry(project_pos, " positive", "l");
      legend->AddEntry(project_neg, " negative", "l");
    }
    latex->DrawLatexNDC(0.3, 0.2, trigger_title);
    latex->DrawLatexNDC(
        0.15, 0.85, Form(" %.1f< p_{t, track}< %.1f GeV/c,", pt_min, pt_max));
    legend->Draw();

    c1->SaveAs(trigger + "comparison.pdf");
  }
  legend->Clear();

  c1->Clear();

  for (int i = 0; i < nSize; i++) {
    c1->Clear();
    QA_problematic[i]->SetTitle("");
    QA_problematic[i]->GetYaxis()->SetTitle("Normalized Counts");
    QA_problematic[i]->GetYaxis()->SetTitleOffset(0.9);

    QA_problematic[i]->GetYaxis()->SetTitleSize(0.07);
    QA_problematic[i]->GetXaxis()->SetLabelSize(0.07);
    QA_problematic[i]->DrawNormalized();
    QA_good[i]->DrawNormalized("same");
    // latex->DrawLatexNDC(0.1, 0.75, "Jets with single track constituent");
    latex->DrawLatexNDC(0.2, 0.92, trigger_title);
    legend->AddEntry(QA_problematic[i], "Problematic(p_{t}>22)", "l");
    legend->AddEntry(QA_good[i], "Good(p_{t}<22)", "l");
    legend->Draw();

    c1->SaveAs(trigger + "comparison.pdf");
    legend->Clear();
  }

  c1->SaveAs(trigger + "comparison.pdf]");
}

void compareBadTracks() {
  vector<TString> triggers = {"MB_dca1_noVPD", "MB_dca3_noVPD",
                              "MB_dca3_withVPD", "MB_dca1_withVPD"};
  // vector<TString> trigger_title = {"DCA<1, no VPD cut", "DCA<3, no VPD cut",
  //                                  "DCA<3, with VPD cut",
  //                                  "DCA<1, with VPD cut"};

  vector<TString> trigger_title = {"+ DCA < 1cm", "Standard",
                                   "+ VPD vertex cut",
                                   "+ DCA < 1cm & VPD vertex cut"};

  for (size_t i = 0; i < triggers.size(); i++) {
    draw(triggers[i], trigger_title[i]);
  }
}