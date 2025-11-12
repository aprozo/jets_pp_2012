#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "config.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace CrossSectionConfig;

// Helper function to load run mapping from file
std::map<int, int> loadRunMap(const std::string &filename,
                              std::map<int, int> &inverseMap) {
  std::map<int, int> runMap;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open run map file: " << filename << std::endl;
    return runMap;
  }

  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    int runid, index;
    if (iss >> runid >> index) {
      runMap[runid] = index;
      inverseMap[index] = runid;
    }
  }
  file.close();
  std::cout << "  Loaded " << runMap.size() << " run mappings" << std::endl;
  return runMap;
}

// Helper function to divide 2D histogram by 1D histogram (per run)
TH2D *divide2DBy1D(TH2D *num, TH1D *denom) {
  TH2D *result = (TH2D *)num->Clone();
  for (int j = 1; j <= result->GetNbinsY(); j++) {
    for (int i = 1; i <= result->GetNbinsX(); i++) {
      double denomContent = denom->GetBinContent(i);
      if (denomContent > 0) {
        double scale = 1.0 / denomContent;
        result->SetBinContent(i, j, result->GetBinContent(i, j) * scale);
        result->SetBinError(i, j, result->GetBinError(i, j) * scale);
      } else {
        result->SetBinContent(i, j, 0);
        result->SetBinError(i, j, 0);
      }
    }
  }
  result->SetDirectory(0);
  return result;
}

// Helper function to multiply 2D histogram by 1D histogram (per run)
TH2D *multiply2DBy1D(TH2D *num, TH1D *factor) {
  TH2D *result = (TH2D *)num->Clone();
  for (int j = 1; j <= result->GetNbinsY(); j++) {
    for (int i = 1; i <= result->GetNbinsX(); i++) {
      double scale = factor->GetBinContent(i);
      result->SetBinContent(i, j, result->GetBinContent(i, j) * scale);
      result->SetBinError(i, j, result->GetBinError(i, j) * scale);
    }
  }
  result->SetDirectory(0);
  return result;
}

// Make run projection with bad run filtering and averaging
TH1D *makeRunProjection(TH2D *hist, const std::vector<int> &badRuns,
                        const std::map<int, int> &inverseRunMap) {
  std::string name = std::string(hist->GetName()) + "_projection";
  TH1D *output = hist->ProjectionY(name.c_str(), 1, 1);
  output->Reset();

  int counterRuns = 0;
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    // Get run ID from bin label
    const char *label = hist->GetXaxis()->GetBinLabel(i);
    if (!label || strlen(label) == 0)
      continue;

    int runid = std::atoi(label);

    // Check if run is in bad runs list
    if (std::find(badRuns.begin(), badRuns.end(), runid) != badRuns.end()) {
      continue;
    }

    // Get projection for this run
    TH1D *projection =
        hist->ProjectionY(Form("%s_run%d", name.c_str(), i), i, i);
    if (projection->Integral() == 0) {
      delete projection;
      continue;
    }

    counterRuns++;
    output->Add(projection);
    delete projection;
  }

  if (counterRuns > 0) {
    output->Scale(1.0 / counterRuns);
  }
  output->SetDirectory(0);

  std::cout << "  Averaged over " << counterRuns << " good runs" << std::endl;
  return output;
}

// Helper function to divide by bin width
void divideByBinWidth(TH1D *hist) {
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    double binWidth = hist->GetBinWidth(i);
    if (binWidth > 0) {
      hist->SetBinContent(i, hist->GetBinContent(i) / binWidth);
      hist->SetBinError(i, hist->GetBinError(i) / binWidth);
    }
  }
}

// Helper function to divide by efficiency with error propagation
void divideByEfficiency(TH1D *data, TH1D *efficiency) {
  for (int i = 1; i <= data->GetNbinsX(); i++) {
    double binContent = data->GetBinContent(i);
    double binError = data->GetBinError(i);
    double binCenter = data->GetBinCenter(i);

    // Find corresponding efficiency bin
    int effBin = efficiency->FindBin(binCenter);
    double eff = efficiency->GetBinContent(effBin);
    double effError = efficiency->GetBinError(effBin);

    if (eff > 0) {
      // Correct by efficiency
      data->SetBinContent(i, binContent / eff);

      // Propagate uncertainties
      if (binContent > 0) {
        double relError = std::sqrt(std::pow(binError / binContent, 2) +
                                    std::pow(effError / eff, 2));
        data->SetBinError(i, (binContent / eff) * relError);
      }
    }
  }
}

// Calculate vertex acceptance ratio (event_ratio) and save QA plot
double calculateVertexRatio(const std::string &filename,
                            const std::string &outputPdf = "") {
  TFile *file = TFile::Open(filename.c_str());
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file for vertex ratio: " << filename
              << std::endl;
    return 1.0;
  }

  TH1D *vertexHist = (TH1D *)file->Get("QA_histograms/vz");
  if (!vertexHist) {
    std::cerr << "Error: Cannot find QA_histograms/vz in " << filename
              << std::endl;
    file->Close();
    return 1.0;
  }
  vertexHist->SetDirectory(0);

  // Rebin
  vertexHist->Rebin(2);

  // Fit gaussian from -30 to 30
  TF1 *gausFit = new TF1("gausFit", "gaus", -30, 30);
  gausFit->SetLineColor(kPink);
  gausFit->SetLineWidth(2);
  vertexHist->Fit(gausFit, "REMQ");

  // Calculate ratio
  double binWidth = vertexHist->GetBinWidth(2);
  int bin1 = vertexHist->FindBin(-30);
  int bin2 = vertexHist->FindBin(30);
  double integral = vertexHist->Integral(bin1, bin2) * binWidth;

  double totalIntegral =
      std::sqrt(2 * M_PI) * gausFit->GetParameter(0) * gausFit->GetParameter(2);
  double ratio = integral / totalIntegral;

  std::cout << "  Vertex acceptance ratio: " << ratio << " (" << 100 * ratio
            << "%)" << std::endl;

  // Create QA plot
  if (!outputPdf.empty()) {
    TCanvas *canvas = new TCanvas("cVertex", "Vertex QA", 800, 600);
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.05);

    vertexHist->SetLineColor(kBlack);
    vertexHist->SetLineWidth(2);
    vertexHist->SetTitle("Vertex z Distribution;v_{z} [cm];Counts");
    vertexHist->Draw("HIST");
    gausFit->Draw("SAME");

    // Add text with fit results
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.035);
    latex->DrawLatex(0.15, 0.85, "Gaussian fit: [-30, 30] cm");
    latex->DrawLatex(0.15, 0.80,
                     Form("#mu = %.2f cm", gausFit->GetParameter(1)));
    latex->DrawLatex(0.15, 0.75,
                     Form("#sigma = %.2f cm", gausFit->GetParameter(2)));
    latex->DrawLatex(0.15, 0.70, Form("Acceptance ratio = %.4f", ratio));

    // Draw acceptance region
    TLine *line = new TLine();
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    double ymax = vertexHist->GetMaximum();
    line->DrawLine(-30, 0, -30, ymax);
    line->DrawLine(30, 0, 30, ymax);

    canvas->SaveAs(outputPdf.c_str());
    delete canvas;
    delete latex;
    delete line;
  }

  delete gausFit;
  delete vertexHist;
  file->Close();
  return ratio;
}

// QA plot for 2D pt vs run histogram
void plotPtVsRun(TH2D *hist, const std::string &outputPdf) {
  TCanvas *canvas = new TCanvas("cPtVsRun", "pT vs Run", 1200, 600);
  canvas->SetLeftMargin(0.08);
  canvas->SetRightMargin(0.15);
  canvas->SetLogz();

  hist->SetTitle(";Run bin;p_{T} [GeV/c];Counts");
  hist->Draw("COLZ");

  canvas->SaveAs(outputPdf.c_str());
  delete canvas;
}

// QA plot for luminosity data
void plotLuminosityData(TH1D *prescale, TH1D *livetime, TH1D *maxnevent,
                        TH1D *luminosity, TH1D *analyzedEvents,
                        const std::string &outputPdf) {
  TCanvas *canvas = new TCanvas("cLumi", "Luminosity Data", 1200, 1000);
  canvas->Divide(2, 3);

  canvas->cd(1);
  gPad->SetLeftMargin(0.12);
  analyzedEvents->SetLineColor(kBlack);
  analyzedEvents->SetLineWidth(2);
  analyzedEvents->SetTitle("Analyzed Events;Run bin;Events");
  analyzedEvents->Draw("HIST");

  canvas->cd(2);
  gPad->SetLeftMargin(0.12);
  prescale->SetLineColor(kBlue);
  prescale->SetLineWidth(2);
  prescale->SetTitle("Prescale;Run bin;Prescale");
  prescale->Draw("HIST");

  canvas->cd(3);
  gPad->SetLeftMargin(0.12);
  livetime->SetLineColor(kGreen + 2);
  livetime->SetLineWidth(2);
  livetime->SetTitle("Live Time;Run bin;Live Time");
  livetime->Draw("HIST");

  canvas->cd(4);
  gPad->SetLeftMargin(0.12);
  maxnevent->SetLineColor(kRed);
  maxnevent->SetLineWidth(2);
  maxnevent->SetTitle("Max N Events;Run bin;Events");
  maxnevent->Draw("HIST");

  canvas->cd(5);
  gPad->SetLeftMargin(0.12);
  luminosity->SetLineColor(kMagenta);
  luminosity->SetLineWidth(2);
  luminosity->SetTitle("Luminosity;Run bin;L [pb^{-1}]");
  luminosity->Draw("HIST");

  canvas->SaveAs(outputPdf.c_str());
  delete canvas;
}

// QA plot comparing spectrum before/after correction
void plotBeforeAfterCorrection(TH1D *before, TH1D *after,
                               const std::string &title,
                               const std::string &correctionName,
                               const std::string &outputPdf,
                               TH1D *correction = nullptr) {
  TCanvas *canvas = new TCanvas("cCorrection", title.c_str(), 1000, 800);
  canvas->Divide(1, 2);

  // Upper pad - spectra
  TPad *pad1 = (TPad *)canvas->cd(1);
  pad1->SetPad(0.0, 0.4, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLogy();
  pad1->SetLeftMargin(0.12);

  TH1D *hBefore = (TH1D *)before->Clone("hBefore");
  TH1D *hAfter = (TH1D *)after->Clone("hAfter");

  hBefore->SetLineColor(kBlue);
  hBefore->SetMarkerColor(kBlue);
  hBefore->SetMarkerStyle(20);
  hBefore->SetTitle(Form("%s;p_{T} [GeV/c];dN/dp_{T}", title.c_str()));
  hBefore->Draw("E1");

  hAfter->SetLineColor(kRed);
  hAfter->SetMarkerColor(kRed);
  hAfter->SetMarkerStyle(21);
  hAfter->Draw("E1 SAME");

  TLegend *leg1 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg1->AddEntry(hBefore, "Before", "lep");
  leg1->AddEntry(hAfter, "After", "lep");
  leg1->Draw();

  // Lower pad - ratio or correction
  TPad *pad2 = (TPad *)canvas->cd(2);
  pad2->SetPad(0.0, 0.0, 1.0, 0.4);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.12);

  if (correction) {
    TH1D *hCorr = (TH1D *)correction->Clone("hCorr");
    hCorr->SetLineColor(kBlack);
    hCorr->SetMarkerColor(kBlack);
    hCorr->SetMarkerStyle(22);
    hCorr->SetTitle(Form(";p_{T} [GeV/c];%s", correctionName.c_str()));
    hCorr->GetYaxis()->SetTitleSize(0.06);
    hCorr->GetYaxis()->SetLabelSize(0.05);
    hCorr->GetXaxis()->SetTitleSize(0.06);
    hCorr->GetXaxis()->SetLabelSize(0.05);
    hCorr->Draw("E1");
  } else {
    TH1D *ratio = (TH1D *)hAfter->Clone("ratio");
    ratio->Divide(hBefore);
    ratio->SetTitle(";p_{T} [GeV/c];After/Before");
    ratio->GetYaxis()->SetTitleSize(0.06);
    ratio->GetYaxis()->SetLabelSize(0.05);
    ratio->GetXaxis()->SetTitleSize(0.06);
    ratio->GetXaxis()->SetLabelSize(0.05);
    ratio->Draw("E1");
  }

  canvas->SaveAs(outputPdf.c_str());
  delete canvas;
  delete hBefore;
  delete hAfter;
  delete leg1;
}

// Main cross-section calculation function
void calculateCrossSection(const std::string &trigger = "JP2",
                           const std::string &jetR = "0.6") {

  std::cout << "\n========================================" << std::endl;
  std::cout << "Calculating cross-section for:" << std::endl;
  std::cout << "  Trigger: " << trigger << std::endl;
  std::cout << "  Jet R: " << jetR << std::endl;
  std::cout << "========================================\n" << std::endl;

  AnalysisConfig cfg(jetR);
  gStyle->SetOptStat(0);

  TString outputDir = Form("%s/results_%s_R%s", cfg.workdir.c_str(),
                           trigger.c_str(), jetR.c_str());
  gSystem->mkdir(outputDir.Data(), kTRUE);

  // =======================================================================
  // Step 1: Load run mapping
  // =======================================================================
  std::cout << "Step 1: Loading run mapping..." << std::endl;
  std::map<int, int> inverseRunMap;
  std::map<int, int> runMap =
      loadRunMap("/home/prozorov/dev/star/jets_pp_2012/macros/ana/run_map.txt",
                 inverseRunMap);

  // =======================================================================
  // Step 2: Calculate vertex acceptance ratio
  // =======================================================================
  std::cout << "\nStep 2: Calculating vertex acceptance..." << std::endl;
  TString dataFileName =
      Form("%s/merged_data_%s_R%s.root", cfg.datapath.c_str(), trigger.c_str(),
           jetR.c_str());
  std::cout << "  Data file: " << dataFileName << std::endl;
  TString vertexQaFile =
      Form("%s/QA_01_vertex_acceptance.pdf", outputDir.Data());
  double eventRatio =
      calculateVertexRatio(dataFileName.Data(), vertexQaFile.Data());

  // =======================================================================
  // Step 3: Create 2D histograms (pt vs run) from trees - trigger matched
  // =======================================================================
  std::cout << "\nStep 3: Creating 2D histograms from ROOT trees..."
            << std::endl;

  // Get number of run bins from analyzed events histogram
  TFile *dataFileTemp = TFile::Open(dataFileName);
  if (!dataFileTemp || dataFileTemp->IsZombie()) {
    std::cerr << "Error: Cannot open data file: " << dataFileName << std::endl;
    return;
  }
  TH1D *hEventsRunTemp = (TH1D *)dataFileTemp->Get("hEventsRun");
  if (!hEventsRunTemp) {
    std::cerr << "Error: Cannot find hEventsRun" << std::endl;
    dataFileTemp->Close();
    return;
  }
  int nRunBins = hEventsRunTemp->GetNbinsX();
  dataFileTemp->Close();

  std::cout << "  Reading from: " << dataFileName << std::endl;
  std::cout << "  Number of run bins: " << nRunBins << std::endl;

  // Create RDataFrame from ROOT file
  ROOT::RDataFrame df("ResultTree", dataFileName.Data());

  // Map runid to bin index
  auto dfMapped = df.Define("run_mapped",
                            [&runMap](int runid) {
                              auto it = runMap.find(runid);
                              return (it != runMap.end()) ? it->second : -1;
                            },
                            {"runid1"});

  // Define pt bins for histogram
  std::vector<double> ptBins(pt_reco_bins.begin(), pt_reco_bins.end());

  TH2D *ptVsRunTrigger = nullptr;

  if (trigger == "JP2" || trigger == "HT2") {
    // For triggered jets - use trigger matching
    std::string triggerBranch =
        "trigger_match_" + std::string(trigger == "JP2" ? "JP2" : "HT2");
    std::cout << "  Applying trigger matching: " << triggerBranch << std::endl;

    auto dfTriggered = dfMapped.Define("triggered_mask", triggerBranch)
                           .Define("triggered_pt", "pt[triggered_mask]")
                           .Filter("triggered_pt.size() > 0");

    // Create 2D histogram
    auto histPtr = dfTriggered.Histo2D(
        {Form("%s_pt_vs_run_trigger", trigger.c_str()),
         Form("%s pt vs run;run bin;p_{T} [GeV/c];counts", trigger.c_str()),
         nRunBins, 0.0, static_cast<double>(nRunBins), (int)ptBins.size() - 1,
         ptBins.data()},
        "run_mapped", "triggered_pt");

    ptVsRunTrigger = (TH2D *)histPtr->Clone();
    ptVsRunTrigger->SetDirectory(0);

  } else {
    // For MB - no trigger matching needed
    std::cout << "  No trigger matching (MB trigger)" << std::endl;

    auto histPtr = dfMapped.Histo2D(
        {Form("%s_pt_vs_run", trigger.c_str()),
         Form("%s pt vs run;run bin;p_{T} [GeV/c];counts", trigger.c_str()),
         nRunBins, 0.0, static_cast<double>(nRunBins), (int)ptBins.size() - 1,
         ptBins.data()},
        "run_mapped", "pt");

    ptVsRunTrigger = (TH2D *)histPtr->Clone();
    ptVsRunTrigger->SetDirectory(0);
  }

  // Set bin labels from inverse run map
  for (int i = 1; i <= nRunBins; i++) {
    auto it = inverseRunMap.find(i - 1);
    if (it != inverseRunMap.end()) {
      ptVsRunTrigger->GetXaxis()->SetBinLabel(i, Form("%d", it->second));
    }
  }

  std::cout << "  Created histogram with " << ptVsRunTrigger->GetEntries()
            << " entries" << std::endl;

  // Create QA plot for pt vs run
  std::cout << "  Creating QA plot..." << std::endl;
  TString ptVsRunQaFile = Form("%s/QA_02_pt_vs_run.pdf", outputDir.Data());
  plotPtVsRun(ptVsRunTrigger, ptVsRunQaFile.Data());

  // =======================================================================
  // Step 4: Load luminosity data
  // =======================================================================
  std::cout << "\nStep 4: Loading luminosity data..." << std::endl;
  TString lumiFileName =
      "/home/prozorov/dev/star/jets_pp_2012/macros/ana/lumi.root";
  TFile *lumiFile = TFile::Open(lumiFileName);
  if (!lumiFile || lumiFile->IsZombie()) {
    std::cerr << "Error: Cannot open lumi file: " << lumiFileName << std::endl;
    return;
  }

  TH1D *prescale = (TH1D *)lumiFile->Get(Form("prescale_%s", trigger.c_str()));
  TH1D *livetime = (TH1D *)lumiFile->Get(Form("livetime_%s", trigger.c_str()));
  TH1D *maxnevent = (TH1D *)lumiFile->Get(Form("nevents_%s", trigger.c_str()));
  TH1D *luminosity =
      (TH1D *)lumiFile->Get(Form("luminosity_%s", trigger.c_str()));

  if (!prescale || !livetime || !maxnevent || !luminosity) {
    std::cerr << "Error: Cannot load luminosity histograms for " << trigger
              << std::endl;
    lumiFile->Close();
    return;
  }

  prescale->SetDirectory(0);
  livetime->SetDirectory(0);
  maxnevent->SetDirectory(0);
  luminosity->SetDirectory(0);
  std::cout << "  Loaded prescale, livetime, nevents, luminosity" << std::endl;
  lumiFile->Close();

  // =======================================================================
  // Step 5: Load analyzed events histogram
  // =======================================================================
  std::cout << "\nStep 5: Loading analyzed events..." << std::endl;
  TFile *dataFile = TFile::Open(dataFileName);
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "Error: Cannot open data file: " << dataFileName << std::endl;
    return;
  }

  TH1D *analyzedEvents = (TH1D *)dataFile->Get("hEventsRun");
  if (!analyzedEvents) {
    std::cerr << "Error: Cannot find hEventsRun in " << dataFileName
              << std::endl;
    dataFile->Close();
    return;
  }
  analyzedEvents->SetDirectory(0);
  std::cout << "  Loaded analyzed events histogram" << std::endl;
  dataFile->Close();

  // Create QA plot for luminosity data
  std::cout << "  Creating luminosity QA plot..." << std::endl;
  TString lumiQaFile = Form("%s/QA_03_luminosity_data.pdf", outputDir.Data());
  plotLuminosityData(prescale, livetime, maxnevent, luminosity, analyzedEvents,
                     lumiQaFile.Data());

  // =======================================================================
  // Step 6: Apply per-run normalization
  // =======================================================================
  std::cout << "\nStep 6: Applying per-run normalization..." << std::endl;

  // Normalize by analyzed events
  std::cout << "  Dividing by analyzed events..." << std::endl;
  TH2D *temp = divide2DBy1D(ptVsRunTrigger, analyzedEvents);

  // Multiply by maxnevent (prescale correction)
  std::cout << "  Multiplying by maxnevent (prescale)..." << std::endl;
  TH2D *temp2 = multiply2DBy1D(temp, maxnevent);
  delete temp;
  temp = temp2;

  // Divide by luminosity
  std::cout << "  Dividing by luminosity..." << std::endl;
  temp2 = divide2DBy1D(temp, luminosity);
  delete temp;
  temp = temp2;

  // Scale by vertex acceptance ratio
  std::cout << "  Scaling by vertex acceptance (1/" << eventRatio << ")..."
            << std::endl;
  // temp->Scale(1.0 / eventRatio);

  // =======================================================================
  // Step 7: Make run projection (average over good runs)
  // =======================================================================
  std::cout << "\nStep 7: Making run projection..." << std::endl;
  TH1D *spectrum = makeRunProjection(temp, cfg.badRunsJP2, inverseRunMap);
  delete temp;

  // Scale by eta acceptance
  std::cout << "  Scaling by eta acceptance (2*(" << (1 - std::stof(jetR))
            << ") = " << cfg.etaAcceptance << ")" << std::endl;
  spectrum->Scale(1.0 / (cfg.etaAcceptance));

  // =======================================================================
  // Step 8: Load trigger and reconstruction efficiencies
  // =======================================================================
  std::cout << "\nStep 8: Loading efficiencies..." << std::endl;

  // Trigger efficiency (reco level)
  TH1D *triggerEfficiency = nullptr;
  if (trigger == "JP2" || trigger == "HT2") {
    TString trigEffFileName =
        Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/matching/"
             "plots/embedding_root_%s_R%s/"
             "%s_over_all.root",
             trigger.c_str(), jetR.c_str(), trigger.c_str());
    TFile *trigEffFile = TFile::Open(trigEffFileName);
    if (trigEffFile && !trigEffFile->IsZombie()) {
      triggerEfficiency = (TH1D *)trigEffFile->Get("h_ratio");
      if (triggerEfficiency) {
        triggerEfficiency->SetDirectory(0);
        std::cout << "  Loaded trigger efficiency (reco)" << std::endl;
      }
      trigEffFile->Close();
    }
  }

  // fakerate efficiency (reco level)
  TH1D *fakeRateEfficiency = nullptr;
  if (trigger == "JP2" || trigger == "HT2") {
    TString fakeRateFileName =
        Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/matching/"
             "plots/embedding_root_%s_R%s/"
             "fakereco_over_all.root",
             trigger.c_str(), jetR.c_str());
    TFile *fakeRateFile = TFile::Open(fakeRateFileName);
    if (fakeRateFile && !fakeRateFile->IsZombie()) {
      fakeRateEfficiency = (TH1D *)fakeRateFile->Get("h_ratio");
      if (fakeRateEfficiency) {
        fakeRateEfficiency->SetDirectory(0);
        std::cout << "  Loaded fake rate efficiency (reco)" << std::endl;
      }
      fakeRateFile->Close();
    }
  }

  // Reconstruction efficiency (MC level)
  TString recoEffFileName =
      Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/matching/"
           "plots/embedding_root_%s_R%s/"
           "matched_over_mc.root",
           trigger.c_str(), jetR.c_str());
  TFile *recoEffFile = TFile::Open(recoEffFileName);
  TH1D *recoEfficiency = nullptr;
  if (recoEffFile && !recoEffFile->IsZombie()) {
    recoEfficiency = (TH1D *)recoEffFile->Get("h_ratio");
    if (recoEfficiency) {
      recoEfficiency->SetDirectory(0);
      std::cout << "  Loaded reconstruction efficiency (MC)" << std::endl;
    }
    recoEffFile->Close();
  }

  // =======================================================================
  // Step 9: Apply trigger efficiency correction
  // =======================================================================
  std::cout << "\nStep 9: Applying trigger efficiency correction..."
            << std::endl;

  // Save spectrum before trigger efficiency correction for QA
  TH1D *spectrumBeforeTrigEff =
      (TH1D *)spectrum->Clone("spectrumBeforeTrigEff");

  // Apply trigger efficiency correction
  if (triggerEfficiency && (trigger == "JP2" || trigger == "HT2")) {
    std::cout << "  Dividing by trigger efficiency..." << std::endl;
    divideByEfficiency(spectrum, triggerEfficiency);

    // Create QA plot
    std::cout << "  Creating trigger efficiency QA plot..." << std::endl;
    TString trigEffQaFile =
        Form("%s/QA_04_trigger_efficiency.pdf", outputDir.Data());
    plotBeforeAfterCorrection(
        spectrumBeforeTrigEff, spectrum, "Trigger Efficiency Correction",
        "Trigger Efficiency", trigEffQaFile.Data(), triggerEfficiency);
  }
  delete spectrumBeforeTrigEff;

  TH1D *spectrumBeforeFakeRate =
      (TH1D *)spectrum->Clone("spectrumBeforeFakeRate");

  // apply fake rate subtraction
  if (trigger == "JP2" || trigger == "HT2") {
    std::cout << "  Applying fake rate subtraction..." << std::endl;
    for (int i = 1; i <= spectrum->GetNbinsX(); i++) {
      double binCenter = spectrum->GetBinCenter(i);
      int effBin = fakeRateEfficiency->FindBin(binCenter);
      double fakeRate = fakeRateEfficiency->GetBinContent(effBin);
      double binContent = spectrum->GetBinContent(i);
      double binError = spectrum->GetBinError(i);

      // Apply fake rate subtraction
      double correctedContent = binContent * (1 - fakeRate);
      double correctedError = binError * (1 - fakeRate);

      spectrum->SetBinContent(i, correctedContent);
      spectrum->SetBinError(i, correctedError);
    }
    // Create QA plot
    std::cout << "  Creating fakerate efficiency QA plot..." << std::endl;
    TString fakeRateQaFile =
        Form("%s/QA_04_fakerate_efficiency.pdf", outputDir.Data());
    plotBeforeAfterCorrection(
        spectrumBeforeFakeRate, spectrum, "Fakerate Efficiency Correction",
        "1 - Fakerate Efficiency", fakeRateQaFile.Data(), fakeRateEfficiency);
  }

  // =======================================================================
  // Step 10: Apply unfolding
  // =======================================================================
  std::cout << "\nStep 10: Loading response matrix and performing unfolding..."
            << std::endl;

  // Save spectrum before unfolding for QA
  TH1D *spectrumBeforeUnfold = (TH1D *)spectrum->Clone("spectrumBeforeUnfold");

  // Load response matrix
  TString responseFileName =
      Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/unfolding/"
           "response_%s_R%s.root",
           trigger.c_str(), jetR.c_str());
  std::cout << "  Response file: " << responseFileName << std::endl;

  TFile *responseFile = TFile::Open(responseFileName);
  if (!responseFile || responseFile->IsZombie()) {
    std::cerr << "Error: Cannot open response file: " << responseFileName
              << std::endl;
    std::cerr << "  Skipping unfolding step" << std::endl;
    // Continue without unfolding
  } else {
    RooUnfoldResponse *response =
        (RooUnfoldResponse *)responseFile->Get("my_response");
    if (!response) {
      std::cerr << "Error: Cannot find response in file" << std::endl;
      std::cerr << "  Skipping unfolding step" << std::endl;
      responseFile->Close();
    } else {
      std::cout << "  Loaded response matrix" << std::endl;
      std::cout << "  Performing Bayesian unfolding with " << cfg.nIterations
                << " iterations..." << std::endl;

      // Perform unfolding
      RooUnfoldBayes unfold(response, spectrum, cfg.nIterations);
      TH1D *unfolded = (TH1D *)unfold.Hunfold();
      unfolded->SetName("spectrum_unfolded");

      // Get covariance matrix
      TMatrixD covMatrix = unfold.Eunfold(RooUnfold::kCovariance);

      std::cout << "  Unfolding complete" << std::endl;
      std::cout << "  Entries before unfolding: " << spectrum->GetEntries()
                << std::endl;
      std::cout << "  Entries after unfolding: " << unfolded->GetEntries()
                << std::endl;

      // Replace spectrum with unfolded
      delete spectrum;
      spectrum = (TH1D *)unfolded->Clone("spectrum_after_unfolding");
      spectrum->SetDirectory(0);

      // Create QA plot
      std::cout << "  Creating unfolding QA plot..." << std::endl;
      TString unfoldQaFile = Form("%s/QA_05_unfolding.pdf", outputDir.Data());
      plotBeforeAfterCorrection(spectrumBeforeUnfold, spectrum,
                                "Unfolding (Bayesian)", "Unfolding Ratio",
                                unfoldQaFile.Data());

      responseFile->Close();
    }
  }
  delete spectrumBeforeUnfold;

  // =======================================================================
  // Step 11: Apply reconstruction efficiency correction
  // =======================================================================
  std::cout << "\nStep 11: Applying reconstruction efficiency correction..."
            << std::endl;

  // Save spectrum before reconstruction efficiency correction for QA
  TH1D *spectrumBeforeRecoEff =
      (TH1D *)spectrum->Clone("spectrumBeforeRecoEff");

  if (recoEfficiency) {
    std::cout << "  Dividing by reconstruction efficiency..." << std::endl;
    divideByEfficiency(spectrum, recoEfficiency);

    // Create QA plot
    std::cout << "  Creating reconstruction efficiency QA plot..." << std::endl;
    TString recoEffQaFile =
        Form("%s/QA_06_reconstruction_efficiency.pdf", outputDir.Data());
    plotBeforeAfterCorrection(
        spectrumBeforeRecoEff, spectrum, "Reconstruction Efficiency Correction",
        "Reconstruction Efficiency", recoEffQaFile.Data(), recoEfficiency);
  }
  delete spectrumBeforeRecoEff;

  // =======================================================================
  // Step 12: Divide by bin width
  // =======================================================================
  std::cout << "\nStep 12: Dividing by bin width..." << std::endl;
  TH1D *hCrossSection = (TH1D *)spectrum->Clone("hCrossSection");
  divideByBinWidth(hCrossSection);

  // Create final QA plot
  std::cout << "  Creating final cross-section QA plot..." << std::endl;
  TString finalQaFile =
      Form("%s/QA_07_final_cross_section.pdf", outputDir.Data());
  TCanvas *cFinal = new TCanvas("cFinal", "Final Cross Section", 800, 600);
  cFinal->SetLogy();
  cFinal->SetLeftMargin(0.12);

  hCrossSection->SetLineColor(kBlack);
  hCrossSection->SetMarkerColor(kBlack);
  hCrossSection->SetMarkerStyle(20);
  hCrossSection->SetTitle(Form("%s R=%s Cross Section;p_{T} "
                               "[GeV/c];d^{2}#sigma/dp_{T}d#eta [pb/(GeV/c)]",
                               trigger.c_str(), jetR.c_str()));
  hCrossSection->Draw("E1");

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->DrawLatex(0.15, 0.85, Form("Trigger: %s", trigger.c_str()));
  latex->DrawLatex(0.15, 0.80, Form("Jet R = %s", jetR.c_str()));
  latex->DrawLatex(0.15, 0.75, "pp #sqrt{s} = 200 GeV");

  cFinal->SaveAs(finalQaFile.Data());
  delete cFinal;
  delete latex;

  // =======================================================================
  // Step 13: Save results
  // =======================================================================
  std::cout << "\nStep 13: Saving results..." << std::endl;

  TString outFileName =
      Form("%s/cross_section_%s_R%s.root", cfg.workdir.c_str(), trigger.c_str(),
           jetR.c_str());
  TFile *outFile = TFile::Open(outFileName, "RECREATE");

  // Save all relevant histograms
  spectrum->Write("spectrum_before_binwidth");
  hCrossSection->Write("cross_section");
  ptVsRunTrigger->Write("pt_vs_run_trigger");
  analyzedEvents->Write("analyzed_events");
  prescale->Write("prescale");
  livetime->Write("livetime");
  maxnevent->Write("maxnevent");
  luminosity->Write("luminosity");

  if (triggerEfficiency)
    triggerEfficiency->Write("trigger_efficiency");
  if (recoEfficiency)
    recoEfficiency->Write("reconstruction_efficiency");

  outFile->Close();
  std::cout << "  Results saved to: " << outFileName << std::endl;

  // Cleanup
  delete spectrum;
  delete hCrossSection;
  if (triggerEfficiency)
    delete triggerEfficiency;
  if (recoEfficiency)
    delete recoEfficiency;

  std::cout << "\n========================================" << std::endl;
  std::cout << "Cross-section calculation complete!" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "\nQA plots created:" << std::endl;
  std::cout << "  1. " << vertexQaFile << std::endl;
  std::cout << "  2. " << ptVsRunQaFile << std::endl;
  std::cout << "  3. " << lumiQaFile << std::endl;
  if (trigger == "JP2" || trigger == "HT2") {
    TString trigEffQaFile =
        Form("%s/QA_04_trigger_efficiency.pdf", outputDir.Data());
    std::cout << "  4. " << trigEffQaFile << std::endl;
  }
  TString unfoldQaFile = Form("%s/QA_05_unfolding.pdf", outputDir.Data());
  std::cout << "  5. " << unfoldQaFile << std::endl;
  TString recoEffQaFile =
      Form("%s/QA_06_reconstruction_efficiency.pdf", outputDir.Data());
  std::cout << "  6. " << recoEffQaFile << std::endl;
  std::cout << "  7. " << finalQaFile << std::endl;
  std::cout << "========================================\n" << std::endl;
}

// Helper function to divide histograms with different binning (using
// interpolation)
TH1D *divideDifferentBins(TH1D *num, TH1D *denom) {
  TH1D *result = (TH1D *)num->Clone();
  for (int i = 1; i <= result->GetNbinsX(); i++) {
    double binCenter = result->GetBinCenter(i);
    double denomValue = denom->Interpolate(binCenter);
    if (denomValue > 0) {
      result->SetBinContent(i, result->GetBinContent(i) / denomValue);
      result->SetBinError(i, result->GetBinError(i) / denomValue);
    } else {
      result->SetBinContent(i, 0);
      result->SetBinError(i, 0);
    }
  }
  result->SetDirectory(0);
  return result;
}

// Compare with Dmitriy's reference cross-section
void compareCrossSections(const std::string &jetR = "0.6") {
  std::cout << "\n========================================" << std::endl;
  std::cout << "Comparing cross-sections with Dmitriy's result" << std::endl;
  std::cout << "  Jet R: " << jetR << std::endl;
  std::cout << "========================================\n" << std::endl;

  AnalysisConfig cfg(jetR);

  // Load Dmitriy's reference
  TString refFileName =
      Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/"
           "jet_cross_section_dmitriyR%s.root",
           jetR.c_str());
  TFile *refFile = TFile::Open(refFileName);
  if (!refFile || refFile->IsZombie()) {
    std::cerr << "Error: Cannot open reference file: " << refFileName
              << std::endl;
    return;
  }

  TH1D *refCrossSection = (TH1D *)refFile->Get("crossSection_systematic");
  if (!refCrossSection) {
    std::cerr << "Error: Cannot find crossSection_systematic" << std::endl;
    refFile->Close();
    return;
  }
  refCrossSection->SetDirectory(0);
  refFile->Close();

  // Create canvas for comparison
  TCanvas *canvas =
      new TCanvas("canvasCS", "Cross-Section Comparison", 800, 800);
  canvas->Divide(1, 2);

  // Upper pad - cross-sections
  TPad *pad1 = (TPad *)canvas->cd(1);
  pad1->SetPad(0.0, 0.5, 1.0, 1.0);
  pad1->SetTopMargin(0.1);
  pad1->SetBottomMargin(0.0);
  pad1->SetLogy();

  refCrossSection->GetYaxis()->SetRangeUser(1e0 + 0.2, 1e7);
  refCrossSection->SetMarkerColor(kViolet);
  refCrossSection->SetLineColor(kViolet);
  refCrossSection->SetMarkerStyle(29);
  refCrossSection->SetMarkerSize(2);
  refCrossSection->Draw("E1");

  TLegend *legend = new TLegend(0.6, 0.5, 0.88, 0.88);
  legend->AddEntry(refCrossSection, "Dmitriy", "lep");

  // Load and plot our results
  std::vector<std::string> triggers = {"JP2", "HT2", "MB"};
  std::vector<int> colors = {kAzure - 1, kRed + 1, kOrange - 1};
  std::vector<int> markers = {20, 21, 22};
  std::map<std::string, TH1D *> myCrossSections;
  std::map<std::string, TH1D *> ratios;

  for (size_t i = 0; i < triggers.size(); i++) {
    std::string trigger = triggers[i];
    TString fileName = Form("%s/cross_section_%s_R%s.root", cfg.workdir.c_str(),
                            trigger.c_str(), jetR.c_str());
    TFile *file = TFile::Open(fileName);
    if (file && !file->IsZombie()) {
      TH1D *cs = (TH1D *)file->Get("cross_section");
      if (cs) {
        cs->SetDirectory(0);
        cs->SetLineColor(colors[i]);
        cs->SetMarkerColor(colors[i]);
        cs->SetMarkerStyle(markers[i]);
        cs->Draw("E1 same");
        myCrossSections[trigger] = cs;
        legend->AddEntry(cs, trigger.c_str(), "lep");

        // Calculate ratio
        ratios[trigger] = divideDifferentBins(cs, refCrossSection);
      }
      file->Close();
    }
  }

  legend->Draw();

  // Lower pad - ratios
  TPad *pad2 = (TPad *)canvas->cd(2);
  pad2->SetPad(0.0, 0.0, 1.0, 0.5);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);

  bool firstRatio = true;
  for (size_t i = 0; i < triggers.size(); i++) {
    std::string trigger = triggers[i];
    if (ratios.find(trigger) != ratios.end()) {
      TH1D *ratio = ratios[trigger];
      ratio->GetYaxis()->SetRangeUser(0, 2.5);
      ratio->GetXaxis()->SetRangeUser(refCrossSection->GetXaxis()->GetXmin() +
                                          0.1,
                                      refCrossSection->GetXaxis()->GetXmax());
      ratio->GetYaxis()->SetTitle("Ratio to Dmitriy");
      ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      ratio->SetLineColor(colors[i]);
      ratio->SetMarkerColor(colors[i]);
      ratio->SetMarkerStyle(markers[i]);

      if (firstRatio) {
        ratio->Draw("E1");
        firstRatio = false;
      } else {
        ratio->Draw("E1 same");
      }
    }
  }

  // Draw line at 1
  TLine *line = new TLine();
  line->SetLineColor(kViolet);
  line->DrawLine(refCrossSection->GetXaxis()->GetXmin(), 1,
                 refCrossSection->GetXaxis()->GetXmax(), 1);

  // Save canvas
  TString outFileName = Form("%s/comparison_with_dmitriy_R%s.pdf",
                             cfg.workdir.c_str(), jetR.c_str());
  canvas->SaveAs(outFileName);
  std::cout << "  Comparison plot saved to: " << outFileName << std::endl;

  // Cleanup
  delete canvas;
  delete refCrossSection;
  for (auto &pair : myCrossSections) {
    delete pair.second;
  }
  for (auto &pair : ratios) {
    delete pair.second;
  }
}

// Process all triggers and jet radii
void cross_section() {
  AnalysisConfig cfg;

  std::vector<std::string> triggers = {"JP2", "HT2"};
  std::vector<std::string> jetRs = {"0.5", "0.6"};

  for (const auto &jetR : jetRs) {
    for (const auto &trigger : triggers) {
      calculateCrossSection(trigger, jetR);
    }
    // Create comparison plot after all triggers for this jetR
    compareCrossSections(jetR);
  }

  std::cout << "\n\nAll cross-sections calculated successfully!" << std::endl;
}
