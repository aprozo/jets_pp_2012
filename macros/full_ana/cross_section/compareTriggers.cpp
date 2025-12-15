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
#include <memory>
#include <set>
#include <sstream>
#include <vector>

using namespace CrossSectionConfig;

// Helper function to load run mapping from file
std::map<std::string, int> loadRunMap(const std::string &filename,
                                      std::map<int, std::string> &inverseMap) {
  std::map<std::string, int> runMap;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open run map file: " << filename << std::endl;
    return runMap;
  }
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string runid;
    int index;
    if (iss >> runid >> index) {
      runMap[runid] = index;
      inverseMap[index] = runid;
    }
  }
  file.close();
  return runMap;
}

void filterBadRuns(TH1D *hist, const std::vector<int> &badRuns) {
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    std::string run = hist->GetXaxis()->GetBinLabel(i);
    int runNumber = std::stoi(run);
    if (std::find(badRuns.begin(), badRuns.end(), runNumber) != badRuns.end()) {
      hist->SetBinContent(i, 0);
      hist->SetBinError(i, 0);
    }
  }
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

set<int> getCurrentRuns(TH1D *hEventsRun,
                        const std::map<std::string, int> &runMap) {
  set<int> runs;
  for (int i = 1; i <= hEventsRun->GetNbinsX(); i++) {
    if (hEventsRun->GetBinContent(i) <= 0)
      continue;
    std::string run = hEventsRun->GetXaxis()->GetBinLabel(i);
    auto it = runMap.find(run);
    if (it == runMap.end())
      continue; // Skip runs not present in the map
    int runIndex = it->second;
    runs.insert(runIndex);
  }
  return runs;
}

// void setSelectedRuns(TH1D *hEventsRun, set<int> &selectedRuns,
//                      const std::map<std::string, int> &runMap) {
//   for (int i = 1; i <= hEventsRun->GetNbinsX(); i++) {
//     std::string run = hEventsRun->GetXaxis()->GetBinLabel(i);
//     auto it = runMap.find(run);
//     if (it == runMap.end())
//       continue; // Skip runs not present in the map
//     int runIndex = it->second;
//     if (selectedRuns.find(runIndex) == selectedRuns.end()) {
//       hEventsRun->SetBinContent(i, 0);
//       hEventsRun->SetBinError(i, 0);
//     }
//   }
// }

// Helper function to divide by efficiency with error propagation
void divideByEfficiency(TH1D *data, TH1D *efficiency) {
  if (!efficiency) {
    std::cerr << "Error: Efficiency histogram is null, skipping 1D correction"
              << std::endl;
    return;
  }

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

void divideByEfficiency(TH2D *data, TH1D *efficiency) {
  if (!efficiency) {
    std::cerr << "Error: Efficiency histogram is null, skipping 2D correction"
              << std::endl;
    return;
  }

  for (int j = 1; j <= data->GetNbinsY(); j++) {
    for (int i = 1; i <= data->GetNbinsX(); i++) {
      double binContent = data->GetBinContent(i, j);
      double binError = data->GetBinError(i, j);
      double binCenter = data->GetXaxis()->GetBinCenter(i);

      // Find corresponding efficiency bin
      int effBin = efficiency->FindBin(binCenter);
      double eff = efficiency->GetBinContent(effBin);
      double effError = efficiency->GetBinError(effBin);

      if (eff > 0) {
        // Correct by efficiency
        data->SetBinContent(i, j, binContent / eff);

        // Propagate uncertainties
        if (binContent > 0) {
          double relError = std::sqrt(std::pow(binError / binContent, 2) +
                                      std::pow(effError / eff, 2));
          data->SetBinError(i, j, (binContent / eff) * relError);
        }
      }
    }
  }
}

void divideByBinWidth(TH1D *hist) {
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    double binWidth = hist->GetBinWidth(i);
    if (binWidth > 0) {
      hist->SetBinContent(i, hist->GetBinContent(i) / binWidth);
      hist->SetBinError(i, hist->GetBinError(i) / binWidth);
    }
  }
}

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

// h2: TH2D* with X = runid (or run index), Y = pT
// result: hRatio = h2 / <pT-spectrum>, same binning as h2

TH2D *makeRunVsPtRatio(TH2D *h2) {
  const int nx = h2->GetNbinsX(); // runs
  const int ny = h2->GetNbinsY(); // pT

  // --- 1. Average pT shape over runs (equal weight per non-empty run) ---

  // use a temporary projection to clone binning
  TH1D *hTmp = h2->ProjectionY("hTmp");
  TH1D *hAvg = (TH1D *)hTmp->Clone("hAvgPtShape");
  hAvg->Reset();

  int nGoodRuns = 0;

  for (int ix = 1; ix <= nx; ++ix) {
    // pT spectrum for a single run (x-bin ix)
    TH1D *hRun = h2->ProjectionY(Form("hPt_run_%d", ix), ix, ix);

    double integral = hRun->Integral();
    if (integral <= 0) {
      delete hRun; // empty run, skip it from the average
      continue;
    }

    // accumulate shapes
    hAvg->Add(hRun);
    ++nGoodRuns;

    delete hRun;
  }

  if (nGoodRuns > 0)
    hAvg->Scale(1.0 / nGoodRuns); // average shape over runs

  // --- 2. Build 2D ratio: data / average shape ---

  TH2D *hRatio = (TH2D *)h2->Clone("hPtVsRun_over_AvgPtShape");
  hRatio->Reset();

  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      double val = h2->GetBinContent(ix, iy);
      double avg = hAvg->GetBinContent(iy); // same pT bin
      double val_err = h2->GetBinError(ix, iy);
      double avg_err = hAvg->GetBinError(iy);

      if (avg > 0.0) {
        hRatio->SetBinContent(ix, iy, val / avg);
        // error propagation
        double rel_err =
            std::sqrt(std::pow(val_err / val, 2) + std::pow(avg_err / avg, 2));
        hRatio->SetBinError(ix, iy, (val / avg) * rel_err);

      } else
        hRatio->SetBinContent(ix, iy, 0.0); // or some sentinel
    }
  }

  // hAvg is your average pT spectrum
  // hRatio is the normalized TH2 for spotting outliers
  return hRatio;
}

// Main cross-section calculation function
void process(const std::string &trigger = "JP2",
             const std::string &jetR = "0.6") {

  std::cout << "\n========================================" << std::endl;
  std::cout << "Calculating cross-section for:" << std::endl;
  std::cout << "  Jet R: " << jetR << std::endl;
  std::cout << "========================================\n" << std::endl;

  AnalysisConfig cfg(jetR);
  gStyle->SetOptStat(0);

  std::string outputDir =
      Form("%s/results_R%s", cfg.workdir.c_str(), jetR.c_str());
  gSystem->mkdir(outputDir.c_str(), kTRUE);

  // find TFile from global list
  TFile *outFile = TFile::Open("compare_triggers.root", "UPDATE");

  if (!outFile || outFile->IsZombie()) {
    std::cerr << "Error: Cannot open output file (compare_triggers.root) for "
                 "writing"
              << std::endl;
    return;
  }

  // =======================================================================
  // Step 1: Calculate event weight factor
  // =======================================================================

  std::cout << "Loading run mapping..." << std::endl;
  std::map<int, std::string> inverseRunMap;
  std::map<std::string, int> runMap =
      loadRunMap(Form("%s/run_map.txt", cfg.workdir.c_str()), inverseRunMap);
  if (runMap.empty()) {
    std::cerr << "Error: Run map is empty, aborting processing for " << trigger
              << std::endl;
    return;
  }

  std::string dataFileName =
      Form("%s/merged_data_%s_R%s.root", cfg.datapath.c_str(), trigger.c_str(),
           jetR.c_str());

  // Get number of run bins from analyzed events histogram
  TFile *dataFile = TFile::Open(dataFileName.c_str());
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "Error: Cannot open data file: " << dataFileName << std::endl;
    return;
  }
  TH1D *analyzedEvents = (TH1D *)dataFile->Get("hEventsRun");
  if (!analyzedEvents) {
    std::cerr << "Error: Cannot find hEventsRun" << std::endl;
    dataFile->Close();
    return;
  }
  set<int> currentRuns = getCurrentRuns(analyzedEvents, runMap);
  if (currentRuns.empty()) {
    std::cerr << "Error: No valid runs found in hEventsRun for " << trigger
              << std::endl;
    dataFile->Close();
    return;
  }

  // =======================================================================
  // Step 1.2: Load luminosity data
  // =======================================================================
  std::cout << "\nStep 1.2: Loading luminosity data..." << std::endl;
  TString lumiFileName =
      "/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/lumi.root";
  TFile *lumiFile = TFile::Open(lumiFileName);
  if (!lumiFile || lumiFile->IsZombie()) {
    std::cerr << "Error: Cannot open lumi file: " << lumiFileName << std::endl;
    return;
  }

  TH1D *maxnevent = (TH1D *)lumiFile->Get(Form("nevents_%s", trigger.c_str()));
  TH1D *luminosity =
      (TH1D *)lumiFile->Get(Form("luminosity_%s", trigger.c_str()));
  if (!maxnevent || !luminosity) {
    std::cerr << "Error: Cannot load luminosity histograms for " << trigger
              << std::endl;
    lumiFile->Close();
    return;
  }

  outFile->cd();
  TDirectory *triggerDir =
      outFile->mkdir(Form("%s_R%s", trigger.c_str(), jetR.c_str()));
  triggerDir->cd();

  TH1D *eventWeight = (TH1D *)maxnevent->Clone("AnalyzedEvents");
  eventWeight->SetTitle(
      Form("Analyzed events %s;runid;counts", trigger.c_str()));
  eventWeight->SetDirectory(0);
  eventWeight->Write();

  filterBadRuns(eventWeight, cfg.badRunsJP2);

  eventWeight->Divide(analyzedEvents);
  eventWeight->SetName("AnalyzedEvents_over_MaxNEvent");
  eventWeight->SetTitle(
      Form("Analyzed events / MaxNEvent %s", trigger.c_str()));
  eventWeight->Write();

  // eventWeight->Divide(luminosity);
  // eventWeight->SetName("EventCrossSection");
  // eventWeight->SetTitle(
  //     Form("Event CrossSection %s;runid;pb", trigger.c_str()));

  // add more bad runs based on eventWeight histogram
  // for (int i = 1; i <= eventWeight->GetNbinsX(); i++) {
  //   double content = eventWeight->GetBinContent(i);
  //   if (content <= 50)
  //     continue;
  //   eventWeight->SetBinContent(i, 0);
  //   eventWeight->SetBinError(i, 0);

  //   std::string run = eventWeight->GetXaxis()->GetBinLabel(i);
  //   int runNumber = std::stoi(run);
  //   if (std::find(cfg.badRunsJP2.begin(), cfg.badRunsJP2.end(), runNumber) ==
  //       cfg.badRunsJP2.end()) {
  //     cfg.badRunsJP2.push_back(runNumber);
  //   }
  // }

  // cout all bad runs
  std::cout << "Bad runs for " << trigger << ": ";
  for (const auto &run : cfg.badRunsJP2) {
    std::cout << run << " ";
  }
  std::cout << std::endl;
  eventWeight->Write();

  double totalWeight = eventWeight->Integral();
  int nCurrentRuns = currentRuns.size();
  double normalizedWeight = totalWeight / nCurrentRuns;

  // create raw jet spectrum histogram
  ROOT::RDataFrame df("ResultTree", dataFileName.c_str());

  // apply runMap to "runid1" branch
  auto dfFiltered =
      df.Define("runIndex",
                [&runMap](int runid1) {
                  auto it = runMap.find(std::to_string(runid1));
                  return it == runMap.end() ? -1 : it->second;
                },
                {"runid1"})
          .Filter(
              [&cfg](int runIndex) {
                // Skip runs that are unknown or marked as bad
                return runIndex >= 0 &&
                       std::find(cfg.badRunsJP2.begin(), cfg.badRunsJP2.end(),
                                 runIndex) == cfg.badRunsJP2.end();
              },
              {"runIndex"});

  std::vector<double> ptBins(pt_reco_bins.begin(), pt_reco_bins.end());
  TH1D *spectrum = nullptr;
  TH2D *spectrumVsRun = nullptr;
  if (trigger == "JP2" || trigger == "HT2") {
    // For triggered jets - use trigger matching
    std::string triggerBranch =
        "trigger_match_" + std::string(trigger == "JP2" ? "JP2" : "HT2");

    auto dfTriggered = dfFiltered.Define("triggered_mask", triggerBranch)
                           .Define("triggered_pt", "pt[triggered_mask]")
                           .Filter("triggered_pt.size() > 0");

    auto hist = dfTriggered.Histo1D(
        {Form("%s_jet_spectrum_R%s", trigger.c_str(), jetR.c_str()),
         Form("%s Jet Spectrum R=%s;p_{T} [GeV/c];counts", trigger.c_str(),
              jetR.c_str()),
         (int)ptBins.size() - 1, ptBins.data()},
        "triggered_pt");
    auto histJetPtvsRun = dfTriggered.Histo2D(
        {Form("%s_jet_spectrum_vs_run_R%s", trigger.c_str(), jetR.c_str()),
         Form("%s Jet Spectrum vs Run R=%s;runid;p_{T} [GeV/c]",
              trigger.c_str(), jetR.c_str()),
         analyzedEvents->GetNbinsX(), analyzedEvents->GetXaxis()->GetXmin(),
         analyzedEvents->GetXaxis()->GetXmax(), (int)ptBins.size() - 1,
         ptBins.data()},
        "runIndex", "triggered_pt");

    if (!hist) {
      std::cerr << "Error: Cannot create histogram for " << trigger
                << std::endl;
      dataFile->Close();
      return;
    }
    spectrum = (TH1D *)hist.GetPtr()->Clone("RawSpectrum");
    spectrum->SetDirectory(0);

    spectrumVsRun = histJetPtvsRun.GetPtr();
    spectrumVsRun = multiply2DBy1D(spectrumVsRun, eventWeight);
    spectrumVsRun->SetName("RawSpectrumVsRun");
    spectrumVsRun->SetDirectory(0);
  } else {
  }
  // } else {
  //   // For MB - no trigger matching needed
  //   auto histPtr =
  //       df.Histo1D({Form("%s_jet_spectrum_R%s", trigger.c_str(),
  //       jetR.c_str()),
  //                   Form("%s Jet Spectrum R=%s;p_{T} [GeV/c];counts",
  //                        trigger.c_str(), jetR.c_str()),
  //                   (int)ptBins.size() - 1, ptBins.data()},
  //                  "pt");
  //   spectrum = (TH1D *)histPtr->Clone("RawSpectrum");
  //   spectrum->SetDirectory(0);
  // }

  // spectrum->Write("RawSpectrumBeforeCorrections");

  // // Trigger efficiency (reco level)
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
  /// event weighted spectrum
  triggerDir->cd();
  spectrum->Write();
  spectrumVsRun->Write();

  TH2D *spectrumVsRunRatio = makeRunVsPtRatio(spectrumVsRun);
  spectrumVsRunRatio->SetName("RawSpectrumVsRun_Ratio");
  spectrumVsRunRatio->GetZaxis()->SetRangeUser(0.4, 2.0);
  spectrumVsRunRatio->Write();
  ////// Efficiency correction
  if (!triggerEfficiency) {
    std::cerr << "Warning: Trigger efficiency histogram missing for " << trigger
              << ", skipping efficiency correction" << std::endl;
  } else {
    divideByEfficiency(spectrum, triggerEfficiency);
    divideByEfficiency(spectrumVsRun, triggerEfficiency);
  }

  spectrum->SetName("CorrectedSpectrum");
  spectrum->Write();
  spectrumVsRun->SetName("CorrectedSpectrumVsRun");
  spectrumVsRun->Write();
  TH2D *correctedVsRunRatio = makeRunVsPtRatio(spectrumVsRun);
  correctedVsRunRatio->SetName("CorrectedSpectrumVsRun_Ratio");
  correctedVsRunRatio->GetZaxis()->SetRangeUser(0.4, 2.0);
  correctedVsRunRatio->Write();

  ////// Normalization

  spectrum->Scale(normalizedWeight);
  spectrum->SetName("NormalizedCorrectedSpectrum");
  spectrum->Write();

  spectrumVsRun->Scale(normalizedWeight);
  spectrumVsRun->SetName("NormalizedCorrectedSpectrumVsRun");
  spectrumVsRun->Write();

  // apply fake rate subtraction

  outFile->Save();
  outFile->Close();
}

void performUnfolding(const std::string &trigger = "JP2",
                      const std::string &jetR = "0.6") {
  TFile *inFile = TFile::Open("compare_triggers.root", "READ");
  if (!inFile || inFile->IsZombie()) {
    std::cerr << "Error: Cannot open compare_triggers.root for unfolding"
              << std::endl;
    return;
  }

  // Get Spectrum histogram
  TH1D *spectrum = (TH1D *)inFile->Get(Form(
      "%s_R%s/NormalizedCorrectedSpectrum", trigger.c_str(), jetR.c_str()));
  if (!spectrum) {
    std::cerr << "Error: Missing normalized spectrum for " << trigger
              << " R=" << jetR
              << ". Run process() first to fill compare_triggers.root."
              << std::endl;
    inFile->Close();
    return;
  }

  // fakerate efficiency (reco level)
  TH1D *fakeRateEfficiency = nullptr;
  TString fakeRateFileName =
      Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/matching/"
           "plots/embedding_root_%s_R%s/"
           "fakereco_over_all.root",
           trigger.c_str(), jetR.c_str());
  TFile *fakeRateFile = TFile::Open(fakeRateFileName);
  if (fakeRateFile && !fakeRateFile->IsZombie()) {
    fakeRateEfficiency = (TH1D *)fakeRateFile->Get("h_ratio");
    fakeRateEfficiency->SetDirectory(0);
    std::cout << "  Loaded fake rate efficiency (reco)" << std::endl;
    fakeRateFile->Close();
  }

  TFile *outFile = TFile::Open("unfolded_spectra.root", "UPDATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "Error: Cannot open unfolded_spectra.root for writing"
              << std::endl;
    inFile->Close();
    return;
  }
  outFile->cd();

  // Apply fake rate correction
  if (!fakeRateEfficiency) {
    std::cerr << "Warning: Fake rate efficiency is missing for " << trigger
              << ", skipping fake rate correction" << std::endl;
  } else {
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
  }

  // Apply unfolding
  // Load response matrix
  TString responseFileName =
      Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/unfolding/"
           "response_%s_R%s.root",
           trigger.c_str(), jetR.c_str());
  TFile *responseFile = TFile::Open(responseFileName);
  if (!responseFile || responseFile->IsZombie()) {
    std::cerr << "Error: Cannot open response file: " << responseFileName
              << std::endl;
    return;
  }

  AnalysisConfig cfg(jetR);
  RooUnfoldResponse *response =
      (RooUnfoldResponse *)responseFile->Get("my_response");
  if (!response) {
    std::cerr << "Error: Response matrix my_response is missing in "
              << responseFileName << std::endl;
    inFile->Close();
    outFile->Close();
    return;
  }
  std::cout << "  Loaded response matrix" << std::endl;
  std::cout << "  Performing Bayesian unfolding with " << cfg.nIterations
            << " iterations..." << std::endl;
  // Perform unfolding
  RooUnfoldBayes unfold(response, spectrum, cfg.nIterations);
  TH1D *unfolded = (TH1D *)unfold.Hunfold();
  unfolded->SetName("spectrum_unfolded");

  //  Apply reconstruction efficiency correction

  // Reconstruction efficiency (MC level)
  TString recoEffFileName =
      Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/matching/"
           "plots/embedding_root_%s_R%s/"
           "matched_over_mc.root",
           trigger.c_str(), jetR.c_str());
  TFile *recoEffFile = TFile::Open(recoEffFileName);
  TH1D *matchingEfficiency = nullptr;
  if (recoEffFile && !recoEffFile->IsZombie()) {
    matchingEfficiency = (TH1D *)recoEffFile->Get("h_ratio");
    if (matchingEfficiency) {
      matchingEfficiency->SetDirectory(0);
      std::cout << "  Loaded reconstruction efficiency (MC)" << std::endl;
    }
  }

  divideByEfficiency(unfolded, matchingEfficiency);

  // divide by bin width pt
  divideByBinWidth(unfolded);
  // normalize by eta acceptance
  unfolded->Scale(1.0 / cfg.etaAcceptance);

  outFile->cd();
  TDirectory *triggerDir =
      outFile->mkdir(Form("%s_R%s", trigger.c_str(), jetR.c_str()));
  triggerDir->cd();
  unfolded->Write("UnfoldedSpectrum");

  outFile->Save();
  outFile->Close();
  inFile->Close();
}

// Plot HT2 vs JP2 (or any triggers listed in config) on the same canvas
// and append the pages to a single PDF.
void plotTriggerComparison(const std::string &jetR, const std::string &pdfOut,
                           TCanvas *c) {
  AnalysisConfig cfg(jetR);
  std::unique_ptr<TFile> inFile(TFile::Open("compare_triggers.root", "READ"));
  if (!inFile || inFile->IsZombie()) {
    std::cerr << "Error: Cannot open compare_triggers.root for plotting"
              << std::endl;
    return;
  }

  struct PlotSpec {
    std::string histName;
    std::string yTitle;
    bool logY;
  };

  std::vector<PlotSpec> plots = {
      {"AnalyzedEvents", "counts", false},
      {"AnalyzedEvents_over_MaxNEvent", "counts", false},
      {"EventCrossSection", "pb", true},
      {"RawSpectrum", "counts", true},
      {"CorrectedSpectrum", "counts", true},
      {"NormalizedCorrectedSpectrum", "pb", true}};

  // Reusable legend and text boxes per plot
  TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.85);
  leg->SetBorderSize(0);
  // Latex removed from plots per request

  for (const auto &plot : plots) {
    c->Clear();
    leg->Clear();
    c->SetLogy(plot.logY);

    bool firstDraw = true;
    size_t colorIdx = 0;
    double maxY = 0.0;
    std::vector<std::unique_ptr<TH1D>> ownedHists;

    for (const auto &trig : cfg.triggers) {
      std::string dirName = Form("%s_R%s", trig.c_str(), jetR.c_str());
      TDirectory *dir = (TDirectory *)inFile->Get(dirName.c_str());
      if (!dir)
        continue;
      dir->cd();
      TH1D *h = (TH1D *)dir->Get(plot.histName.c_str()); // owned by file
      if (!h)
        continue;

      auto hClone = std::unique_ptr<TH1D>(
          (TH1D *)h->Clone(Form("%s_%s_clone_R%s", trig.c_str(),
                                plot.histName.c_str(), jetR.c_str())));
      hClone->SetDirectory(0);
      hClone->SetMarkerStyle(markers[colorIdx % markers.size()]);
      hClone->SetMarkerColor(colors[colorIdx % colors.size()]);
      hClone->SetLineColor(colors[colorIdx % colors.size()]);
      hClone->SetTitle(Form("%s R=%s", plot.histName.c_str(), jetR.c_str()));
      maxY = std::max(maxY, hClone->GetMaximum());

      if (firstDraw) {
        hClone->Draw("PE");
        firstDraw = false;
      } else {
        hClone->Draw("PE SAME");
      }
      leg->AddEntry(hClone.get(), trig.c_str(), "lp");
      ownedHists.push_back(std::move(hClone));
      colorIdx++;
    }

    if (!ownedHists.empty()) {
      ownedHists.front()->GetYaxis()->SetTitle(plot.yTitle.c_str());
      if (!plot.logY) {
        ownedHists.front()->SetMaximum(maxY * 1.2);
      } else {
        ownedHists.front()->SetMaximum(maxY * 10);
      }
      leg->Draw();              // keep overlay legend
      c->Print(pdfOut.c_str()); // overlay page
      std::cout << "Plotted overlay for " << plot.histName << " R=" << jetR
                << std::endl;
    } else {
      // Skip plotting if histogram is missing
      std::cout << "Skipping plot (missing hist): " << plot.histName
                << " R=" << jetR << std::endl;
      continue;
    }

    // Ratio page (separate)
    if (ownedHists.size() > 1) {
      c->Clear();
      c->SetLogy(false);
      c->SetGridy();
      leg->Clear(); // remove legend on ratio page
      double ratioMax = 0.0;
      double ratioMin = 1e9;
      bool hasRatioContent = false;
      TH1D *refHist = ownedHists.front().get();
      bool firstRatio = true;
      const double eps = 1e-8;
      std::vector<std::unique_ptr<TH1D>> ratioHists;

      for (size_t i = 1; i < ownedHists.size(); ++i) {
        auto ratio = std::unique_ptr<TH1D>((TH1D *)ownedHists[i]->Clone(
            Form("ratio_%zu_%s_R%s", i, plot.histName.c_str(), jetR.c_str())));
        for (int b = 1; b <= ratio->GetNbinsX(); ++b) {
          double refVal = refHist->GetBinContent(b);
          double val = ownedHists[i]->GetBinContent(b);
          double refErr = refHist->GetBinError(b);
          double valErr = ownedHists[i]->GetBinError(b);
          if (refVal > eps) {
            hasRatioContent = true;
            double r = val / refVal;
            ratio->SetBinContent(b, r);
            double relErr2 =
                (val > 0 ? std::pow(valErr / std::max(val, eps), 2) : 0.0) +
                std::pow(refErr / refVal, 2);
            ratio->SetBinError(b, r * std::sqrt(relErr2));
            ratioMax = std::max(ratioMax, r + ratio->GetBinError(b));
            ratioMin = std::min(ratioMin, r - ratio->GetBinError(b));
          } else {
            ratio->SetBinContent(b, 0);
            ratio->SetBinError(b, 0);
          }
        }
        ratio->SetMarkerStyle(ownedHists[i]->GetMarkerStyle());
        ratio->SetMarkerColor(ownedHists[i]->GetMarkerColor());
        ratio->SetLineColor(ownedHists[i]->GetLineColor());
        ratio->GetYaxis()->SetTitle("HT2/JP2");
        ratio->GetYaxis()->SetTitleSize(0.05);
        ratio->GetYaxis()->SetTitleOffset(0.9);
        ratio->GetYaxis()->SetLabelSize(0.045);
        ratio->GetXaxis()->SetLabelSize(0.045);
        ratio->GetXaxis()->SetTitleSize(0.05);
        ratio->GetXaxis()->SetTitleOffset(1.0);
        ratio->GetXaxis()->SetTitle(ownedHists[i]->GetXaxis()->GetTitle());

        double yMin = hasRatioContent ? std::max(0.1, ratioMin - 0.2) : 0.5;
        double yMax = 1.5;

        ratio->SetMinimum(yMin);
        ratio->SetMaximum(yMax);

        if (firstRatio) {
          ratio->Draw("PE");
          firstRatio = false;
        } else {
          ratio->Draw("PE SAME");
        }
        ratioHists.push_back(std::move(ratio));
      }

      if (firstRatio) {
        // Nothing drawable; show an empty frame so the page is not blank
        auto frame =
            std::unique_ptr<TH1D>((TH1D *)refHist->Clone("ratio_frame"));
        frame->Reset("ICES");
        frame->SetMinimum(0.5);
        frame->SetMaximum(1.5);
        frame->GetYaxis()->SetTitle("HT2/JP2");
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitleOffset(0.9);
        frame->GetYaxis()->SetLabelSize(0.045);
        frame->GetXaxis()->SetLabelSize(0.045);
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetTitleOffset(1.0);
        frame->GetXaxis()->SetTitle(refHist->GetXaxis()->GetTitle());
        frame->Draw("AXIS");
        firstRatio = false;
      }

      // Draw unity line
      double xmin = ownedHists.front()->GetXaxis()->GetXmin();
      double xmax = ownedHists.front()->GetXaxis()->GetXmax();
      TLine *line = new TLine(xmin, 1.0, xmax, 1.0);
      line->SetLineStyle(2);
      line->Draw();
      // No legend or text on ratio page
      c->Print(pdfOut.c_str()); // ratio page
      std::cout << "Plotted ratio for " << plot.histName << " R=" << jetR
                << std::endl;
    }
  }
}

void plotEfficiencyComparison(const std::string &jetR,
                              const std::string &pdfOut, TCanvas *c) {

  // c->Clear();
  c->SetLogy(false);
  c->SetGridy(true);

  TH1D *triggerEfficiencyHT2 = nullptr;
  TH1D *triggerEfficiencyJP2 = nullptr;

  TString trigEffFileNameHT2 =
      Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/matching/"
           "plots/embedding_root_HT2_R%s/"
           "HT2_over_all.root",
           jetR.c_str());
  TFile *trigEffFileHT2 = TFile::Open(trigEffFileNameHT2);

  TH1D *h = (TH1D *)trigEffFileHT2->Get("h_ratio");
  triggerEfficiencyHT2 = (TH1D *)h->Clone("TriggerEfficiencyHT2");

  TString trigEffFileNameJP2 =
      Form("/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/matching/"
           "plots/embedding_root_JP2_R%s/"
           "JP2_over_all.root",
           jetR.c_str());
  TFile *trigEffFileJP2 = TFile::Open(trigEffFileNameJP2);

  h = (TH1D *)trigEffFileJP2->Get("h_ratio");
  triggerEfficiencyJP2 = (TH1D *)h->Clone("TriggerEfficiencyJP2");

  if (!triggerEfficiencyHT2 || !triggerEfficiencyJP2) {
    std::cerr << "Error: Missing trigger efficiency histograms for plotting"
              << std::endl;
    return;
  }

  triggerEfficiencyJP2->SetMarkerStyle(20);
  triggerEfficiencyJP2->SetMarkerColor(kBlue);
  triggerEfficiencyJP2->SetLineColor(kBlue);
  triggerEfficiencyJP2->SetTitle(Form("Trigger Efficiency R=%s", jetR.c_str()));
  triggerEfficiencyJP2->GetYaxis()->SetTitle("Efficiency");
  triggerEfficiencyJP2->SetMaximum(1.2);
  triggerEfficiencyJP2->Draw("PE");

  triggerEfficiencyHT2->SetMarkerStyle(21);
  triggerEfficiencyHT2->SetMarkerColor(kRed);
  triggerEfficiencyHT2->SetLineColor(kRed);
  triggerEfficiencyHT2->Draw("PE SAME");

  TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.85);
  leg->SetBorderSize(0);
  leg->AddEntry("JP2", "JP2 Trigger", "lp");
  leg->AddEntry("HT2", "HT2 Trigger", "lp");
  leg->Draw();
  c->Print(pdfOut.c_str());

  // draw ratio page
  c->Clear();
  c->SetLogy(false);
  c->SetGridy();
  TH1D *ratio = (TH1D *)triggerEfficiencyHT2->Clone(
      Form("TriggerEfficiencyRatio_R%s", jetR.c_str()));
  ratio->Divide(triggerEfficiencyJP2);
  ratio->SetTitle(Form("Trigger Efficiency Ratio HT2/JP2 R=%s", jetR.c_str()));
  ratio->GetYaxis()->SetTitle("HT2 / JP2");
  ratio->GetYaxis()->SetTitleSize(0.05);
  ratio->GetYaxis()->SetTitleOffset(0.9);
  ratio->GetYaxis()->SetLabelSize(0.045);
  ratio->GetXaxis()->SetLabelSize(0.045);
  ratio->GetXaxis()->SetTitleSize(0.05);
  ratio->GetXaxis()->SetTitleOffset(1.0);
  ratio->GetXaxis()->SetTitle(triggerEfficiencyHT2->GetXaxis()->GetTitle());
  ratio->GetYaxis()->SetRangeUser(0, 1.5);
  ratio->Draw("PE");

  // Draw unity line
  double xmin = triggerEfficiencyHT2->GetXaxis()->GetXmin();
  double xmax = triggerEfficiencyHT2->GetXaxis()->GetXmax();
  TLine *line = new TLine(xmin, 1.0, xmax, 1.0);
  line->SetLineStyle(2);
  line->Draw();

  c->Print(pdfOut.c_str()); // ratio page
  std::cout << "Plotted trigger efficiency comparison R=" << jetR << std::endl;
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
  std::vector<std::string> triggers = {"JP2", "HT2"};
  std::vector<int> colors = {kAzure - 1, kRed + 1, kOrange - 1};
  std::vector<int> markers = {20, 21, 22};
  std::map<std::string, TH1D *> myCrossSections;
  std::map<std::string, TH1D *> ratios;

  for (size_t i = 0; i < triggers.size(); i++) {
    std::string trigger = triggers[i];
    TString fileName = Form("%s/unfolded_spectra.root", cfg.workdir.c_str(),
                            trigger.c_str(), jetR.c_str());
    TFile *file = TFile::Open(fileName);
    if (file && !file->IsZombie()) {
      TH1D *cs = (TH1D *)file->Get(
          Form("%s_R%s/UnfoldedSpectrum", trigger.c_str(), jetR.c_str()));
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
void compareTriggers() {
  AnalysisConfig cfg;
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  std::string pdfOut = "check.pdf";
  can->Print((pdfOut + "[").c_str()); // start PDF

  std::vector<std::string> triggers = {"HT2", "JP2"};
  std::vector<std::string> jetRs = {"0.5"};

  for (const auto &jetR : jetRs) {
    for (const auto &trigger : triggers) {
      process(trigger, jetR);
      performUnfolding(trigger, jetR);
    }
    compareCrossSections(jetR);
  }

  // Append comparison plots for each radius to the same PDF

  for (const auto &jetR : jetRs) {
    plotTriggerComparison(jetR, pdfOut, can);
    plotEfficiencyComparison(jetR, pdfOut, can);
  }

  can->Print((pdfOut + "]").c_str()); // end PDF
}
