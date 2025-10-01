#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>

using namespace std;

void createCrossSectionHistogram(void) {
  // Enable scientific notation on y-axis
  gStyle->SetOptStat(0);

  // Create binning based on the table data
  const int nBins = 12;
  double binEdges[nBins + 1] = {6.9,  8.2,  9.7,  11.5, 13.6, 16.1, 19.0,
                                22.5, 26.6, 31.4, 37.2, 44.0, 52.0};

  // Create histogram with variable bin width
  TH1D *hist = new TH1D("crossSection_statistic",
                        "Differential Cross Section;p_{T} "
                        "[GeV];#frac{#sigma}{#Deltap_{T}#Delta#eta} [pb]",
                        nBins, binEdges);

  //   // Cross section values R 0.6
  //   double crossSection[nBins] =
  //   {6.34e6, 2.99e6, 8.47e5, 3.27e5, 1.19e5, 3.68e4,
  //                                 1.26e4, 3.36e3, 1.03e3, 2.30e2, 4.22e1, 6.30e0};

  //   // Statistical uncertainties
  //   double statErrors[nBins] = {0.02e6, 0.01e6, 0.05e5, 0.02e5, 0.01e5,
  //   0.03e4,
  //                               0.01e4, 0.05e3, 0.02e3, 0.10e2,
  //                               0.37e1, 1.03e0};

  //   // Sum of systematic uncertainties (quadrature sum of syst, simu, JES,
  //   UE) double systErrors[nBins]; double simuErrors[nBins] = {0.30e6, 0.21e6,
  //   0.90e5, 0.30e5, 0.07e5, 0.19e4,
  //                               0.05e4, 0.16e3, 0.03e3, 0.10e2, 0.41e1,
  //                               0.43e0};
  //   double jesErrors[nBins] = {0.76e6, 0.37e6, 0.67e5, 0.39e5, 0.07e5,
  //   0.46e4,
  //                              0.12e4, 0.26e3, 0.21e3, 0.22e2,
  //                              0.90e1, 1.19e0};
  //   double ueErrors[nBins] = {0.46e6, 0.15e6, 0.43e5, 0.16e5, 0.05e5, 0.11e4,
  //                             0.05e4, 0.09e3, 0.04e3, 0.06e2, 0.17e1,
  //                             0.18e0};

  //   // Calculate total systematic error (quadrature sum)
  //   for (int i = 0; i < nBins; i++) {
  //     systErrors[i] =
  //         sqrt(simuErrors[i] * simuErrors[i] + jesErrors[i] * jesErrors[i] +
  //              ueErrors[i] * ueErrors[i]);
  //   }

  // Cross section values R 0.5
  double crossSection[nBins] = {5.46e6, 1.91e6, 7.20e5, 2.40e5, 9.76e4, 3.03e4,
                                1.02e4, 2.95e3, 8.32e2, 1.98e2, 3.81e1, 5.28e0};
  double statErrors[nBins] = {2.00e4, 1.00e4, 3.00e3, 1.00e3, 5.00e2, 2.00e2,
                              1.00e2, 4.00e1, 1.80e1, 8.00e0, 3.00e0, 8.60e-1};
  double systErrors[nBins] = {5.00e5, 1.70e5, 8.10e4, 2.10e4, 9.10e3, 3.20e3,
                              1.00e3, 1.70e2, 1.74e2, 1.10e1, 7.70e0, 1.05e0};

  // Fill histogram with cross section values
  for (int i = 0; i < nBins; i++) {
    hist->SetBinContent(i + 1, crossSection[i]);

    // Set statistical error
    hist->SetBinError(i + 1, statErrors[i]);
  }
  // Draw systematic error bands
  TH1D *systBand = (TH1D *)hist->Clone("crossSection_systematic");
  systBand->SetFillColorAlpha(kRed, 0.35);

  for (int i = 0; i < nBins; i++) {
    systBand->SetBinError(i + 1, systErrors[i]);
  }

  TFile *outputFile =
      new TFile("jet_cross_section_dmitriyR0.5.root", "RECREATE");
  hist->Write();
  systBand->Write();
  outputFile->Close();

  cout << "Histogram created successfully!" << endl;
}