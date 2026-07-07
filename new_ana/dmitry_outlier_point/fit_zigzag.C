// fit_zigzag.C ---------------------------------------------------------------
// The [13.6,16.1) (14.85 GeV) point of Dmitry's Table III is an upward outlier.
// An inclusive-jet spectrum is smooth, so fit a curved power law (log-log cubic)
// to the 10.6-48 GeV points EXCLUDING 14.85, then quote its residual in units of
// its own statistical error.
// -----------------------------------------------------------------------------
#include <TBox.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TStyle.h>

void fit_zigzag()
{
   gStyle->SetOptStat(0);

   TFile *file = TFile::Open("/gpfs01/star/pwg/prozorov/study_pp2012/fable5/jets_pp_2012/"
                             "new_ana/jet_cross_section_dmitriyR0.5.root");
   TH1D *spectrum        = (TH1D *)file->Get("crossSection_statistic"); // the cross section

   const int numberOfBins = spectrum->GetNbinsX();
   const int testBin      = 5;     // [13.6,16.1) = 14.85 GeV, the point under test
   const double fitLowPt  = 5.0;  // fit only the physics range (center >= 10.6)
   const double fitHighPt = 52.0;

   // --- fit ln(y) = pol3(ln pT), all physics points except 14.85 ---
   TGraphErrors logLogGraph;
   int numberOfFitPoints = 0;
   for (int bin = 1; bin <= numberOfBins; ++bin) {
      double jetPt = spectrum->GetBinCenter(bin), crossSection = spectrum->GetBinContent(bin);
      if (crossSection <= 0 || jetPt < fitLowPt || jetPt > fitHighPt || bin == testBin) continue;
      logLogGraph.SetPoint(numberOfFitPoints, TMath::Log(jetPt), TMath::Log(crossSection));
      logLogGraph.SetPointError(numberOfFitPoints, 0, spectrum->GetBinError(bin) / crossSection);
      ++numberOfFitPoints;
   }
   TF1 smoothFit("smoothFit", "pol3", TMath::Log(fitLowPt), TMath::Log(fitHighPt));
   logLogGraph.Fit(&smoothFit, "QN");
   auto smoothModel = [&](double jetPt) {
      double logPt = TMath::Log(jetPt);
      return TMath::Exp(smoothFit.GetParameter(0) + smoothFit.GetParameter(1) * logPt +
                        smoothFit.GetParameter(2) * logPt * logPt +
                        smoothFit.GetParameter(3) * logPt * logPt * logPt);
   };
   const double chiSquarePerNdf = smoothFit.GetChisquare() / smoothFit.GetNDF();

   // --- residuals ---
   printf("\nsmooth log-log cubic fit (14.85 excluded): chi2/ndf = %.2f\n", chiSquarePerNdf);
   printf("  pT      y(data)    fit        (d-f)/f     sigma\n");
   double testResidualPercent = 0, testResidualSigma = 0, coreResidualRms = 0;
   int numberOfCorePoints = 0;
   for (int bin = 1; bin <= numberOfBins; ++bin) {
      double jetPt = spectrum->GetBinCenter(bin), crossSection = spectrum->GetBinContent(bin);
      if (crossSection <= 0 || jetPt < fitLowPt) continue;
      double fitValue      = smoothModel(jetPt);
      double residualPercent = 100 * (crossSection - fitValue) / fitValue;
      double residualSigma   = (crossSection - fitValue) / spectrum->GetBinError(bin);
      printf("  %5.1f  %9.3g  %9.3g  %+7.2f%%  %+6.1f  %s\n", jetPt, crossSection, fitValue,
             residualPercent, residualSigma, bin == testBin ? "<-- TEST" : "");
      if (bin == testBin) {
         testResidualPercent = residualPercent;
         testResidualSigma   = residualSigma;
      }
   }
   // --- figure: spectrum + fit (top), residuals (bottom) ---
   TCanvas *canvas = new TCanvas("canvas", "zigzag", 900, 1000);
   TPad *topPad    = new TPad("topPad", "", 0, 0.36, 1, 1);
   TPad *bottomPad = new TPad("bottomPad", "", 0, 0, 1, 0.36);
   topPad->SetBottomMargin(0.02); topPad->SetLeftMargin(0.13); topPad->SetLogy();
   bottomPad->SetTopMargin(0.02); bottomPad->SetBottomMargin(0.30); bottomPad->SetLeftMargin(0.13);
   topPad->Draw(); bottomPad->Draw();

   TGraphErrors fittedPoints, testPoint, turnOnPoints;
   int numberFitted = 0, numberTest = 0, numberTurnOn = 0;
   for (int bin = 1; bin <= numberOfBins; ++bin) {
      double jetPt = spectrum->GetBinCenter(bin), crossSection = spectrum->GetBinContent(bin);
      double binError = spectrum->GetBinError(bin);
      if (crossSection <= 0) continue;
      if (bin == testBin) {
         testPoint.SetPoint(numberTest, jetPt, crossSection);
         testPoint.SetPointError(numberTest, 0, binError); ++numberTest;
      } else {
         fittedPoints.SetPoint(numberFitted, jetPt, crossSection);
         fittedPoints.SetPointError(numberFitted, 0, binError); ++numberFitted;
      }
   }

   topPad->cd();
   fittedPoints.SetTitle(""); fittedPoints.SetMarkerStyle(20); fittedPoints.SetMarkerSize(1.1);
   fittedPoints.GetXaxis()->SetLimits(fitLowPt, 52); fittedPoints.GetXaxis()->SetLabelSize(0);
   fittedPoints.GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}d#eta  [pb/(GeV/c)]");
   fittedPoints.GetYaxis()->SetTitleOffset(1.25);
   fittedPoints.SetMinimum(2); fittedPoints.SetMaximum(3e7);
   fittedPoints.Draw("AP");
   TF1 fitCurve("fitCurve", "exp([0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3))", fitLowPt, fitHighPt);
   fitCurve.SetParameters(smoothFit.GetParameter(0), smoothFit.GetParameter(1),
                          smoothFit.GetParameter(2), smoothFit.GetParameter(3));
   fitCurve.SetLineColor(kBlue + 1); fitCurve.SetLineWidth(2); fitCurve.Draw("SAME");
   testPoint.SetMarkerStyle(20); testPoint.SetMarkerColor(kRed); testPoint.SetLineColor(kRed);
   testPoint.SetMarkerSize(1.5); testPoint.Draw("P SAME");

  TLegend legend(0.2, 0.2, 0.4, 0.4); legend.SetBorderSize(0); legend.SetTextSize(0.038);
  legend.AddEntry(&fittedPoints, "Dmitry paper", "p");
  legend.AddEntry(&testPoint, "Outlier", "p");
  legend.AddEntry(&fitCurve, "Fit", "l");
  legend.Draw();

   TLatex topLabel; topLabel.SetNDC(); topLabel.SetTextSize(0.038);

   topLabel.DrawLatex(0.3, 0.8, Form ("d#sigma = exp(a_{0}+a_{1}ln p_{T}+a_{2}(ln p_{T})^{2}+a_{3}(ln p_{T})^{3})"));
   topLabel.DrawLatex(0.5, 0.7, Form("#chi^{2}/ndf = %.2f", chiSquarePerNdf));
   topLabel.DrawLatex(0.4, 0.6, Form("14.85 point exlcluded"));
  //  $\ln y = p_0 + p_1 L + p_2 L^2 + p_3 L^3$ ($L=\ln p_T$)

   bottomPad->cd();
   TGraphErrors residualFitted, residualTest;
   int numberResidualFitted = 0, numberResidualTest = 0;
   for (int bin = 1; bin <= numberOfBins; ++bin) {
      double jetPt = spectrum->GetBinCenter(bin), crossSection = spectrum->GetBinContent(bin);
      if (crossSection <= 0 || jetPt < fitLowPt) continue;
      double fitValue        = smoothModel(jetPt);
      double residualPercent = 100 * (crossSection - fitValue) / fitValue;
      double errorPercent    = 100 * spectrum->GetBinError(bin) / fitValue;
      if (bin == testBin) {
         residualTest.SetPoint(numberResidualTest, jetPt, residualPercent);
         residualTest.SetPointError(numberResidualTest, 0, errorPercent); ++numberResidualTest;
      } else {
         residualFitted.SetPoint(numberResidualFitted, jetPt, residualPercent);
         residualFitted.SetPointError(numberResidualFitted, 0, errorPercent); ++numberResidualFitted;
      }
   }
   residualFitted.SetTitle(""); residualFitted.SetMarkerStyle(20); residualFitted.SetMarkerSize(1.1);
   residualFitted.GetXaxis()->SetLimits(fitLowPt, 52);
   residualFitted.GetXaxis()->SetTitle("jet p_{T} - UE  [GeV/c]");
   residualFitted.GetYaxis()->SetTitle("(data-fit)/fit [%]");
   residualFitted.GetXaxis()->SetTitleSize(0.11); residualFitted.GetXaxis()->SetLabelSize(0.09);
   residualFitted.GetYaxis()->SetTitleSize(0.095); residualFitted.GetYaxis()->SetTitleOffset(0.52);
   residualFitted.GetYaxis()->SetLabelSize(0.085); residualFitted.GetYaxis()->SetNdivisions(506);
   residualFitted.SetMinimum(-15); residualFitted.SetMaximum(15);
   residualFitted.Draw("AP");

   TLine *zeroLine = new TLine(fitLowPt, 0, 52, 0); zeroLine->SetLineStyle(2); zeroLine->SetLineColor(kGray + 2);
   zeroLine->Draw();
   residualFitted.Draw("P SAME");
   residualTest.SetMarkerStyle(20); residualTest.SetMarkerColor(kRed); residualTest.SetLineColor(kRed);
   residualTest.SetMarkerSize(1.6); residualTest.Draw("P SAME");
   TLatex bottomLabel; bottomLabel.SetNDC(); bottomLabel.SetTextSize(0.085); bottomLabel.SetTextColor(kRed);
   bottomLabel.DrawLatex(0.28, 0.85, Form("14.85 GeV: %+.1f%%  (%+.1f#sigma)", testResidualPercent, testResidualSigma));

   canvas->SaveAs("/gpfs01/star/pwg/prozorov/study_pp2012/scratch/zigzag/zigzag_fit.pdf");
}
