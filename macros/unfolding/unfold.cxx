#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TLatex.h"
#include "TColor.h"
#include "TRandom3.h"
#include "TEntryList.h"
#include "TMath.h"

const vector<double> pt_reco_bins = {9.7, 11.5, 13.6, 16.1, 19.0, 22.5, 26.6, 31.4, 37.2, 44.0, 52.0, 70.0};
const vector<double> pt_mc_bins = {5, 9.7, 11.5, 13.6, 16.1, 19.0, 22.5, 26.6, 31.4, 37.2, 44.0, 52.0, 70.0};
// const vector<double> pt_reco_bins = {6.9, 8.2, 9.7, 11.5, 13.6, 16.1, 19.0, 22.5, 26.6, 31.4, 37.2, 44.0, 52.0, 70.0};
// const vector<double> pt_mc_bins = {6.9, 8.2, 9.7, 11.5, 13.6, 16.1, 19.0, 22.5, 26.6, 31.4, 37.2, 44.0, 52.0, 70.0};

// get pearson coeffs from covariance matrix
TH2D *getPearsonCoeffs(const TMatrixD &covMatrix)
{
    Int_t nrows = covMatrix.GetNrows();
    Int_t ncols = covMatrix.GetNcols();

    TH2D *PearsonCoeffs = new TH2D((TString) "PearsonCoeffs" + covMatrix.GetName(), "Pearson Coefficients;", nrows, 0, nrows, ncols, 0, ncols);
    for (Int_t row = 0; row < nrows; row++)
    {
        for (Int_t col = 0; col < ncols; col++)
        {
            Double_t pearson = 0.;
            if (covMatrix(row, row) != 0. && covMatrix(col, col) != 0.)
                pearson = covMatrix(row, col) / TMath::Sqrt(covMatrix(row, row) * covMatrix(col, col));
            PearsonCoeffs->SetBinContent(row + 1, col + 1, pearson);
        }
    }

    return PearsonCoeffs;
}

TH2D *remakeMatrix(TH2D *hist, TString xTitle, TString yTitle)
{
    TH2D *newHist;
    // xaxis for reco
    vector<double> new_x_axis;
    if (xTitle == "reco")
        new_x_axis = pt_reco_bins;
    else if (xTitle == "unfolded")
        new_x_axis = pt_mc_bins;

    vector<double> new_y_axis;
    if (yTitle == "mc")
        new_y_axis = pt_reco_bins;
    else if (yTitle == "unfolded")
        new_y_axis = pt_reco_bins;

    newHist = new TH2D(hist->GetName() + (TString) "remake", hist->GetTitle(), new_x_axis.size() - 1, new_x_axis.data(), new_y_axis.size() - 1, new_y_axis.data());
    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        for (int j = 1; j <= hist->GetNbinsY(); j++)
        {
            newHist->SetBinContent(i, j, hist->GetBinContent(i, j));
        }
    }
    newHist->GetXaxis()->SetTitle(xTitle + " p_{T}, GeV/c");
    newHist->GetYaxis()->SetTitle(yTitle + " p_{T}, GeV/c");
    newHist->GetZaxis()->SetTitle(hist->GetZaxis()->GetTitle());
    return newHist;
}

void plotIterations(TCanvas *can, TString outPdf, RooUnfoldResponse *response, TH1D *hTruth, TH1D *hMeasured)
{
    const vector<Int_t> plotIterations = {1, 2, 3, 4, 5};
    const Int_t nIter = plotIterations.size();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);
    can->Clear();
    can->Divide(2, 1);

    TLegend *legend = new TLegend(0.30, 0.62, 0.44, 0.88);
    legend->AddEntry(hTruth, "Truth", "l");
    legend->AddEntry(hMeasured, "Measured", "l");

    can->cd(1);
    gPad->SetLogy();
    hTruth->SetMarkerStyle(20);
    hTruth->SetLineColor(kViolet);
    hTruth->Draw("hist");
    hMeasured->Draw("hist same");

    hMeasured->SetMarkerStyle(21);
    hMeasured->SetLineColor(kTeal - 1);

    TLine *line = new TLine();
    line->SetLineStyle(2);
    line->SetLineColor(kGray + 2);

    for (int iter = 0; iter < nIter; iter++)
    {
        cout << "Unfolding Iteration: " << plotIterations[iter] << endl;
        can->cd(1);
        RooUnfoldBayes unfolding(response, hMeasured, plotIterations[iter]);
        TH1D *hUnfolded = (TH1D *)unfolding.Hunfold();
        hUnfolded->SetName(Form("hUnfolded%i", iter));
        hUnfolded->SetLineColor(iter + 2000);
        hUnfolded->SetMarkerStyle(20);
        hUnfolded->SetMarkerColor(iter + 2000);
        legend->AddEntry(hUnfolded, Form("Iter%i", plotIterations[iter]), "pel");
        hUnfolded->Draw("PE same");

        can->cd(2);
        TH1D *xRatio = (TH1D *)hUnfolded->Clone(Form("xRatio%i", iter));

        xRatio->SetLineColor(iter + 2000);
        xRatio->SetMarkerStyle(20);
        xRatio->SetMarkerColor(iter + 2000);
        xRatio->Divide(hTruth);
        xRatio->GetYaxis()->SetTitle("Unfolded/Truth");
        xRatio->GetYaxis()->SetTitleOffset(1.1);
        xRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
        // Draw +-5% band with lines

        xRatio->Draw(iter == 0 ? "PE" : "PE same");
        line->DrawLine(5, 1.05, 80, 1.05);
        line->DrawLine(5, 0.95, 80, 0.95);
    }
    can->cd();
    legend->Draw();
    can->SaveAs(outPdf);
    can->Clear();
    can->SetLogz(1);

    RooUnfoldBayes lastUnfolding(response, hMeasured, 4);
    TH1D *hUnfolded = (TH1D *)lastUnfolding.Hunfold();

    TH2D *hUnfoldingMatrix = new TH2D(lastUnfolding.UnfoldingMatrix());
    hUnfoldingMatrix->SetTitle("Bin Migration probability");
    hUnfoldingMatrix->SetName("UnfoldingMatrix");
    can->cd();
    hUnfoldingMatrix = remakeMatrix(hUnfoldingMatrix, "reco", "mc");
    hUnfoldingMatrix->GetZaxis()->SetRangeUser(0.00001, 1);
    hUnfoldingMatrix->Draw("colz");

    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis *)hUnfoldingMatrix->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.90);
    palette->SetX2NDC(0.92);

    tex->DrawLatex(0.5, 0.35, "Unfolding Matrix");
    can->SaveAs(outPdf);
    delete hUnfoldingMatrix;

    can->Clear();
    can->cd();
    can->SetLogz(0);

    TH2D *fPearsonCoeffs[nIter];

    for (int iter = 0; iter < nIter; iter++)
    {
        fPearsonCoeffs[iter] = getPearsonCoeffs(lastUnfolding.Eunfold(RooUnfold::kCovariance));
        fPearsonCoeffs[iter]->SetTitle((TString) "Pearson Coefficients;bin pt_{unfolded};bin pt_{unfolded}");
        fPearsonCoeffs[iter]->SetName(Form("PearsonCoeffsIter%i", plotIterations[iter]));
        fPearsonCoeffs[iter] = remakeMatrix(fPearsonCoeffs[iter], "unfolded", "unfolded");
        fPearsonCoeffs[iter]->GetZaxis()->SetRangeUser(-1, 1);
        fPearsonCoeffs[iter]->Draw("colz");
        gPad->Update();
        tex->DrawLatex(0.2, 0.35, Form("Pearson Coefficients Iter %i", plotIterations[iter]));

        TPaletteAxis *palette = (TPaletteAxis *)fPearsonCoeffs[iter]->GetListOfFunctions()->FindObject("palette");
        palette->SetX1NDC(0.91);
        palette->SetX2NDC(0.93);

        can->SaveAs(outPdf);
        delete fPearsonCoeffs[iter];
    }
    can->Clear();
}

void unfold()
{
    gStyle->SetOptStat(0);
    gSystem->Load("~/install/RooUnfold/build/libRooUnfold");
    ///////////////////////////////////////////////////////////////////////////
    const float testFraction = 0.2;

    TH1D *hMeasured;
    TH1D *hTruth;
    TH1D *hMeasuredTest;
    TH1D *hTruthTest;
    TH2D *hResponseMatrix;
    TH2D *detectorResolution;
    // Create RooUnfoldResponse object
    RooUnfoldResponse *response;

    bool isTraining = gSystem->AccessPathName("response.root");

    TFile *responseFile;

    // check if file "response.root" exists
    if (isTraining)
    {
        responseFile = new TFile("response.root", "RECREATE");
        // response = new RooUnfoldResponse(nBinsMeasured, measured_pt_min, measured_pt_max, nBinsTruth, truth_pt_min, truth_pt_max);
        // hResponseMatrix = new TH2D("hResponseMatrix", "Response Matrix; Measured; Truth", nBinsMeasured, measured_pt_min, measured_pt_max, nBinsTruth, truth_pt_min, truth_pt_max);
        // hMeasured = new TH1D("Measured", ";p_{t}, GeV/c;", nBinsMeasured, measured_pt_min, measured_pt_max);
        // hTruth = new TH1D("Truth", ";p_{t}, GeV/c;", nBinsTruth, truth_pt_min, truth_pt_max);
        hResponseMatrix = new TH2D("hResponseMatrix", "Response Matrix; Measured; Truth", pt_reco_bins.size() - 1, pt_reco_bins.data(), pt_mc_bins.size() - 1, pt_mc_bins.data());

        hMeasured = new TH1D("Measured", ";p_{t}, GeV/c; dN/dp_{t}", pt_reco_bins.size() - 1, pt_reco_bins.data());
        hMeasuredTest = (TH1D *)hMeasured->Clone("MeasuredTest");

        hTruth = new TH1D("Truth", ";p_{t}, GeV/c;dN/dp_{t}", pt_mc_bins.size() - 1, pt_mc_bins.data());
        hTruthTest = (TH1D *)hTruth->Clone("TruthTest");

        detectorResolution = new TH2D("detectorResolution", "Detector Resolution; p_{T}^{mc}, GeV/c; p_{T}^{reco} - p_{T}^{mc}, GeV/c", 1000, 0, 100, 1000, -100, 100);
        // response = new RooUnfoldResponse("my_response", "my_response");

        // response->Setup(hMeasured, hTruth);
        // response = new RooUnfoldResponse(nBinsMeasured, measured_pt_min, measured_pt_max, nBinsTruth, truth_pt_min, truth_pt_max);

        // TFile *treeFile = new TFile("~/dev/star/jets_pp_2012/output/matched_jets.root", "READ");
        TFile *treeFile = new TFile("~/matched_jets.root", "READ");
        if (!treeFile || treeFile->IsZombie())
        {
            cout << "Error: cannot open matched_jets.root" << endl;
            return;
        }

        double pt_mc;
        double pt_reco;
        double weight;

        TTree *matchedTree = (TTree *)treeFile->Get("MatchedTree");
        matchedTree->SetBranchAddress("mc_pt", &pt_mc);
        matchedTree->SetBranchAddress("reco_pt", &pt_reco);
        matchedTree->SetBranchAddress("mc_weight", &weight);

        Double_t nEntries = matchedTree->GetEntries();
        TRandom3 rand(0);

        for (Int_t iEntry = 0; iEntry < nEntries; iEntry++)
        {
            Float_t progress = 0.;
            progress = (Float_t)iEntry / (nEntries);
            if (iEntry % 10000 == 0)
            {
                cout << "Training: \r (" << (progress * 100.0) << "%)" << std::flush;
            }
            matchedTree->GetEntry(iEntry);
            if (pt_mc == -9 || pt_reco == -9)
                continue;

            double deltaPt = pt_reco - pt_mc;
            detectorResolution->Fill(pt_mc, deltaPt, weight);

            // if (pt_mc < pt_mc_bins[0] || pt_reco < pt_reco_bins[0])
            //     continue;

            if (rand.Uniform() > testFraction)
            {
                // response->Fill(pt_reco, pt_mc, weight);
                hResponseMatrix->Fill(pt_reco, pt_mc, weight);
                hMeasured->Fill(pt_reco, weight);
                hTruth->Fill(pt_mc, weight);
            }
            else
            {
                hTruthTest->Fill(pt_mc, weight);
                hMeasuredTest->Fill(pt_reco, weight);
            }

        } // end of loop over train entries
        response = new RooUnfoldResponse(hMeasured, hTruth, hResponseMatrix);
        response->UseOverflow();
        response->SetName("my_response");

        detectorResolution->FitSlicesY(0, 0, -1, 10, "QNR G2");

        TH1D *hMean = (TH1D *)gDirectory->Get("detectorResolution_1");
        TH1D *hSigma = (TH1D *)gDirectory->Get("detectorResolution_2");

        TH1D *stdDev = new TH1D("stdDev", "stdDev", hMean->GetNbinsX(), hMean->GetXaxis()->GetXmin(), hMean->GetXaxis()->GetXmax());
        
        for (int i = 1; i <= detectorResolution->GetNbinsX(); i++)
        {
            TH1D *slice = detectorResolution->ProjectionY("slice", i, i);
            stdDev->SetBinContent(i, slice->GetStdDev());
        }

        responseFile->cd();
        detectorResolution->Write();
        hMean->Write();
        hSigma->Write();
        stdDev->Write();

        hMeasured->Write();
        hTruth->Write();
        hMeasuredTest->Write();
        hTruthTest->Write();
        hResponseMatrix->Write();
        response->Write();

        responseFile->Save();
    }
    else
    {
        responseFile = new TFile("response.root", "READ");
        response = (RooUnfoldResponse *)responseFile->Get("my_response");
        hMeasured = (TH1D *)responseFile->Get("Measured");
        hTruth = (TH1D *)responseFile->Get("Truth");
        hMeasuredTest = (TH1D *)responseFile->Get("MeasuredTest");
        hTruthTest = (TH1D *)responseFile->Get("TruthTest");
        hResponseMatrix = (TH2D *)responseFile->Get("hResponseMatrix");
    }

    // ========================================================================================================
    // Draw closure test check
    TCanvas *can = new TCanvas("can", "Closure Test Check", 1400, 600);
    TString outPdf = "closure_check.pdf";
    can->SaveAs(outPdf + "[");
    // plotIterations(can, outPdf, response, hTruth, hMeasured);
    plotIterations(can, outPdf, response, hTruthTest, hMeasuredTest);
    can->SaveAs(outPdf + "]");

    responseFile->Close();
    delete can;
    // Save response and test histograms
}

//==============================================================================
