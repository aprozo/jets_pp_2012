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

// // get pearson coeffs from covariance matrix
// TH2D *getPearsonCoeffs(const TMatrixD &covMatrix)
// {
//     Int_t nrows = covMatrix.GetNrows();
//     Int_t ncols = covMatrix.GetNcols();

//     TH2D *PearsonCoeffs = new TH2D((TString) "PearsonCoeffs" + covMatrix.GetName(), "Pearson Coefficients;", nrows, 0, nrows, ncols, 0, ncols);
//     for (Int_t row = 0; row < nrows; row++)
//     {
//         for (Int_t col = 0; col < ncols; col++)
//         {
//             Double_t pearson = 0.;
//             if (covMatrix(row, row) != 0. && covMatrix(col, col) != 0.)
//                 pearson = covMatrix(row, col) / TMath::Sqrt(covMatrix(row, row) * covMatrix(col, col));
//             PearsonCoeffs->SetBinContent(row + 1, col + 1, pearson);
//         }
//     }

//     PearsonCoeffs->GetZaxis()->SetRangeUser(-1, 1);
//     return PearsonCoeffs;
// }

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
    legend->AddEntry(hTruth, "True", "l");
    legend->AddEntry(hMeasured, "Measured", "l");

    can->cd(1);
    gPad->SetLogy();

    hTruth->SetMarkerStyle(20);
    hTruth->SetLineColor(kViolet);
    hTruth->GetYaxis()->SetTitle("dN/dp_{t}");
    hTruth->GetYaxis()->SetTitleOffset(1.);

    hTruth->Draw("hist");
    hMeasured->Draw("hist same");

    hMeasured->SetMarkerStyle(21);
    hMeasured->SetLineColor(kTeal - 1);

    // TH2D *fPearsonCoeffs[nIter];
    for (int iter = 0; iter < nIter; iter++)
    {
        can->cd(1);
        RooUnfoldBayes unfolding(response, hMeasured, plotIterations[iter]);
        TH1D *hUnfolded = (TH1D *)unfolding.Hunfold();
        hUnfolded->SetName(Form("hUnfolded%i", iter));

        // fPearsonCoeffs[iter] = getPearsonCoeffs(unfolding.Eunfold(RooUnfold::kCovariance));
        hUnfolded->SetLineColor(iter + 2000);
        hUnfolded->SetMarkerStyle(20);
        hUnfolded->SetMarkerColor(iter + 2000);
        legend->AddEntry(hUnfolded, Form("Iter%i", plotIterations[iter]), "pel");
        hUnfolded->Draw("PE same");

        can->cd(2);
        TH1D *xRatio = (TH1D *)hUnfolded->Clone(Form("xRatio%i", iter));
        xRatio->GetYaxis()->SetTitle("dN/dp_{t}");
        xRatio->SetLineColor(iter + 2000);
        xRatio->SetMarkerStyle(20);
        xRatio->SetMarkerColor(iter + 2000);
        xRatio->Divide(hTruth);
        xRatio->GetYaxis()->SetTitle("Unfolded/True");
        xRatio->Draw(iter == 0 ? "PE" : "PE same");
    }
    can->cd();
    legend->Draw();
    can->SaveAs(outPdf);
    can->Clear();
    can->SetLogz(1);

    // RooUnfoldBayes lastUnfolding(response, hMeasured, 4);
    // TH1D *hUnfolded = (TH1D *)lastUnfolding.Hunfold();

    // TH2D *hUnfoldingMatrix = new TH2D(lastUnfolding.UnfoldingMatrix());
    // hUnfoldingMatrix->SetTitle("Bin Migration probability");
    // hUnfoldingMatrix->SetName("UnfoldingMatrix");
    // hUnfoldingMatrix->GetXaxis()->SetTitle("bin=N_{reco}+pt");
    // hUnfoldingMatrix->GetYaxis()->SetTitle("bin=N_{unfold}+pt");
    // can->cd();
    // hUnfoldingMatrix->GetZaxis()->SetTitleOffset(0.7);

    // hUnfoldingMatrix->Draw("colz");
    // gPad->Update();
    // TPaletteAxis *palette = (TPaletteAxis *)hUnfoldingMatrix->GetListOfFunctions()->FindObject("palette");
    // palette->SetX1NDC(0.90);
    // palette->SetX2NDC(0.92);

    // tex->DrawLatex(0.5, 0.35, "Unfolding Matrix");
    // can->SaveAs(outPdf);
    // delete hUnfoldingMatrix;

    // can->Clear();
    // can->cd();
    // can->SetLogz(0);
    // for (int iter = 0; iter < nIter; iter++)
    // {
    //     fPearsonCoeffs[iter]->SetTitle((TString) "Pearson Coefficients;bin =N_{bins}+pt;bin =N_{bins}+pt");
    //     fPearsonCoeffs[iter]->Draw("colz");
    //     gPad->Update();
    //     tex->DrawLatex(0.2, 0.35, Form("Pearson Coefficients Iter %i", plotIterations[iter]));

    //     TPaletteAxis *palette = (TPaletteAxis *)fPearsonCoeffs[iter]->GetListOfFunctions()->FindObject("palette");
    //     palette->SetX1NDC(0.91);
    //     palette->SetX2NDC(0.93);
    //     can->SaveAs(outPdf);
    //     delete fPearsonCoeffs[iter];
    // }
    // can->Clear();
}

void unfold()
{
    gStyle->SetOptStat(0);

    vector<double> ptRecoBinsVec = {5, 7, 10, 13, 50};
    vector<double> ptMcBinsVec = {5, 7, 10, 13, 50};
    new TColor(2000, (255. / 255.), (89. / 255.), (74. / 255.));   // red-ish
    new TColor(2001, (25. / 255.), (170. / 255.), (25. / 255.));   // green-ish
    new TColor(2002, (66. / 255.), (98. / 255.), (255. / 255.));   // blue-ish
    new TColor(2003, (153. / 255.), (0. / 255.), (153. / 255.));   // magenta-ish
    new TColor(2004, (255. / 255.), (166. / 255.), (33. / 255.));  // yellow-ish
    new TColor(2005, (0. / 255.), (170. / 255.), (255. / 255.));   // azur-ish
    new TColor(2006, (204. / 255.), (153. / 255.), (255. / 255.)); // violet-ish
    new TColor(2007, (107. / 255.), (142. / 255.), (35. / 255.));  // olive
    new TColor(2008, (100. / 255.), (149. / 255.), (237. / 255.)); // corn flower blue
    new TColor(2009, (255. / 255.), (69. / 255.), (0. / 255.));    // orange red
    new TColor(2010, (0. / 255.), (128. / 255.), (128. / 255.));   // teal
    new TColor(2011, (176. / 255.), (196. / 255.), (222. / 255.)); // light steel blue
    new TColor(2012, (255. / 255.), (215. / 255.), (0. / 255.));   // gold-ish
    gSystem->Load("~/install/RooUnfold/build/libRooUnfold");

    ///////////////////////////////////////////////////////////////////////////
    // Create/read response RooUnfold
    TFile *responseFile = new TFile("response.root", "RECREATE");
    responseFile->cd();

    TH1D *hMeasured;
    TH1D *hTruth;
    TH1D *hMeasuredTest;
    TH1D *hTruthTest;
    // Create RooUnfoldResponse object
    RooUnfoldResponse *response;

    hMeasured = new TH1D("Meas", ";p_{t}, GeV/c;", 10, 3, 50);
    hTruth = new TH1D("True", ";p_{t}, GeV/c;", 10, 3, 50);
    hMeasuredTest = new TH1D("MeasTest", ";p_{t}, GeV/c;", 10, 3, 50);
    hTruthTest = new TH1D("TrueTest", ";p_{t}, GeV/c;", 10, 3, 50);
    //  RooUnfoldResponse (const TH1* measured, const TH1* truth, const char* name= 0, const char* title= 0);
    response = new RooUnfoldResponse(10, 3, 50, 10, 3, 50, "response", "response");

    TH2D *hResponseMatrix = new TH2D("hResponseMatrix", "Response Matrix; Measured; Truth", 10, 3, 50, 10, 3, 50);

    TFile *treeFile; // Open the file containing the tree.
    treeFile = new TFile("matched_jets.root", "READ");
    if (!treeFile || treeFile->IsZombie())
    {
        return;
    }

    double ptMc;
    double ptReco;
    double weight;

    TTree *matchedTree = (TTree *)treeFile->Get("MatchTree");
    matchedTree->SetBranchAddress("ptMc", &ptMc);
    matchedTree->SetBranchAddress("ptReco", &ptReco);
    matchedTree->SetBranchAddress("weight", &weight);

    Double_t nEntries = matchedTree->GetEntries();

    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        Float_t progress = 0.;
        progress = (Float_t)iEntry / (nEntries);
        if (iEntry % 10000 == 0)
        {
            cout << "Training: \r (" << (progress * 100.0) << "%)" << std::flush;
        }
        matchedTree->GetEntry(iEntry);
        response->Fill(ptReco, ptMc, weight);
        hResponseMatrix->Fill(ptReco, ptMc);

        // cout << "ptMc = " << ptMc << " ptReco = " << ptReco << " weight = " << weight << endl;

        hMeasured->Fill(ptReco, weight);
        hTruth->Fill(ptMc, weight);

    } // end of loop over train entries

    // TTree *missTree = (TTree *)treeFile->Get("MissTree");
    // double ptMcMiss;
    // double weightMiss;
    // missTree->SetBranchAddress("pt", &ptMcMiss);
    // missTree->SetBranchAddress("weight", &weightMiss);

    // for (Int_t iEntry = 0; iEntry < missTree->GetEntries(); iEntry++)
    // {
    //     missTree->GetEntry(iEntry);
    //     response->Miss(ptMcMiss, weightMiss);
    // }

    // TTree *fakeTree = (TTree *)treeFile->Get("FakeTree");
    // double ptRecoFake;
    // double weightFake;
    // fakeTree->SetBranchAddress("pt", &ptRecoFake);
    // fakeTree->SetBranchAddress("weight", &weightFake);
    // for (Int_t iEntry = 0; iEntry < fakeTree->GetEntries(); iEntry++)
    // {
    //     fakeTree->GetEntry(iEntry);
    //     response->Fake(ptRecoFake, weightFake);
    // }
    // ========================================================================================================

    // Take a random sample of the data to test the unfolding
    int nTestEntries = matchedTree->GetEntries() * 0.1;
    TRandom3 rand(0);
    int counter = 0;

    // while (counter < nTestEntries)
    // {
    //     int iEntry = rand.Uniform(0, matchedTree->GetEntries());
    //     matchedTree->GetEntry(iEntry);
    //     hTruthTest->Fill(ptMc, weight);
    //     hMeasuredTest->Fill(ptReco, weight);

    //     float progress = 0.;
    //     progress = (Float_t)counter / (nTestEntries);
    //     if (counter % 10000 == 0)
    //     {
    //         cout << "Testing: \r (" << (progress * 100.0) << "%)" << std::flush;
    //     }

    //     // missTree->GetEntry(rand.Uniform(0, missTree->GetEntries()));
    //     // hTruthTest->Fill(ptMcMiss, weightMiss);

    //     // fakeTree->GetEntry(rand.Uniform(0, fakeTree->GetEntries()));
    //     // hMeasuredTest->Fill(ptRecoFake, weightFake);
    //     counter++;
    // } // end of loop over test entries

    // for (Int_t iEntry = 0; iEntry < missTree->GetEntries(); iEntry++)
    // {
    //     missTree->GetEntry(iEntry);
    //     hTruthTest->Fill(ptMcMiss, weightMiss);
    // }

    // for (Int_t iEntry = 0; iEntry < fakeTree->GetEntries() * (1 - TrainToTestRatio); iEntry++)
    // {
    //     fakeTree->GetEntry(iEntry);
    //     hMeasuredTest->Fill(ptRecoFake, weightFake);
    // }
    // ========================================================================================================
    // Draw closure test check
    TCanvas *can = new TCanvas("can", "Closure Test Check", 1200, 600);
    TString outPdf = "closure_check.pdf";
    can->SaveAs(outPdf + "[");

    plotIterations(can, outPdf, response, hTruth, hMeasured);

    // can->cd();
    // gPad->SetLogy();
    // hMeasured->SetMarkerStyle(20);
    // hMeasured->SetLineColor(kViolet);
    // hMeasured->Draw("hist");
    // hMeasured->GetYaxis()->SetTitle("dN/dp_{t}");
    // hMeasured->GetYaxis()->SetRangeUser(100, 1e12);
    // hTruth->SetMarkerStyle(21);
    // hTruth->SetLineColor(kTeal - 1);
    // hTruth->Draw("hist same");

    // TH1D *hUnfolded;
    // RooUnfoldBayes unfolding(response, hMeasured, 4);
    // hUnfolded = (TH1D *)unfolding.Hunfold();
    // hUnfolded->SetLineColor(1);
    // hUnfolded->SetMarkerStyle(20);
    // hUnfolded->SetMarkerColor(1);
    // hUnfolded->Draw("PE same");
    // can->SaveAs(outPdf);

    // Save response and test histograms
    responseFile->cd();
    hResponseMatrix->Write();
    hMeasuredTest->Write();
    hTruthTest->Write();
    response->Write();

    can->SaveAs(outPdf + "]");

    responseFile->Save();
    // responseFile->Close();
}

//==============================================================================
