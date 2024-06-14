//! HAS to be compiled,
//! root -l macros/PrepUnfolding.cxx+
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>

#include <iostream>
#include <vector>

using namespace std;

int plotHists()
{
    TString data_type[4] = {"MB", "HT2", "mc", "geant"};

    //{"hJetPt"};
    TString hist_names[1] = {"hJetPt_Fine"};

    TFile *dataFile[4];
    TString inputFiles[4];
    TCanvas *can = new TCanvas("c1", "c1", 800, 600);
    TLatex *tex = new TLatex();

    can->SetLogy();
    can->SaveAs("plots.pdf[");

    for (int i = 0; i < 4; i++)
    {
        inputFiles[i] = "output/" + data_type[i] + "_hists.root";
        dataFile[i] = new TFile(inputFiles[i], "READ");
        for (int j = 0; j < 1; j++)
        {
            TH1D *hist = (TH1D *)dataFile[i]->Get(hist_names[j]);
            hist->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            hist->GetYaxis()->SetTitle("dN/dp_{T} [1/GeV/c]");
            hist->Rebin(100);
            hist->GetXaxis()->SetRangeUser(0, 50);
            hist->Draw();
            tex->DrawLatex(0.5, 0.5, data_type[i]);
            can->SaveAs("output/plots/all_pt.pdf");
        }
    }
    can->SaveAs("plots.pdf]");
    return 0;
}
