
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"
#include "TTree.h"

#include "TStyle.h"

using namespace std;

//! ----------------------------------------------------
struct InputTreeEntry {
   int runid;
   int runid1;
   int eventid;
   double weight;
   double refmult;
   int njets;
   double vz;
   int mult;
   float event_sum_pt;
   bool is_rejected;
   double neutral_fraction[1000];
   bool trigger_match[1000];
   double pt[1000];
   int n_constituents[1000];
   int index[1000];
};
void checkExistence(TFile *file, TString name)
{
   if (!file->Get(name)) {
      cout << "Error: " << name << " not found!" << endl;
      exit(1);
   }
}

int ana_trees_simple(TString OutFile = "hists2.root")
{

   gStyle->SetOptStat(0);
   TH1::SetDefaultSumw2();
   float RCut = 0.6;
   float EtaCut = 1.0 - RCut;

   //================================================================================================

   TFile *eventsFile = new TFile("events.root", "READ");
   // TH1D *weights = (TH1D *)eventsFile->Get("hist_JP2_weight");
   TH1D *nEvents_MB = (TH1D *)eventsFile->Get("hEventsRun_MB");

   TH1D *nEvents_JP2 = (TH1D *)eventsFile->Get("hEventsRun_JP2");
   TH1D *hist_sampled_lumi_VPDMB = (TH1D *)eventsFile->Get("hist_sampled_lumi_VPDMB");
   TH1D *hist_prescale_VPDMB = (TH1D *)eventsFile->Get("hist_prescale_VPDMB");
   TH1D *hist_livetime_VPDMB = (TH1D *)eventsFile->Get("hist_livetime_VPDMB");
   TH1D *hist_max_nevents_VPDMB = (TH1D *)eventsFile->Get("hist_max_nevents_VPDMB");
   TH1D *hist_sampled_lumi_JP2 = (TH1D *)eventsFile->Get("hist_sampled_lumi_JP2");
   TH1D *hist_prescale_JP2 = (TH1D *)eventsFile->Get("hist_prescale_JP2");
   TH1D *hist_livetime_JP2 = (TH1D *)eventsFile->Get("hist_livetime_JP2");
   TH1D *hist_max_nevents_JP2 = (TH1D *)eventsFile->Get("hist_max_nevents_JP2");
   TH1D *hist_JP2_weight = (TH1D *)eventsFile->Get("hist_JP2_weight");

   //================================================================================================

   TFile *outputFile = new TFile(OutFile, "RECREATE");
   TH1D *JP2 = new TH1D("JP2", "; p_{T} (GeV/c); Entries (per event)", 100, 0, 100);
   TH1D *JP2_Livetime = (TH1D *)JP2->Clone("JP2_Livetime");
   TH1D *JP2_Prescale = (TH1D *)JP2->Clone("JP2_Prescale");
   TH1D *JP2_SampledLumi = (TH1D *)JP2->Clone("JP2_SampledLumi");
   TH1D *JP2_MaxNEvents = (TH1D *)JP2->Clone("JP2_MaxNEvents");
   TH1D *JP2_MatchTrigger = (TH1D *)JP2->Clone("JP2_MatchTrigger");

   int nRunBins = nEvents_JP2->GetNbinsX();
   TH2D *JP2_run = new TH2D("JP2_run", "; p_{T} (GeV/c); Entries (per event)", 100, 0, 100, nRunBins, 0, nRunBins);

   vector<TString> runBins;
   // set names of bins to run numbers
   for (unsigned int i = 1; i <= nRunBins; ++i) {
      TString runBin = nEvents_JP2->GetXaxis()->GetBinLabel(i);
      runBins.push_back(runBin);
      JP2_run->GetYaxis()->SetBinLabel(i, runBin);
   }
   TH2D *JP2_run_Livetime = (TH2D *)JP2_run->Clone("JP2_run_Livetime");
   TH2D *JP2_run_Prescale = (TH2D *)JP2_run->Clone("JP2_run_Prescale");
   TH2D *JP2_run_SampledLumi = (TH2D *)JP2_run->Clone("JP2_run_SampledLumi");
   TH2D *JP2_run_MaxNEvents = (TH2D *)JP2_run->Clone("JP2_run_MaxNEvents");
   TH2D *JP2_run_MatchTrigger = (TH2D *)JP2_run->Clone("JP2_run_MatchTrigger");

   //================================================================================================

   TFile *mb_hist_file = new TFile("/home/prozorov/dev/star/jets_pp_2012/output/jets_MB.root", "READ");

   TTree *mb_tree = (TTree *)mb_hist_file->Get("ResultTree");

   TH1D *hPtMb = new TH1D("hPtMb", "; p_{T} (GeV/c); Entries (per event)", 100, 0, 100);

   InputTreeEntry mb;
   mb_tree->SetBranchAddress("runid1", &mb.runid1);
   mb_tree->SetBranchAddress("pt", mb.pt);

   for (int i = 0; i < mb_tree->GetEntries(); ++i) {
      mb_tree->GetEntry(i);
      //   show progress
      int run = mb.runid1;
      TString runstring = Form("%i", run);
      int binRun = nEvents_MB->GetXaxis()->FindBin(Form("%i", run));
      double events_in_run = nEvents_MB->GetBinContent(binRun);

      double total_weight = events_in_run;

      for (int j = 0; j < mb.njets; ++j) {
         hPtMb->Fill(mb.pt[j], 1. / total_weight);
      }
   }

   TFile *jp2_jets = new TFile("/home/prozorov/dev/star/jets_pp_2012/output/jets_JP2.root", "READ");
   if (!jp2_jets || jp2_jets->IsZombie()) {
      cout << "Error: JP2 tree file not found!" << endl;
      return -1;
   }
   TTree *jp2_tree = (TTree *)jp2_jets->Get("ResultTree");

   InputTreeEntry jp2;
   jp2_tree->SetBranchAddress("runid1", &jp2.runid1);
   jp2_tree->SetBranchAddress("vz", &jp2.vz);
   jp2_tree->SetBranchAddress("neutral_fraction", jp2.neutral_fraction);
   jp2_tree->SetBranchAddress("trigger_match_jp", jp2.trigger_match);
   jp2_tree->SetBranchAddress("pt", jp2.pt);
   jp2_tree->SetBranchAddress("njets", &jp2.njets);

   for (int i = 0; i < jp2_tree->GetEntries(); ++i) {
      jp2_tree->GetEntry(i);
      //   show progress
      if (i % 100000 == 0) {
         cout << "Processing JP2 entry: " << i << " / " << jp2_tree->GetEntries() << endl;
      }
      if (jp2.njets < 1)
         continue; // skip events with no jets
      if (jp2.vz < -30 || jp2.vz > 30)
         continue; // vz cut
      int run = jp2.runid1;
      if (run == 13049007 || run == 13048092 || run == 13049006 || run == 13051099 || run == 13051074 ||
          run == 13064067 || run == 13070061 || run == 13069023 || run == 13048093 || run == 13068060 ||
          run == 13052061 || run == 13048019 || run == 13061035) // Livetime JP2 very poor!!!!
         continue;

      TString runstring = Form("%i", run);
      int binRun = nEvents_JP2->GetXaxis()->FindBin(Form("%i", run));

      double events_in_run = nEvents_JP2->GetBinContent(binRun);

      double weight_prescale = hist_prescale_VPDMB->GetBinContent(binRun) / hist_prescale_JP2->GetBinContent(binRun);

      double weight_livetime = hist_livetime_JP2->GetBinContent(binRun) / hist_livetime_VPDMB->GetBinContent(binRun);

      double weight_max_nevents =
         hist_max_nevents_VPDMB->GetBinContent(binRun) / hist_max_nevents_JP2->GetBinContent(binRun);

      double total_weight = events_in_run;
      // cout << "======================================" << endl;

      // fill the histogram with the event weight
      for (int j = 0; j < jp2.njets; ++j) {
         total_weight *= 1;
         JP2->Fill(jp2.pt[j], 1. / total_weight);
         JP2_run->Fill(jp2.pt[j], runstring, 1. / total_weight);

         if (!jp2.trigger_match[j])
            continue;
         JP2_MatchTrigger->Fill(jp2.pt[j], 1. / total_weight);
         JP2_run_MatchTrigger->Fill(jp2.pt[j], runstring, 1. / total_weight);

         total_weight *= weight_livetime;
         JP2_Livetime->Fill(jp2.pt[j], 1. / total_weight);
         JP2_run_Livetime->Fill(jp2.pt[j], runstring, 1. / total_weight);

         total_weight *= weight_max_nevents;
         JP2_MaxNEvents->Fill(jp2.pt[j], 1. / total_weight);
         JP2_run_MaxNEvents->Fill(jp2.pt[j], runstring, 1. / total_weight);

         total_weight *= weight_prescale;
         JP2_Prescale->Fill(jp2.pt[j], 1. / total_weight);
         JP2_run_Prescale->Fill(jp2.pt[j], runstring, 1. / total_weight);
      }
   }

   TCanvas *can = new TCanvas("can", "can", 800, 600);
   can->cd();
   TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
   leg->AddEntry(hPtMb, "MB", "l");
   leg->AddEntry(JP2, "JP2", "pel");
   leg->AddEntry(JP2_MatchTrigger, "+MatchTrigger", "l");
   leg->AddEntry(JP2_Livetime, "+Livetime", "l");
   leg->AddEntry(JP2_Prescale, "+Livetime+Prescale", "l");
   leg->AddEntry(JP2_MaxNEvents, "+Livetime+Prescale+MaxEvents", "l");

   JP2->SetLineColor(2009);
   JP2_Livetime->SetLineColor(2001);
   JP2_Prescale->SetLineColor(2002);
   JP2_MaxNEvents->SetLineColor(2003);
   hPtMb->SetLineColor(2004);
   JP2_MatchTrigger->SetLineColor(2006);

   gPad->SetLogy();
   JP2->GetYaxis()->SetRangeUser(1e-8, 1e3);
   JP2->Draw("pe ");

   hPtMb->Draw("same hist");
   leg->Draw("same hist");
   JP2_Livetime->Draw("same hist");
   JP2_MaxNEvents->Draw("same hist");
   JP2_Prescale->Draw("same hist");
   JP2_MatchTrigger->Draw("same hist");
   can->SaveAs("JP2.pdf");

   TCanvas *runRatio = new TCanvas("runRatio", "runRatio", 800, 600);
   runRatio->cd();
   gPad->SetLogy();
   TH1D *weight_prescale_run = (TH1D *)hist_prescale_VPDMB->Clone("weight_prescale_run");
   weight_prescale_run->Divide(hist_prescale_JP2);
   weight_prescale_run->SetLineColor(2002);
   weight_prescale_run->GetYaxis()->SetRangeUser(0.5, 1e3);
   weight_prescale_run->Draw();

   TH1D *weight_livetime_run = (TH1D *)hist_livetime_JP2->Clone("weight_livetime_run");
   weight_livetime_run->Divide(hist_livetime_VPDMB);
   weight_livetime_run->SetLineColor(2001);
   weight_livetime_run->Draw("same");

   TH1D *weight_max_nevents_run = (TH1D *)hist_max_nevents_VPDMB->Clone("weight_max_nevents_run");
   weight_max_nevents_run->Divide(hist_max_nevents_JP2);
   weight_max_nevents_run->SetLineColor(2003);
   weight_max_nevents_run->Draw("same");

   // Ratios to VPD MB
   TH1D *JP2_Livetime_Ratio = (TH1D *)JP2_Livetime->Clone("JP2_Livetime_Ratio");
   JP2_Livetime_Ratio->Divide(hPtMb);
   TH1D *JP2_Prescale_Ratio = (TH1D *)JP2_Prescale->Clone("JP2_Prescale_Ratio");
   JP2_Prescale_Ratio->Divide(hPtMb);
   TH1D *JP2_MaxNEvents_Ratio = (TH1D *)JP2_MaxNEvents->Clone("JP2_MaxNEvents_Ratio");
   JP2_MaxNEvents_Ratio->Divide(hPtMb);

   TH1D *JP2_TrigMatch_Ratio = (TH1D *)JP2_MatchTrigger->Clone("JP2_TrigMatch_Ratio");
   JP2_TrigMatch_Ratio->Divide(hPtMb);

   TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
   leg2->AddEntry(weight_prescale_run, "Prescale MB/JP2", "l");
   leg2->AddEntry(weight_livetime_run, "Livetime JP2/MB", "l");
   leg2->AddEntry(weight_max_nevents_run, "MaxNEvents MB/JP2", "l");
   leg2->Draw();

   // can->cd();
   // JP2_Livetime_Ratio->Draw();
   // JP2_Prescale_Ratio->Draw("same hist");
   // JP2_MaxNEvents_Ratio->Draw("same hist");
   // leg->Draw();
   // can->SaveAs("Ratio.pdf");

   // Write histograms to the output file
   outputFile->cd();
   JP2->Write();
   hPtMb->Write();
   JP2_Livetime->Write();
   JP2_Prescale->Write();
   JP2_MaxNEvents->Write();
   JP2_MatchTrigger->Write();
   JP2_run->Write();
   JP2_run_Livetime->Write();
   JP2_run_Prescale->Write();
   JP2_run_MaxNEvents->Write();
   JP2_run_MatchTrigger->Write();

   can->Write();

   JP2_Livetime_Ratio->Write();
   JP2_Prescale_Ratio->Write();
   JP2_MaxNEvents_Ratio->Write();
   JP2_TrigMatch_Ratio->Write();

   // outputFile->Close();

   return 0;
}