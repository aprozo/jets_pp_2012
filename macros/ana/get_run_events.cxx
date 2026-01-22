#include <TFile.h>
#include <TH1.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

// Define a struct to hold all the run information
struct RunInfo {
   double sampled_lumi_VPDMB;
   double prescale_VPDMB;
   double livetime_VPDMB;
   int max_nevents_VPDMB;
   double sampled_lumi_JP2;
   double prescale_JP2;
   double livetime_JP2;
   int max_nevents_JP2;
   double JP2_weight;
};

std::map<int, RunInfo> CreateRunInfoMap(const char *filename = "../lists/files_lumi_jaime/lumi_VPDMB_JP2.csv")
{
   // Open the file
   std::ifstream file(filename);
   if (!file.is_open()) {
      std::cerr << "Error: Could not open file " << filename << std::endl;
   }

   std::map<int, RunInfo> runInfoMap;

   // Read the header line and discard
   std::string header;
   std::getline(file, header);

   // Read data line by line
   std::string line;
   while (std::getline(file, line)) {
      try {
         std::stringstream ss(line);
         std::string token;

         // Extract run number
         std::getline(ss, token, ',');
         int run = std::stoi(token);

         // Create a new RunInfo object
         RunInfo info;

         // Parse remaining columns
         std::getline(ss, token, ',');
         info.sampled_lumi_VPDMB = std::stod(token);

         std::getline(ss, token, ',');
         info.prescale_VPDMB = std::stod(token);

         std::getline(ss, token, ',');
         info.livetime_VPDMB = std::stod(token);

         std::getline(ss, token, ',');
         info.max_nevents_VPDMB = std::stoi(token);

         std::getline(ss, token, ',');
         info.sampled_lumi_JP2 = std::stod(token);

         std::getline(ss, token, ',');
         info.prescale_JP2 = std::stod(token);

         std::getline(ss, token, ',');
         info.livetime_JP2 = std::stod(token);

         std::getline(ss, token, ',');
         info.max_nevents_JP2 = std::stoi(token);

         std::getline(ss, token, ',');
         info.JP2_weight = std::stod(
            token); // JP2_weight =
                    // (prescale_VPDMB/prescale_JP2)*(livetime_JP2/livetime_VPDMB)*(max_nevents_VPDMB/max_nevents_JP2)

         // Add to map
         runInfoMap[run] = info;
      } catch (const std::exception &e) {
         std::cerr << "Error parsing line: " << line << std::endl;
         std::cerr << "Exception: " << e.what() << std::endl;
         // Continue to next line
         continue;
      }
   }
   // Close the file
   file.close();

   return runInfoMap;
}

int get_run_events()
{

   // Create a map to store run number and associated information
   std::map<int, RunInfo> runInfoMap = CreateRunInfoMap("../lists/files_lumi_james_dunlop/lumi_VPDMB_JP2.csv");

   TString mb_file_name = "../output/MB/tree_jets.root";
   TFile *mb_file = new TFile(mb_file_name, "READ");

   TString jp_file_name = "../output/JP2/tree_jets.root";
   TFile *jp_file = new TFile(jp_file_name, "READ");

   TH1D *hEventsRun_MB = (TH1D *)mb_file->Get("hEventsRun");
   TH1D *hEventsRun_JP2 = (TH1D *)jp_file->Get("hEventsRun");

   hEventsRun_MB->SetName("hEventsRun_MB");
   hEventsRun_JP2->SetName("hEventsRun_JP2");

   TH1D *hRatio_MB_JP2 = (TH1D *)hEventsRun_MB->Clone("hRatio_MB_JP2");
   hRatio_MB_JP2->Divide(hEventsRun_JP2);
   hRatio_MB_JP2->SetTitle("Ratio of MB to JP2 events; Run Number; #frac{MB}{JP2}");

   TH1D *temp = (TH1D *)hRatio_MB_JP2->Clone("temp");
   temp->SetTitle(";Run Number;");
   // clear the bin content and errors
   for (int i = 1; i <= temp->GetNbinsX(); ++i) {
      temp->SetBinContent(i, 0);
      temp->SetBinError(i, 0);
   }
   TH1D *hist_sampled_lumi_VPDMB = (TH1D *)temp->Clone("hist_sampled_lumi_VPDMB");
   TH1D *hist_prescale_VPDMB = (TH1D *)temp->Clone("hist_prescale_VPDMB");
   TH1D *hist_livetime_VPDMB = (TH1D *)temp->Clone("hist_livetime_VPDMB");
   TH1D *hist_max_nevents_VPDMB = (TH1D *)temp->Clone("hist_max_nevents_VPDMB");
   TH1D *hist_sampled_lumi_JP2 = (TH1D *)temp->Clone("hist_sampled_lumi_JP2");
   TH1D *hist_prescale_JP2 = (TH1D *)temp->Clone("hist_prescale_JP2");
   TH1D *hist_livetime_JP2 = (TH1D *)temp->Clone("hist_livetime_JP2");
   TH1D *hist_max_nevents_JP2 = (TH1D *)temp->Clone("hist_max_nevents_JP2");
   TH1D *hist_JP2_weight = (TH1D *)temp->Clone("hist_JP2_weight");

   // Fill the histograms with the corresponding values from the map
   for (const auto &entry : runInfoMap) {
      int run = entry.first;
      const RunInfo &info = entry.second;
      int bin = hist_sampled_lumi_VPDMB->GetXaxis()->FindBin(Form("%i", run));

      hist_sampled_lumi_VPDMB->SetBinContent(bin, info.sampled_lumi_VPDMB);
      hist_prescale_VPDMB->SetBinContent(bin, info.prescale_VPDMB);
      hist_livetime_VPDMB->SetBinContent(bin, info.livetime_VPDMB);
      hist_max_nevents_VPDMB->SetBinContent(bin, info.max_nevents_VPDMB);
      hist_sampled_lumi_JP2->SetBinContent(bin, info.sampled_lumi_JP2);
      hist_prescale_JP2->SetBinContent(bin, info.prescale_JP2);
      hist_livetime_JP2->SetBinContent(bin, info.livetime_JP2);
      hist_max_nevents_JP2->SetBinContent(bin, info.max_nevents_JP2);
      hist_JP2_weight->SetBinContent(bin, info.JP2_weight);

      // JP2_weight = (prescale_VPDMB/prescale_JP2)*(livetime_JP2/livetime_VPDMB)*(max_nevents_VPDMB/max_nevents_JP2)
   }

   hist_sampled_lumi_VPDMB->SetTitle("Sampled Luminosity VPDMB; Run Number; sampled_lumi_VPDMB (nb^{-1})");
   hist_prescale_VPDMB->SetTitle("Prescale VPDMB; Run Number; prescale_VPDMB");
   hist_livetime_VPDMB->SetTitle("Livetime VPDMB; Run Number; livetime_VPDMB");
   hist_max_nevents_VPDMB->SetTitle("Max Events VPDMB; Run Number; max_nevents_VPDMB");
   hist_sampled_lumi_JP2->SetTitle("Sampled Luminosity JP2; Run Number; sampled_lumi_JP2 (nb^{-1})");
   hist_prescale_JP2->SetTitle("Prescale JP2; Run Number; prescale_JP2");
   hist_livetime_JP2->SetTitle("Livetime JP2; Run Number; livetime_JP2");
   hist_max_nevents_JP2->SetTitle("Max Events JP2; Run Number; max_nevents_JP2");

   hist_JP2_weight->SetTitle(
      "#frac{PS_{MB}}/{PS_{JP2}} #frac{livetime_{JP2}}/{livetime_{MB}} #frac{nevents_{MB}}/{nevents_{JP2}}; Run "
      "Number; (prescale_MB/prescale_JP2)*(livetime_JP2/livetime_MB)*(max_nevents_MB/max_nevents_JP2)");

   //   self-consistency check on max nevents
   TH1D *ratio_max_nevents_VPDMB = (TH1D *)hEventsRun_MB->Clone("ratio_max_nevents_VPDMB");
   ratio_max_nevents_VPDMB->Divide(hist_max_nevents_VPDMB);

   TH1D *ratio_max_nevents_JP2 = (TH1D *)hEventsRun_JP2->Clone("ratio_max_nevents_JP2");
   ratio_max_nevents_JP2->Divide(hist_max_nevents_JP2);

   TFile *fout = new TFile("events.root", "RECREATE");
   fout->cd();
   hEventsRun_MB->Write();
   hEventsRun_JP2->Write();
   hRatio_MB_JP2->Write();
   hist_sampled_lumi_VPDMB->Write();
   hist_prescale_VPDMB->Write();
   hist_max_nevents_VPDMB->Write();
   hist_livetime_VPDMB->Write();
   hist_sampled_lumi_JP2->Write();
   hist_prescale_JP2->Write();
   hist_livetime_JP2->Write();
   hist_max_nevents_JP2->Write();
   hist_JP2_weight->Write();
   ratio_max_nevents_VPDMB->Write();
   ratio_max_nevents_JP2->Write();

   fout->Close();

   return 0;
}
