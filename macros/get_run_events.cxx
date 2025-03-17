#include <TH1.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

// Define a struct to hold all the run information
struct RunInfo {
    double sampled_lumi_VPDMB;
    double prescale_VPDMB;
    double lifetime_VPDMB;
    int nevents_VPDMB;
    double sampled_lumi_JP2;
    double prescale_JP2;
    double lifetime_JP2;
    int nevents_JP2;
    double JP2_weight;
};

std::map<int, RunInfo> CreateRunInfoMap(const char* filename = "../lists/files_lumi_jaime/lumi_VPDMB_JP2.csv") {
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
             info.lifetime_VPDMB = std::stod(token);
             
             std::getline(ss, token, ',');
             info.nevents_VPDMB = std::stoi(token);
             
             std::getline(ss, token, ',');
             info.sampled_lumi_JP2 = std::stod(token);
             
             std::getline(ss, token, ',');
             info.prescale_JP2 = std::stod(token);
             
             std::getline(ss, token, ',');
             info.lifetime_JP2 = std::stod(token);
             
             std::getline(ss, token, ',');
             info.nevents_JP2 = std::stoi(token);
             
             std::getline(ss, token, ',');
             info.JP2_weight = std::stod(token);
             
             // Add to map
             runInfoMap[run] = info;
         }
         catch (const std::exception& e) {
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




int get_run_events() {


  // Create a map to store run number and associated information
  std::map<int, RunInfo> runInfoMap = CreateRunInfoMap("../lists/files_lumi_jaime/lumi_VPDMB_JP2.csv");

  TString mb_file_name= "../output/MB/tree_jets.root";
  TFile *mb_file = new TFile(mb_file_name, "READ");

  TString jp_file_name= "../output/JP2/tree_jets.root";
  TFile *jp_file = new TFile(jp_file_name, "READ");

  TH1D *mb_events_run = (TH1D *)mb_file->Get("hEventsRun");
  TH1D *jp_events_run = (TH1D *)jp_file->Get("hEventsRun");

  mb_events_run->SetName("hEventsRun_MB");
  jp_events_run->SetName("hEventsRun_JP2");

  TH1D * ratio = (TH1D *)mb_events_run->Clone("hRatio_MB_JP2");
  ratio->Divide(jp_events_run);
  ratio->SetTitle("Ratio of MB to JP2 events; Run Number; #frac{MB}{JP2}");

  TH1D * temp = (TH1D *)ratio->Clone("temp");
  temp->SetTitle(";Run Number;");
  // clear the bin content and errors
  for (int i = 1; i <= temp->GetNbinsX(); ++i) {
      temp->SetBinContent(i, 0);
      temp->SetBinError(i, 0);
  }
  TH1D * hist_sampled_lumi_VPDMB =  (TH1D *)temp->Clone("hist_sampled_lumi_VPDMB");
  TH1D * hist_prescale_VPDMB =  (TH1D *)temp->Clone("hist_prescale_VPDMB");
  TH1D * hist_lifetime_VPDMB =  (TH1D *)temp->Clone("hist_lifetime_VPDMB");
  TH1D * hist_nevents_VPDMB =  (TH1D *)temp->Clone("hist_max_nevents_VPDMB");
  TH1D * hist_sampled_lumi_JP2 =  (TH1D *)temp->Clone("hist_sampled_lumi_JP2");
  TH1D * hist_prescale_JP2 =  (TH1D *)temp->Clone("hist_prescale_JP2");
  TH1D * hist_lifetime_JP2 =  (TH1D *)temp->Clone("hist_lifetime_JP2");
  TH1D * hist_nevents_JP2 =  (TH1D *)temp->Clone("hist_max_nevents_JP2");
  TH1D * hist_JP2_weight =  (TH1D *)temp->Clone("hist_JP2_weight");

  // Fill the histograms with the corresponding values from the map
  for (const auto& entry : runInfoMap) {
      int run = entry.first;
      const RunInfo& info = entry.second;
      int bin = hist_sampled_lumi_VPDMB->GetXaxis()->FindBin(Form("%i", run));

      hist_sampled_lumi_VPDMB->SetBinContent(bin, info.sampled_lumi_VPDMB);
      hist_prescale_VPDMB->SetBinContent(bin, info.prescale_VPDMB);
      hist_lifetime_VPDMB->SetBinContent(bin, info.lifetime_VPDMB);
      hist_nevents_VPDMB->SetBinContent(bin, info.nevents_VPDMB);
      hist_sampled_lumi_JP2->SetBinContent(bin, info.sampled_lumi_JP2);
      hist_prescale_JP2->SetBinContent(bin, info.prescale_JP2);
      hist_lifetime_JP2->SetBinContent(bin, info.lifetime_JP2);
      hist_nevents_JP2->SetBinContent(bin, info.nevents_JP2);
      hist_JP2_weight->SetBinContent(bin, info.JP2_weight);
  }


  TFile *fout = new TFile("events.root", "RECREATE");
  fout->cd();
  mb_events_run->Write();
  jp_events_run->Write();
  ratio->Write();
  hist_sampled_lumi_VPDMB->Write();
  hist_prescale_VPDMB->Write();
  hist_lifetime_VPDMB->Write();
  hist_nevents_VPDMB->Write();
  hist_sampled_lumi_JP2->Write();
  hist_prescale_JP2->Write();
  hist_lifetime_JP2->Write();
  hist_nevents_JP2->Write();
  hist_JP2_weight->Write();


  fout->Close();



  return 0;
}
