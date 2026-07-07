// hw_ratio.C -----------------------------------------------------------------
// Measure the hardware-vs-simulator trigger-match correction R(pt) for one STAR
// jet-patch trigger (JP0/JP1/JP2), pp200 Run-12 inclusive jets, anti-kT R=0.5.
//
// WHY R EXISTS.  The data and the unfolding response are gated by two slightly
// different trigger "rulers":
//   * the DATA fire on the real HARDWARE decision -- the bit-7 kOnline jet-patch
//     objects stored in the picos (verified 513/513 identical to Dmitry's skim
//     values, i.e. effectively the DSM/L0 decision).  ppAnalysis::match_jp
//     gates on them natively, so trigger_match_<trig> IS the hardware decision.
//   * the RESPONSE (embedding) can only ever be gated by the SOFTWARE trigger
//     SIMULATOR (native isJP<trig>() = patch GetADC over the DSM register),
//     because simulated events have no hardware to fire.
// Near threshold the two rulers disagree at the ~1% level: R measures that
// difference in situ, per reco-pT bin:
//     R(pt) = N(hardware-matched) / N(software-simulator-matched).
//
// R CORRECTS THE DATA, NOT THE EMBEDDING.  Dividing the hardware-gated data
// spectrum by R re-expresses it as "what the software-simulator gate would have
// counted", so the DATA then speak the same trigger language as the response.
//
// ONE DATASET, TWO STAGE-1 TREES (not two data productions).  The picos carry
// BOTH decisions (an offline emulator instance and a kOnline StarTrigSimuOnline
// instance in a single Stage-0 pass), but each Stage-1 ResultTree bakes ONE
// trigger_match + one jp_patch_adc.  R here pairs the HARDWARE trigger_match of
// THIS repo's production (numerator) with the EMULATOR gate jp_patch_adc>=reg+1
// of a sim-native tree (denominator, paperrepro).  Two cheap analysis passes,
// not two productions.  R (and That) are measured once and frozen into
// config.h::TrigEffMeas (C = R x That).
//
// Output: hw_ratio_<trig>.root (hHw, hSim, hR) and hw_ratio_<trig>.pdf.
//
// Usage (ROOT>=6), all arguments optional:
//   root -l -b -q 'hw_ratio.C("JP1")'
//   root -l -b -q 'hw_ratio.C("JP2", 37, "hardware.root", "simulator.root")'
// adcRegisterPlus1 = isJP DSM register + 1 (JP1: 28+1=29, JP2: 36+1=37); 0 auto-sets it.
// -----------------------------------------------------------------------------
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <string>
#include <vector>

void hw_ratio(const char *trig = "JP2", int adcRegisterPlus1 = 0,
              const char *hwFile = "", const char *simFile = "")
{
   const std::string T = trig;
   if (adcRegisterPlus1 == 0)
      adcRegisterPlus1 = (T == "JP1") ? 29 : (T == "JP0") ? 21 : 37; // isJP register + 1

   // Numerator (hardware): THIS repo's unified data merge (all triggers, select
   // internally via trigger_match_<T>).  Denominator (simulator): a sim-native
   // tree gated by the emulator patch ADC.
   const std::string hwPath = std::string(hwFile).empty()
      ? "/gpfs01/star/pwg/prozorov/study_pp2012/fable5/jets_pp_2012/output/merged_data_R0.5.root"
      : hwFile;
   const std::string simPath = std::string(simFile).empty()
      ? "/gpfs01/star/pwg/prozorov/study_pp2012/me/jets_pp_2012_paperrepro"
        "/output/" + T + "/merged_data_" + T + "_R0.5.root"
      : simFile;

   const std::vector<double> QB = {6.9, 8.2, 9.7, 11.5, 13.6, 16.1, 19, 22.5, 26.6, 31.4, 37.2, 44, 52};
   const int nb = (int)QB.size() - 1;

   // Same analysis jet selection on both sides.  Numerator: the HARDWARE match
   // branch.  Denominator: the software-simulator gate (jp_patch_adc>=reg+1;
   // already implies a matched patch, so no separate trigger_match term).
   const std::string jetCut = "abs(det_eta)<0.5 && neutral_fraction<0.95";
   const std::string selHw  = "trigger_match_" + T + " && " + jetCut;
   const std::string selSim = "jp_patch_adc>=" + std::to_string(adcRegisterPlus1) + " && " + jetCut;

   ROOT::RDataFrame dHw("ResultTree", hwPath);
   ROOT::RDataFrame dSim("ResultTree", simPath);
   auto hHwPtr  = dHw.Define("m", selHw).Define("ptSel", "pt_corrected[m]")
                     .Histo1D({"hHw", "", nb, QB.data()}, "ptSel");
   auto hSimPtr = dSim.Define("m", selSim).Define("ptSel", "pt_corrected[m]")
                      .Histo1D({"hSim", "", nb, QB.data()}, "ptSel");

   // R(pt) = hardware / simulator.  The two gates overlap but are NOT nested (R
   // even exceeds 1 in the top bin), so use ordinary error propagation -- NOT
   // the binomial ("B") option, which assumes num is a subset.
   TH1D *hHw  = (TH1D *)hHwPtr->Clone("hHw");
   TH1D *hSim = (TH1D *)hSimPtr->Clone("hSim");
   TH1D *hR   = (TH1D *)hHw->Clone("hR");
   hR->Divide(hSim);

   TFile fout(("hw_ratio_" + T + ".root").c_str(), "RECREATE");
   hHw->Write();
   hSim->Write();
   hR->Write();
   fout.Close();

   printf("\n%s: R(pt) = N(hardware)/N(simulator, ADC>=%d)\n", trig, adcRegisterPlus1);
   for (int i = 1; i <= nb; ++i)
      printf("  [%5.1f,%5.1f)  R=%.4f\n", QB[i - 1], QB[i], hR->GetBinContent(i));

   gStyle->SetOptStat(0);
   TCanvas c("c", "", 700, 800);
   c.Divide(1, 2, 0, 0);

   c.cd(1);
   gPad->SetLogy();
   gPad->SetGrid();
   hHw->SetTitle(Form("%s hardware-vs-simulator trigger match;;jets / bin", trig));
   hHw->SetMarkerStyle(20);
   hHw->SetMarkerColor(kBlue + 1);
   hHw->SetLineColor(kBlue + 1);
   hHw->Draw("E1");
   hSim->SetMarkerStyle(24);
   hSim->SetMarkerColor(kRed + 1);
   hSim->SetLineColor(kRed + 1);
   hSim->Draw("E1 SAME");
   TLegend lg(0.42, 0.72, 0.88, 0.88);
   lg.SetBorderSize(0);
   lg.AddEntry(hHw, "DATA hardware (kOnline)", "lp");
   lg.AddEntry(hSim, Form("software simulator (isJP, ADC#geq%d)", adcRegisterPlus1), "lp");
   lg.Draw();

   c.cd(2);
   gPad->SetGrid();
   hR->SetTitle(";jet p_{T} [GeV];R = N_{hardware} / N_{simulator}");
   hR->GetYaxis()->SetRangeUser(0.90, 1.10);
   hR->SetMarkerStyle(21);
   hR->SetMarkerColor(kBlack);
   hR->SetLineColor(kBlack);
   hR->Draw("E1");
   TLine one(QB.front(), 1, QB.back(), 1);
   one.SetLineStyle(2);
   one.Draw();

   c.SaveAs(("hw_ratio_" + T + ".pdf").c_str());
}
