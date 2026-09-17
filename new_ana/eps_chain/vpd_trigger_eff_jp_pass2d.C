// vpd_trigger_eff_jp_pass2d.C — pass 2 of the min-bias chain-efficiency measurement:
//
//   eps_chain(pT) = P(the VPDMB-nobsmd hardware bit fired AND the VPD vertex matched | a jet at pT)
//
// It is the probability that a jet event would have entered the min-bias stream at all, and it corrects
// the raw min-bias spectrum before unfolding (stage2/unfold.C). The numerator is weighted by the per-run
// VPDMB-nobsmd prescale, so a prescaled-away event still counts as "would have fired".
//
// Measured in THREE jet-patch samples (JP0, JP1, JP2) so that a genuine pT dependence can be told from a
// trigger-sample bias. The jets carry the same NEF window and zLead cut as the min-bias analysis, so the
// efficiency matches the selection it corrects.
//
// Input : <bitsdir>/bits4_chunk_*.root (TrigBits2, from pass 1), the merged Stage-1 data tree of the
//         radius, and the per-run VPDMB-nobsmd luminosity/prescale table.
// Output: <bitsdir>/eps_chain3_R<R>.root with den{0,1,2} (jets in the sample), chain{0,1,2} (prescale-
//         weighted numerator) and raw{0,1,2} (unweighted MB-bit count), index = {JP0, JP1, JP2}.
// Run   : inside star_star.simg:
//         root -l -b -q 'vpd_trigger_eff_jp_pass2d.C("bits4","0.5")'
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "../config.h"
using namespace CrossSectionConfig;

// Jet selection, identical to the min-bias analysis (see CLAUDE.md, "MB jet-quality cuts").
static const double kDetEtaMax = 0.5;  // |eta_det| acceptance of the published measurement
static const double kJetPtMin = 5.0;   // lowest pT the efficiency is measured at
static const double kNefMin = 0.05;    // neutral energy fraction window: reject all-charged jets ...
static const double kNefMax = 0.90;    // ... and tower-only (hot-tower) jets
static const double kZLeadMax = 0.60;  // ptLead / pt on the raw jet: reject single-track jets
static const int kMaxJetsPerEvent = 2000;

// Bits packed per event from the pass-1 tree.
static const unsigned char kBitMbFired = 1;   // VPDMB-nobsmd in the event header
static const unsigned char kBitJp1 = 2;
static const unsigned char kBitJp0 = 4;
static const unsigned char kBitJp2 = 8;
static const unsigned char kBitVpdMatch = 32; // |vz_vpd - vz_tpc| < 6 cm

// The three jet-patch samples the efficiency is measured in.
static const int kNSamples = 3;
static const char *kSampleName[kNSamples] = {"JP0 sample", "JP1 sample", "JP2 sample"};

// pT bins of the efficiency (finer than the analysis bins below 10 GeV, where the turn-on lives).
static const std::vector<double> kPtBins = {5.0, 5.9, 6.9, 8.2, 10, 12, 14, 17, 20, 24, 29, 35, 45};

// The event join key. Stage-1 ran the data pass with -geantnum 1, so ResultTree's (runid, eventid) is
// (adjusted hash of the pico basename, entry index) and the hash COLLIDES across files; the real run
// number (runid1 in ResultTree, realrun in TrigBits2) is mixed in to break those collisions.
static ULong64_t JoinKey(int hashrun, int entryIndex, int realRun)
{
   const ULong64_t packed = (ULong64_t)(unsigned)hashrun << 32 | (unsigned)entryIndex;
   return packed ^ ((ULong64_t)(unsigned)realRun * 0x9E3779B97F4A7C15ULL);
}

// Per-event trigger information carried from pass 1 into the join.
struct EventBits {
   unsigned char bits{0}; // kBitMbFired | kBitJp0 | kBitJp1 | kBitJp2 | kBitVpdMatch
   float prescale{-1};    // VPDMB-nobsmd prescale of the run, -1 when the run has no entry
};

// Per-run VPDMB-nobsmd prescale: columns are run, start, stop, fill, luminosity, prescale.
static std::unordered_map<int, double> ReadPrescales()
{
   std::unordered_map<int, double> prescale;
   std::ifstream table(kWorkDir + "../lists/lum_perrun_VPDMB-nobsmd.txt");
   std::string line;
   while (std::getline(table, line)) {
      std::istringstream row(line);
      long run, tStart, tStop, fill;
      double luminosity, ps;
      if (row >> run >> tStart >> tStop >> fill >> luminosity >> ps) prescale[(int)run] = ps;
   }
   printf("[pass2d] prescales for %zu runs\n", prescale.size());
   return prescale;
}

// The pass-1 trees: one entry per JP0/JP1/JP2-fired event, indexed by the join key.
static std::unordered_map<ULong64_t, EventBits> ReadTriggerBits(const char *bitsdir,
                                                                const std::unordered_map<int, double> &prescale)
{
   std::unordered_map<ULong64_t, EventBits> byEvent;
   TChain chain("TrigBits2");
   chain.Add(Form("%s/bits4_chunk_*.root", bitsdir));
   Int_t hashrun, fidx, realrun;
   Bool_t mbFired, jp0, jp1, jp2, vpdMatch;
   chain.SetBranchAddress("hashrun", &hashrun);
   chain.SetBranchAddress("fidx", &fidx);
   chain.SetBranchAddress("realrun", &realrun);
   chain.SetBranchAddress("mb11", &mbFired);
   chain.SetBranchAddress("jp0", &jp0);
   chain.SetBranchAddress("jp1", &jp1);
   chain.SetBranchAddress("jp2", &jp2);
   chain.SetBranchAddress("match", &vpdMatch);

   const Long64_t nEntries = chain.GetEntries();
   for (Long64_t i = 0; i < nEntries; ++i) {
      chain.GetEntry(i);
      EventBits event;
      if (mbFired) event.bits |= kBitMbFired;
      if (jp1) event.bits |= kBitJp1;
      if (jp0) event.bits |= kBitJp0;
      if (jp2) event.bits |= kBitJp2;
      if (vpdMatch) event.bits |= kBitVpdMatch;
      const auto it = prescale.find(realrun);
      event.prescale = (it != prescale.end()) ? (float)it->second : -1.f;
      byEvent[JoinKey(hashrun, fidx, realrun)] = event;
   }
   printf("[pass2d] TrigBits2 (jp0||jp1||jp2): %lld events loaded\n", nEntries);
   return byEvent;
}

// Denominator, prescale-weighted numerator and unweighted MB count per sample, from the Stage-1 jets
// joined to the pass-1 trigger bits.
static void FillSamples(const char *jetR, const std::unordered_map<ULong64_t, EventBits> &byEvent, TH1D *den[],
                        TH1D *chainNum[], TH1D *raw[])
{
   TChain jets("ResultTree");
   jets.Add((kDataPath + "merged_data_R" + jetR + ".root").c_str());
   Int_t runid, runid1, eventid, njets;
   Bool_t firedJp0, firedJp1, firedJp2;
   double ptCorr[kMaxJetsPerEvent], detEta[kMaxJetsPerEvent], nef[kMaxJetsPerEvent];
   double ptRaw[kMaxJetsPerEvent], ptLead[kMaxJetsPerEvent];
   Bool_t matchJp0[kMaxJetsPerEvent], matchJp1[kMaxJetsPerEvent], matchJp2[kMaxJetsPerEvent];
   jets.SetBranchStatus("*", 0);
   for (auto branch : {"runid", "runid1", "eventid", "njets", "fired_JP0", "fired_JP1", "fired_JP2", "pt_corrected",
                       "det_eta", "neutral_fraction", "trigger_match_JP0", "trigger_match_JP1", "trigger_match_JP2",
                       "pt", "ptLead"})
      jets.SetBranchStatus(branch, 1);
   jets.SetBranchAddress("runid", &runid);
   jets.SetBranchAddress("runid1", &runid1);
   jets.SetBranchAddress("eventid", &eventid);
   jets.SetBranchAddress("njets", &njets);
   jets.SetBranchAddress("fired_JP0", &firedJp0);
   jets.SetBranchAddress("fired_JP1", &firedJp1);
   jets.SetBranchAddress("fired_JP2", &firedJp2);
   jets.SetBranchAddress("pt_corrected", ptCorr);
   jets.SetBranchAddress("det_eta", detEta);
   jets.SetBranchAddress("neutral_fraction", nef);
   jets.SetBranchAddress("trigger_match_JP0", matchJp0);
   jets.SetBranchAddress("trigger_match_JP1", matchJp1);
   jets.SetBranchAddress("trigger_match_JP2", matchJp2);
   jets.SetBranchAddress("pt", ptRaw);
   jets.SetBranchAddress("ptLead", ptLead);

   Long64_t nJoined = 0, nSelected = 0, nJp0Mismatch = 0;
   const Long64_t nEntries = jets.GetEntries();
   for (Long64_t i = 0; i < nEntries; ++i) {
      jets.GetEntry(i);
      if (!firedJp0 && !firedJp1 && !firedJp2) continue;
      ++nSelected;
      const auto it = byEvent.find(JoinKey(runid, eventid, runid1));
      if (it == byEvent.end()) continue;
      ++nJoined;

      const EventBits &event = it->second;
      // sanity check of the join: the JP0 bit must agree between the pico header and Stage-1
      if (((event.bits & kBitJp0) != 0) != (bool)firedJp0) ++nJp0Mismatch;
      const bool mbFired = event.bits & kBitMbFired;
      const bool vpdMatch = event.bits & kBitVpdMatch;

      for (int j = 0; j < njets && j < kMaxJetsPerEvent; ++j) {
         if (std::fabs(detEta[j]) >= kDetEtaMax || ptCorr[j] < kJetPtMin) continue;
         if (nef[j] <= kNefMin || nef[j] >= kNefMax) continue;
         if (ptRaw[j] > 0 && ptLead[j] / ptRaw[j] >= kZLeadMax) continue;
         for (int s = 0; s < kNSamples; ++s) {
            // the jet must sit in a fired patch of that trigger, in an event that trigger recorded
            const bool inSample = (s == 0)   ? (firedJp0 && matchJp0[j])
                                  : (s == 1) ? (firedJp1 && matchJp1[j])
                                             : (firedJp2 && matchJp2[j]);
            if (!inSample) continue;
            den[s]->Fill(ptCorr[j]);
            if (mbFired) raw[s]->Fill(ptCorr[j]);
            if (mbFired && vpdMatch && event.prescale > 0) chainNum[s]->Fill(ptCorr[j], event.prescale);
         }
      }
   }
   printf("[pass2d] R=%s events=%lld joined=%lld  jp0-bit mismatch=%lld\n", jetR, nSelected, nJoined, nJp0Mismatch);
}

// The efficiency per sample and pT bin, as it is read in the log.
static void PrintEfficiency(TH1D *den[], TH1D *chainNum[], TH1D *raw[])
{
   const int nBins = (int)kPtBins.size() - 1;
   for (int s = 0; s < kNSamples; ++s) {
      printf("\n  [%s]  pT bin      eps_chain_hw       raw11  Njets\n", kSampleName[s]);
      for (int b = 1; b <= nBins; ++b) {
         const double nJets = den[s]->GetBinContent(b);
         if (nJets <= 0) continue;
         printf("%5.1f-%-5.1f  %.4f +- %.4f   %5.0f  %9.0f\n", kPtBins[b - 1], kPtBins[b],
                chainNum[s]->GetBinContent(b) / nJets, chainNum[s]->GetBinError(b) / nJets, raw[s]->GetBinContent(b),
                nJets);
      }
   }
}

void vpd_trigger_eff_jp_pass2d(const char *bitsdir, const char *R = "0.5")
{
   // 1. the per-run VPDMB-nobsmd prescales, the weight of the numerator
   const std::unordered_map<int, double> prescale = ReadPrescales();

   // 2. the trigger bits of every JP-fired event, from pass 1
   const std::unordered_map<ULong64_t, EventBits> byEvent = ReadTriggerBits(bitsdir, prescale);

   // 3. the three samples, filled from the Stage-1 jets joined to those bits
   const int nBins = (int)kPtBins.size() - 1;
   TH1D *den[kNSamples], *chainNum[kNSamples], *raw[kNSamples];
   for (int s = 0; s < kNSamples; ++s) {
      den[s] = new TH1D(Form("den%d", s), "", nBins, kPtBins.data());
      chainNum[s] = new TH1D(Form("chain%d", s), "", nBins, kPtBins.data());
      raw[s] = new TH1D(Form("raw%d", s), "", nBins, kPtBins.data());
      chainNum[s]->Sumw2(); // the numerator carries prescale weights
   }
   FillSamples(R, byEvent, den, chainNum, raw);

   // 4. report and write
   PrintEfficiency(den, chainNum, raw);
   TFile out(Form("%s/eps_chain3_R%s.root", bitsdir, R), "RECREATE");
   for (int s = 0; s < kNSamples; ++s) {
      den[s]->Write();
      chainNum[s]->Write();
      raw[s]->Write();
   }
   out.Close();
}
