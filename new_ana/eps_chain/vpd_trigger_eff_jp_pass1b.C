// vpd_trigger_eff_jp_pass1b.C — pass 1 of the min-bias chain-efficiency measurement.
//
// Reads the jet-patch picos listed in <chunk.list> and writes one TrigBits2 row per event that fired
// JP0, JP1 or JP2: the join key into the Stage-1 ResultTree, the real run id (for the prescale lookup),
// which min-bias and jet-patch hardware trigger ids were in the event header, and the offline VPD flags
//     valid = the VPD reconstructed a vertex
//     match = valid && |vz_vpd - vz_tpc| < 6 cm   (the min-bias vertex requirement)
// Only the event header is read, so the pass is I/O-light although it touches every pico.
//
// Input : a text file with one pico path per line.
// Output: <out.root> holding the TrigBits2 tree, the input of vpd_trigger_eff_jp_pass2d.C.
// Run   : inside star_star.simg with lib/eventStructuredAu bound over /usr/local/eventStructuredAu
//         (see run_bits.sh):
//         root -l -b -q 'vpd_trigger_eff_jp_pass1b.C("chunk.list","out.root")'
R__LOAD_LIBRARY(/usr/local/eventStructuredAu/libTStarJetPico.so)
#include "/usr/local/eventStructuredAu/TStarJetPicoEvent.h"
#include "/usr/local/eventStructuredAu/TStarJetPicoEventHeader.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include <climits>
#include <cmath>
#include <fstream>
#include <string>

// STAR trigger ids of Run-12 pp 200 GeV.
static const int kTrigVpdMbNoBsmd = 370011; // VPDMB-nobsmd, the min-bias trigger of the measurement
static const int kTrigVpdMb = 370001;       // VPDMB (with BSMD in readout)
static const int kTrigJP0 = 370601;
static const int kTrigJP1 = 370611;
static const int kTrigJP2 = 370621;

// The min-bias analysis accepts an event when the VPD vertex agrees with the TPC one to this distance.
static const double kVpdMatchWindowCm = 6.0;

// TStarJetPicoEventHeader returns this sentinel when the VPD reconstructed no vertex.
static const double kVpdVzInvalid = -500.0;

// The Stage-1 data pass ran with -geantnum 1, so its ResultTree runid is not the run number but the
// adjusted TString::Hash of the pico BASENAME (ppAnalysis.cxx, GetNOfCurrentEvent). Reproduce it here so
// that pass 2 can join the two trees event by event.
static Int_t AdjustedFileNameHash(const char *path)
{
   TString base = gSystem->BaseName(path);
   UInt_t hash = base.Hash();
   while (hash > (UInt_t)(INT_MAX - 100000)) hash -= INT_MAX / 4;
   if (hash < 1000000) hash += 1000001;
   return (Int_t)hash;
}

void vpd_trigger_eff_jp_pass1b(const char *listfile, const char *outfile)
{
   // 1. the output tree: the join key, the run, the hardware bits and the VPD flags
   TFile fout(outfile, "RECREATE");
   TTree bits("TrigBits2", "hardware MB bits + offline VPD flags, JP-fired events");
   Int_t hashrun, fidx, realrun;
   Bool_t mb11, mb01, jp0, jp1, jp2;
   Bool_t valid, match;
   bits.Branch("hashrun", &hashrun, "hashrun/I"); // ResultTree runid  = adjusted basename hash
   bits.Branch("fidx", &fidx, "fidx/I");          // ResultTree eventid = entry index in the pico
   bits.Branch("realrun", &realrun, "realrun/I"); // ResultTree runid1  = the real run number
   bits.Branch("mb11", &mb11, "mb11/O");
   bits.Branch("mb01", &mb01, "mb01/O");
   bits.Branch("jp0", &jp0, "jp0/O");
   bits.Branch("jp1", &jp1, "jp1/O");
   bits.Branch("jp2", &jp2, "jp2/O");
   bits.Branch("valid", &valid, "valid/O");
   bits.Branch("match", &match, "match/O");

   std::ifstream list(listfile);
   std::string line;
   Long64_t nSeen = 0, nKept = 0;
   int nFiles = 0;

   // 2. one pico at a time
   while (std::getline(list, line)) {
      if (line.empty()) continue;
      TFile *pico = TFile::Open(line.c_str());
      if (!pico || pico->IsZombie()) {
         printf("[pass1b] BAD FILE %s\n", line.c_str());
         if (pico) delete pico;
         continue;
      }
      ++nFiles;
      TTree *jetTree = (TTree *)pico->Get("JetTree");
      TStarJetPicoEvent *event = nullptr;
      jetTree->SetBranchAddress("PicoJetTree", &event);
      // only the header is needed: switch off the heavy per-event payloads
      for (const char *branch : {"fPrimaryTracks*", "fFtpcPrimaryTracks*", "fTowers*", "fV0s*", "fTrigObjs*"})
         jetTree->SetBranchStatus(branch, 0);
      hashrun = AdjustedFileNameHash(line.c_str());

      const Long64_t nEntries = jetTree->GetEntries();
      nSeen += nEntries;
      for (Long64_t i = 0; i < nEntries; ++i) {
         jetTree->GetEntry(i);
         TStarJetPicoEventHeader *header = event->GetHeader();

         // which hardware triggers the event carries
         mb11 = mb01 = jp0 = jp1 = jp2 = false;
         for (int k = 0; k < header->GetNOfTriggerIds(); ++k) {
            switch (header->GetTriggerId(k)) {
            case kTrigVpdMbNoBsmd: mb11 = true; break;
            case kTrigVpdMb: mb01 = true; break;
            case kTrigJP0: jp0 = true; break;
            case kTrigJP1: jp1 = true; break;
            case kTrigJP2: jp2 = true; break;
            }
         }
         if (!jp2 && !jp1 && !jp0) continue; // the measurement lives in the jet-patch samples only

         fidx = (Int_t)i;
         realrun = header->GetRunId();
         const double vpdVz = header->GetVpdVz();
         valid = (vpdVz > kVpdVzInvalid);
         match = valid && std::fabs(vpdVz - header->GetPrimaryVertexZ()) < kVpdMatchWindowCm;
         bits.Fill();
         ++nKept;
      }
      delete pico;
   }

   // 3. opening and deleting the input files moved gDirectory; write into the output file
   fout.cd();
   bits.Write();
   fout.Close();
   printf("[pass1b] %s: files=%d events=%lld jp2kept=%lld\n", outfile, nFiles, nSeen, nKept);
}
