// build_data.C — detector-level jet spectra of the exclusive trigger-level partition, per level and
// physics-eta block, from the merged Stage-1 data tree (ResultTree: one entry per event, per-jet arrays).
//
// Level gates (fired = hardware decision, should = trigger-simulator decision on the same event):
//   l0: fired_JP0                             && should_JP0 && !should_JP1
//   l1: (fired_JP0 || fired_JP1)              && should_JP1 && !should_JP2
//   l2: (fired_JP0 || fired_JP1 || fired_JP2) && should_JP2
//   jp1i, jp2i: the inclusive single triggers (should_JP1, should_JP2 regardless of the higher patch)
//   ht2: fired_HT2 && should_HT2, jets holding the firing tower (no raw floor)
//   mb : VPDMB-nobsmd hardware bit, mode "mb": the merged min-bias tree, no jet-patch match, no floor,
//        min-bias luminosity. Without a patch requirement the min-bias stream carries junk jets the
//        embedding does not model (single constituents at high pT), so its jets take the min-bias
//        quality cuts below, on data and response alike.
// Jet selection per level: trigger match of that level, |eta_det| < 0.9, |eta| < 0.9 (physics eta),
// R_T <= 0.95, raw-pT floor (l1 >= 6.0, l2 >= 8.4 GeV); the UE-subtracted pT is filled with unit weight.
// Eta blocks by physics eta: 00_05 (|eta| < 0.5) and 05_09 (0.5 <= |eta| < 0.9);
// eta = asinh(sinh(eta_det) - vz/225.405), the inverse of the Stage-1 eta_det definition.
// Events: |vz| < 60 cm, run not in the bad-run lists and with a luminosity entry (JP2 / min-bias table).
//
// Systematic variants of the data (variants.h), all filled in the same pass into their own files: the
// filled pT is pT_corr + (1 - f) x rho x area (UE factor f = 0.86 / 1.18) and pT x (1 - 0.01 (1 - R_T))
// for the 1 % track removal; the selection stays on the nominal quantities.
//
// Input : <output>/merged_data[MB]_R<R>.root
// Output: <results>/data_blocks_R<R>[_mb][_<variant>].root with d550_<level>_<blk> (550 x [5,60] GeV, counts)
// root -l -b -q 'build_data.C+("0.5", "jp")'
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TH1D.h>

#include <cstdio>
#include <set>
#include <string>
#include <vector>

#include "common.h"
#include "variants.h"

using namespace CrossSectionConfig;
using namespace Stage2;

// ---- selection constants ---------------------------------------------------------------------------

// the vertex window, the acceptance and the jet-quality cuts are shared with the response (common.h)
static const double kEtaDetRadius = 225.405; // the radius Stage-1 projects a track to for eta_det [cm]

// the three eta blocks written per level: the whole acceptance and the two physics-eta blocks
static const int kNBlocks = 3;
static const char *kBlockName[kNBlocks] = {"all", "00_05", "05_09"};
static const char *kBlockCut[kNBlocks] = {"true", "abs(jeta) < 0.5", "abs(jeta) >= 0.5"};

// one exclusive or inclusive trigger level of the data
struct DataLevel {
   const char *name;     // histogram tag of the level
   const char *gate;     // event-level selection (hardware bit and simulator decision)
   const char *match;    // per-jet trigger match required of the level ("" for min-bias)
   double rawPtFloor;    // lowest raw jet pT the level accepts [GeV]
};

// the levels of the jet-patch / high-tower stream, and the single level of the min-bias stream
static std::vector<DataLevel> LevelsOfMode(bool minBias)
{
   if (minBias) return {{"mb", "isTriggerEvent", "", 0.0}};
   return {{"jp0", "l0", "trigger_match_JP0", 0.0},
           {"jp1", "l1", "trigger_match_JP1", kRawPtFloorJp1},
           {"jp2", "l2", "trigger_match_JP2", kRawPtFloorJp2},
           {"jp1i", "(fired_JP0 || fired_JP1) && should_JP1", "trigger_match_JP1", kRawPtFloorJp1},
           {"jp2i", "(fired_JP0 || fired_JP1 || fired_JP2) && should_JP2", "trigger_match_JP2", kRawPtFloorJp2},
           {"ht2", "fired_HT2 && should_HT2", "trigger_match_HT2", 0.0}};
}

// the per-jet selection of one level inside one eta block, as an RDataFrame expression
static std::string JetSelection(const DataLevel &level, int block, bool minBias)
{
   if (minBias)
      return Form("abs(det_eta) < %g && abs(jeta) < %g && neutral_fraction > %g && neutral_fraction < %g && "
                  "ptLead < %g * pt_corrected && (%s)",
                  kEtaMax, kEtaMax, kMbNeutralFractionMin, kMbNeutralFractionMax, kMbLeadingFractionMax,
                  kBlockCut[block]);
   return Form("%s && abs(det_eta) < %g && abs(jeta) < %g && neutral_fraction <= %g && pt >= %g && (%s)",
               level.match, kEtaMax, kEtaMax, kNeutralFractionMax, level.rawPtFloor, kBlockCut[block]);
}

// ---- the build --------------------------------------------------------------------------------------

// mode: "jp" for the jet-patch / high-tower stream, "mb" for the min-bias stream
void build_data(const char *jetR = "0.5", const char *mode = "jp")
{
   ROOT::EnableImplicitMT(kImtThreads);
   TH1::SetDefaultSumw2();
   AnalysisConfig cfg;
   const bool minBias = std::string(mode) == "mb";

   // 1. the events: good runs with a sampled luminosity, vertex inside the window
   const std::vector<int> badRuns = AllBadRuns(cfg);
   const std::set<int> lumiRuns = LumiRuns(cfg, minBias);
   const std::string input = cfg.datapath + (minBias ? "merged_dataMB_R" : "merged_data_R") + jetR + ".root";
   ROOT::RDataFrame frame("ResultTree", input);
   auto events =
      frame
         .Filter([badRuns, lumiRuns](int run) { return run >= 0 && !IsBadRun(badRuns, run) && lumiRuns.count(run) > 0; },
                 {"runid1"})
         .Filter(Form("vz > -%g && vz < %g", kEventVzMax, kEventVzMax))
         .Define("l0", "fired_JP0 && should_JP0 && !should_JP1")
         .Define("l1", "(fired_JP0 || fired_JP1) && should_JP1 && !should_JP2")
         .Define("l2", "(fired_JP0 || fired_JP1 || fired_JP2) && should_JP2")
         .Define("jeta", Form("ROOT::VecOps::RVec<double> e(det_eta.size()); for (size_t i = 0; i < det_eta.size(); ++i) "
                              "e[i] = asinh(sinh(det_eta[i]) - vz/%g); return e;",
                              kEtaDetRadius));

   // 2. one filled-pT column per variant: the nominal pT plus the UE shift and the track removal
   const std::vector<std::string> variants = Syst::DataBuilds();
   for (size_t iv = 0; iv < variants.size(); ++iv) {
      const Syst::Variant variant = Syst::Parse(variants[iv]);
      events = events.Define(Form("ptfill%zu", iv),
                             Form("pt_corrected + %.6f*jet_area*bg_density - %.6f*(1-neutral_fraction)*pt",
                                  1.0 - variant.ueData, variant.trkThin));
   }

   // 3. book one spectrum per (level, eta block, variant) and one event count per level
   const std::vector<DataLevel> levels = LevelsOfMode(minBias);
   const int nLevels = (int)levels.size();
   const size_t nVariants = variants.size();
   std::vector<ROOT::RDF::RResultPtr<TH1D>> spectra;
   std::vector<ROOT::RDF::RResultPtr<ULong64_t>> nEvents;
   for (int l = 0; l < nLevels; ++l) {
      auto levelEvents = events.Filter(levels[l].gate);
      nEvents.push_back(levelEvents.Count());
      for (int block = 0; block < kNBlocks; ++block) {
         auto selected = levelEvents.Define("sel", JetSelection(levels[l], block, minBias));
         for (size_t iv = 0; iv < nVariants; ++iv) {
            auto filled = selected.Define("ptsel", Form("ptfill%zu[sel]", iv));
            spectra.push_back(filled.Histo1D({Form("d550_%s_%s__%zu", levels[l].name, kBlockName[block], iv),
                                              ";p_{T} [GeV];jets", kNDetFine, kDetFineMin, kDetFineMax},
                                             "ptsel"));
         }
      }
   }

   // 4. run the event loop before any output file is opened (the nominal file is read by other macros)
   *nEvents[0];

   // 5. one output file per variant
   for (size_t iv = 0; iv < nVariants; ++iv) {
      const std::string tag = variants[iv].empty() ? "" : "_" + variants[iv];
      const std::string output = RadiusDir(jetR) + "data_blocks_R" + jetR + (minBias ? "_mb" : "") + tag + ".root";
      TFile out(output.c_str(), "RECREATE");
      for (int l = 0; l < nLevels; ++l) {
         for (int block = 0; block < kNBlocks; ++block) {
            TH1D *spectrum = spectra[(l * kNBlocks + block) * nVariants + iv].GetPtr();
            spectrum->Write(Form("d550_%s_%s", levels[l].name, kBlockName[block]));
         }
      }
      out.Close();
      printf("[data] wrote %s\n", output.c_str());
   }
   for (int l = 0; l < nLevels; ++l)
      printf("[data] level %s: %llu events, %.0f jets (all eta)\n", levels[l].name,
             (unsigned long long)*nEvents[l], spectra[(l * kNBlocks) * nVariants]->GetEntries());
}
