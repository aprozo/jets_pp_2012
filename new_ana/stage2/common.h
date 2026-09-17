#ifndef STAGE2_COMMON_H
#define STAGE2_COMMON_H
//
// common.h — the pieces the Stage-2 engine shares: the pt-hat sample table of the embedding request,
// the run selection of the data, the trigger-level table with its detector windows and raw-pT floors,
// the fine binning of the response, the reader of a Stage-1 matched tree and the jet pairing.
//
// Used by stage2/build_data.C, stage2/build_resp.C, stage2/unfold.C, stage2/systematics.C and
// published/embedding_diff.C. Header only: plain structs and static/inline functions, no state.
//
// Physics content in one place:
//   samples   sigma(pt-hat bin) of the Pythia-6 embedding request and the correction of the three
//             softest bins, so every consumer normalises the embedding the same way;
//   runs      the runs the cross section is measured on: not in a bad-run list and with a sampled
//             luminosity in the stream's own table;
//   levels    the exclusive jet-patch partition jp0/jp1/jp2, min-bias, the inclusive single triggers
//             jp1i/jp2i and the high tower ht2, with the published detector-pT window of each;
//   matching  one-to-one detector-to-particle jet pairing by global minimum dR.
//
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../config.h"

namespace Stage2 {

using CrossSectionConfig::AnalysisConfig;

// ---- pt-hat samples of the embedding request -------------------------------------------------------

// Pythia-6 cross section of each generated pt-hat bin, in mb.
static const std::map<std::string, double> kSigmaMb = {
   {"pt2_3", 9.0012},        {"pt3_4", 1.46253},       {"pt4_5", 0.354566},      {"pt5_7", 0.151622},
   {"pt7_9", 0.0249062},     {"pt9_11", 0.00584527},   {"pt11_15", 0.00230158},  {"pt15_20", 0.000342755},
   {"pt20_25", 4.57002e-05}, {"pt25_35", 9.72535e-06}, {"pt35_45", 4.69889e-07}, {"pt45_55", 2.69202e-08},
   {"pt55_-1", 1.43453e-09}};

// The generator over-reports the cross section of the three softest samples by 22.8 %, 5.1 % and 1.4 %;
// each sample weight carries the inverse as a correction.
static const std::map<std::string, double> kSoftSampleCorr = {
   {"pt2_3", 1. / 1.228}, {"pt3_4", 1. / 1.051}, {"pt4_5", 1. / 1.014}};

static const double kMbToPb = 1e9; // 1 mb in pb, so that every response histogram comes out in pb

// true when the sample name is one of the generated pt-hat bins
inline bool IsKnownSample(const std::string &sample)
{
   return kSigmaMb.count(sample) > 0;
}

// soft-sample correction of a named pt-hat sample (1 for every sample above pt-hat 5 GeV)
inline double SoftSampleCorrection(const std::string &sample)
{
   return kSoftSampleCorr.count(sample) ? kSoftSampleCorr.at(sample) : 1.0;
}

// soft-sample correction from the midpoint of the generated pt-hat bin, for the merged trees, whose
// samples are interleaved and told apart by pthat_mid only
inline double SoftSampleCorrection(double pthatMid)
{
   if (pthatMid < 3) return 1. / 1.228;
   if (pthatMid < 4) return 1. / 1.051;
   if (pthatMid < 5) return 1. / 1.014;
   return 1.0;
}

// weight of one generated event of a per-run sample, in pb: sigma x correction / N_generated
inline double SampleWeightFactor(const std::string &sample, long long nGenerated)
{
   if (nGenerated <= 0) return 0.0;
   return kSigmaMb.at(sample) * kMbToPb * SoftSampleCorrection(sample) / (double)nGenerated;
}

// weight of one generated event of a merged tree, in pb: the per-row mc_weight is sigma / N_generated [mb]
inline double MergedWeightFactor(double mcWeightMb, double pthatMid)
{
   return mcWeightMb * kMbToPb * SoftSampleCorrection(pthatMid);
}

// generated events per (sample, run), from lists/emb2021_events_per_run.txt (counted on the picos)
inline std::map<std::string, std::map<int, long long>> GeneratedEventsPerRun(const AnalysisConfig &cfg)
{
   std::map<std::string, std::map<int, long long>> nGenerated;
   std::ifstream file(cfg.workdir + "../lists/emb2021_events_per_run.txt");
   std::string line;
   while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream fields(line);
      std::string sample;
      int run = 0;
      long long nEvents = 0;
      if (fields >> sample >> run >> nEvents && nEvents > 0) nGenerated[sample][run] = nEvents;
   }
   return nGenerated;
}

// sample and run of a Stage-1 matched tree: matched_e21_<sample>_<run>_<id>_R<R>.root (jet-patch pass)
// or matchedMB_e21_... (min-bias pass). Returns false when the name does not have that shape.
inline bool ParseMatchedFileName(const std::string &path, std::string &sample, int &run)
{
   std::string base = gSystem->BaseName(path.c_str());
   const std::string prefix = base.find("matchedMB_e21_") != std::string::npos ? "matchedMB_e21_" : "matched_e21_";
   const size_t start = base.find(prefix);
   const size_t radius = base.rfind("_R");
   if (start == std::string::npos || radius == std::string::npos) return false;
   base = base.substr(start + prefix.size(), radius - start - prefix.size()); // <sample>_<run>_<id>
   const size_t beforeId = base.rfind('_');
   const size_t beforeRun = base.rfind('_', beforeId - 1);
   if (beforeRun == std::string::npos) return false;
   sample = base.substr(0, beforeRun);
   run = std::atoi(base.substr(beforeRun + 1, beforeId - beforeRun - 1).c_str());
   return true;
}

// ---- run selection ---------------------------------------------------------------------------------

// the runs excluded from the measurement: the detector-quality and dead-time runs of config.h plus
// lists/badrun_extras.list. The luminosity sums drop exactly the same runs.
inline std::vector<int> AllBadRuns(const AnalysisConfig &cfg)
{
   std::vector<int> bad = cfg.badRuns;
   std::ifstream file(cfg.workdir + "../lists/badrun_extras.list");
   std::string line;
   while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#') continue;
      bad.push_back(std::stoi(line));
   }
   return bad;
}

inline bool IsBadRun(const std::vector<int> &bad, int run)
{
   return std::find(bad.begin(), bad.end(), run) != bad.end();
}

// the runs that have a sampled luminosity in the stream's own table: the JP2 table for the jet-patch
// (and high-tower) stream, the VPDMB table for min-bias. A run without an entry cannot be normalised
// and is dropped from data and response alike.
inline std::set<int> LumiRuns(const AnalysisConfig &cfg, bool minBias)
{
   std::set<int> runs;
   const std::string file = cfg.workdir + (minBias ? "inputs/lumi_VPDMB_true.root" : "inputs/lumi_zilong_full.root");
   TFile lumiFile(file.c_str(), "READ");
   auto *lumi = (TH1D *)lumiFile.Get(minBias ? "luminosity_MBtrue" : "luminosity_JP2");
   for (int bin = 1; lumi && bin <= lumi->GetNbinsX(); ++bin) {
      const char *label = lumi->GetXaxis()->GetBinLabel(bin);
      if (label && *label && lumi->GetBinContent(bin) > 0) runs.insert(std::stoi(label));
   }
   return runs;
}

// the runs the jet-patch cross section is measured on, sorted. The embedding is restricted to them so
// that the response carries the data's mix of trigger thresholds.
inline std::vector<int> KeptRuns(const AnalysisConfig &cfg)
{
   const std::vector<int> bad = AllBadRuns(cfg);
   std::vector<int> kept;
   TFile lumiFile((cfg.workdir + "inputs/lumi_zilong_full.root").c_str(), "READ");
   auto *lumi = (TH1D *)lumiFile.Get("luminosity_JP2");
   for (int bin = 1; lumi && bin <= lumi->GetNbinsX(); ++bin) {
      const char *label = lumi->GetXaxis()->GetBinLabel(bin);
      if (!label || !*label || lumi->GetBinContent(bin) <= 0) continue;
      const int run = std::stoi(label);
      if (!IsBadRun(bad, run)) kept.push_back(run);
   }
   std::sort(kept.begin(), kept.end());
   return kept;
}

// ---- trigger levels --------------------------------------------------------------------------------

// levels 0-2: the exclusive jet-patch partition; 3: min-bias; 4-5: the INCLUSIVE single triggers
// (should_JP1 / should_JP2 regardless of the higher patches); 6: the high tower (should_HT2)
static const int kNLevels = 7;
static const char *kLevelName[kNLevels] = {"jp0", "jp1", "jp2", "mb", "jp1i", "jp2i", "ht2"};
// the luminosity table each level is normalised by ("" = the min-bias table)
static const char *kLumiName[kNLevels] = {"JP0", "JP1", "JP2", "", "JP1", "JP2", "HT2"};

static const int kLevelMinBias = 3;
static const int kLevelHighTower = 6;

// jet-patch index of a level: 4/5 are the inclusive JP1/JP2 triggers, the exclusive levels are their own
// index, and 3 (min-bias) / 6 (high tower) keep theirs, which is outside 0-2 and means "no patch".
inline int JetPatchOfLevel(int level)
{
   if (level == 4) return 1;
   if (level == 5) return 2;
   return level;
}

// raw-pT floors of the jet-patch triggers: a JP1 (JP2) jet must carry at least the energy the patch
// threshold implies, else the turn-on of the trigger leaks into the lowest bins
static const double kRawPtFloorJp1 = 6.0;
static const double kRawPtFloorJp2 = 8.4;

inline double RawPtFloorOfPatch(int jetPatch)
{
   if (jetPatch == 1) return kRawPtFloorJp1;
   if (jetPatch == 2) return kRawPtFloorJp2;
   return 0.0;
}

// The published detector-pT windows of the levels, on the fine axis: JP0 below 22.5 GeV, JP1 above 8.2,
// JP2 above 9.7 GeV (this keeps the JP2 jets sitting at their 8.4 GeV raw floor out of the 8.2-9.7 GeV
// bin, which the inversion otherwise see-saws), the high tower above 11.5 GeV. The inclusive and the
// min-bias levels are unrestricted.
static const double kDetPtMaxJp0 = 22.5;
static const double kDetPtMinJp1 = 8.2;
static const double kDetPtMinJp2 = 9.7;
static const double kDetPtMinHt2 = 11.5;

inline bool KeepDet(int level, double detPt)
{
   if (level == 0) return detPt < kDetPtMaxJp0;
   if (level == 1) return detPt > kDetPtMinJp1;
   if (level == 2) return detPt > kDetPtMinJp2;
   if (level == kLevelHighTower) return detPt > kDetPtMinHt2;
   return true;
}

// ---- jet quality cuts shared by the data and the response -------------------------------------------

static const double kEventVzMax = 60.0;           // TPC vertex window of the measurement [cm]
static const double kNeutralFractionMax = 0.95;   // R_T ceiling of a jet-patch / high-tower jet

// The min-bias stream has no patch requirement, so it carries jets the embedding does not model (single
// towers, single tracks). Its jets must be neither a bare tower cluster nor a bare track cluster, and
// must not be carried by one constituent. Applied to the data and to the response alike.
static const double kMbNeutralFractionMin = 0.05;
static const double kMbNeutralFractionMax = 0.90;
static const double kMbLeadingFractionMax = 0.60;

// ---- fine binning and eta blocks -------------------------------------------------------------------

// Everything is accumulated on 0.1 GeV cells and only rebinned to the analysis bins at the very end,
// so that the outlier filter and the level windows act on the true detector pT.
static const int kNDetFine = 550;      // detector axis: 550 cells
static const double kDetFineMin = 5.0; // from 5 GeV (the lowest detector pT any level keeps)
static const double kDetFineMax = 60.0;
static const int kNParFine = 600; // particle axis: 600 cells
static const double kParFineMin = 0.0;
static const double kParFineMax = 60.0;
static const double kFineCell = 0.1; // GeV per cell on both axes

// eta blocks, by physics eta: the jets are split so that the response can migrate between them
static const int kNEtaBlocks = 2;
static const char *kEtaBlockName[kNEtaBlocks] = {"00_05", "05_09"};
static const double kEtaBlockSplit = 0.5; // |eta| below this is block 0, above it block 1
static const double kEtaMax = 0.9;        // physics and detector acceptance of the measurement

inline int EtaBlockOf(double eta)
{
   return std::fabs(eta) < kEtaBlockSplit ? 0 : 1;
}

static const double kFineCellsPerGeV = 10.0; // 1 / kFineCell

// 0-based fine-cell index of a coarse bin edge on an axis of 0.1 GeV cells starting at axisMin
inline int FineIndexOfEdge(double edge, double axisMin)
{
   return (int)std::lround((edge - axisMin) * kFineCellsPerGeV);
}

// ---- matched-tree event reader ---------------------------------------------------------------------

// one reconstructed (detector-level) jet of a matched-tree row
struct DetJet {
   double pt = 0;              // raw jet pT
   double ptCorrected = 0;     // UE-subtracted jet pT (the analysis pT)
   double eta = 0;             // physics eta
   double detEta = 0;          // detector eta (eta at the nominal vertex)
   double phi = 0;             // azimuth
   double neutralFraction = 0; // R_T = tower energy / jet energy
   double ptLead = 0;          // pT of the leading constituent
   double area = 0;            // fastjet area of the jet
   double bgDensity = 0;       // UE density rho of the event
   double nConstituents = 0;   // constituents including the area ghosts
   bool matchJetPatch[3] = {false, false, false}; // jet holds a firing JP0 / JP1 / JP2 patch
   bool matchHighTower = false;                   // jet holds the firing HT2 tower
   int patchAdc = -1;                             // DSM ADC of the jet's patch (-1: not stored)
   int towerAdc = -1;                             // highest single-tower DSM ADC in the jet (-1: not stored)
};

// one generated (particle-level) jet of a matched-tree row
struct ParJet {
   double pt = 0;              // raw jet pT
   double ptCorrected = 0;     // UE-subtracted jet pT
   double eta = 0;             // physics eta
   double phi = 0;             // azimuth
   double neutralFraction = 0; // neutral energy fraction of the generated jet
   double ptLead = 0;          // pT of the leading generated constituent
   double bgDensity = 0;       // particle-level UE density
   double nConstituents = 0;   // generated constituents including the area ghosts
};

// one generated event: the rows of an event are consecutive in the tree, so the jets are re-collected
// from them before anything is filled
struct Event {
   int id = -1;            // eventid of the Stage-1 pass
   double thrownVz = -999; // generated vertex z
   double recoVz = -999;   // reconstructed vertex z (-999: no vertex)
   double pthat = -1;      // generator pt-hat of the event (-1: not stored)
   double pthatMid = -1;   // midpoint of the generated pt-hat bin, the fallback for pthat
   double multiplicity = 0;                    // reconstructed track multiplicity
   bool shouldJetPatch[3] = {false, false, false}; // simulator decision for JP0 / JP1 / JP2
   bool shouldHighTower = false;               // simulator decision for HT2
   bool hasReco = false;                       // the event has a reconstructed vertex
   int patchAdcMax = -1;                       // highest patch ADC of the event (-1: not stored)
   int towerAdcMax = -1;                       // highest tower ADC of the event (-1: not stored)
   int patchThreshold[3] = {0, 0, 0};          // DSM thresholds of JP0 / JP1 / JP2 in this run
   int towerThreshold = 0;                     // DSM threshold of HT2 in this run
   std::vector<DetJet> detJets;                // reconstructed jets of the event
   std::vector<ParJet> parJets;                // generated jets of the event
};

// Row reader of a MatchedTree. The trigger-simulator branches (patch and tower ADCs, thresholds, the
// high-tower decision and match) are optional, so trees of an earlier Stage-1 pass still build every
// variant that does not move a threshold. The fragmentation branches are read only on request, so that
// the response build does not pay their I/O over the 23 GB merged tree.
struct MatchedTreeRows {
   TTreeReader reader;
   TTreeReaderValue<double> mcPt, mcPtCorrected, mcEta, mcPhi;
   TTreeReaderValue<double> recoPt, recoPtCorrected, recoEta, recoDetEta, recoPhi, recoNeutralFraction, recoPtLead,
      recoArea, recoBgDensity;
   TTreeReaderValue<double> thrownVz, recoVz, pthat, pthatMid;
   TTreeReaderValue<bool> matchJp0, matchJp1, matchJp2, shouldJp0, shouldJp1, shouldJp2;
   TTreeReaderValue<int> eventId;
   std::unique_ptr<TTreeReaderValue<double>> mcWeight;                       // sigma_sample / N_sample [mb]
   std::unique_ptr<TTreeReaderValue<bool>> matchHt2, shouldHt2;              // high-tower decision and match
   std::unique_ptr<TTreeReaderValue<int>> patchAdc, towerAdc, patchAdcMax, towerAdcMax;
   std::unique_ptr<TTreeReaderArray<int>> patchThreshold, towerThreshold;    // DSM thresholds of the run
   std::unique_ptr<TTreeReaderValue<double>> mcNeutralFraction, mcPtLead, mcBgDensity; // fragmentation only
   std::unique_ptr<TTreeReaderValue<int>> mcNConstituents, recoNConstituents;
   std::unique_ptr<TTreeReaderValue<double>> recoMultiplicity;

   // needAdc: fail rather than silently fall back when a threshold variant needs the simulator ADCs.
   // needFragmentation: also read the jet-shape and event branches the embedding comparison plots.
   MatchedTreeRows(TTree *tree, bool needAdc, bool needFragmentation = false)
      : reader(tree), mcPt(reader, "mc_pt"), mcPtCorrected(reader, "mc_pt_corrected"), mcEta(reader, "mc_eta"),
        mcPhi(reader, "mc_phi"), recoPt(reader, "reco_pt"), recoPtCorrected(reader, "reco_pt_corrected"),
        recoEta(reader, "reco_eta"), recoDetEta(reader, "reco_det_eta"), recoPhi(reader, "reco_phi"),
        recoNeutralFraction(reader, "reco_neutral_fraction"), recoPtLead(reader, "reco_ptLead"),
        recoArea(reader, "reco_jet_area"), recoBgDensity(reader, "reco_bg_density"), thrownVz(reader, "event_vz"),
        recoVz(reader, "reco_vz"), pthat(reader, "pthat"), pthatMid(reader, "pthat_mid"),
        matchJp0(reader, "reco_trigger_match_JP0"), matchJp1(reader, "reco_trigger_match_JP1"),
        matchJp2(reader, "reco_trigger_match_JP2"), shouldJp0(reader, "evt_should_JP0"),
        shouldJp1(reader, "evt_should_JP1"), shouldJp2(reader, "evt_should_JP2"), eventId(reader, "eventid")
   {
      auto has = [tree](const char *name) { return tree->GetBranch(name) != nullptr; };
      if (has("mc_weight")) mcWeight.reset(new TTreeReaderValue<double>(reader, "mc_weight"));
      if (has("reco_trigger_match_HT2")) matchHt2.reset(new TTreeReaderValue<bool>(reader, "reco_trigger_match_HT2"));
      if (has("evt_should_HT2")) shouldHt2.reset(new TTreeReaderValue<bool>(reader, "evt_should_HT2"));
      if (has("reco_jp_patch_adc")) {
         patchAdc.reset(new TTreeReaderValue<int>(reader, "reco_jp_patch_adc"));
         towerAdc.reset(new TTreeReaderValue<int>(reader, "reco_ht_adc_max"));
         patchAdcMax.reset(new TTreeReaderValue<int>(reader, "evt_jp_adc_max"));
         towerAdcMax.reset(new TTreeReaderValue<int>(reader, "evt_ht_adc_max"));
         patchThreshold.reset(new TTreeReaderArray<int>(reader, "jp_thr"));
         towerThreshold.reset(new TTreeReaderArray<int>(reader, "ht_thr"));
      } else if (needAdc) {
         throw std::runtime_error("threshold variants need the simulator ADC branches (rebuild the matched trees)");
      }
      if (!needFragmentation) return;
      mcNeutralFraction.reset(new TTreeReaderValue<double>(reader, "mc_neutral_fraction"));
      mcPtLead.reset(new TTreeReaderValue<double>(reader, "mc_ptLead"));
      mcBgDensity.reset(new TTreeReaderValue<double>(reader, "mc_bg_density"));
      mcNConstituents.reset(new TTreeReaderValue<int>(reader, "mc_n_constituents"));
      recoNConstituents.reset(new TTreeReaderValue<int>(reader, "reco_n_constituents"));
      recoMultiplicity.reset(new TTreeReaderValue<double>(reader, "reco_multiplicity"));
   }

   // step to the next row; newEvent says whether it opens a new event (the eventid of the Stage-1 data
   // pass collides across files, so the thrown vertex z is part of the event key)
   bool Next(const Event &current, bool &newEvent)
   {
      if (!reader.Next()) return false;
      newEvent = *eventId != current.id || std::fabs(*thrownVz - current.thrownVz) > 1e-6;
      return true;
   }

   // start a new event on the row just read
   void BeginEvent(Event &event)
   {
      event = Event();
      event.id = *eventId;
      event.thrownVz = *thrownVz;
   }

   // add the jets of the current row to the event
   void AddRow(Event &event)
   {
      event.pthatMid = *pthatMid;
      if (*pthat > 0) event.pthat = *pthat;
      if (*recoPt > 0) {
         event.hasReco = true;
         event.recoVz = *recoVz;
         event.shouldJetPatch[0] = *shouldJp0;
         event.shouldJetPatch[1] = *shouldJp1;
         event.shouldJetPatch[2] = *shouldJp2;
         event.shouldHighTower = shouldHt2 ? **shouldHt2 : false;
         if (recoMultiplicity) event.multiplicity = **recoMultiplicity;
         if (patchAdcMax) {
            event.patchAdcMax = **patchAdcMax;
            event.towerAdcMax = **towerAdcMax;
            for (int patch = 0; patch < 3; ++patch) event.patchThreshold[patch] = (*patchThreshold)[patch];
            event.towerThreshold = (*towerThreshold)[2]; // the HT2 entry of the high-tower threshold array
         }
         DetJet jet;
         jet.pt = *recoPt;
         jet.ptCorrected = *recoPtCorrected;
         jet.eta = *recoEta;
         jet.detEta = *recoDetEta;
         jet.phi = *recoPhi;
         jet.neutralFraction = *recoNeutralFraction;
         jet.ptLead = *recoPtLead;
         jet.area = *recoArea;
         jet.bgDensity = *recoBgDensity;
         if (recoNConstituents) jet.nConstituents = (double)**recoNConstituents;
         jet.matchJetPatch[0] = *matchJp0;
         jet.matchJetPatch[1] = *matchJp1;
         jet.matchJetPatch[2] = *matchJp2;
         jet.matchHighTower = matchHt2 ? **matchHt2 : false;
         jet.patchAdc = patchAdc ? **patchAdc : -1;
         jet.towerAdc = towerAdc ? **towerAdc : -1;
         event.detJets.push_back(jet);
      }
      if (*mcPt > 0) {
         ParJet jet;
         jet.pt = *mcPt;
         jet.ptCorrected = *mcPtCorrected;
         jet.eta = *mcEta;
         jet.phi = *mcPhi;
         if (mcNeutralFraction) jet.neutralFraction = **mcNeutralFraction;
         if (mcPtLead) jet.ptLead = **mcPtLead;
         if (mcBgDensity) jet.bgDensity = **mcBgDensity;
         if (mcNConstituents) jet.nConstituents = (double)**mcNConstituents;
         event.parJets.push_back(jet);
      }
   }
};

// the event is usable on the detector side when it has a reconstructed vertex inside the data window
// and within 5 cm of the thrown one
static const double kRecoVzMax = 60.0;
static const double kRecoVzToThrownMax = 5.0;

inline bool HasGoodRecoVertex(const Event &event)
{
   if (!event.hasReco) return false;
   if (std::fabs(event.recoVz) >= kRecoVzMax) return false;
   return std::fabs(event.recoVz - event.thrownVz) < kRecoVzToThrownMax;
}

// ---- detector-to-particle jet pairing ---------------------------------------------------------------

static const double kMatchDeltaR = 0.2; // a detector and a particle jet are the same jet below this dR

// one matched detector-particle jet pair, by index into Event::detJets / Event::parJets
struct JetPair {
   int detJet; // index into Event::detJets
   int parJet; // index into Event::parJets
   double dr;  // separation of the two jets in (eta, phi)
};

// One-to-one pairing by global minimum dR: all candidate pairs are sorted by dR and taken greedily, a
// jet being used at most once, until the separation reaches kMatchDeltaR. detIdx and parIdx select the
// jets that take part; the returned indices point back into the event.
inline std::vector<JetPair> MatchJets(const Event &event, const std::vector<int> &detIdx, const std::vector<int> &parIdx)
{
   std::vector<JetPair> candidates;
   candidates.reserve(detIdx.size() * parIdx.size());
   for (size_t d = 0; d < detIdx.size(); ++d) {
      for (size_t p = 0; p < parIdx.size(); ++p) {
         const DetJet &detJet = event.detJets[detIdx[d]];
         const ParJet &parJet = event.parJets[parIdx[p]];
         const double dEta = detJet.eta - parJet.eta;
         const double dPhi = TVector2::Phi_mpi_pi(detJet.phi - parJet.phi);
         candidates.push_back({(int)d, (int)p, std::sqrt(dEta * dEta + dPhi * dPhi)});
      }
   }
   std::sort(candidates.begin(), candidates.end(),
             [](const JetPair &a, const JetPair &b) { return a.dr < b.dr; });
   std::vector<bool> detUsed(detIdx.size(), false);
   std::vector<bool> parUsed(parIdx.size(), false);
   std::vector<JetPair> pairs;
   for (const JetPair &candidate : candidates) {
      if (detUsed[candidate.detJet] || parUsed[candidate.parJet]) continue;
      if (candidate.dr >= kMatchDeltaR) break;
      detUsed[candidate.detJet] = true;
      parUsed[candidate.parJet] = true;
      pairs.push_back({detIdx[candidate.detJet], parIdx[candidate.parJet], candidate.dr});
   }
   return pairs;
}

// ---- small helpers ----------------------------------------------------------------------------------

// split a string on a separator: "jp1_e23_tiltUp" -> {"jp1", "e23", "tiltUp"}
inline std::vector<std::string> SplitWords(const std::string &text, char separator = ' ')
{
   std::vector<std::string> words;
   for (size_t begin = 0; begin < text.size();) {
      const size_t sep = text.find(separator, begin);
      words.push_back(text.substr(begin, sep == std::string::npos ? std::string::npos : sep - begin));
      if (sep == std::string::npos) break;
      begin = sep + 1;
   }
   return words;
}

// ratio of two histograms bin by bin, the relative errors of the two added in quadrature
inline TH1D *RatioWithErrors(const TH1D *numerator, const TH1D *denominator, const char *name)
{
   TH1D *ratio = (TH1D *)numerator->Clone(name);
   ratio->SetDirectory(0);
   for (int bin = 1; bin <= numerator->GetNbinsX(); ++bin) {
      const double num = numerator->GetBinContent(bin);
      const double den = denominator->GetBinContent(bin);
      const double numError = numerator->GetBinError(bin);
      const double denError = denominator->GetBinError(bin);
      ratio->SetBinContent(bin, den > 0 ? num / den : 0);
      const bool defined = num > 0 && den > 0;
      ratio->SetBinError(bin, defined ? num / den * std::sqrt(numError * numError / (num * num) +
                                                              denError * denError / (den * den))
                                      : 0);
   }
   return ratio;
}

} // namespace Stage2

#endif // STAGE2_COMMON_H
