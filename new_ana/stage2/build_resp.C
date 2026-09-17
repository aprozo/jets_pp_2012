// build_resp.C — per-level, eta-block response of the embedding from the Stage-1 matched trees, either
// one file per pt-hat sample and run (matched_e21_<sample>_<run>_<id>_R<R>.root, build_resp) or one
// merged tree with the samples interleaved (build_resp_merged).
//
// Per event (the rows of one event are consecutive; the jets are re-collected from the rows):
//   particle jets P: |eta| < 0.9, pT = UE-subtracted particle pT; filled for EVERY generated event
//                    (events without a reconstructed vertex included), weight  F_sample x w_soft
//   detector side  : only if the event has a reconstructed vertex with |vz| < 60 cm within 5 cm of the
//                    thrown one; level l = highest simulated trigger decision (should_JP2 > JP1 > JP0);
//                    jets D: jet-patch match of level l, |eta_det| < 0.9, |eta| < 0.9, R_T <= 0.95,
//                    raw pT >= floor(l) (0 / 6.0 / 8.4 GeV);  weight  F_sample x w_soft x w_vz(vz_reco);
//                    the min-bias level (mb, every event with a vertex, no patch match) takes instead
//                    the min-bias jet-quality cuts (common.h), the high-tower level (ht2: should_HT2,
//                    jets holding the firing tower) has no raw floor
//   pairs          : one-to-one D-P pairing by global minimum dR, dR < 0.2; pairs with
//                    0 <= pT_P < 60 and 5 <= pT_D < 60 fill the response, the rest stay as fake + miss.
// Eta blocks (physics eta): 00_05 (|eta| < 0.5), 05_09 (0.5 <= |eta| < 0.9), for D (rows) and P (columns).
// Runs: the per-run build is restricted to the runs kept in the data, so the response carries the data's
// mix of trigger thresholds; N_sample counts the generated events of the same runs.
// F_sample = sigma_sample [mb] x 1e9 x corr_sample / N_sample (lists/emb2021_events_per_run.txt), so every
// histogram is in pb.  w_soft: soft_reweight.h at the event pt-hat; w_vz: vertex_reweight.h.
//
// Systematic variants (variants.h) act on the detector side at fill time:
//   energy scale:     filled pT = pT_corr + s(R_T) x pT_raw, s = +-sqrt(((1-R_T) 0.011)^2 + (R_T 0.032)^2);
//   thresholds:       the jet-patch / high-tower matches replayed from the stored DSM ADCs at threshold +-1;
//   underlying event: filled pT = pT_corr + (1 - f) x rho x area, f = 0.86 / 1.18.
// The selection (floors, eta, R_T) stays on the nominal quantities.
//
// Output <results>/response_blocks_R<R>[<suffix>][_<variant>].root, all fine-binned (detector 550 x [5,60],
// particle 600 x [0,60]), once per pt-hat sample (tag "_<sample>") and once summed (no tag):
//   x_<pb>                 particle spectrum (all events)         [+ _entries: unit weights]
//   b_<level>_<db>         detector spectrum, all D jets of the level
//   A_<level>_<db>_<pb>    response (detector pT, particle pT) of the in-range pairs
//   md_<level>_<db>, mp_<level>_<pb>   matched detector / particle spectra of the in-range pairs
// root -l -b -q 'build_resp.C+("0.5", "<results>/matching_files.list")'
// root -l -b -q -e '.L build_resp.C+' -e 'build_resp_merged("0.5", "e23")'
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "../soft_reweight.h"
#include "../vertex_reweight.h"
#include "common.h"
#include "variants.h"

using namespace CrossSectionConfig;
using namespace Stage2;

// the variant of this build; nominal until build_resp / build_resp_merged parses its argument
static Syst::Variant gVariant;

// ---- the histogram set of one pt-hat sample ---------------------------------------------------------

// Every quantity comes twice: once weighted (the physics) and once with unit weights (the "_entries"
// copy), which is what the outlier filter and the replica loop of unfold.C count.
struct Hists {
   TH1D *parSpectrum[kNEtaBlocks];                  // x:  particle jets of every generated event
   TH1D *parEntries[kNEtaBlocks];
   TH1D *detSpectrum[kNLevels][kNEtaBlocks];        // b:  detector jets of the level
   TH1D *detEntries[kNLevels][kNEtaBlocks];
   TH1D *matchedDet[kNLevels][kNEtaBlocks];         // md: detector jets of the in-range pairs
   TH1D *matchedDetEntries[kNLevels][kNEtaBlocks];
   TH1D *matchedPar[kNLevels][kNEtaBlocks];         // mp: particle jets of the in-range pairs
   TH1D *matchedParEntries[kNLevels][kNEtaBlocks];
   TH2D *response[kNLevels][kNEtaBlocks][kNEtaBlocks];  // A:  (detector pT, particle pT) of the pairs
   TH2D *responseEntries[kNLevels][kNEtaBlocks][kNEtaBlocks];

   // tag: "" for the sum over the samples, "_<sample>" for one pt-hat sample
   Hists(const char *tagText)
   {
      // names via std::string: Form()'s circular buffer wraps within this constructor (98 names)
      // and would corrupt a tag that itself came from Form()
      const std::string tag = tagText;
      for (int j = 0; j < kNEtaBlocks; ++j) {
         const std::string block = kEtaBlockName[j];
         parSpectrum[j] = NewPar("x_" + block + tag);
         parEntries[j] = NewPar("x_" + block + "_entries" + tag);
      }
      for (int l = 0; l < kNLevels; ++l) {
         for (int i = 0; i < kNEtaBlocks; ++i) {
            const std::string levelBlock = std::string(kLevelName[l]) + "_" + kEtaBlockName[i];
            detSpectrum[l][i] = NewDet("b_" + levelBlock + tag);
            detEntries[l][i] = NewDet("b_" + levelBlock + "_entries" + tag);
            matchedDet[l][i] = NewDet("md_" + levelBlock + tag);
            matchedDetEntries[l][i] = NewDet("md_" + levelBlock + "_entries" + tag);
            matchedPar[l][i] = NewPar("mp_" + levelBlock + tag);
            matchedParEntries[l][i] = NewPar("mp_" + levelBlock + "_entries" + tag);
            for (int j = 0; j < kNEtaBlocks; ++j) {
               const std::string cell = levelBlock + "_" + kEtaBlockName[j];
               response[l][i][j] = NewResponse("A_" + cell + tag);
               responseEntries[l][i][j] = NewResponse("A_" + cell + "_entries" + tag);
            }
         }
      }
   }

   // apply an action to every histogram of the set, always in the same order
   template <class Action>
   void ForEach(Action act)
   {
      for (int j = 0; j < kNEtaBlocks; ++j) {
         act(parSpectrum[j]);
         act(parEntries[j]);
      }
      for (int l = 0; l < kNLevels; ++l) {
         for (int i = 0; i < kNEtaBlocks; ++i) {
            act(detSpectrum[l][i]);
            act(detEntries[l][i]);
            act(matchedDet[l][i]);
            act(matchedDetEntries[l][i]);
            act(matchedPar[l][i]);
            act(matchedParEntries[l][i]);
            for (int j = 0; j < kNEtaBlocks; ++j) {
               act(response[l][i][j]);
               act(responseEntries[l][i][j]);
            }
         }
      }
   }

private:
   static TH1D *NewDet(const std::string &name) { return new TH1D(name.c_str(), "", kNDetFine, kDetFineMin, kDetFineMax); }
   static TH1D *NewPar(const std::string &name) { return new TH1D(name.c_str(), "", kNParFine, kParFineMin, kParFineMax); }
   static TH2D *NewResponse(const std::string &name)
   {
      return new TH2D(name.c_str(), "", kNDetFine, kDetFineMin, kDetFineMax, kNParFine, kParFineMin, kParFineMax);
   }
};

// ---- the filled pT of a jet under the current variant -----------------------------------------------

// particle pT: the UE-subtracted one, or the raw one for the no-UE jet definition
static inline double ParPt(const ParJet &jet)
{
   return gVariant.ueParticle ? jet.ptCorrected : jet.pt;
}

// detector pT: the UE-subtracted one plus the energy-scale shift (a fraction of the raw pT) and the
// change of the UE subtraction the variant asks for
static inline double DetPt(const DetJet &jet)
{
   double pt = jet.ptCorrected;
   if (gVariant.jesSign != 0) pt += Syst::JesShift(jet.neutralFraction, gVariant.jesSign) * jet.pt;
   if (gVariant.ueResp != 1.0) pt += (1.0 - gVariant.ueResp) * jet.bgDensity * jet.area;
   return pt;
}

// ---- DSM threshold members --------------------------------------------------------------------------

// As published, the trigger LEVEL keeps the simulator decision at the nominal threshold and the shift
// acts on the jet-to-patch (jet-to-tower) match only. Should() replays a decision at a shifted threshold
// from the stored maximum ADC of the event (patches are stored down to JP0 - 2, towers down to HT3); it
// is kept for the ADC-scan cross-checks and is called with shift 0 by the members.
static inline bool Should(bool stored, int adcMax, int threshold, int shift)
{
   if (shift == 0) return stored;
   return adcMax > threshold + shift || (stored && !(adcMax > threshold));
}

static inline bool Match(bool stored, int adc, int threshold, int shift)
{
   return shift == 0 ? stored : adc > threshold + shift;
}

// ---- one event ---------------------------------------------------------------------------------------

// does the jet carry the trigger object and the quality the level asks for?
static bool PassesLevelQuality(const Event &event, const DetJet &jet, int jetPatch)
{
   if (jetPatch < 3)
      return Match(jet.matchJetPatch[jetPatch], jet.patchAdc, event.patchThreshold[jetPatch], gVariant.thrShift) &&
             jet.neutralFraction <= kNeutralFractionMax;
   if (jetPatch == kLevelHighTower)
      return Match(jet.matchHighTower, jet.towerAdc, event.towerThreshold, gVariant.thrShift) &&
             jet.neutralFraction <= kNeutralFractionMax;
   // min-bias: no trigger object, the junk-jet protection of build_data.C instead
   return jet.neutralFraction > kMbNeutralFractionMin && jet.neutralFraction < kMbNeutralFractionMax &&
          jet.ptLead < kMbLeadingFractionMax * jet.ptCorrected;
}

// detector side of one event for one level: the detector spectrum and, paired with the particle jets,
// the response. parIdx are the particle jets of the event that are inside the acceptance.
static void FillLevel(const Event &event, Hists &hists, int level, const std::vector<int> &parIdx, double weight)
{
   const int jetPatch = JetPatchOfLevel(level);
   const double rawPtFloor = RawPtFloorOfPatch(jetPatch);

   // 1. the detector jets of this level
   std::vector<int> detIdx;
   for (size_t i = 0; i < event.detJets.size(); ++i) {
      const DetJet &jet = event.detJets[i];
      if (!PassesLevelQuality(event, jet, jetPatch)) continue;
      if (std::fabs(jet.detEta) >= kEtaMax || std::fabs(jet.eta) >= kEtaMax || jet.pt < rawPtFloor) continue;
      detIdx.push_back((int)i);
      const int block = EtaBlockOf(jet.eta);
      hists.detSpectrum[level][block]->Fill(DetPt(jet), weight);
      hists.detEntries[level][block]->Fill(DetPt(jet));
   }

   // 2. the pairs: everything outside the response range stays a fake (detector) or a miss (particle)
   for (const JetPair &pair : MatchJets(event, detIdx, parIdx)) {
      const DetJet &detJet = event.detJets[pair.detJet];
      const ParJet &parJet = event.parJets[pair.parJet];
      const double detPt = DetPt(detJet);
      const double parPt = ParPt(parJet);
      const bool inRange = parPt >= kParFineMin && parPt < kParFineMax && detPt >= kDetFineMin && detPt < kDetFineMax;
      if (!inRange) continue;
      const int detBlock = EtaBlockOf(detJet.eta);
      const int parBlock = EtaBlockOf(parJet.eta);
      hists.response[level][detBlock][parBlock]->Fill(detPt, parPt, weight);
      hists.responseEntries[level][detBlock][parBlock]->Fill(detPt, parPt);
      hists.matchedDet[level][detBlock]->Fill(detPt, weight);
      hists.matchedDetEntries[level][detBlock]->Fill(detPt);
      hists.matchedPar[level][parBlock]->Fill(parPt, weight);
      hists.matchedParEntries[level][parBlock]->Fill(parPt);
   }
}

// one generated event into the histogram set; minBiasOnly skips the trigger levels (the min-bias pass
// of Stage-1 has no jet-patch information worth reading)
static void ProcessEvent(const Event &event, Hists &hists, bool minBiasOnly = false)
{
   const double pthat = event.pthat > 0 ? event.pthat : event.pthatMid;
   const double softWeight = SoftReweight::weight(pthat);

   // 1. the particle level: every generated event, whether it was reconstructed or not
   std::vector<int> parIdx;
   for (size_t i = 0; i < event.parJets.size(); ++i) {
      const ParJet &jet = event.parJets[i];
      if (std::fabs(jet.eta) >= kEtaMax) continue;
      parIdx.push_back((int)i);
      const int block = EtaBlockOf(jet.eta);
      hists.parSpectrum[block]->Fill(ParPt(jet), softWeight);
      hists.parEntries[block]->Fill(ParPt(jet));
   }

   // 2. the detector level: only events with a reconstructed vertex close to the thrown one
   if (!HasGoodRecoVertex(event)) return;
   const double detWeight = softWeight * VertexReweight::weight(event.recoVz);
   bool shouldPatch[3];
   for (int patch = 0; patch < 3; ++patch)
      shouldPatch[patch] = Should(event.shouldJetPatch[patch], event.patchAdcMax, event.patchThreshold[patch], 0);
   const bool shouldHighTower = Should(event.shouldHighTower, event.towerAdcMax, event.towerThreshold, 0);

   // 3. the exclusive level of the event is its highest simulated jet-patch decision
   const int exclusiveLevel = shouldPatch[2] ? 2 : shouldPatch[1] ? 1 : shouldPatch[0] ? 0 : -1;
   if (exclusiveLevel >= 0 && !minBiasOnly) FillLevel(event, hists, exclusiveLevel, parIdx, detWeight);
   FillLevel(event, hists, kLevelMinBias, parIdx, detWeight); // min-bias: no trigger requirement
   if (minBiasOnly) return;
   if (shouldPatch[1]) FillLevel(event, hists, 4, parIdx, detWeight); // inclusive JP1
   if (shouldPatch[2]) FillLevel(event, hists, 5, parIdx, detWeight); // inclusive JP2
   if (shouldHighTower) FillLevel(event, hists, kLevelHighTower, parIdx, detWeight);
}

// ---- accumulation ------------------------------------------------------------------------------------

// scale a sample to pb, add it to the running sum and write it out. The "_entries" copies keep their
// unit weights: they are counts, not a cross section.
static void ScaleAddAndWrite(Hists &sample, Hists &total, double sampleWeight, TFile &out)
{
   std::vector<TH1 *> ofSample;
   std::vector<TH1 *> ofTotal;
   sample.ForEach([&](TH1 *h) { ofSample.push_back(h); });
   total.ForEach([&](TH1 *h) { ofTotal.push_back(h); });
   out.cd();
   for (size_t i = 0; i < ofSample.size(); ++i) {
      if (!TString(ofSample[i]->GetName()).Contains("_entries")) ofSample[i]->Scale(sampleWeight);
      ofTotal[i]->Add(ofSample[i]);
      ofSample[i]->Write();
      delete ofSample[i];
   }
}

// ---- response from the per-run trees of one embedding pass -------------------------------------------

// the files of the list grouped by pt-hat sample, keeping only the runs the data keeps
static std::map<std::string, std::vector<std::string>> GroupFilesBySample(const char *listfile,
                                                                         const std::vector<int> &keptRuns)
{
   std::map<std::string, std::vector<std::string>> files;
   std::ifstream list(listfile);
   std::string line;
   std::string sample;
   int run = 0;
   int nSkipped = 0;
   while (std::getline(list, line)) {
      if (line.size() <= 4 || !ParseMatchedFileName(line, sample, run)) continue;
      if (std::find(keptRuns.begin(), keptRuns.end(), run) == keptRuns.end()) {
         ++nSkipped;
         continue;
      }
      files[sample].push_back(line);
   }
   printf("[resp] %d files of runs outside the kept set skipped\n", nSkipped);
   return files;
}

// read one pt-hat sample into its histogram set; returns the number of events that carried jets and
// reports the generated-event count of the runs actually read
static long long ReadSampleFiles(const std::vector<std::string> &files,
                                 const std::map<int, long long> &nGeneratedPerRun, Hists &hists,
                                 long long &nGenerated, int &nFilesRead, int &nFilesMissing)
{
   long long nEvents = 0;
   nGenerated = 0;
   nFilesRead = 0;
   nFilesMissing = 0;
   std::set<int> runsSeen; // a run split over two picos (memory) counts its generated events once
   for (const std::string &path : files) {
      std::string fileSample;
      int run = 0;
      ParseMatchedFileName(path, fileSample, run);
      auto generated = nGeneratedPerRun.find(run);
      if (generated == nGeneratedPerRun.end()) {
         ++nFilesMissing;
         continue;
      }
      if (runsSeen.insert(run).second) nGenerated += generated->second;
      TFile *file = TFile::Open(path.c_str());
      if (!file || file->IsZombie()) {
         ++nFilesMissing;
         continue;
      }
      MatchedTreeRows rows((TTree *)file->Get("MatchedTree"), gVariant.thrShift != 0);
      ++nFilesRead;
      Event event;
      bool newEvent = false;
      while (rows.Next(event, newEvent)) {
         if (newEvent) {
            if (event.id >= 0) {
               ProcessEvent(event, hists);
               ++nEvents;
            }
            rows.BeginEvent(event);
         }
         rows.AddRow(event);
      }
      if (event.id >= 0) {
         ProcessEvent(event, hists);
         ++nEvents;
      }
      file->Close();
   }
   return nEvents;
}

// suffix: "" for the nominal (jet-patch pass) response, "_mb2021m" for a list of min-bias-pass trees;
// variant: a name of variants.h ("" nominal), appended to the output name
void build_resp(const char *jetR = "0.5", const char *listfile = "matching_files.list", const char *suffix = "",
                const char *variant = "")
{
   TH1::SetDefaultSumw2();
   AnalysisConfig cfg;
   gVariant = Syst::Parse(variant);
   const std::string variantTag = gVariant.name.empty() ? "" : "_" + gVariant.name;

   // 1. the run set and the generated-event bookkeeping
   const std::vector<int> keptRuns = KeptRuns(cfg);
   printf("[resp] %zu kept data runs define the embedding run set; variant '%s'\n", keptRuns.size(), variant);
   std::map<std::string, std::map<int, long long>> nGeneratedPerRun = GeneratedEventsPerRun(cfg);
   const std::map<std::string, std::vector<std::string>> files = GroupFilesBySample(listfile, keptRuns);

   // 2. one pass per sample, each scaled to pb and added to the total
   Hists total("");
   const std::string output = RadiusDir(jetR) + "response_blocks_R" + jetR + suffix + variantTag + ".root";
   TFile out(output.c_str(), "RECREATE");
   long long nEventsAll = 0;
   for (const auto &entry : files) {
      const std::string &sample = entry.first;
      if (!IsKnownSample(sample)) {
         printf("[resp] unknown sample %s, skipped\n", sample.c_str());
         continue;
      }
      Hists hists(Form("_%s", sample.c_str()));
      long long nGenerated = 0;
      int nFilesRead = 0;
      int nFilesMissing = 0;
      const long long nEvents = ReadSampleFiles(entry.second, nGeneratedPerRun[sample], hists, nGenerated,
                                                nFilesRead, nFilesMissing);
      const double sampleWeight = SampleWeightFactor(sample, nGenerated);
      printf("[resp] %-8s %4d files (%d skipped)  events with jets %lld  N_generated %lld  F = %.4e pb/event\n",
             sample.c_str(), nFilesRead, nFilesMissing, nEvents, nGenerated, sampleWeight);
      nEventsAll += nEvents;
      ScaleAddAndWrite(hists, total, sampleWeight, out);
   }

   // 3. the sum over the samples
   out.cd();
   total.ForEach([&](TH1 *h) { h->Write(); });
   out.Close();
   printf("[resp] %lld events with jets in %zu samples -> %s\n", nEventsAll, files.size(), output.c_str());
}

// ---- response from a merged tree ---------------------------------------------------------------------

// Response from a MERGED 2023 matched tree: which = "e23" reads the jet-patch pass
// (merged_matching_R<R>, every level) and "mb2023" the min-bias pass (merged_matchingMB_R<R>, the
// leading-track sDCA veto of the min-bias data, min-bias level only). mc_weight = sigma_sample / N_sample
// [mb] per row; the samples are interleaved in the tree and separated by pthat_mid, which is what the
// outlier filter needs. No run restriction (the merged tree is chunked, not per run).
// Output <results>/response_blocks_R<R>_<which>[_<variant>].root
void build_resp_merged(const char *jetR = "0.5", const char *which = "e23", const char *variant = "")
{
   TH1::SetDefaultSumw2(false); // bin contents only (unfold.C counts on the _entries histograms)
   TH1::AddDirectory(false);
   const bool minBiasOnly = std::string(which) == "mb2023";
   AnalysisConfig cfg;
   gVariant = Syst::Parse(variant);
   const std::string variantTag = gVariant.name.empty() ? "" : "_" + gVariant.name;
   const std::string input = cfg.datapath + (minBiasOnly ? "merged_matchingMB_R" : "merged_matching_R") + jetR + ".root";
   const std::string output = RadiusDir(jetR) + "response_blocks_R" + jetR + "_" + which + variantTag + ".root";

   TFile *file = TFile::Open(input.c_str());
   if (!file || file->IsZombie()) {
      printf("[resp] cannot open %s\n", input.c_str());
      return;
   }
   // the min-bias level has no trigger, so its threshold members equal the nominal
   MatchedTreeRows rows((TTree *)file->Get("MatchedTree"), gVariant.thrShift != 0 && !minBiasOnly);
   if (!rows.mcWeight) throw std::runtime_error("merged tree without mc_weight");
   printf("[resp] %s, variant '%s'\n", input.c_str(), variant);

   // 1. one histogram set per pt-hat sample, the samples told apart by the pt-hat bin midpoint
   std::map<std::string, Hists *> perSample;
   std::map<std::string, double> weightOfSample;
   std::map<std::string, long long> nEventsOfSample;
   Event event;
   std::string sample;
   Hists *current = nullptr;
   long long nRows = 0;
   bool newEvent = false;
   while (rows.Next(event, newEvent)) {
      if (newEvent) {
         if (event.id >= 0) {
            ProcessEvent(event, *current, minBiasOnly);
            ++nEventsOfSample[sample];
         }
         rows.BeginEvent(event);
         sample = Form("pm%.1f", *rows.pthatMid);
         auto it = perSample.find(sample);
         if (it == perSample.end()) {
            weightOfSample[sample] = MergedWeightFactor(**rows.mcWeight, *rows.pthatMid);
            it = perSample.emplace(sample, new Hists(("_" + sample).c_str())).first;
         }
         current = it->second;
      }
      rows.AddRow(event);
      if (++nRows % 20000000 == 0) {
         printf("[resp] %lld rows\n", nRows);
         fflush(stdout);
      }
   }
   if (event.id >= 0) {
      ProcessEvent(event, *current, minBiasOnly);
      ++nEventsOfSample[sample];
   }
   file->Close();

   // 2. scale each sample to pb, add it to the total and write both out
   TFile out(output.c_str(), "RECREATE");
   Hists total("");
   for (auto &entry : perSample) {
      const double sampleWeight = weightOfSample[entry.first];
      printf("[resp] sample %-8s events %lld  F = %.4e pb/event\n", entry.first.c_str(),
             nEventsOfSample[entry.first], sampleWeight);
      ScaleAddAndWrite(*entry.second, total, sampleWeight, out);
      delete entry.second;
   }
   total.ForEach([&](TH1 *h) { h->Write(); });
   out.Close();
   printf("[resp] %zu samples -> %s\n", perSample.size(), output.c_str());
}
