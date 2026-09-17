// embedding_diff.C — the two embeddings compared at fixed particle pT: jet energy scale, resolution,
// efficiency, fragmentation at particle and detector level, vertex and background density. The sample
// weights, the event grouping and the jet pairing are the ones the response build uses (stage2/common.h),
// so a difference seen here is a difference of the samples and not of the treatment.
//
// which = "e21": the per-run matched trees of the 2021 request, weighted sigma x corr / N_generated;
// which = "e23": the merged 2023 tree, weighted by its per-row mc_weight.
// Output <results>/embedding_diff_<e21|e23>.root; figures: embedding_diff_plot.C.
// root -l -b -q 'embedding_diff.C+("e21", "<list of the per-run 2021 matched trees>")'
// root -l -b -q 'embedding_diff.C+("e23", "<output>/merged_matching_R0.5.root")'
#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "../soft_reweight.h"
#include "../stage2/common.h"
#include "../vertex_reweight.h"

using namespace CrossSectionConfig;
using namespace Stage2;

// the comparison runs on the published bins up to 52 GeV, in particle pT
static const std::vector<double> kEdges = {5.0,  6.9,  8.2,  9.7,  11.5, 13.6, 16.1,
                                           19.0, 22.5, 26.6, 31.4, 37.2, 44.0, 52.0};
static const double kComparisonPtMin = 5.0;  // lowest jet pT entering the comparison [GeV]
static const double kComparisonPtMax = 60.0; // highest jet pT entering the comparison [GeV]

// the two jet selections compared side by side
static const int kNSelections = 2;
static const int kSelJetPatch = 0; // the response-like selection: the exclusive jet-patch levels
static const int kSelAll = 1;      // every jet of every event with a vertex, no trigger
static const char *kSelectionName[kNSelections] = {"jp", "all"};

// every distribution of one embedding, per particle-pT bin and |eta| < 0.5 (physics eta)
struct Plots {
   TProfile *scale[kNSelections];     // <pT_det / pT_part>
   TProfile *scaleSquared[kNSelections]; // its second moment, for the resolution
   TProfile *detNeutralFraction[kNSelections];
   TProfile *detNConstituents[kNSelections];
   TProfile *detLeadingFraction[kNSelections];
   TProfile *detBgDensity[kNSelections];
   TH1D *nMatched[kNSelections];      // particle jets matched with 5 <= pT_det < 60
   TH1D *nParticleJets;               // all particle jets
   TProfile *parNeutralFraction;      // particle-level fragmentation, all particle jets
   TProfile *parNConstituents;
   TProfile *parLeadingFraction;
   TProfile *parBgDensity;
   TH1D *thrownVz;                    // event level: every generated event
   TH1D *recoVz;                      // event level: events with a vertex inside the gate
   TH1D *multiplicity;
   TH1D *eventBgDensity;
   TH1D *pthat;
   double nEvents = 0;             // weighted generated events
   double nEventsWithVertex = 0;   // weighted events with a vertex inside the gate

   Plots(const char *tag)
   {
      const int nBins = (int)kEdges.size() - 1;
      for (int s = 0; s < kNSelections; ++s) {
         const char *sel = kSelectionName[s];
         scale[s] = new TProfile(Form("r_%s_%s", sel, tag), "", nBins, kEdges.data());
         scaleSquared[s] = new TProfile(Form("r2_%s_%s", sel, tag), "", nBins, kEdges.data());
         detNeutralFraction[s] = new TProfile(Form("nefD_%s_%s", sel, tag), "", nBins, kEdges.data());
         detNConstituents[s] = new TProfile(Form("nconD_%s_%s", sel, tag), "", nBins, kEdges.data());
         detLeadingFraction[s] = new TProfile(Form("leadD_%s_%s", sel, tag), "", nBins, kEdges.data());
         detBgDensity[s] = new TProfile(Form("bgD_%s_%s", sel, tag), "", nBins, kEdges.data());
         nMatched[s] = new TH1D(Form("nM_%s_%s", sel, tag), "", nBins, kEdges.data());
      }
      nParticleJets = new TH1D(Form("nP_%s", tag), "", nBins, kEdges.data());
      parNeutralFraction = new TProfile(Form("nefP_%s", tag), "", nBins, kEdges.data());
      parNConstituents = new TProfile(Form("nconP_%s", tag), "", nBins, kEdges.data());
      parLeadingFraction = new TProfile(Form("leadP_%s", tag), "", nBins, kEdges.data());
      parBgDensity = new TProfile(Form("bgP_%s", tag), "", nBins, kEdges.data());
      thrownVz = new TH1D(Form("vzThrown_%s", tag), "", 60, -150, 150);
      recoVz = new TH1D(Form("vzReco_%s", tag), "", 60, -150, 150);
      multiplicity = new TH1D(Form("mult_%s", tag), "", 60, 0, 60);
      eventBgDensity = new TH1D(Form("bgEv_%s", tag), "", 60, 0, 3);
      pthat = new TH1D(Form("pthat_%s", tag), "", 70, 0, 70);
   }

   void Write()
   {
      for (int s = 0; s < kNSelections; ++s) {
         scale[s]->Write();
         scaleSquared[s]->Write();
         detNeutralFraction[s]->Write();
         detNConstituents[s]->Write();
         detLeadingFraction[s]->Write();
         detBgDensity[s]->Write();
         nMatched[s]->Write();
      }
      nParticleJets->Write();
      parNeutralFraction->Write();
      parNConstituents->Write();
      parLeadingFraction->Write();
      parBgDensity->Write();
      thrownVz->Write();
      recoVz->Write();
      multiplicity->Write();
      eventBgDensity->Write();
      pthat->Write();
   }
};

// the matched pairs of one event under one selection: the detector quantities plotted against the
// particle pT of the jet they belong to
static void FillPairs(const Event &event, Plots &plots, int selection, const std::vector<int> &parIdx,
                      double weight)
{
   // the exclusive level of the event is its highest simulated jet-patch decision
   const int level = event.shouldJetPatch[2] ? 2 : event.shouldJetPatch[1] ? 1 : event.shouldJetPatch[0] ? 0 : -1;
   if (selection == kSelJetPatch && level < 0) return;
   const double rawPtFloor = selection == kSelJetPatch ? RawPtFloorOfPatch(level) : 0.0;

   std::vector<int> detIdx;
   for (size_t i = 0; i < event.detJets.size(); ++i) {
      const DetJet &jet = event.detJets[i];
      if (selection == kSelJetPatch && !jet.matchJetPatch[level]) continue;
      if (jet.neutralFraction > kNeutralFractionMax) continue;
      if (std::fabs(jet.detEta) >= kEtaMax || std::fabs(jet.eta) >= kEtaMax || jet.pt < rawPtFloor) continue;
      detIdx.push_back((int)i);
   }

   for (const JetPair &pair : MatchJets(event, detIdx, parIdx)) {
      const DetJet &detJet = event.detJets[pair.detJet];
      const ParJet &parJet = event.parJets[pair.parJet];
      const bool inRange = parJet.ptCorrected >= kComparisonPtMin && parJet.ptCorrected < kComparisonPtMax &&
                           detJet.ptCorrected >= kComparisonPtMin && detJet.ptCorrected < kComparisonPtMax;
      if (!inRange || std::fabs(parJet.eta) >= kEtaBlockSplit) continue;
      const double parPt = parJet.ptCorrected;
      const double ratio = detJet.ptCorrected / parJet.ptCorrected;
      plots.scale[selection]->Fill(parPt, ratio, weight);
      plots.scaleSquared[selection]->Fill(parPt, ratio * ratio, weight);
      plots.detNeutralFraction[selection]->Fill(parPt, detJet.neutralFraction, weight);
      plots.detNConstituents[selection]->Fill(parPt, detJet.nConstituents, weight);
      plots.detLeadingFraction[selection]->Fill(parPt, detJet.pt > 0 ? detJet.ptLead / detJet.pt : 0, weight);
      plots.detBgDensity[selection]->Fill(parPt, detJet.bgDensity, weight);
      plots.nMatched[selection]->Fill(parPt, weight);
   }
}

// one generated event; sampleWeight is the pb per generated event of its pt-hat sample
static void ProcessEvent(const Event &event, Plots &plots, double sampleWeight)
{
   const double pthat = event.pthat > 0 ? event.pthat : event.pthatMid;
   const double weight = sampleWeight * SoftReweight::weight(pthat);

   // 1. the particle level: every generated event, whether it was reconstructed or not
   plots.nEvents += weight;
   plots.thrownVz->Fill(event.thrownVz, weight);
   plots.pthat->Fill(pthat, weight);
   std::vector<int> parIdx;
   for (size_t i = 0; i < event.parJets.size(); ++i) {
      const ParJet &jet = event.parJets[i];
      if (std::fabs(jet.eta) >= kEtaMax) continue;
      parIdx.push_back((int)i);
      if (std::fabs(jet.eta) >= kEtaBlockSplit || jet.ptCorrected < kComparisonPtMin) continue;
      plots.nParticleJets->Fill(jet.ptCorrected, weight);
      plots.parNeutralFraction->Fill(jet.ptCorrected, jet.neutralFraction, weight);
      plots.parNConstituents->Fill(jet.ptCorrected, jet.nConstituents, weight);
      plots.parLeadingFraction->Fill(jet.ptCorrected, jet.pt > 0 ? jet.ptLead / jet.pt : 0, weight);
      plots.parBgDensity->Fill(jet.ptCorrected, jet.bgDensity, weight);
   }

   // 2. the detector level: only events with a reconstructed vertex close to the thrown one
   if (!HasGoodRecoVertex(event)) return;
   const double detWeight = weight * VertexReweight::weight(event.recoVz);
   plots.nEventsWithVertex += detWeight;
   plots.recoVz->Fill(event.recoVz, detWeight);
   plots.multiplicity->Fill(event.multiplicity, detWeight);
   if (!event.detJets.empty()) plots.eventBgDensity->Fill(event.detJets[0].bgDensity, detWeight);
   FillPairs(event, plots, kSelJetPatch, parIdx, detWeight);
   FillPairs(event, plots, kSelAll, parIdx, detWeight);
}

// ---- the two inputs ----------------------------------------------------------------------------------

// the per-run 2021 trees: one pass per pt-hat sample, weighted sigma x corr / N_generated over the runs
// of the list
static long long ReadPerRunTrees(const char *listfile, Plots &plots)
{
   AnalysisConfig cfg;
   const std::map<std::string, std::map<int, long long>> nGeneratedPerRun = GeneratedEventsPerRun(cfg);
   std::map<std::string, std::vector<std::string>> files;
   std::map<std::string, long long> nGeneratedOfSample;
   {
      std::ifstream list(listfile);
      std::string line;
      while (std::getline(list, line)) {
         std::string sample;
         int run = 0;
         if (line.size() < 5 || !ParseMatchedFileName(line, sample, run)) continue;
         const auto ofSample = nGeneratedPerRun.find(sample);
         if (ofSample == nGeneratedPerRun.end()) continue;
         const auto ofRun = ofSample->second.find(run);
         if (ofRun == ofSample->second.end()) continue;
         files[sample].push_back(line);
         nGeneratedOfSample[sample] += ofRun->second;
      }
   }
   long long nEvents = 0;
   for (const auto &entry : files) {
      const std::string &sample = entry.first;
      const double sampleWeight = SampleWeightFactor(sample, nGeneratedOfSample[sample]);
      for (const std::string &path : entry.second) {
         TFile *file = TFile::Open(path.c_str());
         if (!file || file->IsZombie()) continue;
         MatchedTreeRows rows((TTree *)file->Get("MatchedTree"), false, true);
         Event event;
         bool newEvent = false;
         while (rows.Next(event, newEvent)) {
            if (newEvent) {
               if (event.id >= 0) {
                  ProcessEvent(event, plots, sampleWeight);
                  ++nEvents;
               }
               rows.BeginEvent(event);
            }
            rows.AddRow(event);
         }
         if (event.id >= 0) {
            ProcessEvent(event, plots, sampleWeight);
            ++nEvents;
         }
         file->Close();
      }
      printf("[diag] %-8s %zu files done\n", sample.c_str(), entry.second.size());
      fflush(stdout);
   }
   return nEvents;
}

// the merged 2023 tree: the samples are interleaved and each row carries its own mc_weight
static long long ReadMergedTree(const char *input, Plots &plots)
{
   TFile *file = TFile::Open(input);
   MatchedTreeRows rows((TTree *)file->Get("MatchedTree"), false, true);
   Event event;
   double sampleWeight = 0;
   long long nEvents = 0;
   long long nRows = 0;
   bool newEvent = false;
   while (rows.Next(event, newEvent)) {
      if (newEvent) {
         if (event.id >= 0) {
            ProcessEvent(event, plots, sampleWeight);
            ++nEvents;
         }
         rows.BeginEvent(event);
         sampleWeight = MergedWeightFactor(**rows.mcWeight, *rows.pthatMid);
      }
      rows.AddRow(event);
      if (++nRows % 20000000 == 0) {
         printf("[diag] %lld rows\n", nRows);
         fflush(stdout);
      }
   }
   if (event.id >= 0) {
      ProcessEvent(event, plots, sampleWeight);
      ++nEvents;
   }
   return nEvents;
}

void embedding_diff(const char *which, const char *input)
{
   const std::string output = RadiusDir("0.5") + "embedding_diff_" + which + ".root";
   Plots plots(which);
   const long long nEvents =
      std::string(which) == "e21" ? ReadPerRunTrees(input, plots) : ReadMergedTree(input, plots);

   TFile out(output.c_str(), "RECREATE");
   plots.Write();
   TH1D meta("meta", "", 2, 0, 2); // the weighted event counts the figures normalise by
   meta.SetBinContent(1, plots.nEvents);
   meta.SetBinContent(2, plots.nEventsWithVertex);
   meta.Write();
   out.Close();
   printf("[diag] %s: %lld events, weighted %.4g, with vertex in gate %.4g (%.3f) -> %s\n", which, nEvents,
          plots.nEvents, plots.nEventsWithVertex, plots.nEvents > 0 ? plots.nEventsWithVertex / plots.nEvents : 0,
          output.c_str());
}
