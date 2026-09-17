// embedding_particles.C — composition of the generated particles of an embedding pico: the fraction of
// particles and of generated pT carried by each PDG species inside |eta| < 1 above the 0.2 GeV
// constituent cut, plus the share of the generated pT that falls BELOW that cut (in total, and for
// photons and pi0 separately).
//
// This is how the two embedding requests are told apart: the 2021 sample has pi0/eta/Sigma0 already
// decayed in the stored generator record, the 2023 one keeps them undecayed, which moves pT between the
// pi0 and the photon rows and changes the soft fraction below the constituent cut.
// The maker stores the PDG code of a generated track in its dEdx field.
//
// Input : one embedding pico (JetTreeMc). Read ONE pico per process — the pico TClonesArrays are shared
//         statics, and a second file clobbers the first.
// Output: the table on stdout.
// Run   : root -l -b -q -e 'gSystem->Load("/usr/local/eventStructuredAu/libTStarJetPico.so");
//         gSystem->AddIncludePath("-I/usr/local/eventStructuredAu");' 'embedding_particles.C+("<pico>")'
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoPrimaryTrack.h"
#include <cstdio>
#include <map>
#include <cmath>

// Acceptance of the composition table.
static const double kEtaMax = 1.0;
// The Stage-1 jet-finding constituent cut; particles below it never reach a jet.
static const double kPtConstituentMin = 0.2;

static const int kPdgPhoton = 22;
static const int kPdgPi0 = 111;

// The species the table lists, in the order they are printed.
static const int kSpeciesId[] = {211, 321, 2212, 22, 111, 221, 130, 310, 3122, 3212, 2112, 11, 13, 3112, 3222, 3312, 3322};
static const char *kSpeciesName[] = {"pi+-", "K+-", "p",      "gamma",  "pi0",    "eta",   "K0L", "K0S", "Lambda",
                                     "Sigma0", "n", "e",      "mu",     "Sigma-", "Sigma+", "Xi-", "Xi0"};

void embedding_particles(const char *file, long long nmax = 30000)
{
   TFile pico(file);
   TTree *tree = (TTree *)pico.Get("JetTreeMc");
   TStarJetPicoEvent *event = nullptr;
   tree->SetBranchAddress("PicoJetTree", &event);

   // counts and pT sums per PDG species, above the constituent cut
   std::map<int, double> countById, ptById;
   double nAbove = 0, ptAbove = 0, nCharged = 0;
   // generated pT in total and below the constituent cut
   double ptAll = 0, ptBelow = 0, ptBelowPhoton = 0, ptBelowPi0 = 0;
   long long nEvents = 0;

   // 1. sum over the generated tracks of every event
   for (long long i = 0; i < tree->GetEntries() && i < nmax; ++i) {
      tree->GetEntry(i);
      ++nEvents;
      for (int k = 0; k < event->GetPrimaryTracks()->GetEntries(); ++k) {
         auto *track = (TStarJetPicoPrimaryTrack *)event->GetPrimaryTracks()->At(k);
         if (std::fabs(track->GetEta()) > kEtaMax) continue;
         const int pdg = std::abs(TMath::Nint(track->GetdEdx()));
         const double pt = track->GetPt();
         ptAll += pt;
         if (pt < kPtConstituentMin) {
            ptBelow += pt;
            if (pdg == kPdgPhoton) ptBelowPhoton += pt;
            if (pdg == kPdgPi0) ptBelowPi0 += pt;
            continue;
         }
         countById[pdg] += 1;
         ptById[pdg] += pt;
         nAbove += 1;
         ptAbove += pt;
         if (std::fabs(track->GetCharge()) > 0) nCharged += 1;
      }
   }

   // 2. the table
   printf("%s\n  events %lld, particles per event (|eta|<1, pT>0.2) %.2f, charged fraction %.3f\n", file, nEvents,
          nAbove / nEvents, nCharged / nAbove);
   printf("  pT below 0.2 GeV: all %.4f, photons %.4f, pi0 %.4f of the generated pT\n", ptBelow / ptAll,
          ptBelowPhoton / ptAll, ptBelowPi0 / ptAll);
   printf("  %-8s %9s %9s\n", "species", "N frac", "pT frac");
   for (size_t j = 0; j < sizeof(kSpeciesId) / sizeof(int); ++j) {
      const int pdg = kSpeciesId[j];
      if (!countById.count(pdg)) continue;
      printf("  %-8s %9.4f %9.4f\n", kSpeciesName[j], countById[pdg] / nAbove, ptById[pdg] / ptAbove);
   }
}
