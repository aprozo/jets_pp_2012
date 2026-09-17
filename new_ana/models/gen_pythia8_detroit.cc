// gen_pythia8_detroit.cc — Pythia 8 with the RHIC Detroit tune, in one hard-process pT-hat bin.
//
// Every physics setting is the list kSettings below, in Pythia 8's own "Key = value" language:
// the beams, the hard process, the Detroit tune on top of Monash 2013 (Aguilar et al., PRD 105
// (2022) 016011, Table III: the underlying-event parameters fitted to RHIC data), the 13 particles
// left undecayed (the list of the embedding request), the print level. Edit the list and rerun
// `run_models.sh gen pythia8detroit`; the driver recompiles what changed. gen_pythia8_monash.cc is
// the same program without the tune block.
//
// The pT-hat bin, the number of events and the seed come from the command line.
// Output: the particle tree of tree.h; the cross section of the bin is Info::sigmaGen after the
//         last event.
// Usage : gen_pythia8_detroit <pthat_min> <pthat_max|-1> <nevents> <seed> <out.root>
#include "tree.h"
#include "Pythia8/Pythia.h"
#include <cstdio>
#include <cstdlib>
#include <string>

// ---- settings -------------------------------------------------------------------------------------
static const char *kTuneName = "Detroit tune (RHIC underlying-event tune on Monash 2013)";
static const char *kSettings[] = {
   // beams: protons at the energy of the 2012 RHIC run
   "Beams:idA = 2212",
   "Beams:idB = 2212",
   "Beams:eCM = 200.",
   // hard process: inclusive QCD 2 -> 2
   "HardQCD:all = on",
   // the Detroit tune, PRD 105 (2022) 016011 Table III: multiparton interactions and colour reconnection
   "MultipartonInteractions:ecmRef = 200.",     // reference energy of the MPI cut-off, GeV
   "MultipartonInteractions:pT0Ref = 1.40",     // MPI cut-off at the reference energy, GeV
   "MultipartonInteractions:ecmPow = 0.135",    // energy dependence of the cut-off
   "MultipartonInteractions:bProfile = 2",      // double-Gaussian overlap of the two protons
   "MultipartonInteractions:coreRadius = 0.56", // radius of the dense core
   "MultipartonInteractions:coreFraction = 0.78", // fraction of matter in the core
   "ColourReconnection:range = 5.4",            // colour-reconnection strength
   // the 13 species left undecayed, as in the embedding request (remove a line to let Pythia decay it)
   "211:mayDecay = off",  // pi+-
   "321:mayDecay = off",  // K+-
   "310:mayDecay = off",  // K0S
   "130:mayDecay = off",  // K0L
   "3122:mayDecay = off", // Lambda
   "3112:mayDecay = off", // Sigma-
   "3222:mayDecay = off", // Sigma+
   "3312:mayDecay = off", // Xi-
   "3322:mayDecay = off", // Xi0
   "3334:mayDecay = off", // Omega-
   "111:mayDecay = off",  // pi0
   "221:mayDecay = off",  // eta
   "3212:mayDecay = off", // Sigma0
   // printout: the initialisation and the final cross-section table only
   "Next:numberCount = 0",
   "Next:numberShowEvent = 0",
   "Next:numberShowInfo = 0",
   "Next:numberShowProcess = 0",
};
// ---------------------------------------------------------------------------------------------------

// mb -> pb, the unit the particle trees carry the cross section in.
static const double kMbToPb = 1e9;

int main(int argc, char **argv)
{
   if (argc != 6) {
      fprintf(stderr, "usage: gen_pythia8_detroit <pthat_min> <pthat_max|-1> <nevents> <seed> <out.root>\n");
      return 1;
   }
   const double lo = atof(argv[1]);
   const double hi = atof(argv[2]);
   const long nev = atol(argv[3]);
   const int seed = atoi(argv[4]);

   // 1. the settings above, then the bin and the seed of this job
   Pythia8::Pythia pythia;
   for (const char *setting : kSettings) pythia.readString(setting);
   pythia.readString(Form("PhaseSpace:pTHatMin = %g", lo));
   if (hi > 0) pythia.readString(Form("PhaseSpace:pTHatMax = %g", hi));
   pythia.readString("Random:setSeed = on");
   pythia.readString(Form("Random:seed = %d", seed));
   if (!pythia.init()) return 2;

   // 2. generate and store the final state
   Models::Writer particleTree(
      argv[5], Form("Pythia 8.313 %s, HardQCD:all, pTHat %g-%g, pp 200 GeV, 13 particles undecayed, seed %d",
                    kTuneName, lo, hi, seed));
   for (long i = 0; i < nev; ++i) {
      if (!pythia.next()) continue;
      const Pythia8::Event &event = pythia.event;
      for (int j = 0; j < event.size(); ++j) {
         if (event[j].isFinal())
            particleTree.ev.Add(event[j].px(), event[j].py(), event[j].pz(), event[j].e(), event[j].id());
      }
      particleTree.ev.pthat = pythia.info.pTHat();
      particleTree.Fill();
   }

   // 3. the cross section of the bin is only final after the last event
   pythia.stat();
   const double sigma_pb = pythia.info.sigmaGen() * kMbToPb;
   const double sigma_err_pb = pythia.info.sigmaErr() * kMbToPb;
   particleTree.Close(sigma_pb, sigma_err_pb);
   printf("gen_pythia8_detroit: %ld events, pthat %g-%g GeV, sigma = %g +- %g pb -> %s\n", nev, lo, hi, sigma_pb,
          sigma_err_pb, argv[5]);
   return 0;
}
