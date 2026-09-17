// gen_pythia6.cc — Pythia 6.4 (the generator of the embedding) in one hard-process pT-hat bin.
//
// Every physics setting is in the block "settings" below: the tune of the embedding request
// (Perugia 2012, PYTUNE 370, CTEQ6L1 through LHAPDF), PARP(90) = 0.213, QCD 2->2 (MSEL 1), pp at
// 200 GeV, and the 13 particles the request leaves for GEANT to decay. Edit the block and rerun
// `run_models.sh gen pythia6`; the driver recompiles what changed.
//
// [tune] on the command line selects another PYTUNE index: the Perugia 2012 variants 371 radHi,
// 372 radLo, 373 IBK, 376 FL, 377 FT of the hadronisation-correction uncertainty (PARP(90) stays
// 0.213). With [partons.root] every event is written twice — the record before hadronisation
// (MSTP(111) = 0, i.e. the partons after the shower) to that file and the particles after PYEXEC to
// <out.root> — the same events with the same cross section, which is what C_had is built from.
//
// Output: the particle tree of tree.h; the cross section of the bin is PYINT5 XSEC(0,3) after the last
//         event. Pythia 6.4.29 of the LCG view (the embedding itself ran 6.4.28).
// Usage : gen_pythia6 <pthat_min> <pthat_max|-1> <nevents> <seed> <out.root> [tune] [partons.root]
#include "tree.h"
#include <cstdio>
#include <cstdlib>
#include <string>

extern "C" {
void pytune_(int *itune);
void pyinit_(const char *frame, const char *beam, const char *target, double *win, int, int, int);
void pyevnt_();
void pyexec_();
void pystat_(int *mstat);
int pycomp_(int *kf);
extern struct { int n, npad, k[5][4000]; double p[5][4000], v[5][4000]; } pyjets_;
extern struct { int mdcy[3][500], mdme[2][8000]; double brat[8000]; int kfdp[5][8000]; } pydat3_;
extern struct { int mstp[200]; double parp[200]; int msti[200]; double pari[200]; } pypars_;
extern struct { int msel, mselpd, msub[500], kfin[81][2]; double ckin[200]; } pysubs_;
extern struct { int ngenpd, ngen[3][501]; double xsec[3][501]; } pyint5_;
extern struct { int mrpy[6]; double rrpy[100]; } pydatr_;
}

// ---- settings -------------------------------------------------------------------------------------
// Perugia 2012, the tune of the embedding request (CTEQ6L1 parton distributions through LHAPDF).
static const int kDefaultTune = 370;

// PARP(90): the exponent of the energy dependence of the multiparton-interaction cut-off; the
// embedding request uses 0.213 instead of the 0.24 of the tune.
static const double kParp90 = 0.213;

// MSEL 1: inclusive QCD 2 -> 2. The pT-hat bin goes into CKIN(3), CKIN(4) from the command line.
static const int kProcessSelection = 1;

// pp centre-of-mass energy of the 2012 RHIC run, GeV.
static const double kSqrtS = 200.0;

// The 13 species the embedding request keeps stable in the generator, GEANT decaying them later.
// Remove a line to let Pythia decay that species.
static const int kStableIds[] = {
   211,  // pi+-
   321,  // K+-
   310,  // K0S
   130,  // K0L
   3122, // Lambda
   3112, // Sigma-
   3222, // Sigma+
   3312, // Xi-
   3322, // Xi0
   3334, // Omega-
   111,  // pi0
   221,  // eta
   3212, // Sigma0
};
// ---------------------------------------------------------------------------------------------------

// mb -> pb, the unit the particle trees carry the cross section in.
static const double kMbToPb = 1e9;

// Copy the current PYJETS record — the entries with a stable/decayed-particle status code — into the
// tree and write the event out.
static void RecordEvent(Models::Writer &out)
{
   for (int j = 0; j < pyjets_.n; ++j) {
      const int status = pyjets_.k[0][j];
      if (status < 1 || status > 10) continue;
      out.ev.Add(pyjets_.p[0][j], pyjets_.p[1][j], pyjets_.p[2][j], pyjets_.p[3][j], pyjets_.k[1][j]);
   }
   out.ev.pthat = pypars_.pari[16];
   out.Fill();
}

int main(int argc, char **argv)
{
   if (argc < 6 || argc > 8) {
      fprintf(stderr, "usage: gen_pythia6 <pthat_min> <pthat_max|-1> <nevents> <seed> <out.root> [tune] [partons.root]\n");
      return 1;
   }
   const double lo = atof(argv[1]);
   const double hi = atof(argv[2]);
   const long nev = atol(argv[3]);
   const int seed = atoi(argv[4]);
   const bool partons = argc > 7;

   // 1. the generator settings: the tune first, then everything the tune does not fix
   pydatr_.mrpy[0] = seed;
   int tune = argc > 6 ? atoi(argv[6]) : kDefaultTune;
   pytune_(&tune);
   if (partons) pypars_.mstp[110] = 0; // MSTP(111) = 0: stop before hadronisation
   pypars_.parp[89] = kParp90;
   for (int id : kStableIds) {
      int kf = id;
      pydat3_.mdcy[0][pycomp_(&kf) - 1] = 0; // MDCY(KC,1) = 0: the species does not decay
   }
   pysubs_.msel = kProcessSelection;
   pysubs_.ckin[2] = lo; // CKIN(3), CKIN(4): the pT-hat bin
   pysubs_.ckin[3] = hi;
   double win = kSqrtS;
   pyinit_("CMS", "p", "p", &win, 3, 1, 1);

   // 2. the output trees, stamped with the settings
   char settings[256];
   snprintf(settings, sizeof settings,
            "Pythia 6.4.29, PYTUNE %d (Perugia 2012%s, CTEQ6L1), PARP(90)=%g, MSEL %d, "
            "CKIN(3)=%g CKIN(4)=%g, pp %g GeV, %zu particles undecayed, seed %d",
            tune, tune == kDefaultTune ? "" : " variant", kParp90, kProcessSelection, lo, hi, kSqrtS,
            sizeof kStableIds / sizeof kStableIds[0], seed);
   Models::Writer particleTree(argv[5], settings);
   Models::Writer *partonTree =
      partons ? new Models::Writer(argv[7], std::string(settings) + ", parton level (before hadronisation)") : nullptr;

   // 3. generate; with a parton file the same event is recorded before and after hadronisation
   for (long i = 0; i < nev; ++i) {
      pyevnt_();
      if (partonTree) {
         RecordEvent(*partonTree);
         pyexec_();
      }
      RecordEvent(particleTree);
   }

   // 4. the cross section of the bin is only final after the last event
   int one = 1;
   pystat_(&one);
   const double sigma_pb = pyint5_.xsec[2][0] * kMbToPb;
   particleTree.Close(sigma_pb, 0.0);
   if (partonTree) {
      partonTree->Close(sigma_pb, 0.0);
      delete partonTree;
   }
   printf("gen_pythia6: %ld events, tune %d, pthat %g-%g GeV, sigma = %g pb -> %s%s%s\n", nev, tune, lo, hi, sigma_pb,
          argv[5], partons ? " + " : "", partons ? argv[7] : "");
   return 0;
}
