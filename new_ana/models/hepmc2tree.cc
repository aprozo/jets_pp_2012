// hepmc2tree.cc — convert the HepMC output of Herwig 7 into the particle tree of tree.h, so that the
// Herwig sample is read by jetspec.cc exactly like the Pythia ones.
//
// Final state = HepMC status 1; pthat = the scale of the PDF record, which for QCD 2->2 is the hard-
// process kT; the cross section of the bin = the GenCrossSection of the last event (the generator's
// running estimate, already in pb).
//
// Input : <in.hepmc> written by the HepMCFile analysis handler of herwig.in.
// Output: <out.root>, the particle tree, stamped with the <settings> string given on the command line.
// Usage : hepmc2tree <in.hepmc> <out.root> <settings>
#include "tree.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenPdfInfo.h"
#include "HepMC3/ReaderFactory.h"
#include <cstdio>

// HepMC status code of a final-state particle.
static const int kStatusFinal = 1;

int main(int argc, char **argv)
{
   if (argc != 4) {
      fprintf(stderr, "usage: hepmc2tree <in.hepmc> <out.root> <settings>\n");
      return 1;
   }
   auto reader = HepMC3::deduce_reader(argv[1]);
   if (!reader) {
      fprintf(stderr, "hepmc2tree: cannot read %s\n", argv[1]);
      return 2;
   }

   Models::Writer particleTree(argv[2], argv[3]);
   double sigma_pb = 0, sigma_err_pb = 0;
   long nev = 0;
   while (!reader->failed()) {
      HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
      reader->read_event(event);
      if (reader->failed()) break;
      event.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

      for (const auto &particle : event.particles()) {
         if (particle->status() != kStatusFinal) continue;
         const auto &p = particle->momentum();
         particleTree.ev.Add(p.px(), p.py(), p.pz(), p.e(), particle->pid());
      }
      const auto pdf = event.pdf_info();
      particleTree.ev.pthat = pdf ? pdf->scale : -1;

      // the cross section is cumulative: the value of the last event is the one for the whole file
      const auto crossSection = event.cross_section();
      if (crossSection) {
         sigma_pb = crossSection->xsec();
         sigma_err_pb = crossSection->xsec_err();
      }
      particleTree.Fill();
      ++nev;
   }
   reader->close();
   particleTree.Close(sigma_pb, sigma_err_pb);
   printf("hepmc2tree: %ld events, sigma = %g +- %g pb -> %s\n", nev, sigma_pb, sigma_err_pb, argv[2]);
   return 0;
}
