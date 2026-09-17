// tree.h — the particle tree shared by the generator drivers (gen_pythia6.cc, gen_pythia8.cc,
// hepmc2tree.cc) and the jet analysis (jetspec.cc).
//
// Per event: the generator final state within |eta| < 3 (px, py, pz, E in GeV and the PDG id of each
// particle) and the hard-process pT (pthat, GeV; for Herwig the scale of the PDF record).
// Per file: the cross section of the generated pT-hat bin (TParameter sigma_pb, sigma_err_pb) and the
// generator settings (TNamed settings); the number of generated events is the number of tree entries.
// The weight of one event in a spectrum is therefore sigma_pb / N_generated.
//
// Writer: construct with the output file and the settings string, fill ev, call Fill() per event and
// Close(sigma, err) once. Reader: construct with the file, loop while Next(), read ev / sigma_pb / nev.
#ifndef MODELS_TREE_H
#define MODELS_TREE_H
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <cmath>
#include <stdexcept>
#include <string>

namespace Models {

// Array size of one event; the softest pT-hat bins of pp at 200 GeV stay far below it.
const int kMaxPart = 4000;
// Only particles inside this pseudorapidity are stored: wide enough for jets up to R = 0.5 at
// |eta_jet| < 0.5 with the off-axis UE cones, narrow enough to keep the trees small.
const double kEtaStore = 3.0;

struct Event {
   Int_t n = 0;              // number of stored particles
   Float_t px[kMaxPart];     // momentum components in GeV
   Float_t py[kMaxPart];
   Float_t pz[kMaxPart];
   Float_t e[kMaxPart];      // energy in GeV
   Int_t pid[kMaxPart];      // PDG id
   Double_t pthat = -1;      // hard-process pT of the event in GeV

   void Clear()
   {
      n = 0;
      pthat = -1;
   }

   // Store one final-state particle, unless it is outside |eta| < kEtaStore or the event is full.
   void Add(double x, double y, double z, double E, int id)
   {
      const double p = std::sqrt(x * x + y * y + z * z);
      if (p <= 0 || std::fabs(z) >= p || n >= kMaxPart) return;
      const double eta = 0.5 * std::log((p + z) / (p - z));
      if (std::fabs(eta) > kEtaStore) return;
      px[n] = x;
      py[n] = y;
      pz[n] = z;
      e[n] = E;
      pid[n] = id;
      ++n;
   }
};

class Writer {
public:
   Writer(const std::string &file, const std::string &settings)
      : f_(file.c_str(), "RECREATE"), t_("particles", "generator final state, |eta| < 3")
   {
      t_.Branch("n", &ev.n, "n/I");
      t_.Branch("px", ev.px, "px[n]/F");
      t_.Branch("py", ev.py, "py[n]/F");
      t_.Branch("pz", ev.pz, "pz[n]/F");
      t_.Branch("e", ev.e, "e[n]/F");
      t_.Branch("pid", ev.pid, "pid[n]/I");
      t_.Branch("pthat", &ev.pthat, "pthat/D");
      TNamed("settings", settings.c_str()).Write();
   }

   void Fill()
   {
      t_.Fill();
      ev.Clear();
   }

   // The cross section of the pT-hat bin is known only after the last event.
   void Close(double sigma_pb, double sigma_err_pb)
   {
      f_.cd();
      t_.Write();
      TParameter<double>("sigma_pb", sigma_pb).Write();
      TParameter<double>("sigma_err_pb", sigma_err_pb).Write();
      f_.Close();
   }

   Event ev;

private:
   TFile f_;
   TTree t_;
};

class Reader {
public:
   explicit Reader(const std::string &file) : f_(file.c_str())
   {
      if (f_.IsZombie()) throw std::runtime_error("cannot open " + file);
      t_ = (TTree *)f_.Get("particles");
      auto *sigma = (TParameter<double> *)f_.Get("sigma_pb");
      if (!t_ || !sigma) throw std::runtime_error("no particle tree / cross section in " + file);
      sigma_pb = sigma->GetVal();
      nev = t_->GetEntries();
      t_->SetBranchAddress("n", &ev.n);
      t_->SetBranchAddress("px", ev.px);
      t_->SetBranchAddress("py", ev.py);
      t_->SetBranchAddress("pz", ev.pz);
      t_->SetBranchAddress("e", ev.e);
      t_->SetBranchAddress("pid", ev.pid);
      t_->SetBranchAddress("pthat", &ev.pthat);
   }

   bool Next() { return i_ < nev && t_->GetEntry(i_++) > 0; }

   double sigma_pb = 0; // cross section of the generated pT-hat bin, pb
   long long nev = 0;   // generated events in the file
   Event ev;

private:
   TFile f_;
   TTree *t_ = nullptr;
   long long i_ = 0;
};

} // namespace Models
#endif
