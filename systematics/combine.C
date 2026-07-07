// combine.C — build the systematic band per trigger from the variation xsec
// files that the pipeline produced (one real re-run per config.h::Systematics()
// entry, each stamped with its provenance).
//
// For each trigger it reads the nominal xsec_<T>_R0.5.root and scans the workdir
// for xsec_<T>_R0.5_<name>.root, takes the per-bin up/down envelope (quadrature)
// of every variation around the nominal, and writes xsec_<T>_R0.5_systband.root:
//   canonical  the nominal spectrum
//   systematic the nominal with bin errors = the systematic envelope
//   reference  Dmitry's paper band (copied through)
//
// Workdir + trigger list come from config.h (the single source of truth).
// Usage: root -l -b -q combine.C
#include "../new_ana/config.h"

#include <TFile.h>
#include <TH1D.h>
#include <TList.h>
#include <TNamed.h>
#include <TString.h>
#include <TSystemDirectory.h>
#include <cmath>
#include <string>
#include <vector>

using namespace CrossSectionConfig;
AnalysisConfig cfg;

static void combine_one(const std::string &trig)
{
   const std::string dir = cfg.workdir;
   TFile *fnom = TFile::Open(Form("%sxsec_%s_R0.5.root", dir.c_str(), trig.c_str()));
   if (!fnom || fnom->IsZombie()) { printf("no nominal xsec for %s\n", trig.c_str()); return; }
   TH1D *nom = (TH1D *)fnom->Get("canonical");
   TH1D *ref = (TH1D *)fnom->Get("reference");
   if (!nom) { printf("no canonical for %s\n", trig.c_str()); return; }
   nom->SetDirectory(0);
   if (ref) ref->SetDirectory(0);

   const int nb = nom->GetNbinsX();
   std::vector<double> up(nb + 1, 0.0), dn(nb + 1, 0.0);

   const TString pref = Form("xsec_%s_R0.5_", trig.c_str());  // nominal ("...R0.5.root") does NOT match
   const TString band = Form("xsec_%s_R0.5_systband.root", trig.c_str());
   TSystemDirectory sdir(dir.c_str(), dir.c_str());
   TList *files = sdir.GetListOfFiles();
   if (files) {
      TIter next(files);
      while (TObject *o = next()) {
         TString fn = o->GetName();
         if (!fn.BeginsWith(pref) || !fn.EndsWith(".root") || fn == band) continue;
         TFile *fv = TFile::Open(Form("%s%s", dir.c_str(), fn.Data()));
         if (!fv || fv->IsZombie()) continue;
         TH1D *hv = (TH1D *)fv->Get("canonical");
         TNamed *prov = (TNamed *)fv->Get("variation");
         if (hv) {
            for (int b = 1; b <= nb; ++b) {
               const double d = hv->GetBinContent(b) - nom->GetBinContent(b);
               if (d > 0) up[b] = std::sqrt(up[b] * up[b] + d * d);
               else       dn[b] = std::sqrt(dn[b] * dn[b] + d * d);
            }
            printf("  + %s (variation=%s)\n", fn.Data(), prov ? prov->GetTitle() : "?");
         }
         fv->Close();
      }
   }

   TH1D *syst = (TH1D *)nom->Clone("systematic");
   syst->SetDirectory(0);
   for (int b = 1; b <= nb; ++b)
      syst->SetBinError(b, 0.5 * (up[b] + dn[b])); // symmetrized envelope

   TFile fout(Form("%sxsec_%s_R0.5_systband.root", dir.c_str(), trig.c_str()), "RECREATE");
   nom->Write("canonical");
   syst->Write("systematic");
   if (ref) ref->Write("reference");
   fout.Close();
   printf("wrote %sxsec_%s_R0.5_systband.root\n", dir.c_str(), trig.c_str());
}

void combine()
{
   for (const auto &t : cfg.triggers)
      combine_one(t);
}
