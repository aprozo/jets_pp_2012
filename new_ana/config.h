#ifndef CROSS_SECTION_CONFIG_H
#define CROSS_SECTION_CONFIG_H
//
// Shared configuration for the Stage-2 physics pipeline. NO environment flags —
// every choice is hardcoded to its validated physical value and edited here.
//   unfolding/unfold.cxx      builds the Miss/Fake response per trigger
//   cross_section.cpp         unfolds (RooUnfoldBayes), normalizes, writes xsec
//   cross_section_inverse/    the same, unfolded by Dmitry's matrix inversion
//   plot_alltriggers.C        overlays every trigger against Dmitry's Table III
//
// Triggers are ANALYSIS FILTERS, not a production split: one Stage-1 production
// (merged_data_R<R>.root, merged_matching_R<R>.root) carries every trigger's
// per-jet trigger_match_<T> / per-event fired_<T> bit, selected here at Stage-2.
// The data-side trigger correction is the measured HYBRID C(pt) = T-hat (turn-on)
// then R (plateau ruler) — both physical, applied before unfolding. Each trigger
// is quoted only in its efficient window (QuoteLo); below it the standalone
// spectrum is turn-on extrapolation, not a measurement, and is dropped.

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <TColor.h>
#include <TNamed.h>
#include <TString.h>

namespace CrossSectionConfig {

// ---- fixed paths (edit here if the repo moves) -----------------------------
const std::string kWorkDir  = "/gpfs01/star/pwg/prozorov/study_pp2012/fable5/jets_pp_2012/new_ana/";
const std::string kDataPath = "/gpfs01/star/pwg/prozorov/study_pp2012/fable5/jets_pp_2012/output/";

// ---- fixed analysis choices ------------------------------------------------
const int kNIter      = 2; // RooUnfoldBayes iterations (validated: no rising tail)
const int kImtThreads = 4; // bounded IMT (unbounded IMT over gpfs on >4GB trees segfaults)

// Allocate the custom ROOT color slots used throughout the analysis.
inline void DefineCustomColors()
{
   static bool defined = false;
   if (defined)
      return;
   defined = true;

   new TColor(2000, 255 / 255., 89 / 255., 74 / 255.);
   new TColor(2001, 25 / 255., 170 / 255., 25 / 255.);
   new TColor(2002, 66 / 255., 98 / 255., 255 / 255.);
   new TColor(2003, 153 / 255., 0 / 255., 153 / 255.);
   new TColor(2004, 255 / 255., 166 / 255., 33 / 255.);
   new TColor(2005, 0 / 255., 170 / 255., 255 / 255.);
   new TColor(2006, 204 / 255., 153 / 255., 255 / 255.);
   new TColor(2007, 107 / 255., 142 / 255., 35 / 255.);
   new TColor(2008, 100 / 255., 149 / 255., 237 / 255.);
   new TColor(2009, 255 / 255., 69 / 255., 0 / 255.);
   new TColor(2010, 0 / 255., 128 / 255., 128 / 255.);
   new TColor(2011, 176 / 255., 196 / 255., 222 / 255.);
   new TColor(2012, 255 / 255., 215 / 255., 0 / 255.);
}

// Per-trigger reco-pT validity floor (GeV, on the UE-subtracted reco pT): the
// start of the trigger's efficient region. Shared by the data selection
// (cross_section.cpp::Raw) and the response reco side (unfold.cxx).
inline double TrigPtFloor(const std::string &trigger)
{
   if (trigger == "JP2") return 9.7;
   if (trigger == "JP1") return 8.2;
   if (trigger == "JP0") return 6.9;
   if (trigger == "HT2") return 11.5;
   return 0.0;
}

// Per-trigger QUOTE-window low edge (GeV): the lowest pT the trigger is quoted
// at. Below it the standalone spectrum is unfolding extrapolation below the
// turn-on, so it is dropped from the quoted spectrum and the comparison.
inline double QuoteLo(const std::string &trigger)
{
   if (trigger == "JP2") return 13.6;
   if (trigger == "JP1") return 8.2;
   if (trigger == "JP0") return 13.6;
   if (trigger == "HT2") return 11.5;
   return 0.0;
}

// PT binning for reconstructed jets (reco-level): fine 1-GeV mid-pT, coarser
// high-pT. Bayes regularization works better when reco bins are finer than truth.
const std::vector<double> pt_reco_bins = {6.9, 8.2, 9.7, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,
                                          19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5,
                                          30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 41.0,
                                          43.0, 45.0, 47.0, 49.0, 51.0, 53.0, 55.0, 57.0, 60.0, 68.0, 80.0};
// MC truth bins (Dmitry-aligned). The trailing 52->86 cell is a FEED-DOWN
// BUFFER (never quoted): it catches matched pairs whose truth mc>52 but whose
// reco (~0.8x mc) still lands in the measured 40-52 region, so the last QUOTED
// bin (44-52) passes train/test closure.
const std::vector<double> pt_mc_bins = {6.9,  8.2,  9.7,  11.5, 13.6, 16.1, 19.0, 22.5,
                                        26.6, 31.4, 37.2, 44.0, 52.0, 86.0};

inline const std::vector<double> &McBins()   { return pt_mc_bins; }
inline const std::vector<double> &RecoBins() { return pt_reco_bins; }

// =====================================================================
// MEASURED trigger correction C(pt), applied to the JP1/JP2 data before
// unfolding (cross_section.cpp::Raw). HYBRID: T-hat (turn-on, pt<16.1) then R
// (plateau ruler, pt>=16.1).
// =====================================================================
// T-hat = [P(emulator-match|pt) data(JP0-fired base) / embedding] — the residual
// data-vs-embedding turn-on SHAPE (the embedding calorimeter fires "hotter" than
// the real one). R = the hardware-vs-software match ratio, applied ONLY on the
// plateau (below 16.1 T-hat already carries the hardware/software mismatch, so
// multiplying R in there double-counts it). JP0 / HT2 get no correction.
inline double TrigEffMeas(const std::string &trigger, double pt)
{
   if (trigger == "JP1") {
      if      (pt <  8.2) return 0.8339; // T-hat turn-on
      else if (pt <  9.7) return 0.8957;
      else if (pt < 11.5) return 0.9378;
      else if (pt < 13.6) return 0.9707;
      else if (pt < 16.1) return 0.9929;
      else if (pt < 19.0) return 0.9725; // R plateau (ruler)
      else if (pt < 22.5) return 0.9741;
      else if (pt < 26.6) return 0.9783;
      else if (pt < 31.4) return 0.9877;
      else if (pt < 37.2) return 0.9744;
      else                return 1.0;
   }
   if (trigger == "JP2") {
      if      (pt <  8.2) return 0.5266; // T-hat turn-on
      else if (pt <  9.7) return 0.6630;
      else if (pt < 11.5) return 0.7737;
      else if (pt < 13.6) return 0.8607;
      else if (pt < 16.1) return 0.9575;
      else if (pt < 19.0) return 0.9571; // R plateau (ruler)
      else if (pt < 22.5) return 0.9814;
      else if (pt < 26.6) return 0.9954;
      else if (pt < 31.4) return 1.0000;
      else if (pt < 37.2) return 0.9722;
      else                return 1.0;
   }
   return 1.0; // JP0 / HT2: no data-side trigger correction
}

std::unordered_map<int, int> LoadRunToBinMap(const std::string &path)
{
   std::unordered_map<int, int> m;
   std::ifstream in(path.c_str());
   if (!in)
      throw std::runtime_error("Cannot open: " + path);

   std::string line;
   int run = 0, bin = 0;
   while (std::getline(in, line)) {
      bin++; // histograms start from 1
      if (line.empty() || line[0] == '#')
         continue;
      std::istringstream ss(line);
      if (!(ss >> run))
         continue;
      m[run] = bin;
   }
   return m;
}

struct AnalysisConfig {
   std::string workdir  = kWorkDir;
   std::string datapath = kDataPath;

   std::unordered_map<int, int> runMap = LoadRunToBinMap(workdir + "run_map.txt");

   // Triggers processed by every run (edit the list to run a subset). Jet
   // radius fixed at R=0.5.
   std::vector<std::string> triggers = {"JP0", "JP1", "JP2", "HT2"};
   std::vector<std::string> jetRs = {"0.5"};

   // Dmitry's bad runs (8 dmitry-bad + 7 deadtime). The runtime Leff in
   // cross_section.cpp follows this list automatically.
   std::vector<int> badRuns = {
      13050011, 13059087, 13055015, 13069004, 13066101,
      13066102, 13066104, 13066109,
      13048092, 13049006, 13049007, 13051074, 13052061, 13069023, 13070061
   };

   int nIterations = kNIter;
};

// ====================== SYSTEMATIC VARIATIONS (config-as-code) ======================
// A systematic is a named DELTA from the nominal — a typed, committed object, not
// an env var or a text row. The pipeline macros take a variation NAME as an
// argument (default "nominal" = identity) and look up its Systematic here.
//   jesShift/jerSmear  shift the reconstructed energy in the RESPONSE -> SHAPE
//                       systematics; the response is rebuilt (needsResponse()).
//   lumiScale          luminosity scale (the spectrum is divided by it).
//   trigEffScale       flat scale on C(pt)  (the spectrum is divided by it).
//   nIter              RooUnfoldBayes iterations (unfolding systematic; reuses the
//                       nominal response).
// The variation name is stamped into every output ROOT file (provenance), so the
// config travels with the data.
// Calorimeter / tracking scale uncertainties (Dmitry's values, star-jet
// default.nix:262-274): BEMC tower scale 3.2%, TPC track scale 1.1%, track
// efficiency 1%. The jet-level shift is weighted by the jet's OWN neutral
// fraction rt: d(pt)/pt = sqrt(((1-rt)*kTrackScale)^2 + (rt*kTowerScale)^2).
const double kTowerScaleUnc = 0.032;
const double kTrackScaleUnc = 0.011;
const double kTrackEffUnc   = 0.010;

struct Systematic {
   std::string name;
   double jesShift = 0.0;     // flat fractional reco-pT shift (generic studies)
   int    emcSign = 0;        // +-1: per-jet EMC scale shift, rt-weighted
                              // sqrt(((1-rt)*0.011)^2 + (rt*0.032)^2) (response)
   int    trkSign = 0;        // +-1: track-efficiency equivalent, per-jet
                              // 0.01*(1-rt) reco-pT shift (response; == data
                              // thinned in the opposite direction)
   double jerSmear = 0.0;     // fractional extra Gaussian reco smear (response)
   double ueFraction = 1.0;   // detector-side UE-subtraction fraction (data;
                              // nominal 1.0, variations 0.86 / 1.18)
   double lumiScale = 1.0;    // luminosity scale (NOT in the band by default —
                              // the 10% lumi normalization is quoted separately)
   double trigEffScale = 1.0; // flat scale on C(pt)
   int    nIter = kNIter;     // Bayes iterations
   double jpxLambda = -1.0;   // JPX Tikhonov damping (<0 = promotion.h default)
   int    embStatSign = 0;    // +-1: JPX embedding-statistics toys, nominal
                              // +- 1 sigma(toys) (response-statistics term)
   Systematic() = default;
   Systematic(std::string n) : name(n) {} // for {"nominal"} and the factories below
   bool needsResponse() const { return jesShift != 0.0 || jerSmear != 0.0 || emcSign != 0 || trkSign != 0; }
   // The per-jet reco-pT scale factor of the response-side shape variations.
   double RecoShift(double rt) const
   {
      double s = jesShift;
      if (emcSign != 0)
         s += emcSign * std::sqrt(std::pow((1.0 - rt) * kTrackScaleUnc, 2) +
                                  std::pow(rt * kTowerScaleUnc, 2));
      if (trkSign != 0) s += trkSign * kTrackEffUnc * (1.0 - rt);
      return s;
   }
};

// Named constructors (C++17: no designated initializers) — each makes the intent
// of a variation obvious at the definition site.
inline Systematic SystJES(const std::string &n, double jes)  { Systematic s; s.name = n; s.jesShift = jes; return s; }
inline Systematic SystEMC(const std::string &n, int sign)    { Systematic s; s.name = n; s.emcSign = sign; return s; }
inline Systematic SystTRK(const std::string &n, int sign)    { Systematic s; s.name = n; s.trkSign = sign; return s; }
inline Systematic SystJER(const std::string &n, double jer)  { Systematic s; s.name = n; s.jerSmear = jer; return s; }
inline Systematic SystUE(const std::string &n, double f)     { Systematic s; s.name = n; s.ueFraction = f; return s; }
inline Systematic SystNorm(const std::string &n, double lumi, double trig) { Systematic s; s.name = n; s.lumiScale = lumi; s.trigEffScale = trig; return s; }
inline Systematic SystIter(const std::string &n, int ni)     { Systematic s; s.name = n; s.nIter = ni; return s; }
inline Systematic SystJpxL(const std::string &n, double l)   { Systematic s; s.name = n; s.jpxLambda = l; return s; }
inline Systematic SystEmbS(const std::string &n, int sign)   { Systematic s; s.name = n; s.embStatSign = sign; return s; }

// The committed systematic matrix — MIRRORS Dmitry's published composition
// (star-jet default.nix:1213-1245: quadrature of EMC scale, track efficiency,
// embedding statistics, UE fraction; luminosity is a separate 10%
// normalization statement, "not shown", NOT folded into the per-bin band):
//   emcUp/Down    response reco-pT shifted per jet by the rt-weighted
//                 tower(3.2%)/track(1.1%) scale uncertainty
//   trkEffUp/Down 1% track-efficiency equivalent, per-jet 0.01*(1-rt)
//   ueUp/Down     detector-side UE-subtraction fraction 1.18 / 0.86 (data)
//   trigEffUp/Down +-1.5% on the measured C(pt) (its own measurement
//                 precision: plateau fit +-0.8%, turn-on bins +-1-2%)
//   unfoldReg     Bayes nIter 2->3 (per-trigger pipelines)
//   jpxDamp       JPX Tikhonov 0->0.030 (unregularized -> damped)
//   embStatUp/Down JPX embedding-statistics toys (+-1 sigma of 200 Poisson
//                 resamplings of the response ingredients)
inline std::vector<Systematic> Systematics()
{
   return {
      {"nominal"},
      SystEMC("emcUp", +1), SystEMC("emcDown", -1),
      SystTRK("trkEffUp", +1), SystTRK("trkEffDown", -1),
      SystUE("ueUp", 1.18), SystUE("ueDown", 0.86),
      SystNorm("trigEffUp", 1.0, 1.015), SystNorm("trigEffDown", 1.0, 0.985),
      SystIter("unfoldReg", 3),
      SystJpxL("jpxDamp", 0.030),
      SystEmbS("embStatUp", +1), SystEmbS("embStatDown", -1),
   };
}

inline Systematic FindSystematic(const std::string &name)
{
   for (const auto &s : Systematics())
      if (s.name == name)
         return s;
   throw std::runtime_error("unknown systematic: " + name);
}

// Output-filename tag (""=nominal).
inline std::string SystTag(const Systematic &s)
{
   return (s.name == "nominal" || s.name.empty()) ? std::string("") : ("_" + s.name);
}

// Provenance: write the full variation config as TNamed keys into the currently
// open TFile (call after cd-ing into it), so "what produced this histogram" is
// answerable from the file itself.
inline void StampProvenance(const Systematic &s)
{
   TNamed("variation", s.name.c_str()).Write();
   TNamed("jesShift", Form("%.4f", s.jesShift)).Write();
   TNamed("emcSign", Form("%d", s.emcSign)).Write();
   TNamed("trkSign", Form("%d", s.trkSign)).Write();
   TNamed("jerSmear", Form("%.4f", s.jerSmear)).Write();
   TNamed("ueFraction", Form("%.4f", s.ueFraction)).Write();
   TNamed("lumiScale", Form("%.4f", s.lumiScale)).Write();
   TNamed("trigEffScale", Form("%.4f", s.trigEffScale)).Write();
   TNamed("nIter", Form("%d", s.nIter)).Write();
   TNamed("jpxLambda", Form("%.4f", s.jpxLambda)).Write();
   TNamed("embStatSign", Form("%d", s.embStatSign)).Write();
}

const std::vector<int> colors = {2000, 2002, 2003, 2004, 2005, 2006, 2007, 2008};
const std::vector<int> markers = {20, 21, 22, 23, 33, 34};

} // namespace CrossSectionConfig

#endif // CROSS_SECTION_CONFIG_H
