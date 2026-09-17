#ifndef CROSS_SECTION_CONFIG_H
#define CROSS_SECTION_CONFIG_H
//
// Shared configuration of Stage-2 (stage2/build_data.C, build_resp.C,
// unfold.C): paths, the published truth binning and the bad-run list. No environment
// flags — every choice is hardcoded and edited here.
// (published/, analysis/ and models/ include it through the stage2 macros or directly.)
// The shared code of the engine — sample weights, run selection, trigger levels, the
// matched-tree reader and the jet pairing — lives in stage2/common.h.

#include <string>
#include <vector>

#include <TSystem.h>

namespace CrossSectionConfig {

const std::string kWorkDir  = "/gpfs01/star/pwg/prozorov/study_pp2012/fable5/jets_pp_2012/new_ana/";
const std::string kDataPath = "/gpfs01/star/pwg/prozorov/study_pp2012/fable5/jets_pp_2012/output/";

// Every per-radius artifact (responses, unfolded cross sections, plots) lives in
// <kWorkDir>results/R<R>/; the shared inputs (luminosities, chain efficiency, the
// published table) stay in <kWorkDir>inputs/.
inline std::string RadiusDir(const std::string &jetR)
{
   const std::string dir = kWorkDir + "results/R" + jetR + "/";
   gSystem->mkdir(dir.c_str(), /*recursive=*/true);
   return dir;
}

const int kImtThreads = 4; // bounded IMT (unbounded IMT over gpfs on >4GB trees segfaults)

// MC truth bins (the published binning). The trailing 52->86 cell is a FEED-DOWN
// BUFFER (never quoted): it catches matched pairs whose truth mc>52 but whose
// reco (~0.8x mc) still lands in the measured 40-52 region, so the last QUOTED
// bin (44-52) passes train/test closure.
const std::vector<double> pt_mc_bins = {6.9,  8.2,  9.7,  11.5, 13.6, 16.1, 19.0, 22.5,
                                        26.6, 31.4, 37.2, 44.0, 52.0, 86.0};

inline const std::vector<double> &McBins() { return pt_mc_bins; }

struct AnalysisConfig {
   std::string workdir  = kWorkDir;
   std::string datapath = kDataPath;

   // Bad runs (8 detector-quality + 7 dead-time); lists/badrun_extras.list adds
   // more at runtime, and the luminosity sums follow the same list.
   std::vector<int> badRuns = {
      13050011, 13059087, 13055015, 13069004, 13066101,
      13066102, 13066104, 13066109,
      13048092, 13049006, 13049007, 13051074, 13052061, 13069023, 13070061
   };
};

} // namespace CrossSectionConfig

#endif // CROSS_SECTION_CONFIG_H
