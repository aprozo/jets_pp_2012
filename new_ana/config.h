#ifndef CROSS_SECTION_CONFIG_H
#define CROSS_SECTION_CONFIG_H

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace CrossSectionConfig {

// // PT binning for reconstructed jets (reco-level)
const std::vector<double> pt_reco_bins = {
    10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
    21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0,
    32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 42.0, 44.0,
    46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 68.0, 80.0, 90.0};
// PT binning for MC truth jets (particle-level)
const std::vector<double> pt_mc_bins = {5.0,  6.9,  8.2,  9.7,  11.5,
                                        13.6, 16.1, 19.0, 22.5, 26.6,
                                        31.4, 37.2, 44.0, 52.0, 80.0};

std::unordered_map<int, int> LoadRunToBinMap(const std::string &path) {
  std::unordered_map<int, int> m;
  std::ifstream in(path.c_str());
  if (!in)
    throw std::runtime_error("Cannot open: " + path);

  std::string line;
  int run = 0, bin = 0;

  while (std::getline(in, line)) {
    bin++; // histograms start from 1
    // allow comments / blank lines
    if (line.empty())
      continue;
    if (line[0] == '#')
      continue;
    std::istringstream ss(line);
    if (!(ss >> run))
      continue; // or throw if you want strict parsing

    m[run] = bin;
  }
  return m;
}

// Analysis configuration
struct AnalysisConfig {
  std::string workdir = "/home/prozorov/dev/star/jets_pp_2012/new_ana/";
  std::string datapath = "/home/prozorov/dev/star/jets_pp_2012/output/";

  std::unordered_map<int, int> runMap =
      LoadRunToBinMap(workdir + "run_map.txt");

  // Available triggers
  std::vector<std::string> triggers = {"JP2", "HT2"};
  std::vector<std::string> jetRs = {"0.5", "0.6"};
  // Bad runs to exclude (JP2 trigger)
  std::vector<int> badRuns = {
      13070061, 13048019, 13048092, 13048093, 13049006,
      13049007, 13051074, 13051099, 13052061, 13061035,
      13064067, 13068060, 13069023, 13063012

      // // ///////extra strict runs
      // ,
      // 13043029, 13043034, 13043046, 13043062, 13044015, 13045123, 13045139,
      // 13046009, 13046016, 13046022, 13046027, 13046116, 13047021, 13047121,
      // 13048009, 13048010, 13048011, 13048012, 13048013, 13048014, 13048015,
      // 13048016, 13048017, 13048018, 13048030, 13048031, 13048032, 13048040,
      // 13048041, 13048042, 13048043, 13048044, 13048045, 13048049, 13048050,
      // 13048051, 13048052, 13048053, 13048087, 13048088, 13048089, 13048090,
      // 13048091,

      // 13065018, 13050023, 13057019,
      // ////
      // 13063022, 13057009,
      // ////
      // 13071064, 13072002, 13072003, 13072005
  };

  // Number of unfolding iterations (Bayesian)
  int nIterations = 4;
};

// Color palette for plots
const std::vector<int> colors = {2000, 2002, 2003, 2004,
                                 2005, 2006, 2007, 2008};

// Marker styles for plots
const std::vector<int> markers = {20, 21, 22, 23, 33, 34};

} // namespace CrossSectionConfig

#endif // CROSS_SECTION_CONFIG_H
