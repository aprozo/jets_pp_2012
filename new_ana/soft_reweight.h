#ifndef NEW_ANA_SOFT_REWEIGHT_H
#define NEW_ANA_SOFT_REWEIGHT_H
//
// Soft pT reweight of the embedding (the published analysis, Eq. 5-6, §4.4). The Pythia6
// Perugia2012 sample was generated with CKIN(3) = 2 GeV, so the soft regime is artificially
// suppressed; each event is reweighted by the cross-section ratio
//
//     w(pthat) = sigma_[0,inf](pthat) / sigma_[2,inf](pthat)
//              ~ 1 / (1 + (1.22 - 0.33k + 0.17k^2) exp(-0.82k)),   k = (pthat - 2 GeV)/(1 GeV)
//
// the coefficients being a fit to that ratio at sqrt(s) = 200 GeV.
//
// The true partonic pthat is not stored in the pico, so pthat is taken from the MatchedTree
// row: the generator value when present, else `pthat_mid`, the midpoint of the generated
// pt-hat bin (filled in macros/matching_mc_reco.cxx from the input filename).
#include <cmath>

namespace SoftReweight {

inline double weight(double pthat_gev) {
    if (pthat_gev <= 0) return 1.0;   // no info -> no reweight
    const double kappa = pthat_gev - 2.0;
    const double poly  = 1.22 - 0.33 * kappa + 0.17 * kappa * kappa;
    return 1.0 / (1.0 + poly * std::exp(-0.82 * kappa));
}

}  // namespace SoftReweight

#endif  // NEW_ANA_SOFT_REWEIGHT_H
