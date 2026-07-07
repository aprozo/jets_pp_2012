#ifndef NEW_ANA_SOFT_REWEIGHT_H
#define NEW_ANA_SOFT_REWEIGHT_H
//
// Soft pT reweight from Dmitry's analysis note (Eq. 5-6, §4.4).
//
// The Pythia6 Perugia2012 sample was generated with CKIN(3) = 2 GeV, so
// the soft regime is artificially suppressed.  Dmitry restores the
// hard/soft balance by reweighting each event by the ratio
//
//     ω(p̂_T) = σ_[0, ∞](p̂_T) / σ_[2, ∞](p̂_T)
//
// which he parameterises (Eq. 6) as
//
//     ω(p̂_T) ≈ 1 / (1 + (1.22 − 0.33κ + 0.17κ²) exp(−0.82κ))
//              with κ = (p̂_T − 2 GeV) / (1 GeV)
//
// Coefficients are fixed by a fit to the Pythia6 Perugia2012 σ ratio at
// √s = 200 GeV.  Effect peaks at low p̂_T (~0.6 at κ≲1, → 1 at large κ).
//
// User-side note:  the true partonic p̂_T is not stored in the pico, so
// we approximate p̂_T by the midpoint of the Pythia pt-hat bin the event
// was generated in.  Each MatchedTree row carries `pthat_mid` (populated
// in macros/matching_mc_reco.cxx by parsing the input filename).
#include <cmath>

namespace SoftReweight {

inline double weight(double pthat_gev) {
    if (pthat_gev <= 0) return 1.0;   // no info → no reweight
    const double kappa = pthat_gev - 2.0;
    const double poly  = 1.22 - 0.33 * kappa + 0.17 * kappa * kappa;
    return 1.0 / (1.0 + poly * std::exp(-0.82 * kappa));
}

}  // namespace SoftReweight

#endif  // NEW_ANA_SOFT_REWEIGHT_H
