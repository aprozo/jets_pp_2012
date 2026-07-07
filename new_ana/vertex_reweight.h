#ifndef NEW_ANA_VERTEX_REWEIGHT_H
#define NEW_ANA_VERTEX_REWEIGHT_H
//
// User-side analogue of Dmitry's SetVertexReweightingParams() — reshapes the
// embedding TPC primary-vertex-z distribution to match the data distribution.
//
// Formula (mirrors Dmitry/star-jet/StJetPlots/StVertexReweighting.h:80-85):
//
//     w(v_z) = exp[(1/sigma^2 - 1/sigma_target^2) * v_z^2 / 2]
//            * (1 + (v_z/b)^2) / (1 + (v_z/b_target)^2)
//            * norm(sigma, b) / norm(sigma_target, b_target)
//
// where norm(sigma, b=infty) = sigma * sqrt(2*pi), and otherwise
//       norm(sigma, b)       = b * exp(b^2 / (2 sigma^2))
//                              * pi * erfc(b / (sqrt(2)*sigma)).
//
// Defaults match Dmitry's default.nix `process_embed_overrides`:
//     vertex_reweight_sigma         = 45 cm  (embedding TPC vertex width)
//     vertex_reweight_b             = +infty
//     vertex_reweight_sigma_target  = 70 cm  (data TPC vertex width)
//     vertex_reweight_b_target      = 80 cm
//
// Used as an event-level weight on top of mc_weight in the response build
// (bayes/unfold.cxx).
#include <cmath>
#include <limits>

namespace VertexReweight {

struct Params {
    // KNOWN OPEN SYSTEMATIC (2026-06-12): sigma should be THIS embedding's
    // actual truth-vz spread. Dmitry's 45.0 matches HIS v3 embedding
    // (measured 45.2); MY 20235003 embedding measures 37.3. Setting 37.3
    // here was TESTED and is statistically UNSTABLE: the weight to reach the
    // sigma_target=70/b=80 data shape explodes on the sparse high-|vz|
    // truth-only events (x2.5 at vz=60, x5 at 80) -> jagged spectra. Kept at
    // 45 (smooth, the validated configuration). Measured 2026-06-27: a capped
    // source-sigma override lifts the spectrum near-FLAT ~+8% -> carried as a
    // normalization systematic, central stays sigma=45.
    double sigma           = 45.0;
    double b               = std::numeric_limits<double>::infinity();
    double sigma_target    = 70.0;
    double b_target        = 80.0;
};

// Normalisation constant of the modified-Gaussian model used by Dmitry:
//     f(v_z) ∝ exp(-v_z^2 / (2 sigma^2)) * (1 + (v_z/b)^2)^{-1}
// reduces to a plain Gaussian when b -> infty.
inline double norm(double sigma, double b) {
    if (std::isinf(b)) {
        return sigma * std::sqrt(2.0 * M_PI);
    }
    return b * std::exp(b * b / (2.0 * sigma * sigma))
           * M_PI * std::erfc(b / (std::sqrt(2.0) * sigma));
}

inline double weight(double vz, const Params &p = Params{}) {
    const double inv_sigma2_diff =
        1.0 / (p.sigma * p.sigma) - 1.0 / (p.sigma_target * p.sigma_target);
    const double exp_factor      = std::exp(0.5 * inv_sigma2_diff * vz * vz);
    const double poly_num        =
        std::isinf(p.b) ? 1.0 : (1.0 + (vz / p.b) * (vz / p.b));
    const double poly_den        = 1.0 + (vz / p.b_target) * (vz / p.b_target);
    const double poly_factor     = poly_num / poly_den;
    const double norm_factor     = norm(p.sigma, p.b) / norm(p.sigma_target, p.b_target);
    return exp_factor * poly_factor * norm_factor;
}

}  // namespace VertexReweight

#endif  // NEW_ANA_VERTEX_REWEIGHT_H
