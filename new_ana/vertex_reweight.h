#ifndef NEW_ANA_VERTEX_REWEIGHT_H
#define NEW_ANA_VERTEX_REWEIGHT_H
//
// Vertex-z reweighting of the embedding (the published analysis): reshapes the embedding
// TPC primary-vertex-z distribution into the data one. Event-level weight on top of the
// sample weight in the response build (stage2/build_resp.C).
//
//     w(v_z) = exp[(1/sigma^2 - 1/sigma_target^2) * v_z^2 / 2]
//            * (1 + (v_z/b)^2) / (1 + (v_z/b_target)^2)
//            * norm(sigma, b) / norm(sigma_target, b_target)
//
// norm() being the normalisation of the modified-Gaussian model
//     f(v_z) ~ exp(-v_z^2 / (2 sigma^2)) * (1 + (v_z/b)^2)^-1,
// which reduces to a plain Gaussian as b -> infinity.
#include <cmath>
#include <limits>

namespace VertexReweight {

struct Params {
    // sigma is the published 45 cm rather than this embedding's own truth-vz width
    // (37.3 cm): the steeper weight that width implies explodes on the sparse
    // high-|vz| truth-only events and makes the spectra jagged. Its near-flat effect
    // on the spectrum is carried as a normalisation systematic instead.
    double sigma           = 45.0; // embedding TPC vertex width
    double b               = std::numeric_limits<double>::infinity();
    double sigma_target    = 70.0; // data TPC vertex width
    double b_target        = 80.0;
};

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
