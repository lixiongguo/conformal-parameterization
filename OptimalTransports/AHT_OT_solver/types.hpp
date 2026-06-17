#pragma once

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

namespace aht_ot {

/// 2D vector field: .first = x-component, .second = y-component
using Tensor2D = std::pair<Eigen::MatrixXd, Eigen::MatrixXd>;
using ImageMat = Eigen::MatrixXd;
using RGBImage = std::array<ImageMat, 3>;

template <typename T>
inline T sqr(T x) { return x * x; }

constexpr double kEps = 1e-12;

enum class FlowType {
    NonLocal,  ///< Eq. (28): u_t = (1/μ₀) Du ∇⊥ Δ⁻¹ div(P⊥)
    Local      ///< Eq. (29): u_t = -(1/μ₀) Du ∇⊥ div(P⊥)
};

enum class PType {
    PureOT = 0,       ///< P = u
    WithComparison = 1 ///< P with image comparison term (Eq. in §7.4)
};

enum class Equalization {
    MassRatio = 0,
    HistEq = 1
};

struct AHTOptions {
    int square_edge_len = 8;
    PType p_type = PType::PureOT;
    FlowType flow = FlowType::NonLocal;
    Equalization equalization = Equalization::MassRatio;
    int max_iterations = 32;
    double convergence_threshold = 0.08;  ///< mean deformation size
    double pure_omt_ratio = 0.3;          ///< α in comparison term
    double cfl_factor = 0.5;              ///< dt scaling vs CFL bound
};

}  // namespace aht_ot
