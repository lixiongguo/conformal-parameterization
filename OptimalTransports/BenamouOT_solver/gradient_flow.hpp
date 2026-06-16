#pragma once

#include <cstdio>
#include <utility>

#include "field_ops.hpp"
#include "initial_mapping.hpp"
#include "poisson.hpp"

namespace aht_ot {

/// Sec 7.4 comparison term P (density-grid version)
inline Tensor2D compute_P(const Tensor2D& u, const ImageMat& mu0,
                          const ImageMat& mu1, double alpha) {
    const int s1 = static_cast<int>(mu0.rows());
    const int s2 = static_cast<int>(mu0.cols());

    ImageMat mu0x = ImageMat::Zero(s1, s2);
    ImageMat mu0y = ImageMat::Zero(s1, s2);
    for (int i = 1; i < s1 - 1; ++i) {
        for (int j = 1; j < s2 - 1; ++j) {
            mu0x(i, j) = (mu0(i + 1, j) - mu0(i - 1, j)) * 0.5;
            mu0y(i, j) = (mu0(i, j + 1) - mu0(i, j - 1)) * 0.5;
        }
    }

    const ImageMat mu0_safe = safe_density(mu0);
    Tensor2D P;
    P.first = Eigen::MatrixXd::Zero(s1, s2);
    P.second = Eigen::MatrixXd::Zero(s1, s2);

    for (int i = 0; i < s1; ++i) {
        for (int j = 0; j < s2; ++j) {
            const int xl = std::max(
                0, std::min(s1 - 1, static_cast<int>(std::round(u.first(i, j))) - 1));
            const int yl = std::max(
                0, std::min(s2 - 1, static_cast<int>(std::round(u.second(i, j))) - 1));
            const double diff = sqr(mu0_safe(i, j) - mu1(xl, yl));
            const double inv_mu = 1.0 / mu0_safe(i, j);
            const double inv_mu2 = inv_mu * inv_mu;
            P.first(i, j) = inv_mu2 * diff * mu0x(i, j) +
                            2.0 * inv_mu * diff * mu0x(i, j) +
                            2.0 * sqr(alpha) * u.first(i, j);
            P.second(i, j) = inv_mu2 * diff * mu0y(i, j) +
                             2.0 * inv_mu * diff * mu0y(i, j) +
                             2.0 * sqr(alpha) * u.second(i, j);
        }
    }
    return P;
}

/// u_t = (1/μ₀) Du ∇^⊥ φ,  φ from non-local or local flow; returns (u_t, CFL dt bound)
inline std::pair<Tensor2D, double> compute_velocity(
    const Tensor2D& u, const ImageMat& mu0, const Tensor2D& P_in,
    int square_edge_len, FlowType flow) {
    Tensor2D P = rotate_perp(P_in);
    const ImageMat div_P = divergence(P);

    const int n1 = static_cast<int>(div_P.rows());
    const int n2 = static_cast<int>(div_P.cols());

    Tensor2D grad_phi;
    if (flow == FlowType::NonLocal) {
        Eigen::VectorXd rhs(n1 * n2);
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                rhs(i * n2 + j) = -div_P(i, j);
            }
        }
        PoissonSolver5pt poisson(n1, n2, static_cast<double>(square_edge_len));
        const Eigen::VectorXd f_vec = poisson.solve(rhs);
        ImageMat f(n1, n2);
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                f(i, j) = f_vec(i * n2 + j);
            }
        }
        grad_phi = rotate_perp(gradient(f));
    } else {
        grad_phi = rotate_perp(gradient(div_P));
        grad_phi.first *= -1.0;
        grad_phi.second *= -1.0;
    }

    const Tensor2D Du = jacobian(u);
    const ImageMat mu0_safe = safe_density(mu0);

    Tensor2D ut;
    ut.first = Eigen::MatrixXd::Zero(n1, n2);
    ut.second = Eigen::MatrixXd::Zero(n1, n2);
    Tensor2D inv_grad;
    inv_grad.first = Eigen::MatrixXd::Zero(n1, n2);
    inv_grad.second = Eigen::MatrixXd::Zero(n1, n2);

    for (int i = 1; i < n1 - 1; ++i) {
        for (int j = 1; j < n2 - 1; ++j) {
            const double inv_mu = 1.0 / mu0_safe(i, j);
            ut.first(i, j) = inv_mu * Du.first(i, j) * grad_phi.first(i, j);
            ut.second(i, j) = inv_mu * Du.second(i, j) * grad_phi.second(i, j);
            inv_grad.first(i, j) = inv_mu * grad_phi.first(i, j);
            inv_grad.second(i, j) = inv_mu * grad_phi.second(i, j);
        }
    }

    double min_inv = std::numeric_limits<double>::max();
    for (int i = 1; i < n1 - 1; ++i) {
        for (int j = 1; j < n2 - 1; ++j) {
            const double a1 = std::abs(inv_grad.first(i, j));
            const double a2 = std::abs(inv_grad.second(i, j));
            if (a1 > kEps) min_inv = std::min(min_inv, 1.0 / a1);
            if (a2 > kEps) min_inv = std::min(min_inv, 1.0 / a2);
        }
    }
    const double dt_bound =
        (min_inv < std::numeric_limits<double>::max()) ? min_inv : 1.0;
    return {ut, dt_bound};
}

inline double mean_deformation_size(const Tensor2D& a, const Tensor2D& b) {
    const int s1 = static_cast<int>(a.first.rows());
    const int s2 = static_cast<int>(a.first.cols());
    double sum = 0.0;
    int cnt = 0;
    for (int i = 1; i < s1 - 1; ++i) {
        for (int j = 1; j < s2 - 1; ++j) {
            const double dx = a.first(i, j) - b.first(i, j);
            const double dy = a.second(i, j) - b.second(i, j);
            sum += std::sqrt(dx * dx + dy * dy);
            ++cnt;
        }
    }
    return (cnt > 0) ? sum / static_cast<double>(cnt) : 0.0;
}

/// Sec 7.4 — polar decomposition + gradient descent loop
inline Tensor2D gradient_descent(Tensor2D u0_abs, const ImageMat& mu0,
                                 const ImageMat& mu1, const AHTOptions& opt) {
    const int s1 = static_cast<int>(u0_abs.first.rows());
    const int s2 = static_cast<int>(u0_abs.first.cols());

    Tensor2D u;
    u.first = Eigen::MatrixXd(s1, s2);
    u.second = Eigen::MatrixXd(s1, s2);
    for (int i = 0; i < s1; ++i) {
        for (int j = 0; j < s2; ++j) {
            u.first(i, j) = u0_abs.first(i, j) - static_cast<double>(i + 1);
            u.second(i, j) = u0_abs.second(i, j) - static_cast<double>(j + 1);
        }
    }

    Tensor2D last_ut;
    last_ut.first = Eigen::MatrixXd::Zero(s1, s2);
    last_ut.second = Eigen::MatrixXd::Zero(s1, s2);

    for (int iter = 0; iter < opt.max_iterations; ++iter) {
        Tensor2D P = (opt.p_type == PType::PureOT)
                         ? u
                         : compute_P(u, mu0, mu1, opt.pure_omt_ratio);

        const auto [ut, dt_bound] = compute_velocity(
            u, mu0, P, opt.square_edge_len, opt.flow);
        const double alpha = opt.cfl_factor * dt_bound;

        for (int i = 0; i < s1; ++i) {
            for (int j = 0; j < s2; ++j) {
                u.first(i, j) -= alpha * ut.first(i, j);
                u.second(i, j) -= alpha * ut.second(i, j);
            }
        }

        const double deform = mean_deformation_size(ut, last_ut);
        const double curl = mean_absolute_curl(u);
        std::printf("  iter %3d  |ut|=%.6f  mean|curl|=%.6f  dt=%.4f\n", iter,
                    deform, curl, alpha);

        if (deform < opt.convergence_threshold) {
            std::printf("  converged (deformation below threshold)\n");
            break;
        }
        last_ut = ut;
    }

    Tensor2D u_abs;
    u_abs.first = Eigen::MatrixXd(s1, s2);
    u_abs.second = Eigen::MatrixXd(s1, s2);
    for (int i = 0; i < s1; ++i) {
        for (int j = 0; j < s2; ++j) {
            u_abs.first(i, j) = static_cast<double>(i + 1) + u.first(i, j);
            u_abs.second(i, j) = static_cast<double>(j + 1) + u.second(i, j);
        }
    }
    return u_abs;
}

}  // namespace aht_ot
