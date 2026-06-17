#pragma once

// FFT-OT style grid-to-grid optimal transport solver.
//
// The implementation follows Lei & Gu's fixed-point formulation:
//
//   phi_{k+1} = Delta^{-1} {
//       sqrt((phi_xx + 1)^2 + (phi_yy + 1)^2 + 2 phi_xy^2
//            + 2 f / (g o (Id + grad phi))) - 2
//   }
//
// with homogeneous Neumann boundary conditions.  The Poisson equation is solved
// by a separable DCT diagonalization of the cell-centered Neumann Laplacian.
// The DCT kernel below is intentionally dependency-free; it can be replaced by
// FFTW/Eigen::FFT without changing the solver interface.

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fftot {

struct FFTOTOptions {
    int max_iterations = 200;
    double tolerance = 1e-8;
    double min_density = 1e-12;
    bool rescale_target_mass = true;
    bool verbose = false;
};

struct FFTOTResult {
    Eigen::MatrixXd potential;  // Kantorovich potential phi on cell centers.
    Eigen::MatrixXd map_x;      // mapped x-coordinate at cell centers.
    Eigen::MatrixXd map_y;      // mapped y-coordinate at cell centers.
    int iterations = 0;
    bool converged = false;
    double final_delta = std::numeric_limits<double>::infinity();
    double transport_cost = 0.0;
};

struct Domain {
    Eigen::Vector2d min{0.0, 0.0};
    Eigen::Vector2d max{1.0, 1.0};

    Eigen::Vector2d extent() const { return max - min; }
};

namespace detail {

inline double clamp(double v, double lo, double hi) {
    return std::max(lo, std::min(v, hi));
}

inline double sample_clamped(const Eigen::MatrixXd& m, int r, int c) {
    r = std::max(0, std::min(r, static_cast<int>(m.rows()) - 1));
    c = std::max(0, std::min(c, static_cast<int>(m.cols()) - 1));
    return m(r, c);
}

inline double bilinear_sample(const Eigen::MatrixXd& grid, const Domain& domain,
                              const Eigen::Vector2d& p) {
    const int rows = static_cast<int>(grid.rows());
    const int cols = static_cast<int>(grid.cols());
    const Eigen::Vector2d ext = domain.extent();
    if (rows <= 0 || cols <= 0 || ext.x() <= 0.0 || ext.y() <= 0.0) {
        throw std::invalid_argument("bilinear_sample: invalid grid or domain");
    }

    const double ux = clamp((p.x() - domain.min.x()) / ext.x(), 0.0, 1.0);
    const double uy = clamp((p.y() - domain.min.y()) / ext.y(), 0.0, 1.0);
    const double fc = clamp(ux * cols - 0.5, 0.0, static_cast<double>(cols - 1));
    const double fr = clamp(uy * rows - 0.5, 0.0, static_cast<double>(rows - 1));

    const int r0 = static_cast<int>(std::floor(fr));
    const int c0 = static_cast<int>(std::floor(fc));
    const int r1 = std::min(r0 + 1, rows - 1);
    const int c1 = std::min(c0 + 1, cols - 1);
    const double tr = fr - r0;
    const double tc = fc - c0;

    const double a = grid(r0, c0) * (1.0 - tc) + grid(r0, c1) * tc;
    const double b = grid(r1, c0) * (1.0 - tc) + grid(r1, c1) * tc;
    return a * (1.0 - tr) + b * tr;
}

inline Eigen::MatrixXd make_dct_matrix(int n) {
    if (n <= 0) {
        throw std::invalid_argument("make_dct_matrix: n must be positive");
    }

    const double pi = std::acos(-1.0);
    Eigen::MatrixXd c(n, n);
    for (int k = 0; k < n; ++k) {
        const double scale = (k == 0) ? std::sqrt(1.0 / n) : std::sqrt(2.0 / n);
        for (int i = 0; i < n; ++i) {
            c(k, i) = scale * std::cos(pi * (static_cast<double>(i) + 0.5) *
                                            static_cast<double>(k) /
                                            static_cast<double>(n));
        }
    }
    return c;
}

class DctPoissonSolver {
public:
    DctPoissonSolver(int rows, int cols, double hx, double hy)
        : rows_(rows), cols_(cols), hx_(hx), hy_(hy),
          dct_rows_(make_dct_matrix(rows)),
          dct_cols_(make_dct_matrix(cols)) {
        if (rows_ <= 0 || cols_ <= 0 || hx_ <= 0.0 || hy_ <= 0.0) {
            throw std::invalid_argument("DctPoissonSolver: invalid dimensions");
        }
    }

    Eigen::MatrixXd solve(Eigen::MatrixXd rhs) const {
        rhs.array() -= rhs.mean();
        Eigen::MatrixXd freq = dct_rows_ * rhs * dct_cols_.transpose();

        const double pi = std::acos(-1.0);
        for (int r = 0; r < rows_; ++r) {
            for (int c = 0; c < cols_; ++c) {
                if (r == 0 && c == 0) {
                    freq(r, c) = 0.0;
                    continue;
                }
                const double lambda =
                    2.0 * (std::cos(pi * static_cast<double>(c) /
                                     static_cast<double>(cols_)) -
                           1.0) /
                        (hx_ * hx_) +
                    2.0 * (std::cos(pi * static_cast<double>(r) /
                                     static_cast<double>(rows_)) -
                           1.0) /
                        (hy_ * hy_);
                freq(r, c) /= lambda;
            }
        }

        Eigen::MatrixXd out = dct_rows_.transpose() * freq * dct_cols_;
        out.array() -= out.mean();
        return out;
    }

private:
    int rows_ = 0;
    int cols_ = 0;
    double hx_ = 1.0;
    double hy_ = 1.0;
    Eigen::MatrixXd dct_rows_;
    Eigen::MatrixXd dct_cols_;
};

}  // namespace detail

class FFTOTSolver {
public:
    explicit FFTOTSolver(FFTOTOptions options = {}) : options_(options) {}

    FFTOTResult solve(const Eigen::MatrixXd& source_density,
                      const Eigen::MatrixXd& target_density,
                      Domain domain = {}) const {
        if (source_density.rows() != target_density.rows() ||
            source_density.cols() != target_density.cols()) {
            throw std::invalid_argument(
                "FFTOTSolver requires source and target grids with the same size");
        }
        if (source_density.rows() < 2 || source_density.cols() < 2) {
            throw std::invalid_argument(
                "FFTOTSolver requires at least a 2x2 density grid");
        }

        const int rows = static_cast<int>(source_density.rows());
        const int cols = static_cast<int>(source_density.cols());
        const Eigen::Vector2d ext = domain.extent();
        if (ext.x() <= 0.0 || ext.y() <= 0.0) {
            throw std::invalid_argument("FFTOTSolver requires a positive domain");
        }

        Eigen::MatrixXd f = positive_density(source_density);
        Eigen::MatrixXd g = positive_density(target_density);
        if (options_.rescale_target_mass) {
            const double target_sum = g.sum();
            if (target_sum <= 0.0) {
                throw std::invalid_argument("target density has zero mass");
            }
            g *= f.sum() / target_sum;
        }

        const double hx = ext.x() / static_cast<double>(cols);
        const double hy = ext.y() / static_cast<double>(rows);
        detail::DctPoissonSolver poisson(rows, cols, hx, hy);

        Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(rows, cols);
        Eigen::MatrixXd next_phi = phi;
        Eigen::MatrixXd map_x(rows, cols);
        Eigen::MatrixXd map_y(rows, cols);
        double delta = std::numeric_limits<double>::infinity();
        bool converged = false;
        int iter = 0;

        for (; iter < options_.max_iterations; ++iter) {
            Eigen::MatrixXd rhs(rows, cols);

            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    const Derivatives d = derivatives(phi, r, c, hx, hy);
                    const Eigen::Vector2d x = cell_center(r, c, rows, cols, domain);
                    const Eigen::Vector2d y = clamp_to_domain(
                        Eigen::Vector2d(x.x() + d.phi_x, x.y() + d.phi_y),
                        domain);

                    map_x(r, c) = y.x();
                    map_y(r, c) = y.y();

                    const double gy = std::max(
                        detail::bilinear_sample(g, domain, y), options_.min_density);
                    const double radicand =
                        (d.phi_xx + 1.0) * (d.phi_xx + 1.0) +
                        (d.phi_yy + 1.0) * (d.phi_yy + 1.0) +
                        2.0 * d.phi_xy * d.phi_xy +
                        2.0 * f(r, c) / gy;
                    rhs(r, c) = std::sqrt(std::max(0.0, radicand)) - 2.0;
                }
            }

            next_phi = poisson.solve(rhs);
            delta = (next_phi - phi).norm() /
                    std::sqrt(static_cast<double>(rows * cols));
            phi.swap(next_phi);

            if (delta < options_.tolerance) {
                converged = true;
                ++iter;
                break;
            }
        }

        fill_map(phi, rows, cols, hx, hy, domain, map_x, map_y);

        FFTOTResult result;
        result.potential = std::move(phi);
        result.map_x = std::move(map_x);
        result.map_y = std::move(map_y);
        result.iterations = iter;
        result.converged = converged;
        result.final_delta = delta;
        result.transport_cost = estimate_cost(f, result.map_x, result.map_y, domain);
        return result;
    }

private:
    struct Derivatives {
        double phi_x = 0.0;
        double phi_y = 0.0;
        double phi_xx = 0.0;
        double phi_yy = 0.0;
        double phi_xy = 0.0;
    };

    Eigen::MatrixXd positive_density(const Eigen::MatrixXd& density) const {
        if (!density.allFinite()) {
            throw std::invalid_argument("density contains non-finite values");
        }
        Eigen::MatrixXd out = density.cwiseMax(options_.min_density);
        if (out.sum() <= 0.0) {
            throw std::invalid_argument("density has zero mass");
        }
        return out;
    }

    static Eigen::Vector2d cell_center(int r, int c, int rows, int cols,
                                       const Domain& domain) {
        const Eigen::Vector2d ext = domain.extent();
        return Eigen::Vector2d(
            domain.min.x() + (static_cast<double>(c) + 0.5) * ext.x() /
                                 static_cast<double>(cols),
            domain.min.y() + (static_cast<double>(r) + 0.5) * ext.y() /
                                 static_cast<double>(rows));
    }

    static Eigen::Vector2d clamp_to_domain(const Eigen::Vector2d& p,
                                           const Domain& domain) {
        return Eigen::Vector2d(
            detail::clamp(p.x(), domain.min.x(), domain.max.x()),
            detail::clamp(p.y(), domain.min.y(), domain.max.y()));
    }

    static Derivatives derivatives(const Eigen::MatrixXd& phi, int r, int c,
                                   double hx, double hy) {
        const double cc = detail::sample_clamped(phi, r, c);
        const double xp = detail::sample_clamped(phi, r, c + 1);
        const double xm = detail::sample_clamped(phi, r, c - 1);
        const double yp = detail::sample_clamped(phi, r + 1, c);
        const double ym = detail::sample_clamped(phi, r - 1, c);

        Derivatives d;
        d.phi_x = (xp - xm) / (2.0 * hx);
        d.phi_y = (yp - ym) / (2.0 * hy);
        d.phi_xx = (xp + xm - 2.0 * cc) / (hx * hx);
        d.phi_yy = (yp + ym - 2.0 * cc) / (hy * hy);
        d.phi_xy = (detail::sample_clamped(phi, r + 1, c + 1) +
                    detail::sample_clamped(phi, r - 1, c - 1) -
                    detail::sample_clamped(phi, r + 1, c - 1) -
                    detail::sample_clamped(phi, r - 1, c + 1)) /
                   (4.0 * hx * hy);
        return d;
    }

    static void fill_map(const Eigen::MatrixXd& phi, int rows, int cols,
                         double hx, double hy, const Domain& domain,
                         Eigen::MatrixXd& map_x, Eigen::MatrixXd& map_y) {
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                const Derivatives d = derivatives(phi, r, c, hx, hy);
                const Eigen::Vector2d x = cell_center(r, c, rows, cols, domain);
                const Eigen::Vector2d y = clamp_to_domain(
                    Eigen::Vector2d(x.x() + d.phi_x, x.y() + d.phi_y), domain);
                map_x(r, c) = y.x();
                map_y(r, c) = y.y();
            }
        }
    }

    static double estimate_cost(const Eigen::MatrixXd& f,
                                const Eigen::MatrixXd& map_x,
                                const Eigen::MatrixXd& map_y,
                                const Domain& domain) {
        const int rows = static_cast<int>(f.rows());
        const int cols = static_cast<int>(f.cols());
        const double cell_area =
            domain.extent().x() * domain.extent().y() /
            static_cast<double>(rows * cols);
        double cost = 0.0;
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                const Eigen::Vector2d x = cell_center(r, c, rows, cols, domain);
                const double dx = map_x(r, c) - x.x();
                const double dy = map_y(r, c) - x.y();
                cost += f(r, c) * (dx * dx + dy * dy) * cell_area;
            }
        }
        return cost;
    }

    FFTOTOptions options_;
};

}  // namespace fftot
