#pragma once

// Adapter exposing FFT-OT through the shared ot::Solver interface.

#include "../ot_problem.hpp"
#include "fft_ot.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>

namespace fftot {

class FFTOTTransportPlan : public ot::TransportPlan {
public:
    FFTOTTransportPlan(FFTOTResult result, ot::Domain domain)
        : result_(std::move(result)), domain_(domain) {
        rows_ = static_cast<int>(result_.map_x.rows());
        cols_ = static_cast<int>(result_.map_x.cols());
        if (rows_ <= 0 || cols_ <= 0 || result_.map_y.rows() != rows_ ||
            result_.map_y.cols() != cols_) {
            throw std::invalid_argument("FFTOTTransportPlan: invalid map grid");
        }
    }

    Eigen::Vector2d forward(const Eigen::Vector2d& x) const override {
        return Eigen::Vector2d(sample(result_.map_x, x), sample(result_.map_y, x));
    }

    const FFTOTResult& native_result() const { return result_; }

private:
    double sample(const Eigen::MatrixXd& grid, const Eigen::Vector2d& x) const {
        const Eigen::Vector2d ext = domain_.extent();
        const double ux =
            std::clamp((x.x() - domain_.min.x()) / ext.x(), 0.0, 1.0);
        const double uy =
            std::clamp((x.y() - domain_.min.y()) / ext.y(), 0.0, 1.0);
        const double fc =
            std::clamp(ux * cols_ - 0.5, 0.0, static_cast<double>(cols_ - 1));
        const double fr =
            std::clamp(uy * rows_ - 0.5, 0.0, static_cast<double>(rows_ - 1));

        const int r0 = static_cast<int>(std::floor(fr));
        const int c0 = static_cast<int>(std::floor(fc));
        const int r1 = std::min(r0 + 1, rows_ - 1);
        const int c1 = std::min(c0 + 1, cols_ - 1);
        const double tr = fr - r0;
        const double tc = fc - c0;

        const double a = grid(r0, c0) * (1.0 - tc) + grid(r0, c1) * tc;
        const double b = grid(r1, c0) * (1.0 - tc) + grid(r1, c1) * tc;
        return a * (1.0 - tr) + b * tr;
    }

    FFTOTResult result_;
    ot::Domain domain_;
    int rows_ = 0;
    int cols_ = 0;
};

class FFTOTSolverAdapter : public ot::Solver {
public:
    explicit FFTOTSolverAdapter(FFTOTOptions options = {}) : options_(options) {}

    ot::Result solve(const ot::Problem& problem) override {
        if (!problem.source.is_grid() || !problem.target.is_grid()) {
            throw std::invalid_argument(
                "FFTOTSolver requires grid source and target densities");
        }
        if ((problem.source.domain.min - problem.target.domain.min).norm() > 1e-12 ||
            (problem.source.domain.max - problem.target.domain.max).norm() > 1e-12) {
            throw std::invalid_argument(
                "FFTOTSolver requires source and target on the same domain");
        }

        Domain domain;
        domain.min = problem.source.domain.min;
        domain.max = problem.source.domain.max;

        FFTOTSolver solver(options_);
        FFTOTResult native =
            solver.solve(problem.source.grid, problem.target.grid, domain);

        ot::Result out;
        out.iterations = native.iterations;
        out.converged = native.converged;
        out.transport_cost = native.transport_cost;
        out.plan = std::make_shared<FFTOTTransportPlan>(std::move(native),
                                                        problem.source.domain);
        return out;
    }

private:
    FFTOTOptions options_;
};

}  // namespace fftot
