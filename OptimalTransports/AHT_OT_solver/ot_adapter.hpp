#pragma once

// Adapter exposing the AHT (Benamou-style) grid-to-grid solver through the
// shared ot::Solver interface (see ../ot_problem.hpp).

#include "../ot_problem.hpp"
#include "aht_ot.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace aht_ot {

// Transport plan backed by the dense AHT forward map (Tensor2D), which stores
// for each density-grid cell its mapped (row,col) position in 1-based grid
// index coordinates.
class AhtTransportPlan : public ot::TransportPlan {
public:
    AhtTransportPlan(Tensor2D map, ot::Domain domain)
        : map_(std::move(map)), domain_(domain) {
        rows_ = static_cast<int>(map_.first.rows());
        cols_ = static_cast<int>(map_.first.cols());
    }

    Eigen::Vector2d forward(const Eigen::Vector2d& x) const override {
        const Eigen::Vector2d ext = domain_.extent();
        const double ux = (x.x() - domain_.min.x()) / ext.x();
        const double uy = (x.y() - domain_.min.y()) / ext.y();

        // cell-centred continuous grid coordinates
        double fc = std::clamp(ux * cols_ - 0.5, 0.0, static_cast<double>(cols_ - 1));
        double fr = std::clamp(uy * rows_ - 0.5, 0.0, static_cast<double>(rows_ - 1));

        const double mapped_row = sample(map_.first, fr, fc);   // 1-based row
        const double mapped_col = sample(map_.second, fr, fc);  // 1-based col

        const double vy = (mapped_row - 0.5) / rows_;
        const double vx = (mapped_col - 0.5) / cols_;
        return Eigen::Vector2d(domain_.min.x() + vx * ext.x(),
                               domain_.min.y() + vy * ext.y());
    }

private:
    static double sample(const Eigen::MatrixXd& m, double fr, double fc) {
        const int r0 = static_cast<int>(std::floor(fr));
        const int c0 = static_cast<int>(std::floor(fc));
        const int r1 = std::min(r0 + 1, static_cast<int>(m.rows()) - 1);
        const int c1 = std::min(c0 + 1, static_cast<int>(m.cols()) - 1);
        const double tr = fr - r0;
        const double tc = fc - c0;
        const double a = m(r0, c0) * (1 - tc) + m(r0, c1) * tc;
        const double b = m(r1, c0) * (1 - tc) + m(r1, c1) * tc;
        return a * (1 - tr) + b * tr;
    }

    Tensor2D map_;
    ot::Domain domain_;
    int rows_ = 0;
    int cols_ = 0;
};

class AhtOTSolver : public ot::Solver {
public:
    explicit AhtOTSolver(AHTOptions options = {}) : options_(options) {}

    ot::Result solve(const ot::Problem& problem) override {
        if (!problem.source.is_grid() || !problem.target.is_grid()) {
            throw std::invalid_argument(
                "AhtOTSolver requires both source and target to be grid densities");
        }
        const ImageMat img0 = problem.source.grid;
        const ImageMat img1 = problem.target.grid;

        const AHTResult res = aht_ot::solve(img0, img1, options_);

        ot::Result out;
        out.plan = std::make_shared<AhtTransportPlan>(res.optimal_map,
                                                      problem.source.domain);
        out.converged = true;
        return out;
    }

private:
    AHTOptions options_;
};

}  // namespace aht_ot
