#pragma once

// Adapter exposing the semi-discrete (MongeAmpere, geogram-backed) solver
// through the shared ot::Solver interface (see ../ot_problem.hpp).
//
// Source: a grid density (continuous measure). Target: weighted Dirac points.
// The transport map sends a source point to the target site whose Laguerre
// (power) cell contains it, i.e. the nearest power vertex.

#include "../ot_problem.hpp"
#include "../geometry/density_grid.hpp"
#include "optimal_transport.hpp"

#include <Eigen/Dense>
#include <limits>
#include <memory>
#include <stdexcept>

namespace MA {

class SemiTransportPlan : public ot::TransportPlan {
public:
    SemiTransportPlan(Eigen::MatrixX2d sites, Eigen::VectorXd weights)
        : sites_(std::move(sites)), weights_(std::move(weights)) {}

    Eigen::Vector2d forward(const Eigen::Vector2d& x) const override {
        int best = -1;
        double best_val = std::numeric_limits<double>::max();
        for (int i = 0; i < sites_.rows(); ++i) {
            const double dx = x.x() - sites_(i, 0);
            const double dy = x.y() - sites_(i, 1);
            const double val = dx * dx + dy * dy - weights_(i);
            if (val < best_val) {
                best_val = val;
                best = i;
            }
        }
        if (best < 0) return x;
        return Eigen::Vector2d(sites_(best, 0), sites_(best, 1));
    }

private:
    Eigen::MatrixX2d sites_;
    Eigen::VectorXd weights_;
};

class SemiOTSolver : public ot::Solver {
public:
    SemiOTSolver(double eps_g = 1e-7, std::size_t max_iter = 100,
                 bool verbose = false)
        : eps_g_(eps_g), max_iter_(max_iter), verbose_(verbose) {}

    ot::Result solve(const ot::Problem& problem) override {
        if (!problem.source.is_grid()) {
            throw std::invalid_argument(
                "SemiOTSolver requires a grid-density source");
        }
        if (!problem.target.is_points()) {
            throw std::invalid_argument(
                "SemiOTSolver requires weighted-point target");
        }

        const ot::Domain& dom = problem.source.domain;
        const otgeo::Box box{dom.min.x(), dom.min.y(), dom.max.x(), dom.max.y()};
        const otgeo::DensityGrid source(problem.source.grid, box);

        const Eigen::MatrixX2d X = problem.target.points;
        Eigen::VectorXd masses = problem.target.masses;

        // Normalize target masses to the source total mass (mass balance).
        const double src_mass = source.total_mass();
        const double tgt_mass = masses.sum();
        if (tgt_mass > 0.0) masses *= (src_mass / tgt_mass);

        Eigen::VectorXd weights = Eigen::VectorXd::Zero(X.rows());
        Statistics stats{0, 0};
        ot_solve(source, X, masses, weights, eps_g_, max_iter_, verbose_, &stats);

        ot::Result out;
        out.plan = std::make_shared<SemiTransportPlan>(X, weights);
        out.iterations = static_cast<int>(stats.niter);
        out.converged = true;
        return out;
    }

private:
    double eps_g_;
    std::size_t max_iter_;
    bool verbose_;
};

}  // namespace MA
