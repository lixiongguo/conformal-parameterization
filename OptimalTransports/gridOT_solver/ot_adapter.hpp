#pragma once

// Adapter exposing the grid-based solver (otmap) through the shared ot::Solver
// interface (see ../ot_problem.hpp).
//
// The grid solver transports the implicit UNIFORM density on a regular grid to
// a target grid density, so the source distribution must be uniform and the
// target a (square) grid density. The resulting otmap::TransportMap provides
// forward/inverse maps on the unit square, which we rescale to the requested
// domain.

#include "../ot_problem.hpp"
#include "otsolver_2dgrid.h"
#include "transport_map.h"

#include <Eigen/Dense>
#include <memory>
#include <stdexcept>

namespace otmap {

class GridTransportPlan : public ot::TransportPlan {
public:
    GridTransportPlan(std::shared_ptr<TransportMap> tmap, ot::Domain domain)
        : tmap_(std::move(tmap)), domain_(domain) {
        tmap_->init_forward();
        tmap_->init_inverse();
    }

    Eigen::Vector2d forward(const Eigen::Vector2d& x) const override {
        return from_unit(tmap_->fwd(to_unit(x)));
    }

    Eigen::Vector2d inverse(const Eigen::Vector2d& y) const override {
        return from_unit(tmap_->inv(to_unit(y)));
    }

private:
    Eigen::Vector2d to_unit(const Eigen::Vector2d& p) const {
        const Eigen::Vector2d e = domain_.extent();
        return Eigen::Vector2d((p.x() - domain_.min.x()) / e.x(),
                               (p.y() - domain_.min.y()) / e.y());
    }
    Eigen::Vector2d from_unit(const Eigen::Vector2d& u) const {
        const Eigen::Vector2d e = domain_.extent();
        return Eigen::Vector2d(domain_.min.x() + u.x() * e.x(),
                               domain_.min.y() + u.y() * e.y());
    }

    std::shared_ptr<TransportMap> tmap_;
    ot::Domain domain_;
};

class GridOTSolver : public ot::Solver {
public:
    explicit GridOTSolver(SolverOptions options = SolverOptions())
        : options_(options) {}

    ot::Result solve(const ot::Problem& problem) override {
        if (!problem.target.is_grid()) {
            throw std::invalid_argument("GridOTSolver requires a grid target");
        }
        const Eigen::MatrixXd& g = problem.target.grid;
        const int rows = static_cast<int>(g.rows());
        const int cols = static_cast<int>(g.cols());
        if (rows != cols) {
            throw std::invalid_argument(
                "GridOTSolver requires a square target grid");
        }
        const int n = rows;

        // Flatten target into the solver's face ordering: id = j + i*gridSize.
        Eigen::VectorXd density(n * n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) density(j + i * n) = g(i, j);

        GridBasedTransportSolver solver;
        solver.set_verbose_level(0);
        solver.init(n);
        TransportMap tmap = solver.solve(density, options_);

        ot::Result out;
        out.plan = std::make_shared<GridTransportPlan>(
            std::make_shared<TransportMap>(tmap), problem.target.domain);
        out.converged = true;
        return out;
    }

private:
    SolverOptions options_;
};

}  // namespace otmap
