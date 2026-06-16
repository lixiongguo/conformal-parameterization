#pragma once

// Common input/output interface shared by the three optimal-transport solvers
// in this folder:
//
//   * BenamouOT_solver  (namespace aht_ot)  grid density  -> grid density
//   * gridOT_solver     (namespace otmap )  uniform grid  -> grid density
//   * SemiOT_solver     (namespace MA    )  grid density  -> weighted points
//
// Every solver computes an optimal transport between a SOURCE distribution mu0
// and a TARGET distribution mu1. This header factors out the pieces they have
// in common:
//
//   input  : two distributions   -> ot::Problem { source, target }
//   output : a transport plan     -> ot::Result  { plan, cost }
//
// The header only depends on Eigen so it stays light; each solver provides a
// thin adapter (subclass of ot::Solver) in its own directory that converts to
// and from its native types.

#include <Eigen/Dense>

#include <memory>
#include <stdexcept>

namespace ot {

// Axis-aligned rectangular domain that anchors both grids and point sets in a
// shared coordinate system. Defaults to the unit square [0,1] x [0,1].
struct Domain {
    Eigen::Vector2d min{0.0, 0.0};
    Eigen::Vector2d max{1.0, 1.0};

    Eigen::Vector2d extent() const { return max - min; }

    // Map normalized coordinates in [0,1]^2 to world coordinates.
    Eigen::Vector2d from_unit(const Eigen::Vector2d& u) const {
        return min + u.cwiseProduct(extent());
    }
};

// How a distribution is represented. The two forms together cover the inputs of
// all three solvers.
enum class Kind {
    GridDensity,    // dense non-negative density sampled on a regular grid
    WeightedPoints  // discrete atoms: positions + masses (e.g. Dirac mixture)
};

// A probability/mass distribution over a 2D domain.
//
//   * GridDensity   : `grid` holds per-cell density values (rows x cols),
//                     laid out over `domain`. Used as both source and target by
//                     BenamouOT, and as the target density by gridOT (its source
//                     is the implicit uniform grid, see Distribution::uniform).
//   * WeightedPoints: `points` (N x 2) are atom positions in `domain`, `masses`
//                     (N) their weights. Used as the target by SemiOT.
struct Distribution {
    Kind kind = Kind::GridDensity;
    Domain domain{};

    // GridDensity payload.
    Eigen::MatrixXd grid;

    // WeightedPoints payload.
    Eigen::MatrixX2d points;
    Eigen::VectorXd masses;

    static Distribution from_grid(Eigen::MatrixXd density, Domain domain = {}) {
        Distribution d;
        d.kind = Kind::GridDensity;
        d.grid = std::move(density);
        d.domain = domain;
        return d;
    }

    static Distribution from_points(Eigen::MatrixX2d positions,
                                    Eigen::VectorXd weights,
                                    Domain domain = {}) {
        if (positions.rows() != weights.rows()) {
            throw std::invalid_argument(
                "Distribution::from_points: positions and weights size mismatch");
        }
        Distribution d;
        d.kind = Kind::WeightedPoints;
        d.points = std::move(positions);
        d.masses = std::move(weights);
        d.domain = domain;
        return d;
    }

    // Uniform density on a rows x cols grid (the implicit source of gridOT).
    static Distribution uniform(int rows, int cols, Domain domain = {}) {
        return from_grid(Eigen::MatrixXd::Ones(rows, cols), domain);
    }

    bool is_grid() const { return kind == Kind::GridDensity; }
    bool is_points() const { return kind == Kind::WeightedPoints; }

    // Total mass (integral). For a grid this is the plain sum of cell values;
    // for points it is the sum of masses.
    double total_mass() const {
        return is_grid() ? grid.sum() : masses.sum();
    }
};

// Two distributions: this is the common input of every solver.
struct Problem {
    Distribution source;  // mu0
    Distribution target;  // mu1
};

// The common output: a transport plan that maps the source domain onto the
// target domain. `forward` pushes a source point to where its mass lands in the
// target; `inverse` pulls a target point back. Solvers that only produce one
// direction may leave the other unimplemented (throwing).
class TransportPlan {
public:
    virtual ~TransportPlan() = default;

    virtual Eigen::Vector2d forward(const Eigen::Vector2d& x) const = 0;

    virtual Eigen::Vector2d inverse(const Eigen::Vector2d& /*y*/) const {
        throw std::logic_error("TransportPlan::inverse not available");
    }
};

// Result returned by every solver: the plan plus a few common diagnostics.
struct Result {
    std::shared_ptr<TransportPlan> plan;
    double transport_cost = 0.0;   // total (squared) Wasserstein-2 cost, if known
    int iterations = 0;
    bool converged = false;
};

// Abstract solver: takes two distributions, returns a transport plan.
// Each backend (aht_ot / otmap / MA) supplies a subclass adapter.
class Solver {
public:
    virtual ~Solver() = default;
    virtual Result solve(const Problem& problem) = 0;
};

}  // namespace ot
