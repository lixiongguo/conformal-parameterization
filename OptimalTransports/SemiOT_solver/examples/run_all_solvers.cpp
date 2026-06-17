// How to call the three OT solvers from outside.
//
// Every backend implements the same ot::Solver interface (see ot_problem.hpp):
//
//     ot::Problem { ot::Distribution source; ot::Distribution target; }
//          |
//          v
//     ot::Result = solver.solve(problem);
//          |
//          v
//     result.plan->forward(x)   // push a source point to the target
//
// So calling code only ever:
//   1. builds two ot::Distribution (grid density or weighted points),
//   2. picks a concrete solver,
//   3. calls solve() and uses result.plan.
//
// The `run()` helper below takes an ot::Solver& and therefore works with ANY
// backend without knowing which one it is.
//
// Build: this example needs all three backends and their dependencies
// (geogram for SemiOT, stb for AHT, surface_mesh for gridOT). It is gated
// behind the CMake option OT_BUILD_ALL_SOLVERS (default OFF). The geogram-only
// `run_semi_ot` example shows the same pattern with just one backend.

#include "ot_problem.hpp"

#include "BenamouOT_solver/ot_adapter.hpp"  // aht_ot::AhtOTSolver
#include "SemiOT_solver/ot_adapter.hpp"     // MA::SemiOTSolver
#include "gridOT_solver/ot_adapter.hpp"     // otmap::GridOTSolver

#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <random>

// A simple non-negative grid density: a Gaussian bump over the unit square.
static Eigen::MatrixXd gaussian_grid(int rows, int cols, double cx, double cy,
                                     double sigma) {
    Eigen::MatrixXd g(rows, cols);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const double x = (c + 0.5) / cols;
            const double y = (r + 0.5) / rows;
            const double d2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
            g(r, c) = std::exp(-d2 / (2.0 * sigma * sigma)) + 1e-3;
        }
    }
    return g;
}

// Generic driver: identical code path for every backend.
static void run(const char* name, ot::Solver& solver,
                const ot::Problem& problem) {
    std::printf("=== %s ===\n", name);
    const ot::Result result = solver.solve(problem);

    const Eigen::Vector2d samples[] = {
        Eigen::Vector2d(0.25, 0.25),
        Eigen::Vector2d(0.50, 0.50),
        Eigen::Vector2d(0.75, 0.75),
    };
    for (const Eigen::Vector2d& s : samples) {
        const Eigen::Vector2d t = result.plan->forward(s);
        std::printf("  T(%.2f, %.2f) -> (%.4f, %.4f)\n", s.x(), s.y(), t.x(),
                    t.y());
    }
    std::printf("  iterations = %d, converged = %d\n\n", result.iterations,
                static_cast<int>(result.converged));
}

int main() {
    const int R = 64, C = 64;

    // (1) BenamouOT / AHT  --  grid density -> grid density
    {
        ot::Problem problem;
        problem.source = ot::Distribution::from_grid(
            gaussian_grid(R, C, 0.35, 0.35, 0.18));
        problem.target = ot::Distribution::from_grid(
            gaussian_grid(R, C, 0.65, 0.65, 0.18));

        aht_ot::AhtOTSolver solver;  // concrete backend
        run("BenamouOT (AHT):  grid -> grid", solver, problem);
    }

    // (2) gridOT  --  implicit uniform source -> grid density target
    {
        ot::Problem problem;
        problem.source = ot::Distribution::uniform(R, C);
        problem.target = ot::Distribution::from_grid(
            gaussian_grid(R, C, 0.60, 0.40, 0.15));

        otmap::GridOTSolver solver;
        run("gridOT:           uniform -> grid", solver, problem);
    }

    // (3) SemiOT  --  grid density source -> weighted Dirac points target
    {
        ot::Problem problem;
        problem.source = ot::Distribution::from_grid(
            gaussian_grid(R, C, 0.50, 0.50, 0.25));

        const int N = 64;
        Eigen::MatrixX2d pts(N, 2);
        Eigen::VectorXd masses(N);
        std::mt19937 rng(1);
        std::uniform_real_distribution<double> U(0.05, 0.95);
        for (int i = 0; i < N; ++i) {
            pts(i, 0) = U(rng);
            pts(i, 1) = U(rng);
            masses(i) = 1.0;
        }
        problem.target = ot::Distribution::from_points(pts, masses);

        MA::SemiOTSolver solver(/*eps_g=*/1e-6, /*max_iter=*/200,
                                /*verbose=*/false);
        run("SemiOT:           grid -> points", solver, problem);
    }

    return 0;
}
