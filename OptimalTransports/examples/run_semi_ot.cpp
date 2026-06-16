// End-to-end demo of the common OT interface driving the geogram-backed
// semi-discrete solver: a uniform grid source transported to random weighted
// Dirac points.

#include "ot_problem.hpp"
#include "SemiOT_solver/ot_adapter.hpp"

#include <Eigen/Dense>
#include <cstdio>
#include <random>

int main() {
    // Source: uniform density on a grid over the unit square.
    const int R = 64, C = 64;
    Eigen::MatrixXd grid = Eigen::MatrixXd::Ones(R, C);
    ot::Distribution source = ot::Distribution::from_grid(grid);

    // Target: N random weighted Dirac points in the unit square.
    const int N = 50;
    Eigen::MatrixX2d pts(N, 2);
    Eigen::VectorXd masses(N);
    std::mt19937 rng(0);
    std::uniform_real_distribution<double> U(0.05, 0.95);
    for (int i = 0; i < N; ++i) {
        pts(i, 0) = U(rng);
        pts(i, 1) = U(rng);
        masses(i) = 1.0;
    }
    ot::Distribution target = ot::Distribution::from_points(pts, masses);

    ot::Problem problem{source, target};
    MA::SemiOTSolver solver(/*eps_g=*/1e-6, /*max_iter=*/200, /*verbose=*/true);
    ot::Result result = solver.solve(problem);

    const Eigen::Vector2d q(0.5, 0.5);
    const Eigen::Vector2d t = result.plan->forward(q);
    std::printf("iterations = %d ; T(0.5,0.5) -> (%.4f, %.4f)\n",
                result.iterations, t.x(), t.y());
    return 0;
}
