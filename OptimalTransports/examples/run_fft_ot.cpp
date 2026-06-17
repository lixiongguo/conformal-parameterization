// End-to-end demo for FFT-OT: grid density -> grid density on the unit square.

#include "FFTOT_solver/ot_adapter.hpp"
#include "ot_problem.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <memory>

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

int main() {
    const int rows = 32;
    const int cols = 32;

    ot::Problem problem;
    problem.source =
        ot::Distribution::from_grid(gaussian_grid(rows, cols, 0.35, 0.35, 0.16));
    problem.target =
        ot::Distribution::from_grid(gaussian_grid(rows, cols, 0.65, 0.65, 0.16));

    fftot::FFTOTOptions opt;
    opt.max_iterations = 200;
    opt.tolerance = 1e-5;

    fftot::FFTOTSolverAdapter solver(opt);
    const ot::Result result = solver.solve(problem);

    const Eigen::Vector2d samples[] = {
        Eigen::Vector2d(0.25, 0.25),
        Eigen::Vector2d(0.50, 0.50),
        Eigen::Vector2d(0.75, 0.75),
    };
    for (const Eigen::Vector2d& s : samples) {
        const Eigen::Vector2d t = result.plan->forward(s);
        std::printf("T(%.2f, %.2f) -> (%.4f, %.4f)\n", s.x(), s.y(), t.x(),
                    t.y());
    }

    const auto fft_plan =
        std::dynamic_pointer_cast<fftot::FFTOTTransportPlan>(result.plan);
    const double final_delta =
        fft_plan ? fft_plan->native_result().final_delta : 0.0;
    std::printf(
        "iterations = %d ; converged = %d ; final_delta = %.6e ; W2 cost ~= %.8f\n",
        result.iterations, static_cast<int>(result.converged), final_delta,
        result.transport_cost);
    return 0;
}
