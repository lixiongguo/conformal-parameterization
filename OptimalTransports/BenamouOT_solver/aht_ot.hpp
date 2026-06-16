#pragma once

/// AHT optimal mass transport — Haker, Angenent, Tannenbaum (IJCV 2004)
/// Sections 7.2-7.4: initial MP map + gradient descent curl removal

#include <cstddef>
#include <cstdio>
#include <string>

#include "gradient_flow.hpp"
#include "image_ops.hpp"
#include "initial_mapping.hpp"

namespace aht_ot {

struct AHTResult {
    Tensor2D initial_map;
    Tensor2D optimal_map;
    ImageMat density0;
    ImageMat density1;
    double final_mean_curl = 0.0;
};

/// Full pipeline: equalize → density maps → initial MP map → gradient descent
inline AHTResult solve(const ImageMat& image0, const ImageMat& image1,
                       AHTOptions opt = {}) {
    if (image0.rows() != image1.rows() || image0.cols() != image1.cols()) {
        throw std::runtime_error("Images must have the same size");
    }

    ImageMat img0 = image0;
    ImageMat img1 = image1;
    equalize_masses(img0, img1, opt.equalization);

    const ImageMat mu0 = compute_density_map(img0, opt.square_edge_len);
    const ImageMat mu1 = compute_density_map(img1, opt.square_edge_len);

    std::printf("Density grid: %td x %td (block=%d)\n",
                static_cast<ptrdiff_t>(mu0.rows()),
                static_cast<ptrdiff_t>(mu0.cols()), opt.square_edge_len);

    const Tensor2D u0 = compute_initial_mapping(mu0, mu1);
    std::printf("Initial MP mapping (sec 7.2) done.\n");
    std::printf("Gradient descent (sec 7.4, flow=%s)...\n",
                opt.flow == FlowType::NonLocal ? "non-local" : "local");

    const Tensor2D u_opt =
        gradient_descent(u0, mu0, mu1, opt);

    AHTResult result;
    result.initial_map = u0;
    result.optimal_map = u_opt;
    result.density0 = mu0;
    result.density1 = mu1;
    Tensor2D disp;
    disp.first = Eigen::MatrixXd(u_opt.first.rows(), u_opt.first.cols());
    disp.second = Eigen::MatrixXd(u_opt.second.rows(), u_opt.second.cols());
    for (int i = 0; i < disp.first.rows(); ++i) {
        for (int j = 0; j < disp.first.cols(); ++j) {
            disp.first(i, j) = u_opt.first(i, j) - static_cast<double>(i + 1);
            disp.second(i, j) = u_opt.second(i, j) - static_cast<double>(j + 1);
        }
    }
    result.final_mean_curl = mean_absolute_curl(disp);
    return result;
}

inline AHTResult solve_files(const std::string& path0, const std::string& path1,
                             AHTOptions opt = {}) {
    return solve(imread_gray(path0), imread_gray(path1), opt);
}

}  // namespace aht_ot
