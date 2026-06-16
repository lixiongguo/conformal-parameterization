#pragma once

#include "field_ops.hpp"

namespace aht_ot {

/// Sec 7.2 — separable 1D optimal transport along x then y (eq 11-13)
inline Tensor2D compute_initial_mapping(const ImageMat& mu0,
                                        const ImageMat& mu1) {
    const int s1 = static_cast<int>(mu0.rows());
    const int s2 = static_cast<int>(mu0.cols());

    Eigen::VectorXd a(s1);
    double cum0 = 0.0, cum1 = 0.0;
    bool first_correction = false;
    int idx1 = 0;

    for (int idx0 = 0; idx0 < s1; ++idx0) {
        cum0 += row_sum(mu0, idx0);
        const int last = idx1;
        while (cum1 < cum0 && idx1 < s1) {
            cum1 += row_sum(mu1, idx1);
            ++idx1;
        }
        if (idx1 < s1 && (idx1 - last) > 0) {
            --idx1;
            cum1 -= row_sum(mu1, std::max(0, idx1));
        } else {
            if (first_correction)
                cum0 -= row_sum(mu0, idx0);
            else
                first_correction = true;
        }
        a(idx0) = static_cast<double>(idx1 + 1);
    }

    Eigen::VectorXd b(s2);
    cum0 = 0.0;
    cum1 = 0.0;
    first_correction = false;
    idx1 = 0;

    for (int idx0 = 0; idx0 < s2; ++idx0) {
        cum0 += col_sum(mu0, idx0);
        const int last = idx1;
        while (cum1 < cum0 && idx1 < s2) {
            cum1 += col_sum(mu1, idx1);
            ++idx1;
        }
        if (idx1 < s2 && (idx1 - last) > 0) {
            --idx1;
            cum1 -= col_sum(mu1, std::max(0, idx1));
        } else {
            if (first_correction)
                cum0 -= col_sum(mu0, idx0);
            else
                first_correction = true;
        }
        b(idx0) = static_cast<double>(idx1 + 1);
    }

    Tensor2D u0;
    u0.first = Eigen::MatrixXd(s1, s2);
    u0.second = Eigen::MatrixXd(s1, s2);
    for (int i = 0; i < s1; ++i) {
        for (int j = 0; j < s2; ++j) {
            u0.first(i, j) = a(i);
            u0.second(i, j) = b(j);
        }
    }
    return u0;
}

inline ImageMat compute_density_map(const ImageMat& image, int block_size) {
    const int d1 = static_cast<int>(image.rows()) / block_size;
    const int d2 = static_cast<int>(image.cols()) / block_size;
    ImageMat mu(d1, d2);
    for (int i = 0; i < d1; ++i) {
        for (int j = 0; j < d2; ++j) {
            mu(i, j) =
                image
                    .block(i * block_size, j * block_size, block_size, block_size)
                    .sum() /
                static_cast<double>(block_size * block_size);
        }
    }
    return mu;
}

}  // namespace aht_ot
