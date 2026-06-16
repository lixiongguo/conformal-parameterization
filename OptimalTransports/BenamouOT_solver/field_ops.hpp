#pragma once

#include "types.hpp"

namespace aht_ot {

inline double row_sum(const ImageMat& m, int r) {
    double s = 0.0;
    for (int c = 0; c < m.cols(); ++c) s += m(r, c);
    return s;
}

inline double col_sum(const ImageMat& m, int c) {
    double s = 0.0;
    for (int r = 0; r < m.rows(); ++r) s += m(r, c);
    return s;
}

/// perp gradient: (-v_y, v_x), 90 deg CCW rotation
inline Tensor2D rotate_perp(const Tensor2D& v) {
    return {-v.second, v.first};
}

inline ImageMat divergence(const Tensor2D& v) {
    const int s1 = static_cast<int>(v.first.rows());
    const int s2 = static_cast<int>(v.first.cols());
    ImageMat div = ImageMat::Zero(s1, s2);
    for (int i = 1; i < s1 - 1; ++i) {
        for (int j = 1; j < s2 - 1; ++j) {
            const double dv_dx =
                (v.second(i + 1, j) - v.second(i - 1, j)) * 0.5;
            const double dv_dy =
                (v.first(i, j + 1) - v.first(i, j - 1)) * 0.5;
            div(i, j) = dv_dx + dv_dy;
        }
    }
    return div;
}

inline Tensor2D gradient(const ImageMat& f) {
    const int s1 = static_cast<int>(f.rows());
    const int s2 = static_cast<int>(f.cols());
    Tensor2D g;
    g.first = Eigen::MatrixXd::Zero(s1, s2);
    g.second = Eigen::MatrixXd::Zero(s1, s2);
    for (int i = 1; i < s1 - 1; ++i) {
        for (int j = 1; j < s2 - 1; ++j) {
            g.first(i, j) = (f(i + 1, j) - f(i - 1, j)) * 0.5;
            g.second(i, j) = (f(i, j + 1) - f(i, j - 1)) * 0.5;
        }
    }
    return g;
}

/// Jacobian Du (central differences on displacement field)
inline Tensor2D jacobian(const Tensor2D& u) {
    const int s1 = static_cast<int>(u.first.rows());
    const int s2 = static_cast<int>(u.first.cols());
    Tensor2D J;
    J.first = Eigen::MatrixXd::Zero(s1, s2);
    J.second = Eigen::MatrixXd::Zero(s1, s2);
    for (int i = 1; i < s1 - 1; ++i) {
        for (int j = 1; j < s2 - 1; ++j) {
            J.first(i, j) = (u.first(i + 1, j) - u.first(i - 1, j)) * 0.5;
            J.second(i, j) = (u.second(i, j + 1) - u.second(i, j - 1)) * 0.5;
        }
    }
    return J;
}

/// 2D curl: ∂u_y/∂x - ∂u_x/∂y
inline ImageMat curl_scalar(const Tensor2D& u) {
    const int s1 = static_cast<int>(u.first.rows());
    const int s2 = static_cast<int>(u.first.cols());
    ImageMat c = ImageMat::Zero(s1, s2);
    for (int i = 1; i < s1 - 1; ++i) {
        for (int j = 1; j < s2 - 1; ++j) {
            const double duy_dx =
                (u.second(i + 1, j) - u.second(i - 1, j)) * 0.5;
            const double dux_dy =
                (u.first(i, j + 1) - u.first(i, j - 1)) * 0.5;
            c(i, j) = duy_dx - dux_dy;
        }
    }
    return c;
}

inline double mean_absolute_curl(const Tensor2D& u) {
    const ImageMat c = curl_scalar(u);
    const int s1 = static_cast<int>(c.rows());
    const int s2 = static_cast<int>(c.cols());
    double sum = 0.0;
    int cnt = 0;
    for (int i = 1; i < s1 - 1; ++i) {
        for (int j = 1; j < s2 - 1; ++j) {
            sum += std::abs(c(i, j));
            ++cnt;
        }
    }
    return (cnt > 0) ? sum / static_cast<double>(cnt) : 0.0;
}

inline ImageMat safe_density(const ImageMat& mu) {
    ImageMat out = mu;
    for (int i = 0; i < out.rows(); ++i) {
        for (int j = 0; j < out.cols(); ++j) {
            if (out(i, j) < 1.0) out(i, j) = 1.0;
        }
    }
    return out;
}

}  // namespace aht_ot
