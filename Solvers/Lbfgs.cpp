#include "Lbfgs.hpp"

namespace p2p_harmonic {

MatC lbfgsIter(const MatX& invH0, const MatC& g, const MatC& x) {
    const int k = static_cast<int>(x.cols());
    MatC s = x.leftCols(k - 1) - x.rightCols(k - 1);
    MatC y = g.leftCols(k - 1) - g.rightCols(k - 1);

    Eigen::VectorXd rho(k - 1);
    for (int i = 0; i < k - 1; ++i) {
        const double denom = std::real(s.col(i).dot(y.col(i)));
        rho[i] = denom != 0.0 ? 1.0 / denom : 0.0;
    }

    VecC q = -g.col(0);
    Eigen::VectorXd alpha = Eigen::VectorXd::Zero(k - 1);
    for (int i = 0; i < k - 1; ++i) {
        if (rho[i] == 0.0) break;
        alpha[i] = rho[i] * std::real(s.col(i).dot(q));
        q -= alpha[i] * y.col(i);
    }

    VecC r = invH0 * q;
    for (int i = k - 2; i >= 0; --i) {
        if (rho[i] == 0.0) break;
        const double beta = rho[i] * std::real(y.col(i).dot(r));
        r += s.col(i) * (alpha[i] - beta);
    }
    return r;
}

}  // namespace p2p_harmonic
