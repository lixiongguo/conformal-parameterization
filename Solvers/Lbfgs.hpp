#pragma once

#include <Eigen/Core>

#include <complex>

namespace p2p_harmonic {

using Complex = std::complex<double>;
using VecC = Eigen::VectorXcd;
using MatC = Eigen::MatrixXcd;
using MatX = Eigen::MatrixXd;

MatC lbfgsIter(const MatX& invH0, const MatC& g, const MatC& x);

}  // namespace p2p_harmonic
