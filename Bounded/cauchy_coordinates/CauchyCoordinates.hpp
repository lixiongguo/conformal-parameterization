#pragma once

#include <Eigen/Core>

#include <complex>
#include <stdexcept>
#include <vector>

namespace cauchy {

using Complex = std::complex<double>;
using VecC = Eigen::VectorXcd;
using MatC = Eigen::MatrixXcd;

struct Options {
    double boundary_epsilon = 1e-12;
};

struct PrecisionError {
    double constant_precision = 0.0;
    double linear_precision = 0.0;
};

class CauchyError : public std::runtime_error {
public:
    explicit CauchyError(const std::string& message) : std::runtime_error(message) {}
};

MatC computeCauchyCoordinates(
    const VecC& polygon,
    const VecC& query_points,
    const Options& options = Options{});

MatC cauchyCoordinates(
    const VecC& polygon,
    const VecC& query_points,
    const Options& options = Options{});

MatC computeCauchyCoordinateDerivatives(
    const VecC& polygon,
    const VecC& query_points,
    const Options& options = Options{});

MatC derivativesOfCauchyCoord(
    const VecC& polygon,
    const VecC& query_points,
    const Options& options = Options{});

MatC computeCauchyCoordinateSecondDerivatives(
    const VecC& polygon,
    const VecC& query_points,
    const Options& options = Options{});

MatC secondDerivativesOfCauchyCoord(
    const VecC& polygon,
    const VecC& query_points,
    const Options& options = Options{});

VecC evaluateCauchyMap(
    const VecC& polygon,
    const VecC& target_vertices,
    const VecC& query_points,
    const Options& options = Options{});

PrecisionError checkPrecision(
    const VecC& polygon,
    const VecC& query_points,
    const Options& options = Options{});

}  // namespace cauchy
