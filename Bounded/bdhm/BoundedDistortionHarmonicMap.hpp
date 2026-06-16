#pragma once

#include "../cauchy_coordinates/CauchyCoordinates.hpp"

#include <Eigen/Core>

#include <complex>
#include <stdexcept>
#include <vector>

namespace bdhm {

using Complex = std::complex<double>;
using VecC = Eigen::VectorXcd;
using MatC = Eigen::MatrixXcd;

enum class SolverMethod {
    GradientDescent,
    LBFGS
};

struct PointConstraint {
    Complex source = Complex(0.0, 0.0);
    Complex target = Complex(0.0, 0.0);
    double weight = 1.0;
};

struct Options {
    double distortion_bound = 2.0;
    double min_jacobian = 1e-8;
    double distortion_penalty = 100.0;
    double jacobian_penalty = 100.0;
    double position_weight = 1e4;
    double anti_holomorphic_weight = 1e-3;
    double phi_reference_weight = 1e-4;
    double initial_step = 1e-2;
    int max_iterations = 1000;
    double gradient_tolerance = 1e-10;
    SolverMethod solver_method = SolverMethod::GradientDescent;
    int lbfgs_history = 10;
};

struct Stats {
    double max_distortion = 0.0;
    double min_jacobian = 0.0;
    double max_cone_violation = 0.0;
};

struct SolveResult {
    VecC phi;
    VecC psi;
    VecC mapped_energy_samples;
    VecC mapped_constraints;
    Stats stats;
    double final_energy = 0.0;
    int iterations = 0;
};

using NloP2PHarmonicInput = Options;
using NloP2PHarmonicResult = SolveResult;

class BdhmError : public std::runtime_error {
public:
    explicit BdhmError(const std::string& message) : std::runtime_error(message) {}
};

double distortionToKappa(double distortion_bound);

VecC evaluateMap(const MatC& C, const VecC& phi, const VecC& psi);

Stats computeStats(const MatC& D, const VecC& phi, const VecC& psi, double distortion_bound);

SolveResult solve(
    const VecC& boundary_polygon,
    const VecC& energy_samples,
    const std::vector<PointConstraint>& constraints,
    const Options& options = Options{});

NloP2PHarmonicResult nloP2PHarmonic(
    const VecC& boundary_polygon,
    const VecC& energy_samples,
    const std::vector<PointConstraint>& constraints,
    const NloP2PHarmonicInput& options = NloP2PHarmonicInput{});

}  // namespace bdhm
