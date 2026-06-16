#pragma once

#include <Eigen/Core>

#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

namespace bounded_distortion {

using Complex = std::complex<double>;

struct Anchor {
    int vertex = -1;
    Eigen::Vector2d target = Eigen::Vector2d::Zero();
};

struct Options {
    double distortion_bound = 2.0;
    double min_alpha_real = 1e-6;
    double lscm_weight = 0.0;
    double distortion_penalty = 100.0;
    double positivity_penalty = 100.0;
    double reference_weight = 1e-3;
    double smoothness_weight = 1e-2;
    double initial_step = 1e-2;
    int outer_iterations = 8;
    int inner_iterations = 200;
    double gradient_tolerance = 1e-10;
};

struct FaceCoefficients {
    Complex alpha = Complex(0.0, 0.0);
    Complex beta = Complex(0.0, 0.0);
    Complex delta = Complex(0.0, 0.0);
};

struct FaceStats {
    FaceCoefficients coeffs;
    double sigma_max = 0.0;
    double sigma_min = 0.0;
    double distortion = std::numeric_limits<double>::infinity();
    double jacobian = 0.0;
    double cone_violation = 0.0;
};

struct SolveResult {
    Eigen::MatrixXd uv;
    std::vector<double> frame_angles;
    std::vector<FaceStats> faces;
    double max_distortion = 0.0;
    double min_jacobian = 0.0;
    double final_energy = 0.0;
    int iterations = 0;
};

class BoundedDistortionError : public std::runtime_error {
public:
    explicit BoundedDistortionError(const std::string& message) : std::runtime_error(message) {}
};

double distortionToKappa(double distortion_bound);

std::vector<FaceCoefficients> computeFaceCoefficients(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv,
    const std::vector<double>& frame_angles = {});

std::vector<FaceStats> computeFaceStats(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv,
    double distortion_bound,
    const std::vector<double>& frame_angles = {});

std::vector<double> alignedFrameAngles(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv);

SolveResult solveBoundedDistortionMap(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& initial_uv,
    const std::vector<Anchor>& anchors,
    const Options& options = Options{});

SolveResult solveBoundedDistortionLscm(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& initial_uv,
    const std::vector<Anchor>& anchors,
    const Options& options = Options{});

}  // namespace bounded_distortion
