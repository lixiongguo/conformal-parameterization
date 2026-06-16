#include "BoundedDistortionMapping.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace bounded_distortion {
namespace {

constexpr double kEps = 1e-12;

Complex toComplex(const Eigen::MatrixXd& m, int row) {
    return Complex(m(row, 0), m(row, 1));
}

void validateInputs(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, const Eigen::MatrixXd& uv) {
    if (vertices.cols() != 2) {
        throw BoundedDistortionError("vertices must be an n x 2 matrix");
    }
    if (uv.rows() != vertices.rows() || uv.cols() != 2) {
        throw BoundedDistortionError("initial uv must have the same n x 2 shape as vertices");
    }
    if (faces.cols() != 3) {
        throw BoundedDistortionError("faces must be an m x 3 triangle index matrix");
    }
    for (int f = 0; f < faces.rows(); ++f) {
        for (int l = 0; l < 3; ++l) {
            const int v = faces(f, l);
            if (v < 0 || v >= vertices.rows()) {
                throw BoundedDistortionError("face contains an out-of-range vertex index");
            }
        }
    }
}

Eigen::Matrix3cd localSystem(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    int face,
    double frame_angle) {
    const Complex origin = toComplex(vertices, faces(face, 0));
    const Complex rot = std::polar(1.0, frame_angle);

    Eigen::Matrix3cd system;
    for (int l = 0; l < 3; ++l) {
        const Complex z = rot * (toComplex(vertices, faces(face, l)) - origin);
        system(l, 0) = z;
        system(l, 1) = std::conj(z);
        system(l, 2) = Complex(1.0, 0.0);
    }

    if (std::abs(system.determinant()) < kEps) {
        throw BoundedDistortionError("degenerate triangle in local system");
    }
    return system;
}

FaceCoefficients faceCoefficients(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv,
    int face,
    double frame_angle) {
    const Eigen::Matrix3cd system = localSystem(vertices, faces, face, frame_angle);
    Eigen::Vector3cd target;
    for (int l = 0; l < 3; ++l) {
        target[l] = toComplex(uv, faces(face, l));
    }

    const Eigen::Vector3cd coeffs = system.colPivHouseholderQr().solve(target);
    FaceCoefficients out;
    out.alpha = coeffs[0];
    out.beta = coeffs[1];
    out.delta = coeffs[2];
    return out;
}

Eigen::Matrix<double, 4, 6> faceLinearMap(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    int face,
    double frame_angle) {
    const Eigen::Matrix3cd inv = localSystem(vertices, faces, face, frame_angle).inverse();
    Eigen::Matrix<double, 4, 6> map = Eigen::Matrix<double, 4, 6>::Zero();

    for (int l = 0; l < 3; ++l) {
        const Complex ca = inv(0, l);
        const Complex cb = inv(1, l);
        const int x_col = 2 * l;
        const int y_col = 2 * l + 1;

        map(0, x_col) = ca.real();
        map(0, y_col) = -ca.imag();
        map(1, x_col) = ca.imag();
        map(1, y_col) = ca.real();

        map(2, x_col) = cb.real();
        map(2, y_col) = -cb.imag();
        map(3, x_col) = cb.imag();
        map(3, y_col) = cb.real();
    }

    return map;
}

std::vector<std::pair<int, int>> uniqueEdges(const Eigen::MatrixXi& faces) {
    std::set<std::pair<int, int>> edge_set;
    for (int f = 0; f < faces.rows(); ++f) {
        for (int e = 0; e < 3; ++e) {
            int a = faces(f, e);
            int b = faces(f, (e + 1) % 3);
            if (a > b) std::swap(a, b);
            edge_set.emplace(a, b);
        }
    }
    return {edge_set.begin(), edge_set.end()};
}

double faceArea(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, int face) {
    const Eigen::Vector2d a = vertices.row(faces(face, 0));
    const Eigen::Vector2d b = vertices.row(faces(face, 1));
    const Eigen::Vector2d c = vertices.row(faces(face, 2));
    const Eigen::Vector2d ab = b - a;
    const Eigen::Vector2d ac = c - a;
    return 0.5 * std::abs(ab.x() * ac.y() - ab.y() * ac.x());
}

struct EnergyGradient {
    double energy = 0.0;
    Eigen::MatrixXd gradient;
};

EnergyGradient energyAndGradient(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv,
    const Eigen::MatrixXd& reference_uv,
    const std::vector<double>& frame_angles,
    const std::vector<std::pair<int, int>>& edges,
    const Options& options) {
    const double kappa = distortionToKappa(options.distortion_bound);
    EnergyGradient out;
    out.gradient = Eigen::MatrixXd::Zero(uv.rows(), 2);

    for (int f = 0; f < faces.rows(); ++f) {
        const Eigen::Matrix<double, 4, 6> map = faceLinearMap(vertices, faces, f, frame_angles[f]);

        Eigen::Matrix<double, 6, 1> local_uv;
        for (int l = 0; l < 3; ++l) {
            const int v = faces(f, l);
            local_uv[2 * l] = uv(v, 0);
            local_uv[2 * l + 1] = uv(v, 1);
        }

        const Eigen::Vector4d y = map * local_uv;
        const double alpha_re = y[0];
        const double beta_re = y[2];
        const double beta_im = y[3];
        const double beta_norm = std::sqrt(beta_re * beta_re + beta_im * beta_im);

        Eigen::Vector4d dy = Eigen::Vector4d::Zero();

        if (options.lscm_weight > 0.0) {
            const double weight = options.lscm_weight * faceArea(vertices, faces, f);
            out.energy += weight * (beta_re * beta_re + beta_im * beta_im);
            dy[2] += 2.0 * weight * beta_re;
            dy[3] += 2.0 * weight * beta_im;
        }

        const double cone = beta_norm - kappa * alpha_re;
        if (cone > 0.0) {
            out.energy += options.distortion_penalty * cone * cone;
            dy[0] += options.distortion_penalty * (-2.0 * kappa * cone);
            if (beta_norm > kEps) {
                dy[2] += options.distortion_penalty * (2.0 * cone * beta_re / beta_norm);
                dy[3] += options.distortion_penalty * (2.0 * cone * beta_im / beta_norm);
            }
        }

        const double positive = options.min_alpha_real - alpha_re;
        if (positive > 0.0) {
            out.energy += options.positivity_penalty * positive * positive;
            dy[0] += options.positivity_penalty * (-2.0 * positive);
        }

        const Eigen::Matrix<double, 6, 1> local_grad = map.transpose() * dy;
        for (int l = 0; l < 3; ++l) {
            const int v = faces(f, l);
            out.gradient(v, 0) += local_grad[2 * l];
            out.gradient(v, 1) += local_grad[2 * l + 1];
        }
    }

    if (options.reference_weight > 0.0) {
        const Eigen::MatrixXd diff = uv - reference_uv;
        out.energy += options.reference_weight * diff.squaredNorm();
        out.gradient += 2.0 * options.reference_weight * diff;
    }

    if (options.smoothness_weight > 0.0) {
        for (const auto& edge : edges) {
            const int a = edge.first;
            const int b = edge.second;
            const Eigen::RowVector2d current = uv.row(a) - uv.row(b);
            const Eigen::RowVector2d reference = reference_uv.row(a) - reference_uv.row(b);
            const Eigen::RowVector2d diff = current - reference;
            out.energy += options.smoothness_weight * diff.squaredNorm();
            out.gradient.row(a) += 2.0 * options.smoothness_weight * diff;
            out.gradient.row(b) -= 2.0 * options.smoothness_weight * diff;
        }
    }

    return out;
}

void applyAnchors(Eigen::MatrixXd& uv, const std::vector<Anchor>& anchors) {
    for (const Anchor& anchor : anchors) {
        if (anchor.vertex < 0 || anchor.vertex >= uv.rows()) {
            throw BoundedDistortionError("anchor contains an out-of-range vertex index");
        }
        uv.row(anchor.vertex) = anchor.target.transpose();
    }
}

void zeroAnchorGradient(Eigen::MatrixXd& gradient, const std::vector<Anchor>& anchors) {
    for (const Anchor& anchor : anchors) {
        gradient.row(anchor.vertex).setZero();
    }
}

}  // namespace

double distortionToKappa(double distortion_bound) {
    if (distortion_bound < 1.0) {
        throw BoundedDistortionError("distortion_bound must be >= 1");
    }
    return (distortion_bound - 1.0) / (distortion_bound + 1.0);
}

std::vector<FaceCoefficients> computeFaceCoefficients(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv,
    const std::vector<double>& frame_angles) {
    validateInputs(vertices, faces, uv);
    if (!frame_angles.empty() && static_cast<int>(frame_angles.size()) != faces.rows()) {
        throw BoundedDistortionError("frame_angles must be empty or one value per face");
    }

    std::vector<FaceCoefficients> coeffs(static_cast<std::size_t>(faces.rows()));
    for (int f = 0; f < faces.rows(); ++f) {
        const double angle = frame_angles.empty() ? 0.0 : frame_angles[static_cast<std::size_t>(f)];
        coeffs[static_cast<std::size_t>(f)] = faceCoefficients(vertices, faces, uv, f, angle);
    }
    return coeffs;
}

std::vector<FaceStats> computeFaceStats(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv,
    double distortion_bound,
    const std::vector<double>& frame_angles) {
    const double kappa = distortionToKappa(distortion_bound);
    const std::vector<FaceCoefficients> coeffs = computeFaceCoefficients(vertices, faces, uv, frame_angles);

    std::vector<FaceStats> stats(coeffs.size());
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        const double alpha_abs = std::abs(coeffs[i].alpha);
        const double beta_abs = std::abs(coeffs[i].beta);

        FaceStats s;
        s.coeffs = coeffs[i];
        s.sigma_max = alpha_abs + beta_abs;
        s.sigma_min = alpha_abs - beta_abs;
        s.jacobian = alpha_abs * alpha_abs - beta_abs * beta_abs;
        s.distortion = s.sigma_min > kEps ? s.sigma_max / s.sigma_min : std::numeric_limits<double>::infinity();
        s.cone_violation = std::max(0.0, beta_abs - kappa * coeffs[i].alpha.real());
        stats[i] = s;
    }
    return stats;
}

std::vector<double> alignedFrameAngles(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv) {
    const std::vector<FaceCoefficients> coeffs = computeFaceCoefficients(vertices, faces, uv);
    std::vector<double> angles(coeffs.size(), 0.0);
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        angles[i] = std::arg(coeffs[i].alpha);
    }
    return angles;
}

SolveResult solveBoundedDistortionMap(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& initial_uv,
    const std::vector<Anchor>& anchors,
    const Options& options) {
    validateInputs(vertices, faces, initial_uv);

    Eigen::MatrixXd uv = initial_uv;
    applyAnchors(uv, anchors);
    const Eigen::MatrixXd reference_uv = uv;
    const std::vector<std::pair<int, int>> edges = uniqueEdges(faces);

    std::vector<double> frame_angles = alignedFrameAngles(vertices, faces, uv);
    double last_energy = std::numeric_limits<double>::infinity();
    int iterations = 0;

    for (int outer = 0; outer < options.outer_iterations; ++outer) {
        frame_angles = alignedFrameAngles(vertices, faces, uv);

        for (int inner = 0; inner < options.inner_iterations; ++inner) {
            EnergyGradient eg = energyAndGradient(vertices, faces, uv, reference_uv, frame_angles, edges, options);
            zeroAnchorGradient(eg.gradient, anchors);

            const double grad_norm = eg.gradient.norm();
            if (grad_norm < options.gradient_tolerance) {
                last_energy = eg.energy;
                break;
            }

            double step = options.initial_step;
            bool accepted = false;
            for (int trial = 0; trial < 20; ++trial) {
                Eigen::MatrixXd candidate = uv - step * eg.gradient;
                applyAnchors(candidate, anchors);
                const double new_energy =
                    energyAndGradient(vertices, faces, candidate, reference_uv, frame_angles, edges, options).energy;
                if (new_energy <= eg.energy || trial == 19) {
                    uv = candidate;
                    last_energy = new_energy;
                    accepted = true;
                    break;
                }
                step *= 0.5;
            }
            ++iterations;
            if (!accepted) break;
        }
    }

    SolveResult result;
    result.uv = uv;
    result.frame_angles = alignedFrameAngles(vertices, faces, uv);
    result.faces = computeFaceStats(vertices, faces, uv, options.distortion_bound, result.frame_angles);
    result.max_distortion = 0.0;
    result.min_jacobian = std::numeric_limits<double>::infinity();
    for (const FaceStats& stat : result.faces) {
        result.max_distortion = std::max(result.max_distortion, stat.distortion);
        result.min_jacobian = std::min(result.min_jacobian, stat.jacobian);
    }
    result.final_energy = last_energy;
    result.iterations = iterations;
    return result;
}

SolveResult solveBoundedDistortionLscm(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& initial_uv,
    const std::vector<Anchor>& anchors,
    const Options& options) {
    Options lscm_options = options;
    if (lscm_options.lscm_weight <= 0.0) {
        lscm_options.lscm_weight = 1.0;
    }
    lscm_options.reference_weight = 0.0;
    lscm_options.smoothness_weight = 0.0;
    return solveBoundedDistortionMap(vertices, faces, initial_uv, anchors, lscm_options);
}

}  // namespace bounded_distortion
