#include "BoundedDistortionMapping.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <utility>
#include <vector>

namespace bounded_distortion {
namespace {

constexpr double kEps = 1e-12;
constexpr int kConstraintsPerFace = 2;  // 0: Re(α)≥ε    1: |β|≤κ Re(α)
constexpr int kMaxFactorizations = 16;

Complex toComplex(const Eigen::MatrixXd& m, int row) {
    return Complex(m(row, 0), m(row, 1));
}

void validateInputs(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, const Eigen::MatrixXd& uv) {
    if (vertices.cols() != 2 && vertices.cols() != 3) {
        throw BoundedDistortionError("vertices must be an n x 2 or n x 3 matrix");
    }
    if (uv.rows() != vertices.rows() || uv.cols() != 2) {
        throw BoundedDistortionError("initial uv must have the same number of rows as vertices and 2 columns");
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

// 平面网格用全局 xy；曲面网格用与 LSCM 相同的三角形局部正交基。
bool triangleLocalPoints(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    int face,
    Eigen::Vector2d p[3],
    double& area) {
    if (faces(face, 0) == faces(face, 1) || faces(face, 1) == faces(face, 2) ||
        faces(face, 2) == faces(face, 0)) {
        return false;
    }

    if (vertices.cols() == 2) {
        for (int l = 0; l < 3; ++l) {
            p[l] = vertices.row(faces(face, l));
        }
        const Eigen::Vector2d ab = p[1] - p[0];
        const Eigen::Vector2d ac = p[2] - p[0];
        area = 0.5 * std::abs(ab.x() * ac.y() - ab.y() * ac.x());
        return area >= kEps;
    }

    Eigen::Vector3d P[3];
    for (int l = 0; l < 3; ++l) {
        P[l] = vertices.row(faces(face, l));
    }
    const Eigen::Vector3d e1 = P[1] - P[0];
    const Eigen::Vector3d e2 = P[2] - P[0];
    const Eigen::Vector3d n = e1.cross(e2);
    const double nlen = n.norm();
    area = 0.5 * nlen;
    if (area < kEps || e1.squaredNorm() < kEps) {
        return false;
    }
    const Eigen::Vector3d xhat = e1 / e1.norm();
    const Eigen::Vector3d zhat = n / nlen;
    const Eigen::Vector3d yhat = zhat.cross(xhat);
    p[0] = Eigen::Vector2d(0.0, 0.0);
    p[1] = Eigen::Vector2d(e1.norm(), 0.0);
    p[2] = Eigen::Vector2d(e2.dot(xhat), e2.dot(yhat));
    return true;
}

bool localSystem(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    int face,
    double frame_angle,
    Eigen::Matrix3cd& system,
    double& area) {
    Eigen::Vector2d p[3];
    if (!triangleLocalPoints(vertices, faces, face, p, area)) {
        return false;
    }

    const Complex origin(p[0].x(), p[0].y());
    const Complex rot = std::polar(1.0, frame_angle);
    for (int l = 0; l < 3; ++l) {
        const Complex z = rot * (Complex(p[l].x(), p[l].y()) - origin);
        system(l, 0) = z;
        system(l, 1) = std::conj(z);
        system(l, 2) = Complex(1.0, 0.0);
    }
    return std::abs(system.determinant()) >= kEps && system.allFinite();
}

bool faceCoefficients(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv,
    int face,
    double frame_angle,
    FaceCoefficients& out) {
    Eigen::Matrix3cd system;
    double area = 0.0;
    if (!localSystem(vertices, faces, face, frame_angle, system, area)) {
        out.alpha = Complex(1.0, 0.0);
        out.beta = Complex(0.0, 0.0);
        out.delta = Complex(0.0, 0.0);
        return false;
    }

    Eigen::Vector3cd target;
    for (int l = 0; l < 3; ++l) {
        target[l] = toComplex(uv, faces(face, l));
    }
    const Eigen::Vector3cd coeffs = system.colPivHouseholderQr().solve(target);
    if (!coeffs.allFinite()) {
        out.alpha = Complex(1.0, 0.0);
        out.beta = Complex(0.0, 0.0);
        out.delta = Complex(0.0, 0.0);
        return false;
    }
    out.alpha = coeffs[0];
    out.beta = coeffs[1];
    out.delta = coeffs[2];
    return true;
}

struct FaceMap {
    int v[3] = {-1, -1, -1};
    double area = 0.0;
    double M[4][6] = {};
    bool ok = false;
};

bool buildFaceMap(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    int face,
    double frame_angle,
    double kappa,
    FaceMap& out) {
    out = FaceMap();
    out.v[0] = faces(face, 0);
    out.v[1] = faces(face, 1);
    out.v[2] = faces(face, 2);

    Eigen::Matrix3cd system;
    if (!localSystem(vertices, faces, face, frame_angle, system, out.area)) {
        return false;
    }
    const Eigen::Matrix3cd inv = system.inverse();
    if (!inv.allFinite()) {
        return false;
    }

    for (int l = 0; l < 3; ++l) {
        const Complex ca = inv(0, l);
        const Complex cb = inv(1, l);
        const int x_col = 2 * l;
        const int y_col = 2 * l + 1;
        out.M[0][x_col] = ca.real();
        out.M[0][y_col] = -ca.imag();
        out.M[1][x_col] = ca.imag();
        out.M[1][y_col] = ca.real();
        out.M[2][x_col] = cb.real();
        out.M[2][y_col] = -cb.imag();
        out.M[3][x_col] = cb.imag();
        out.M[3][y_col] = cb.real();
    }

    (void)kappa;
    out.ok = true;
    return true;
}

// 局部坐标下两条不等式的取值和梯度。g1 在 |β|≈0 时用 0 次梯度（此时约束通常不活跃）。
void faceConstraint(
    const FaceMap& fm,
    const double loc[6],
    double kappa,
    double eps,
    double g[2],
    double dg0[6],
    double dg1[6]) {
    double alpha_re = 0.0;
    double beta_re = 0.0;
    double beta_im = 0.0;
    for (int k = 0; k < 6; ++k) {
        alpha_re += fm.M[0][k] * loc[k];
        beta_re += fm.M[2][k] * loc[k];
        beta_im += fm.M[3][k] * loc[k];
    }
    const double beta_norm = std::sqrt(beta_re * beta_re + beta_im * beta_im);
    g[0] = eps - alpha_re;
    g[1] = beta_norm - kappa * alpha_re;
    for (int k = 0; k < 6; ++k) {
        dg0[k] = -fm.M[0][k];
        const double dbeta = (beta_norm > 1e-14)
            ? (beta_re * fm.M[2][k] + beta_im * fm.M[3][k]) / beta_norm
            : 0.0;
        dg1[k] = dbeta - kappa * fm.M[0][k];
    }
}

void localUv(const FaceMap& fm, const Eigen::MatrixXd& uv, double loc[6]) {
    for (int l = 0; l < 3; ++l) {
        loc[2 * l] = uv(fm.v[l], 0);
        loc[2 * l + 1] = uv(fm.v[l], 1);
    }
}

std::vector<std::pair<int, int>> uniqueEdges(const Eigen::MatrixXi& faces) {
    std::set<std::pair<int, int>> edge_set;
    for (int f = 0; f < faces.rows(); ++f) {
        for (int e = 0; e < 3; ++e) {
            int a = faces(f, e);
            int b = faces(f, (e + 1) % 3);
            if (a == b) continue;
            if (a > b) std::swap(a, b);
            edge_set.emplace(a, b);
        }
    }
    return {edge_set.begin(), edge_set.end()};
}

struct FreeLayout {
    std::vector<int> index;  // 2 * vertex + comp → free id, -1 if pinned
    std::vector<int> vertex;
    std::vector<int> comp;
    int n = 0;
};

FreeLayout makeFreeLayout(int nV, const std::vector<Anchor>& anchors) {
    std::vector<char> pinned(static_cast<std::size_t>(nV), 0);
    for (const Anchor& anchor : anchors) {
        if (anchor.vertex >= 0 && anchor.vertex < nV) {
            pinned[static_cast<std::size_t>(anchor.vertex)] = 1;
        }
    }
    FreeLayout layout;
    layout.index.assign(static_cast<std::size_t>(2 * nV), -1);
    for (int v = 0; v < nV; ++v) {
        if (pinned[static_cast<std::size_t>(v)]) continue;
        for (int c = 0; c < 2; ++c) {
            layout.index[static_cast<std::size_t>(2 * v + c)] = layout.n++;
            layout.vertex.push_back(v);
            layout.comp.push_back(c);
        }
    }
    return layout;
}

int freeId(const FreeLayout& layout, int vertex, int comp) {
    return layout.index[static_cast<std::size_t>(2 * vertex + comp)];
}

void applyStep(Eigen::MatrixXd& uv, const FreeLayout& layout, const Eigen::VectorXd& dx, double t) {
    for (int i = 0; i < layout.n; ++i) {
        uv(layout.vertex[static_cast<std::size_t>(i)], layout.comp[static_cast<std::size_t>(i)]) += t * dx(i);
    }
}

struct Eval {
    double quad = 0.0;
    double merit = 0.0;
    double violation = 0.0;
    Eigen::VectorXd grad;
    std::vector<double> gval;
};

Eval evaluate(
    const Eigen::MatrixXd& uv,
    const Eigen::MatrixXd& reference_uv,
    const std::vector<FaceMap>& maps,
    const std::vector<std::pair<int, int>>& edges,
    const FreeLayout& layout,
    const std::vector<double>& lambda,
    double rho_cone,
    double rho_pos,
    double kappa,
    double eps,
    const Options& options) {
    Eval out;
    out.grad = Eigen::VectorXd::Zero(layout.n);
    out.gval.assign(maps.size() * static_cast<std::size_t>(kConstraintsPerFace), 0.0);

    if (options.lscm_weight > 0.0) {
        for (const FaceMap& fm : maps) {
            if (!fm.ok) continue;
            double loc[6];
            localUv(fm, uv, loc);
            double beta_re = 0.0;
            double beta_im = 0.0;
            for (int k = 0; k < 6; ++k) {
                beta_re += fm.M[2][k] * loc[k];
                beta_im += fm.M[3][k] * loc[k];
            }
            const double weight = options.lscm_weight * fm.area;
            out.quad += weight * (beta_re * beta_re + beta_im * beta_im);
            for (int k = 0; k < 6; ++k) {
                const int id = freeId(layout, fm.v[k / 2], k % 2);
                if (id < 0) continue;
                out.grad(id) += 2.0 * weight * (beta_re * fm.M[2][k] + beta_im * fm.M[3][k]);
            }
        }
    }

    if (options.reference_weight > 0.0) {
        const double w = options.reference_weight;
        for (int i = 0; i < layout.n; ++i) {
            const int v = layout.vertex[static_cast<std::size_t>(i)];
            const int c = layout.comp[static_cast<std::size_t>(i)];
            const double diff = uv(v, c) - reference_uv(v, c);
            out.quad += w * diff * diff;
            out.grad(i) += 2.0 * w * diff;
        }
    }

    if (options.smoothness_weight > 0.0) {
        const double w = options.smoothness_weight;
        for (const auto& edge : edges) {
            const Eigen::RowVector2d diff =
                (uv.row(edge.first) - uv.row(edge.second)) -
                (reference_uv.row(edge.first) - reference_uv.row(edge.second));
            out.quad += w * diff.squaredNorm();
            for (int c = 0; c < 2; ++c) {
                const int ia = freeId(layout, edge.first, c);
                const int ib = freeId(layout, edge.second, c);
                if (ia >= 0) out.grad(ia) += 2.0 * w * diff(c);
                if (ib >= 0) out.grad(ib) -= 2.0 * w * diff(c);
            }
        }
    }

    out.merit = out.quad;
    if (rho_cone <= 0.0 && rho_pos <= 0.0) {
        return out;
    }

    for (std::size_t f = 0; f < maps.size(); ++f) {
        const FaceMap& fm = maps[f];
        if (!fm.ok) continue;
        double loc[6];
        localUv(fm, uv, loc);
        double g[2], dg0[6], dg1[6];
        faceConstraint(fm, loc, kappa, eps, g, dg0, dg1);
        const double* dg[2] = {dg0, dg1};
        for (int ci = 0; ci < kConstraintsPerFace; ++ci) {
            const double rho = (ci == 0) ? rho_pos : rho_cone;
            if (rho <= 0.0) continue;
            const std::size_t gid = f * static_cast<std::size_t>(kConstraintsPerFace) + static_cast<std::size_t>(ci);
            out.gval[gid] = g[ci];
            out.violation = std::max(out.violation, g[ci]);

            const double lam = lambda[gid];
            const double shift = -lam / rho;
            if (g[ci] < shift) {
                out.merit += -0.5 * lam * lam / rho;
                continue;
            }
            out.merit += lam * g[ci] + 0.5 * rho * g[ci] * g[ci];
            const double coeff = lam + rho * g[ci];
            for (int k = 0; k < 6; ++k) {
                const int id = freeId(layout, fm.v[k / 2], k % 2);
                if (id < 0) continue;
                out.grad(id) += coeff * dg[ci][k];
            }
        }
    }
    return out;
}

void accumulateQuadHessian(
    std::vector<Eigen::Triplet<double>>& trips,
    const std::vector<FaceMap>& maps,
    const std::vector<std::pair<int, int>>& edges,
    const FreeLayout& layout,
    const Options& options) {
    if (options.lscm_weight > 0.0) {
        for (const FaceMap& fm : maps) {
            if (!fm.ok) continue;
            const double weight = options.lscm_weight * fm.area;
            int ids[6];
            for (int k = 0; k < 6; ++k) ids[k] = freeId(layout, fm.v[k / 2], k % 2);
            for (int a = 0; a < 6; ++a) {
                if (ids[a] < 0) continue;
                for (int b = 0; b < 6; ++b) {
                    if (ids[b] < 0) continue;
                    const double h = 2.0 * weight *
                        (fm.M[2][a] * fm.M[2][b] + fm.M[3][a] * fm.M[3][b]);
                    if (h != 0.0) trips.emplace_back(ids[a], ids[b], h);
                }
            }
        }
    }

    if (options.reference_weight > 0.0) {
        const double h = 2.0 * options.reference_weight;
        for (int i = 0; i < layout.n; ++i) trips.emplace_back(i, i, h);
    }

    if (options.smoothness_weight > 0.0) {
        const double h = 2.0 * options.smoothness_weight;
        for (const auto& edge : edges) {
            for (int c = 0; c < 2; ++c) {
                const int ia = freeId(layout, edge.first, c);
                const int ib = freeId(layout, edge.second, c);
                if (ia >= 0) trips.emplace_back(ia, ia, h);
                if (ib >= 0) trips.emplace_back(ib, ib, h);
                if (ia >= 0 && ib >= 0) {
                    trips.emplace_back(ia, ib, -h);
                    trips.emplace_back(ib, ia, -h);
                }
            }
        }
    }

    const double reg = 1e-10;
    for (int i = 0; i < layout.n; ++i) trips.emplace_back(i, i, reg);
}

void accumulatePenaltyHessian(
    std::vector<Eigen::Triplet<double>>& trips,
    const std::vector<FaceMap>& maps,
    const Eigen::MatrixXd& uv,
    const FreeLayout& layout,
    const std::vector<double>& lambda,
    double rho_cone,
    double rho_pos,
    double kappa,
    double eps) {
    for (std::size_t f = 0; f < maps.size(); ++f) {
        const FaceMap& fm = maps[f];
        if (!fm.ok) continue;
        double loc[6];
        localUv(fm, uv, loc);
        double g[2], dg0[6], dg1[6];
        faceConstraint(fm, loc, kappa, eps, g, dg0, dg1);
        const double* dg[2] = {dg0, dg1};
        double beta_re = 0.0, beta_im = 0.0;
        for (int k = 0; k < 6; ++k) {
            beta_re += fm.M[2][k] * loc[k];
            beta_im += fm.M[3][k] * loc[k];
        }
        const double beta_norm = std::sqrt(beta_re * beta_re + beta_im * beta_im);
        int ids[6];
        for (int k = 0; k < 6; ++k) ids[k] = freeId(layout, fm.v[k / 2], k % 2);
        for (int ci = 0; ci < kConstraintsPerFace; ++ci) {
            const double rho = (ci == 0) ? rho_pos : rho_cone;
            if (rho <= 0.0) continue;
            const double lam = lambda[f * static_cast<std::size_t>(kConstraintsPerFace) + static_cast<std::size_t>(ci)];
            if (g[ci] < -lam / rho) continue;
            const double coeff = lam + rho * g[ci];
            for (int a = 0; a < 6; ++a) {
                if (ids[a] < 0) continue;
                for (int b = 0; b < 6; ++b) {
                    if (ids[b] < 0) continue;
                    double h = rho * dg[ci][a] * dg[ci][b];
                    if (ci == 1 && coeff > 0.0 && beta_norm > 1e-14) {
                        const double jtj = fm.M[2][a] * fm.M[2][b] + fm.M[3][a] * fm.M[3][b];
                        const double dir_a = (beta_re * fm.M[2][a] + beta_im * fm.M[3][a]) / beta_norm;
                        const double dir_b = (beta_re * fm.M[2][b] + beta_im * fm.M[3][b]) / beta_norm;
                        h += coeff * (jtj - dir_a * dir_b) / beta_norm;
                    }
                    if (h != 0.0) trips.emplace_back(ids[a], ids[b], h);
                }
            }
        }
    }
}

bool solveNewton(
    const Eigen::SparseMatrix<double>& H,
    const Eigen::VectorXd& grad,
    Eigen::VectorXd& dx) {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(H);
    if (solver.info() != Eigen::Success) return false;
    dx = solver.solve(-grad);
    if (solver.info() != Eigen::Success || !dx.allFinite()) return false;
    if (grad.dot(dx) >= 0.0) return false;
    return true;
}

double orientViolation(const std::vector<FaceMap>& maps, const Eigen::MatrixXd& uv, double kappa) {
    double violation = 0.0;
    for (const FaceMap& fm : maps) {
        if (!fm.ok) continue;
        double loc[6];
        localUv(fm, uv, loc);
        double alpha_re = 0.0, alpha_im = 0.0, beta_re = 0.0, beta_im = 0.0;
        for (int k = 0; k < 6; ++k) {
            alpha_re += fm.M[0][k] * loc[k];
            alpha_im += fm.M[1][k] * loc[k];
            beta_re += fm.M[2][k] * loc[k];
            beta_im += fm.M[3][k] * loc[k];
        }
        const double alpha_abs = std::sqrt(alpha_re * alpha_re + alpha_im * alpha_im);
        const double beta_abs = std::sqrt(beta_re * beta_re + beta_im * beta_im);
        violation = std::max(violation, beta_abs - kappa * alpha_abs);
        if (alpha_abs * alpha_abs - beta_abs * beta_abs <= 0.0) {
            violation = std::max(violation, beta_abs + 1.0);
        }
    }
    return violation;
}

void applyAnchors(Eigen::MatrixXd& uv, const std::vector<Anchor>& anchors) {
    for (const Anchor& anchor : anchors) {
        if (anchor.vertex < 0 || anchor.vertex >= uv.rows()) {
            throw BoundedDistortionError("anchor contains an out-of-range vertex index");
        }
        uv.row(anchor.vertex) = anchor.target.transpose();
    }
}

double medianAbsAlpha(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& uv) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(faces.rows()));
    for (int f = 0; f < faces.rows(); ++f) {
        FaceCoefficients coeff;
        if (!faceCoefficients(vertices, faces, uv, f, 0.0, coeff)) continue;
        const double a = std::abs(coeff.alpha);
        if (a > kEps && std::isfinite(a)) values.push_back(a);
    }
    if (values.empty()) return 1.0;
    const std::size_t mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid), values.end());
    return std::max(values[mid], kEps);
}

}  // namespace

double distortionToKappa(double distortion_bound) {
    if (!(distortion_bound >= 1.0) || !std::isfinite(distortion_bound)) {
        throw BoundedDistortionError("distortion_bound must be a finite number >= 1");
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
        faceCoefficients(vertices, faces, uv, f, angle, coeffs[static_cast<std::size_t>(f)]);
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
        if (std::abs(coeffs[i].alpha) > kEps) {
            angles[i] = std::arg(coeffs[i].alpha);
        }
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
    const FreeLayout layout = makeFreeLayout(static_cast<int>(uv.rows()), anchors);
    const double kappa = distortionToKappa(options.distortion_bound);

    double eps = options.min_alpha_real;
    if (!(eps > 0.0) || !std::isfinite(eps)) eps = 1e-6;
    const double alpha_scale = medianAbsAlpha(vertices, faces, uv);
    eps = std::min(eps, 1e-4 * alpha_scale);

    double rho_cone = std::max(0.0, options.distortion_penalty);
    double rho_pos = std::max(0.0, options.positivity_penalty);
    std::vector<double> lambda(static_cast<std::size_t>(faces.rows()) * kConstraintsPerFace, 0.0);

    Eigen::MatrixXd best_uv = uv;
    double best_violation = std::numeric_limits<double>::infinity();
    double best_energy = std::numeric_limits<double>::infinity();
    int iterations = 0;
    int factorizations = 0;

    const int outer_count = std::max(1, options.outer_iterations);
    const int newton_per_outer = std::max(1, std::min(options.inner_iterations, 8));
    const double viol_tol = 1e-4;
    const double grad_tol = 1e-6;

    const Eigen::RowVector2d uv_lo = reference_uv.colwise().minCoeff();
    const Eigen::RowVector2d uv_hi = reference_uv.colwise().maxCoeff();
    double span = (uv_hi - uv_lo).norm();
    if (!(span > 0.0)) span = 1.0;
    const double max_disp = 0.25 * span;

    auto consider = [&](const Eigen::MatrixXd& cand, double violation, double energy) {
        const bool better_viol = violation < best_violation - 1e-9;
        const bool tie = std::abs(violation - best_violation) <= 1e-9 && energy < best_energy;
        if (better_viol || tie) {
            best_uv = cand;
            best_violation = violation;
            best_energy = energy;
        }
    };

    if (layout.n == 0) {
        SolveResult result;
        result.uv = uv;
        result.frame_angles = alignedFrameAngles(vertices, faces, uv);
        result.faces = computeFaceStats(vertices, faces, uv, options.distortion_bound, result.frame_angles);
        result.max_distortion = 0.0;
        result.min_jacobian = std::numeric_limits<double>::infinity();
        for (const FaceStats& stat : result.faces) {
            if (!std::isfinite(stat.jacobian) || !std::isfinite(stat.distortion)) continue;
            result.max_distortion = std::max(result.max_distortion, stat.distortion);
            result.min_jacobian = std::min(result.min_jacobian, stat.jacobian);
        }
        if (!std::isfinite(result.min_jacobian)) result.min_jacobian = 0.0;
        result.final_energy = 0.0;
        result.iterations = 0;
        return result;
    }

    for (int outer = 0; outer < outer_count; ++outer) {
        if (factorizations >= kMaxFactorizations) break;
        std::fill(lambda.begin(), lambda.end(), 0.0);

        const std::vector<double> frame_angles = alignedFrameAngles(vertices, faces, uv);
        std::vector<FaceMap> maps(static_cast<std::size_t>(faces.rows()));
        for (int f = 0; f < faces.rows(); ++f) {
            buildFaceMap(vertices, faces, f, frame_angles[static_cast<std::size_t>(f)], kappa, maps[static_cast<std::size_t>(f)]);
        }

        std::vector<Eigen::Triplet<double>> base_trips;
        base_trips.reserve(static_cast<std::size_t>(layout.n) + maps.size() * 36);
        accumulateQuadHessian(base_trips, maps, edges, layout, options);

        for (int inner = 0; inner < newton_per_outer; ++inner) {
            if (factorizations >= kMaxFactorizations) break;

            Eval cur = evaluate(uv, reference_uv, maps, edges, layout, lambda, rho_cone, rho_pos, kappa, eps, options);
            if (!std::isfinite(cur.merit) || !cur.grad.allFinite()) break;
            const double cur_orient = orientViolation(maps, uv, kappa);
            consider(uv, cur_orient, cur.quad);

            const double grad_norm = cur.grad.norm();
            if (cur_orient <= viol_tol && grad_norm <= grad_tol * (1.0 + std::sqrt(std::max(cur.quad, 0.0)))) {
                outer = outer_count;
                break;
            }

            std::vector<Eigen::Triplet<double>> trips = base_trips;
            accumulatePenaltyHessian(trips, maps, uv, layout, lambda, rho_cone, rho_pos, kappa, eps);
            Eigen::SparseMatrix<double> H(layout.n, layout.n);
            H.setFromTriplets(trips.begin(), trips.end());
            H.makeCompressed();
            double dmax = 0.0;
            for (int i = 0; i < layout.n; ++i) dmax = std::max(dmax, std::abs(H.coeff(i, i)));
            const double damp = 1e-3 * std::max(dmax, 1.0);
            for (int i = 0; i < layout.n; ++i) H.coeffRef(i, i) += damp;

            Eigen::VectorXd dx;
            ++factorizations;
            if (!solveNewton(H, cur.grad, dx)) break;

            double step = 1.0;
            const double peak = dx.cwiseAbs().maxCoeff();
            if (peak > max_disp) step = max_disp / peak;
            bool accepted = false;
            const double descent = cur.grad.dot(dx);
            for (int trial = 0; trial < 16; ++trial) {
                Eigen::MatrixXd candidate = uv;
                applyStep(candidate, layout, dx, step);
                Eval next = evaluate(candidate, reference_uv, maps, edges, layout, lambda, rho_cone, rho_pos, kappa, eps, options);
                if (!std::isfinite(next.merit)) {
                    step *= 0.5;
                    continue;
                }
                const double next_orient = orientViolation(maps, candidate, kappa);
                const bool merit_ok = next.merit <= cur.merit + 1e-4 * step * descent;
                const bool orient_ok = next_orient <= cur_orient + 1e-5;
                if (merit_ok && orient_ok) {
                    uv = candidate;
                    accepted = true;
                    consider(uv, next_orient, next.quad);
                    for (std::size_t i = 0; i < lambda.size(); ++i) {
                        const double rho = (static_cast<int>(i % kConstraintsPerFace) == 0) ? rho_pos : rho_cone;
                        if (rho <= 0.0) continue;
                        lambda[i] = std::max(0.0, lambda[i] + rho * next.gval[i]);
                    }
                    if (next_orient > 0.8 * std::max(cur_orient, viol_tol)) {
                        if (rho_cone > 0.0) rho_cone = std::min(rho_cone * 2.0, 1e7);
                        if (rho_pos > 0.0) rho_pos = std::min(rho_pos * 2.0, 1e7);
                    }
                    break;
                }
                step *= 0.5;
            }
            ++iterations;
            if (!accepted || step < 1e-8) break;
        }
    }

    uv = best_uv;

    SolveResult result;
    result.uv = uv;
    result.frame_angles = alignedFrameAngles(vertices, faces, uv);
    result.faces = computeFaceStats(vertices, faces, uv, options.distortion_bound, result.frame_angles);
    result.max_distortion = 0.0;
    result.min_jacobian = std::numeric_limits<double>::infinity();
    for (const FaceStats& stat : result.faces) {
        if (std::isfinite(stat.jacobian)) {
            result.min_jacobian = std::min(result.min_jacobian, stat.jacobian);
        }
        if (!std::isfinite(stat.distortion)) {
            result.max_distortion = std::numeric_limits<double>::infinity();
        } else {
            result.max_distortion = std::max(result.max_distortion, stat.distortion);
        }
    }
    if (!std::isfinite(result.min_jacobian)) result.min_jacobian = 0.0;
    result.final_energy = std::isfinite(best_energy) ? best_energy : 0.0;
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
