#include "AbfPlusPlus.h"
#include "LinAbf.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

double angle2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c)
{
    Eigen::Vector2d u = b - a;
    Eigen::Vector2d v = c - a;
    const double nu = u.norm();
    const double nv = v.norm();
    if (nu < 1e-14 || nv < 1e-14) {
        return kPi / 3.0;
    }
    u /= nu;
    v /= nv;
    const double d = std::max(-1.0, std::min(1.0, u.dot(v)));
    return std::acos(d);
}

double angle3D(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
    Eigen::Vector3d u = b - a;
    Eigen::Vector3d v = c - a;
    const double nu = u.norm();
    const double nv = v.norm();
    if (nu < 1e-14 || nv < 1e-14) {
        return kPi / 3.0;
    }
    u /= nu;
    v /= nv;
    const double d = std::max(-1.0, std::min(1.0, u.dot(v)));
    return std::acos(d);
}

double computeAngleResidual(const Mesh& mesh)
{
    double sumResidual = 0.0;
    int count = 0;

    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) {
            continue;
        }

        HalfEdgeCIter he = f->he;
        int v[3];
        for (int i = 0; i < 3; i++) {
            v[i] = he->vertex->index;
            he = he->next;
        }

        const Eigen::Vector3d p0 = mesh.vertices[v[0]].position;
        const Eigen::Vector3d p1 = mesh.vertices[v[1]].position;
        const Eigen::Vector3d p2 = mesh.vertices[v[2]].position;
        const Eigen::Vector2d u0 = mesh.vertices[v[0]].uv;
        const Eigen::Vector2d u1 = mesh.vertices[v[1]].uv;
        const Eigen::Vector2d u2 = mesh.vertices[v[2]].uv;

        const double beta0 = angle3D(p0, p1, p2);
        const double beta1 = angle3D(p1, p2, p0);
        const double beta2 = angle3D(p2, p0, p1);

        const double alpha0 = angle2D(u0, u1, u2);
        const double alpha1 = angle2D(u1, u2, u0);
        const double alpha2 = angle2D(u2, u0, u1);

        const double residual =
            (std::abs(alpha0 - beta0) + std::abs(alpha1 - beta1) + std::abs(alpha2 - beta2)) / 3.0;
        sumResidual += residual;
        count++;
    }

    return count > 0 ? sumResidual / static_cast<double>(count) : 0.0;
}

void dampUvUpdate(Mesh& mesh, const std::vector<Eigen::Vector2d>& prevUv, double blend)
{
    const double t = std::max(0.0, std::min(1.0, blend));
    int idx = 0;
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++, idx++) {
        v->uv = (1.0 - t) * prevUv[idx] + t * v->uv;
    }
}

} // namespace

AbfPlusPlus::AbfPlusPlus(Mesh& mesh0, int maxIterations)
    : Parameterization(mesh0)
    , m_maxIterations(std::max(1, maxIterations))
{
}

void AbfPlusPlus::parameterize()
{
    // geogram-style pipeline: linear ABF initialization followed by nonlinear refinement.
    LinAbf linInit(mesh, std::max(4, m_maxIterations / 2));
    linInit.parameterize();

    double prevResidual = std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < m_maxIterations; iter++) {
        std::vector<Eigen::Vector2d> prevUv;
        prevUv.reserve(mesh.vertices.size());
        for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
            prevUv.push_back(v->uv);
        }

        // Nonlinear stage: rerun one linearized step and blend to stabilize.
        LinAbf linStep(mesh, 1);
        linStep.parameterize();
        dampUvUpdate(mesh, prevUv, 0.78);

        const double residual = computeAngleResidual(mesh);
        if (std::abs(prevResidual - residual) < 5e-7) {
            break;
        }
        prevResidual = residual;
    }

    normalize();
}

