#include "FastHGPSimple.h"
#include "../Parameterization/SimpleParam/LSCM/Lscm.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace {

constexpr double kEps = 1e-12;

void computeFastHgpLocalBasisCoordinates(std::vector<Eigen::Vector2d>& coords,
                                         const std::vector<Eigen::Vector3d>& positions)
{
    Eigen::Vector3d x = positions[1] - positions[0];
    Eigen::Vector3d y = positions[2] - positions[0];

    Eigen::Vector3d xhat = x;
    xhat.normalize();
    Eigen::Vector3d zhat = xhat.cross(y);
    zhat.normalize();
    Eigen::Vector3d yhat = zhat.cross(xhat);
    yhat.normalize();

    coords.push_back(Eigen::Vector2d(0, 0));
    coords.push_back(Eigen::Vector2d(x.norm(), 0));
    coords.push_back(Eigen::Vector2d(y.dot(xhat), y.dot(yhat)));
}

double signedUvArea(const Eigen::Vector2d& a,
                    const Eigen::Vector2d& b,
                    const Eigen::Vector2d& c)
{
    const Eigen::Vector2d e1 = b - a;
    const Eigen::Vector2d e2 = c - a;
    return 0.5 * (e1.x() * e2.y() - e2.x() * e1.y());
}

} // namespace

FastHGPSimple::FastHGPSimple(Mesh& mesh0, int maxIterations0)
    : Parameterization(mesh0),
      maxIterations(std::max(1, maxIterations0))
{
}

void FastHGPSimple::parameterize()
{
    if (mesh.boundaries.empty()) {
        return;
    }

    initializeWithLscm();
    buildTriangleData();
    buildFreeVariables();

    if (triangles.empty() || freeVertices.empty()) {
        normalize();
        return;
    }

    double previousEnergy = std::numeric_limits<double>::infinity();
    int stagnantIters = 0;

    for (int iter = 0; iter < maxIterations; ++iter) {
        Eigen::VectorXd gradient;
        const double energy = computeEnergyAndGradient(gradient);
        if (!std::isfinite(energy) || !gradient.allFinite()) {
            break;
        }

        const double gradNorm = gradient.norm();
        if (gradNorm < 1e-8) {
            break;
        }

        Eigen::VectorXd direction = -gradient;
        const double directionNorm = direction.norm();
        if (directionNorm < 1e-12) {
            break;
        }

        direction /= directionNorm;
        const std::vector<Eigen::Vector2d> uv = currentUvs();
        const double tMax = localInjectivityStep(uv, direction);
        const double t = decreasingEnergyStep(energy, gradient, direction, tMax);

        if (t * directionNorm < 1e-9) {
            break;
        }

        Eigen::VectorXd x = packFreeUvs();
        x += t * direction;
        unpackFreeUvs(x);

        if (std::abs(previousEnergy - energy) < 1e-8 * (energy + 1.0)) {
            ++stagnantIters;
            if (stagnantIters >= 5) {
                break;
            }
        } else {
            stagnantIters = 0;
        }
        previousEnergy = energy;
    }

    normalize();
}

void FastHGPSimple::initializeWithLscm()
{
    Lscm lscm(mesh);
    lscm.parameterize();
}

void FastHGPSimple::buildTriangleData()
{
    triangles.clear();

    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) {
            continue;
        }

        std::vector<int> ids;
        std::vector<Eigen::Vector3d> positions;
        HalfEdgeCIter h = f->he;
        do {
            ids.push_back(h->vertex->index);
            positions.push_back(h->vertex->position);
            h = h->next;
        } while (h != f->he);

        if (ids.size() != 3) {
            continue;
        }

        std::vector<Eigen::Vector2d> local;
        computeFastHgpLocalBasisCoordinates(local, positions);

        Eigen::Matrix2d p;
        p.col(0) = local[1] - local[0];
        p.col(1) = local[2] - local[0];

        const double det = p.determinant();
        if (std::abs(det) < kEps) {
            continue;
        }

        TriData tri;
        tri.v[0] = ids[0];
        tri.v[1] = ids[1];
        tri.v[2] = ids[2];
        tri.invLocal = p.inverse();
        tri.area = 0.5 * std::abs(det);
        triangles.push_back(tri);
    }
}

void FastHGPSimple::buildFreeVariables()
{
    variableIndex.assign(mesh.vertices.size(), -1);
    freeVertices.clear();

    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        if (!v->isBoundary()) {
            variableIndex[v->index] = static_cast<int>(freeVertices.size());
            freeVertices.push_back(v->index);
        }
    }
}

Eigen::VectorXd FastHGPSimple::packFreeUvs() const
{
    Eigen::VectorXd x(2 * static_cast<int>(freeVertices.size()));
    const int n = static_cast<int>(freeVertices.size());
    for (int i = 0; i < n; ++i) {
        const Eigen::Vector2d& uv = mesh.vertices[freeVertices[i]].uv;
        x(i) = uv.x();
        x(i + n) = uv.y();
    }
    return x;
}

void FastHGPSimple::unpackFreeUvs(const Eigen::VectorXd& x)
{
    const int n = static_cast<int>(freeVertices.size());
    for (int i = 0; i < n; ++i) {
        Eigen::Vector2d& uv = mesh.vertices[freeVertices[i]].uv;
        uv.x() = x(i);
        uv.y() = x(i + n);
    }
}

std::vector<Eigen::Vector2d> FastHGPSimple::currentUvs() const
{
    std::vector<Eigen::Vector2d> uv(mesh.vertices.size());
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        uv[v->index] = v->uv;
    }
    return uv;
}

std::vector<Eigen::Vector2d> FastHGPSimple::uvsWithStep(const Eigen::VectorXd& direction, double t) const
{
    std::vector<Eigen::Vector2d> uv = currentUvs();
    const int n = static_cast<int>(freeVertices.size());
    for (int i = 0; i < n; ++i) {
        uv[freeVertices[i]].x() += t * direction(i);
        uv[freeVertices[i]].y() += t * direction(i + n);
    }
    return uv;
}

double FastHGPSimple::symmetricDirichletEnergy(const std::vector<Eigen::Vector2d>& uv) const
{
    double total = 0.0;

    for (const TriData& tri : triangles) {
        Eigen::Matrix2d q;
        q.col(0) = uv[tri.v[1]] - uv[tri.v[0]];
        q.col(1) = uv[tri.v[2]] - uv[tri.v[0]];

        const Eigen::Matrix2d j = q * tri.invLocal;
        const double det = j.determinant();
        if (det <= kEps) {
            return std::numeric_limits<double>::infinity();
        }

        const Eigen::Matrix2d c = j.transpose() * j;
        total += 0.5 * tri.area * (c.trace() + c.inverse().trace());
    }

    return total;
}

double FastHGPSimple::computeEnergyAndGradient(Eigen::VectorXd& gradient) const
{
    const int n = static_cast<int>(freeVertices.size());
    gradient = Eigen::VectorXd::Zero(2 * n);

    double total = 0.0;
    const std::vector<Eigen::Vector2d> uv = currentUvs();

    for (const TriData& tri : triangles) {
        Eigen::Matrix2d q;
        q.col(0) = uv[tri.v[1]] - uv[tri.v[0]];
        q.col(1) = uv[tri.v[2]] - uv[tri.v[0]];

        const Eigen::Matrix2d j = q * tri.invLocal;
        const double det = j.determinant();
        if (det <= kEps) {
            return std::numeric_limits<double>::infinity();
        }

        const Eigen::Matrix2d c = j.transpose() * j;
        const Eigen::Matrix2d cInv = c.inverse();
        const Eigen::Matrix2d cInv2 = cInv * cInv;

        total += 0.5 * tri.area * (c.trace() + cInv.trace());

        const Eigen::Matrix2d gradJ = tri.area * (j - j * cInv2);
        const Eigen::Matrix2d gradQ = gradJ * tri.invLocal.transpose();

        Eigen::Vector2d g[3];
        g[1] = gradQ.col(0);
        g[2] = gradQ.col(1);
        g[0] = -g[1] - g[2];

        for (int k = 0; k < 3; ++k) {
            const int var = variableIndex[tri.v[k]];
            if (var < 0) {
                continue;
            }
            gradient(var) += g[k].x();
            gradient(var + n) += g[k].y();
        }
    }

    return total;
}

double FastHGPSimple::localInjectivityStep(const std::vector<Eigen::Vector2d>& uv,
                                           const Eigen::VectorXd& direction) const
{
    double t = 1.0;
    const int n = static_cast<int>(freeVertices.size());

    for (const TriData& tri : triangles) {
        const Eigen::Vector2d a = uv[tri.v[0]];
        const Eigen::Vector2d b = uv[tri.v[1]];
        const Eigen::Vector2d c = uv[tri.v[2]];

        Eigen::Vector2d da = Eigen::Vector2d::Zero();
        Eigen::Vector2d db = Eigen::Vector2d::Zero();
        Eigen::Vector2d dc = Eigen::Vector2d::Zero();

        const int ia = variableIndex[tri.v[0]];
        const int ib = variableIndex[tri.v[1]];
        const int ic = variableIndex[tri.v[2]];
        if (ia >= 0) da = Eigen::Vector2d(direction(ia), direction(ia + n));
        if (ib >= 0) db = Eigen::Vector2d(direction(ib), direction(ib + n));
        if (ic >= 0) dc = Eigen::Vector2d(direction(ic), direction(ic + n));

        const double area0 = signedUvArea(a, b, c);
        const double area1 = signedUvArea(a + da, b + db, c + dc);
        if (area0 <= kEps) {
            return 0.0;
        }
        if (area1 <= kEps) {
            const double denom = area0 - area1;
            if (denom > kEps) {
                t = std::min(t, 0.9 * area0 / denom);
            }
        }
    }

    return std::max(0.0, std::min(1.0, t));
}

double FastHGPSimple::decreasingEnergyStep(double energy,
                                           const Eigen::VectorXd& gradient,
                                           const Eigen::VectorXd& direction,
                                           double tMax) const
{
    double t = std::max(0.0, std::min(1.0, tMax));
    const double alpha = 1e-4;
    const double beta = 0.5;
    const double gd = gradient.dot(direction);

    if (t <= 0.0 || gd >= 0.0) {
        return 0.0;
    }

    for (int i = 0; i < 30; ++i) {
        const std::vector<Eigen::Vector2d> uv = uvsWithStep(direction, t);
        const double trialEnergy = symmetricDirichletEnergy(uv);
        if (std::isfinite(trialEnergy) && trialEnergy < energy + alpha * t * gd) {
            return t;
        }
        t *= beta;
    }

    return t;
}
