#include "LinAbf.h"
#include "Lscm.h"
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace {

struct FaceData {
    int v[3];
    Eigen::Vector2d p[3];
    double targetAngles[3];
    double area2;
};

constexpr double kPi = 3.14159265358979323846;

double cornerAngleFrom3D(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
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

double cornerAngleFrom2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c)
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

void collectBoundaryLoops(const Mesh& mesh, std::vector<std::vector<int>>& loops)
{
    loops.clear();
    loops.reserve(mesh.boundaries.size());
    for (size_t i = 0; i < mesh.boundaries.size(); i++) {
        HalfEdgeCIter start = mesh.boundaries[i];
        HalfEdgeCIter he = start;
        std::vector<int> loop;
        do {
            loop.push_back(he->vertex->index);
            he = he->next;
        } while (he != start);
        if (loop.size() >= 3) {
            loops.push_back(loop);
        }
    }
}

double loopPerimeter(const Mesh& mesh, const std::vector<int>& loop)
{
    double perimeter = 0.0;
    const int n = static_cast<int>(loop.size());
    for (int i = 0; i < n; i++) {
        const int a = loop[i];
        const int b = loop[(i + 1) % n];
        perimeter += (mesh.vertices[b].position - mesh.vertices[a].position).norm();
    }
    return perimeter;
}

void pinLoopOnCircle(Mesh& mesh, const std::vector<int>& loop, const Eigen::Vector2d& center, double radius)
{
    const int n = static_cast<int>(loop.size());
    if (n < 3) {
        return;
    }

    constexpr double kTwoPi = 6.28318530717958647692;
    for (int i = 0; i < n; i++) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        const double theta = kTwoPi * t;
        mesh.vertices[loop[i]].uv = Eigen::Vector2d(
            center.x() + radius * std::cos(theta),
            center.y() + radius * std::sin(theta));
    }
}

void initializeBoundaryUv(Mesh& mesh)
{
    std::vector<std::vector<int>> loops;
    collectBoundaryLoops(mesh, loops);
    if (loops.empty()) {
        return;
    }

    std::sort(loops.begin(), loops.end(), [&mesh](const std::vector<int>& a, const std::vector<int>& b) {
        return loopPerimeter(mesh, a) > loopPerimeter(mesh, b);
    });

    pinLoopOnCircle(mesh, loops[0], Eigen::Vector2d(0.0, 0.0), 0.95);

    const int holes = static_cast<int>(loops.size()) - 1;
    if (holes <= 0) {
        return;
    }

    constexpr double kTwoPi = 6.28318530717958647692;
    const double centerRing = (holes == 1) ? 0.0 : 0.55;
    const double holeRadius = std::max(0.04, std::min(0.20, 0.35 / static_cast<double>(holes)));

    for (int i = 0; i < holes; i++) {
        const double angle = kTwoPi * static_cast<double>(i) / static_cast<double>(holes);
        const Eigen::Vector2d center(centerRing * std::cos(angle), centerRing * std::sin(angle));
        pinLoopOnCircle(mesh, loops[i + 1], center, holeRadius);
    }
}

void buildFaceData(const Mesh& mesh, std::vector<FaceData>& faces)
{
    faces.clear();
    faces.reserve(mesh.faces.size());

    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) {
            continue;
        }

        FaceData fd{};
        HalfEdgeCIter he = f->he;
        Eigen::Vector3d pos[3];
        for (int k = 0; k < 3; k++) {
            fd.v[k] = he->vertex->index;
            pos[k] = he->vertex->position;
            he = he->next;
        }

        Eigen::Vector3d x3 = pos[1] - pos[0];
        Eigen::Vector3d y3 = pos[2] - pos[0];
        Eigen::Vector3d xhat = x3.normalized();
        Eigen::Vector3d zhat = xhat.cross(y3);
        if (zhat.norm() < 1e-12) {
            zhat = Eigen::Vector3d(0.0, 0.0, 1.0);
        } else {
            zhat.normalize();
        }
        Eigen::Vector3d yhat = zhat.cross(xhat);
        if (yhat.norm() < 1e-12) {
            yhat = Eigen::Vector3d(0.0, 1.0, 0.0);
        } else {
            yhat.normalize();
        }

        fd.p[0] = Eigen::Vector2d(0.0, 0.0);
        fd.p[1] = Eigen::Vector2d(x3.dot(xhat), x3.dot(yhat));
        fd.p[2] = Eigen::Vector2d(y3.dot(xhat), y3.dot(yhat));
        fd.area2 = (fd.p[1] - fd.p[0]).x() * (fd.p[2] - fd.p[0]).y() -
                   (fd.p[1] - fd.p[0]).y() * (fd.p[2] - fd.p[0]).x();

        fd.targetAngles[0] = cornerAngleFrom3D(pos[0], pos[1], pos[2]);
        fd.targetAngles[1] = cornerAngleFrom3D(pos[1], pos[2], pos[0]);
        fd.targetAngles[2] = cornerAngleFrom3D(pos[2], pos[0], pos[1]);

        faces.push_back(fd);
    }
}

bool solveWeightedLinearSystem(
    Mesh& mesh,
    const std::vector<FaceData>& faces,
    const std::unordered_set<int>& boundarySet,
    const std::vector<double>& faceWeights)
{
    const int nVertices = static_cast<int>(mesh.vertices.size());
    std::vector<int> interiorIndex(nVertices, -1);
    int nInterior = 0;
    for (int i = 0; i < nVertices; i++) {
        if (!boundarySet.count(i)) {
            interiorIndex[i] = nInterior++;
        }
    }
    if (nInterior <= 0) {
        return true;
    }

    const int rows = static_cast<int>(faces.size()) * 2;
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(rows * 18);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(rows);

    for (size_t fi = 0; fi < faces.size(); fi++) {
        const FaceData& f = faces[fi];
        const double area2 = f.area2;
        if (std::abs(area2) < 1e-14) {
            continue;
        }

        const double wFace = std::sqrt(std::max(1e-6, faceWeights[fi]));
        const int row0 = static_cast<int>(2 * fi);
        const int row1 = row0 + 1;

        for (int i = 0; i < 3; i++) {
            const int j = (i + 1) % 3;
            const int k = (i + 2) % 3;
            const Eigen::Vector2d e = f.p[k] - f.p[j];
            const double gx = e.y() / area2;
            const double gy = -e.x() / area2;

            const int vi = f.v[i];
            const double c_u_row0 = wFace * gx;
            const double c_v_row0 = -wFace * gy;
            const double c_u_row1 = wFace * gy;
            const double c_v_row1 = wFace * gx;

            if (boundarySet.count(vi)) {
                const double u = mesh.vertices[vi].uv.x();
                const double v = mesh.vertices[vi].uv.y();
                b(row0) -= c_u_row0 * u + c_v_row0 * v;
                b(row1) -= c_u_row1 * u + c_v_row1 * v;
            } else {
                const int base = 2 * interiorIndex[vi];
                triplets.push_back(Eigen::Triplet<double>(row0, base, c_u_row0));
                triplets.push_back(Eigen::Triplet<double>(row0, base + 1, c_v_row0));
                triplets.push_back(Eigen::Triplet<double>(row1, base, c_u_row1));
                triplets.push_back(Eigen::Triplet<double>(row1, base + 1, c_v_row1));
            }
        }
    }

    Eigen::SparseMatrix<double> A(rows, 2 * nInterior);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseMatrix<double> At = A.transpose();
    Eigen::SparseMatrix<double> AtA = At * A;
    Eigen::VectorXd Atb = At * b;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(AtA);
    if (solver.info() != Eigen::Success) {
        return false;
    }
    Eigen::VectorXd x = solver.solve(Atb);
    if (solver.info() != Eigen::Success) {
        return false;
    }

    for (int i = 0; i < nVertices; i++) {
        if (!boundarySet.count(i)) {
            const int base = 2 * interiorIndex[i];
            mesh.vertices[i].uv = Eigen::Vector2d(x(base), x(base + 1));
        }
    }
    return true;
}

double updateWeightsFromAngleResidual(
    const Mesh& mesh,
    const std::vector<FaceData>& faces,
    std::vector<double>& faceWeights)
{
    double avgResidual = 0.0;
    int valid = 0;

    for (size_t fi = 0; fi < faces.size(); fi++) {
        const FaceData& f = faces[fi];
        const Eigen::Vector2d uv0 = mesh.vertices[f.v[0]].uv;
        const Eigen::Vector2d uv1 = mesh.vertices[f.v[1]].uv;
        const Eigen::Vector2d uv2 = mesh.vertices[f.v[2]].uv;

        const double alpha0 = cornerAngleFrom2D(uv0, uv1, uv2);
        const double alpha1 = cornerAngleFrom2D(uv1, uv2, uv0);
        const double alpha2 = cornerAngleFrom2D(uv2, uv0, uv1);

        const double r =
            (std::abs(alpha0 - f.targetAngles[0]) +
             std::abs(alpha1 - f.targetAngles[1]) +
             std::abs(alpha2 - f.targetAngles[2])) / 3.0;

        const double signedArea2 =
            (uv1.x() - uv0.x()) * (uv2.y() - uv0.y()) -
            (uv1.y() - uv0.y()) * (uv2.x() - uv0.x());
        const bool flipped = signedArea2 <= 1e-12;

        double w = 1.0 / (r + 1e-3);
        if (flipped) {
            w *= 0.35;
        }
        faceWeights[fi] = std::max(0.15, std::min(12.0, w));

        avgResidual += r;
        valid++;
    }

    return valid > 0 ? avgResidual / static_cast<double>(valid) : 0.0;
}

} // namespace

LinAbf::LinAbf(Mesh& mesh0, int maxIterations)
    : Parameterization(mesh0)
    , m_maxIterations(std::max(1, maxIterations))
{
}

void LinAbf::parameterize()
{
    if (mesh.boundaries.empty()) {
        // Closed meshes are not handled by this simplified ABF pipeline.
        Lscm fallback(mesh);
        fallback.parameterize();
        return;
    }

    initializeBoundaryUv(mesh);

    std::unordered_set<int> boundarySet;
    for (size_t i = 0; i < mesh.boundaries.size(); i++) {
        HalfEdgeCIter start = mesh.boundaries[i];
        HalfEdgeCIter he = start;
        do {
            boundarySet.insert(he->vertex->index);
            he = he->next;
        } while (he != start);
    }

    std::vector<FaceData> faces;
    buildFaceData(mesh, faces);
    if (faces.empty()) {
        normalize();
        return;
    }

    std::vector<double> weights(faces.size(), 1.0);
    double prevResidual = std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < m_maxIterations; iter++) {
        if (!solveWeightedLinearSystem(mesh, faces, boundarySet, weights)) {
            break;
        }
        const double residual = updateWeightsFromAngleResidual(mesh, faces, weights);
        if (std::abs(prevResidual - residual) < 1e-6) {
            break;
        }
        prevResidual = residual;
    }

    normalize();
}

