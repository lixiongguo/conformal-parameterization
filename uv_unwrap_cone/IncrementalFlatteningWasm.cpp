#include "IncrementalFlatteningWasm.h"

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "Mesh.h"

namespace {

constexpr double kPi = 3.14159265358979323846;

struct FaceMetricData {
    int v0;
    int v1;
    int v2;
    int e01;
    int e12;
    int e20;
};

double clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

bool computeTriangleAngles(double l01, double l12, double l20,
                           double& a0, double& a1, double& a2) {
    const double eps = 1e-12;
    if (l01 <= eps || l12 <= eps || l20 <= eps) return false;
    if (l01 + l12 <= l20 + eps || l12 + l20 <= l01 + eps || l20 + l01 <= l12 + eps) return false;

    const double c0 = clamp((l01 * l01 + l20 * l20 - l12 * l12) / (2.0 * l01 * l20), -1.0, 1.0);
    const double c1 = clamp((l01 * l01 + l12 * l12 - l20 * l20) / (2.0 * l01 * l12), -1.0, 1.0);
    const double c2 = clamp((l12 * l12 + l20 * l20 - l01 * l01) / (2.0 * l12 * l20), -1.0, 1.0);
    a0 = std::acos(c0);
    a1 = std::acos(c1);
    a2 = std::acos(c2);
    return std::isfinite(a0) && std::isfinite(a1) && std::isfinite(a2);
}

double roundQuarterPi(double x) {
    return std::round(2.0 * x / kPi) * (kPi * 0.5);
}

void buildFaceMetricData(const Mesh& mesh, std::vector<FaceMetricData>& outFaces) {
    outFaces.clear();
    outFaces.reserve(mesh.faces.size());
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        HalfEdgeCIter h = f->he;
        FaceMetricData fm;
        fm.v0 = h->vertex->index;
        fm.v1 = h->next->vertex->index;
        fm.v2 = h->next->next->vertex->index;
        fm.e01 = h->edge->index;
        fm.e12 = h->next->edge->index;
        fm.e20 = h->next->next->edge->index;
        outFaces.push_back(fm);
    }
}

void buildVertexAreas(const Mesh& mesh, Eigen::VectorXd& areas) {
    const int n = static_cast<int>(mesh.vertices.size());
    areas = Eigen::VectorXd::Constant(n, 1e-6);
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        const double a = f->area();
        if (a <= 0.0) continue;
        HalfEdgeCIter h = f->he;
        areas[h->vertex->index] += a / 3.0;
        areas[h->next->vertex->index] += a / 3.0;
        areas[h->next->next->vertex->index] += a / 3.0;
    }
}

void buildBoundaryFlags(const Mesh& mesh, std::vector<bool>& isBoundary) {
    const int n = static_cast<int>(mesh.vertices.size());
    isBoundary.assign(n, false);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        isBoundary[v->index] = v->isBoundary();
    }
}

Eigen::SparseMatrix<double> buildCotanLaplacian(const Mesh& mesh) {
    const int n = static_cast<int>(mesh.vertices.size());
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.edges.size() * 4 + n);

    std::vector<double> diag(n, 0.0);
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
        const int i = e->he->vertex->index;
        const int j = e->he->flip->vertex->index;
        double w = e->cotanWeigth();
        if (!std::isfinite(w)) w = 0.0;
        diag[i] += w;
        diag[j] += w;
        triplets.emplace_back(i, j, -w);
        triplets.emplace_back(j, i, -w);
    }
    for (int i = 0; i < n; ++i) {
        triplets.emplace_back(i, i, diag[i] + 1e-10);
    }

    Eigen::SparseMatrix<double> L(n, n);
    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}

void edgeLengthsFromScale(const Mesh& mesh, const Eigen::VectorXd& u,
                          const std::vector<double>& baseLen, std::vector<double>& edgeLen) {
    edgeLen.resize(mesh.edges.size());
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
        const int i = e->he->vertex->index;
        const int j = e->he->flip->vertex->index;
        edgeLen[e->index] = baseLen[e->index] * std::exp(0.5 * (u[i] + u[j]));
    }
}

void computeMetricCurvature(const std::vector<FaceMetricData>& faces,
                            const std::vector<double>& edgeLen,
                            const std::vector<bool>& isBoundary,
                            Eigen::VectorXd& K) {
    const int n = static_cast<int>(isBoundary.size());
    Eigen::VectorXd angleSums = Eigen::VectorXd::Zero(n);

    for (const FaceMetricData& f : faces) {
        const double l01 = edgeLen[f.e01];
        const double l12 = edgeLen[f.e12];
        const double l20 = edgeLen[f.e20];
        double a0 = 0.0, a1 = 0.0, a2 = 0.0;
        if (!computeTriangleAngles(l01, l12, l20, a0, a1, a2)) {
            continue;
        }
        angleSums[f.v0] += a0;
        angleSums[f.v1] += a1;
        angleSums[f.v2] += a2;
    }

    K = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < n; ++i) {
        const double tgt = isBoundary[i] ? kPi : (2.0 * kPi);
        K[i] = tgt - angleSums[i];
    }
}

bool solveConstrainedMinNorm(
    const Eigen::VectorXd& vertexAreas,
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& Lrow,
    const Eigen::VectorXd& Korig,
    const std::vector<std::pair<int, double>>& constraints,
    Eigen::VectorXd& uOut) {
    const int n = static_cast<int>(vertexAreas.size());
    const int m = static_cast<int>(constraints.size());
    if (m == 0) {
        uOut = Eigen::VectorXd::Zero(n);
        return true;
    }

    const int N = n + m;
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(n + 2 * (Lrow.nonZeros() + m));

    for (int i = 0; i < n; ++i) {
        trips.emplace_back(i, i, std::max(vertexAreas[i], 1e-8));
    }

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(N);
    for (int r = 0; r < m; ++r) {
        const int rowIdx = constraints[r].first;
        const double targetK = constraints[r].second;
        rhs[n + r] = targetK - Korig[rowIdx];

        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Lrow, rowIdx); it; ++it) {
            const int c = it.col();
            const double v = it.value();
            if (std::abs(v) < 1e-14) continue;
            trips.emplace_back(c, n + r, v);
            trips.emplace_back(n + r, c, v);
        }
    }

    Eigen::SparseMatrix<double> KKT(N, N);
    KKT.setFromTriplets(trips.begin(), trips.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(KKT);
    solver.factorize(KKT);
    if (solver.info() != Eigen::Success) return false;
    Eigen::VectorXd x = solver.solve(rhs);
    if (solver.info() != Eigen::Success) return false;
    uOut = x.head(n);
    return true;
}

void buildUvFromScale(const Eigen::VectorXd& u, std::vector<double>& uvFlat) {
    const int n = static_cast<int>(u.size());
    uvFlat.clear();
    uvFlat.reserve(static_cast<size_t>(n) * 2);
    if (n == 0) return;

    double minU = std::numeric_limits<double>::infinity();
    double maxU = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
        minU = std::min(minU, u[i]);
        maxU = std::max(maxU, u[i]);
    }
    const double rangeU = std::max(maxU - minU, 1e-12);

    for (int i = 0; i < n; ++i) {
        const double x = (n == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(n - 1);
        const double y = (u[i] - minU) / rangeU;
        uvFlat.push_back(x);
        uvFlat.push_back(y);
    }
}

}  // namespace

bool runIncrementalFlattening(Mesh& mesh, const IncrementalOptions& options, IncrementalState& outState) {
    outState = IncrementalState();
    outState.nVertices = static_cast<int>(mesh.vertices.size());
    outState.nEdges = static_cast<int>(mesh.edges.size());
    outState.nFaces = static_cast<int>(mesh.faces.size());

    std::vector<FaceMetricData> faces;
    buildFaceMetricData(mesh, faces);

    buildBoundaryFlags(mesh, outState.isBoundary);
    outState.baseLen.assign(static_cast<size_t>(outState.nEdges), 0.0);
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
        outState.baseLen[e->index] = e->length();
    }

    Eigen::VectorXd areas;
    buildVertexAreas(mesh, areas);
    Eigen::SparseMatrix<double> L = buildCotanLaplacian(mesh);
    Eigen::SparseMatrix<double, Eigen::RowMajor> Lrow = L;

    outState.u = Eigen::VectorXd::Zero(outState.nVertices);
    outState.edgeLen = outState.baseLen;

    computeMetricCurvature(faces, outState.edgeLen, outState.isBoundary, outState.Korig);

    outState.flattened.assign(static_cast<size_t>(outState.nVertices), false);
    for (int i = 0; i < outState.nVertices; ++i) {
        if (outState.isBoundary[i]) outState.flattened[i] = true;
    }

    double eps = options.eps0;
    for (int iter = 0; iter < options.maxIters; ++iter) {
        Eigen::VectorXd Kcur;
        computeMetricCurvature(faces, outState.edgeLen, outState.isBoundary, Kcur);

        int remainingInterior = 0;
        for (int i = 0; i < outState.nVertices; ++i) {
            if (outState.isBoundary[i]) continue;
            if (!outState.flattened[i]) {
                if (std::abs(Kcur[i]) <= eps) {
                    outState.flattened[i] = true;
                } else {
                    remainingInterior++;
                }
            }
        }

        std::vector<std::pair<int, double>> constraints;
        constraints.reserve(static_cast<size_t>(outState.nVertices));
        for (int i = 0; i < outState.nVertices; ++i) {
            if (outState.flattened[i] && !outState.isBoundary[i]) {
                constraints.emplace_back(i, 0.0);
            }
        }

        if (!constraints.empty()) {
            if (!solveConstrainedMinNorm(areas, Lrow, outState.Korig, constraints, outState.u)) {
                outState.error = "Flattening solve failed";
                return false;
            }
            edgeLengthsFromScale(mesh, outState.u, outState.baseLen, outState.edgeLen);
        }

        outState.flattenIterations = iter + 1;
        if (remainingInterior <= options.maxCones) break;
        if (eps >= options.epsMax) break;
        eps += options.epsStep;
    }

    for (int i = 0; i < outState.nVertices; ++i) {
        if (!outState.isBoundary[i] && !outState.flattened[i]) outState.candidates.insert(i);
    }

    std::set<int> remaining = outState.candidates;
    int roundIter = 0;
    while (!remaining.empty()) {
        Eigen::VectorXd Kcur;
        computeMetricCurvature(faces, outState.edgeLen, outState.isBoundary, Kcur);

        int bestV = -1;
        double bestErr = std::numeric_limits<double>::infinity();
        double bestRounded = 0.0;
        for (int vid : remaining) {
            const double r = roundQuarterPi(Kcur[vid]);
            const double err = std::abs(Kcur[vid] - r);
            if (err < bestErr) {
                bestErr = err;
                bestV = vid;
                bestRounded = r;
            }
        }
        if (bestV < 0) break;

        outState.roundedTargets[bestV] = bestRounded;
        remaining.erase(bestV);

        std::vector<std::pair<int, double>> constraints;
        constraints.reserve(static_cast<size_t>(outState.nVertices));
        for (int i = 0; i < outState.nVertices; ++i) {
            if (outState.flattened[i] && !outState.isBoundary[i]) constraints.emplace_back(i, 0.0);
        }
        for (const auto& kv : outState.roundedTargets) constraints.emplace_back(kv.first, kv.second);

        if (!solveConstrainedMinNorm(areas, Lrow, outState.Korig, constraints, outState.u)) {
            outState.error = "Rounding solve failed at iteration " + std::to_string(roundIter);
            return false;
        }
        edgeLengthsFromScale(mesh, outState.u, outState.baseLen, outState.edgeLen);
        roundIter++;
    }

    computeMetricCurvature(faces, outState.edgeLen, outState.isBoundary, outState.Kfinal);
    buildUvFromScale(outState.u, outState.uvFlat);
    return true;
}
