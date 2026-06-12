#include "NRosyVectorFields.h"
#include "PrincipalCurvatureField.h"

#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <queue>
#include <utility>

using namespace Eigen;

static constexpr double kPi = 3.14159265358979323846;

// ===================================================================
//  Constructor
// ===================================================================

NRosyVectorFields::NRosyVectorFields(Mesh& mesh, int N)
    : mesh_(mesh)
    , N_(std::max(1, N))
{
    nVerts_ = static_cast<int>(mesh_.vertices.size());
    nEdges_ = static_cast<int>(mesh_.edges.size());
    nFaces_ = countInteriorFaces();

    faceField_.resize(nFaces_, 3);
    faceField_.setZero();
    faceTheta_.resize(nFaces_);
    faceTheta_.setZero();
    edgeJumps_.resize(nEdges_);
    edgeJumps_.setZero();
    singIndices_.resize(nVerts_);
    singIndices_.setZero();
}

// ===================================================================
//  Public computation methods
// ===================================================================

bool NRosyVectorFields::computeFromCurvature()
{
    lastAlgorithm_ = "PrincipalCurvature";

    if (nFaces_ < 1) return false;
    const auto t0 = std::chrono::high_resolution_clock::now();

    // Build Eigen-format geometry
    MatrixXd V(nVerts_, 3);
    for (VertexCIter v = mesh_.vertices.begin(); v != mesh_.vertices.end(); ++v)
        V.row(v->index) = v->position;

    // Count interior faces and build F
    int fi = 0;
    VectorXi F(nFaces_ * 3);
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        F[fi * 3]     = f->he->vertex->index;
        F[fi * 3 + 1] = f->he->next->vertex->index;
        F[fi * 3 + 2] = f->he->next->next->vertex->index;
        ++fi;
    }

    // Per-vertex normals (area-weighted)
    MatrixXd vN(nVerts_, 3);
    vN.setZero();
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        Vector3d p0 = f->he->vertex->position;
        Vector3d p1 = f->he->next->vertex->position;
        Vector3d p2 = f->he->next->next->vertex->position;
        Vector3d fn = (p1 - p0).cross(p2 - p0);
        double w = fn.norm();
        fn.normalize();
        vN.row(f->he->vertex->index)           += w * fn;
        vN.row(f->he->next->vertex->index)     += w * fn;
        vN.row(f->he->next->next->vertex->index) += w * fn;
    }
    for (int i = 0; i < nVerts_; ++i) {
        double len = vN.row(i).norm();
        if (len > 1e-12) vN.row(i) /= len;
    }

    // Per-face normals
    MatrixXd fN(nFaces_, 3);
    fi = 0;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        Vector3d p0 = f->he->vertex->position;
        Vector3d p1 = f->he->next->vertex->position;
        Vector3d p2 = f->he->next->next->vertex->position;
        Vector3d fn = (p1 - p0).cross(p2 - p0);
        fn.normalize();
        fN.row(fi) = fn;
        ++fi;
    }

    PrincipalCurvatureField::estimateFromWeingarten(V, F, vN, fN, faceField_);

    // Compute per-face theta
    fi = 0;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        Vector3d e1, e2, n;
        buildLocalFrame(f, e1, e2, n);
        Vector3d d = faceField_.row(fi);
        d -= d.dot(n) * n;
        double len = d.norm();
        if (len < 1e-12) d = e1; else d /= len;
        faceTheta_(fi) = std::atan2(d.dot(e2), d.dot(e1));
        ++fi;
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    lastTimeMs_ = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return true;
}

bool NRosyVectorFields::computeTrivialConnection(const std::vector<int>& singularities)
{
    lastAlgorithm_ = "TrivialConnection/N-" + std::to_string(N_) + "RoSy";

    if (nFaces_ < 1) return false;
    if (static_cast<int>(singularities.size()) != nVerts_) return false;

    const auto t0 = std::chrono::high_resolution_clock::now();

    // 1. Compute vertex curvature
    VectorXd K(nVerts_);
    for (VertexCIter v = mesh_.vertices.begin(); v != mesh_.vertices.end(); ++v)
        K(v->index) = vertexCurvature(v);

    // 2. Store singularities
    for (int i = 0; i < nVerts_; ++i)
        singIndices_(i) = singularities[i];

    // 3. Build incidence B
    SparseMatrix<double> B(nVerts_, nEdges_);
    buildIncidence(B);

    // 4. RHS: 2π·q/N - K
    const double period = 2.0 * kPi / static_cast<double>(N_);
    VectorXd rhs(nVerts_);
    for (int i = 0; i < nVerts_; ++i)
        rhs(i) = period * static_cast<double>(singularities[i]) - K(i);

    // 5. Solve KKT for φ
    VectorXd phi(nEdges_);
    if (!solveKKT(B, rhs, phi)) {
        lastAlgorithm_ += " (solve failed)";
        return false;
    }

    // 6. Round period jumps, compute residual (N-RoSy only, N>1)
    edgeJumps_.resize(nEdges_);
    VectorXd residual(nEdges_);
    if (N_ > 1) {
        for (int ei = 0; ei < nEdges_; ++ei) {
            double p = std::round(phi(ei) / period);
            edgeJumps_(ei) = p;
            residual(ei) = phi(ei) - p * period;
        }
    } else {
        residual = phi;
        edgeJumps_.setZero();
    }

    // 7. Propagate per-face angles
    propagateAngles(residual, period);

    // 8. Compute per-face unit directions
    {
        std::vector<int> globalToLocal(mesh_.faces.size(), -1);
        int li = 0;
        for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f)
            if (!f->isBoundary()) globalToLocal[f->index] = li++;

        for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
            if (f->isBoundary()) continue;
            const int fi = globalToLocal[f->index];

            Vector3d e1, e2, n;
            buildLocalFrame(f, e1, e2, n);
            double a = faceTheta_(fi);
            faceField_.row(fi) = (std::cos(a) * e1 + std::sin(a) * e2).normalized();
        }
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    lastTimeMs_ = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return true;
}

void NRosyVectorFields::autoSingularities()
{
    VectorXd K(nVerts_);
    for (VertexCIter v = mesh_.vertices.begin(); v != mesh_.vertices.end(); ++v)
        K(v->index) = vertexCurvature(v);

    const int chi = nVerts_ - nEdges_ + nFaces_;
    const int Qtarget = N_ * chi;

    singIndices_.setZero();
    if (Qtarget == 0) return;

    // Rank vertices by |curvature|
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(nVerts_);
    for (int i = 0; i < nVerts_; ++i)
        ranked.push_back({std::abs(K(i)), i});
    std::sort(ranked.begin(), ranked.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    int sign = (Qtarget > 0) ? 1 : -1;
    int remaining = std::abs(Qtarget);
    for (size_t i = 0; i < ranked.size() && remaining > 0; ++i) {
        singIndices_(ranked[i].second) += sign;
        --remaining;
    }
}

// ===================================================================
//  Private helpers
// ===================================================================

int NRosyVectorFields::countInteriorFaces() const
{
    int cnt = 0;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f)
        if (!f->isBoundary()) ++cnt;
    return cnt;
}

void NRosyVectorFields::buildLocalFrame(FaceCIter f,
                                        Vector3d& e1,
                                        Vector3d& e2,
                                        Vector3d& n) const
{
    Vector3d p0 = f->he->vertex->position;
    Vector3d p1 = f->he->next->vertex->position;
    Vector3d p2 = f->he->next->next->vertex->position;

    e1 = p1 - p0;
    if (e1.norm() < 1e-12) e1 = Vector3d::UnitX();
    e1.normalize();

    n = (p1 - p0).cross(p2 - p0);
    if (n.norm() < 1e-12) n = Vector3d::UnitZ();
    n.normalize();

    e2 = n.cross(e1);
    if (e2.norm() < 1e-12) e2 = Vector3d::UnitY();
    e2.normalize();
}

double NRosyVectorFields::vertexCurvature(VertexCIter v) const
{
    double sumAngles = 0.0;
    HalfEdgeCIter h = v->he;
    do {
        if (!h->onBoundary)
            sumAngles += h->next->next->angle();
        h = h->flip->next;
    } while (h != v->he);
    return 2.0 * kPi - sumAngles;
}

double NRosyVectorFields::transportRotation(HalfEdgeCIter h) const
{
    Vector3d e1i, e2i, ni;
    Vector3d e1j, e2j, nj;
    buildLocalFrame(h->face,     e1i, e2i, ni);
    buildLocalFrame(h->flip->face, e1j, e2j, nj);

    Vector3d u = h->flip->vertex->position - h->vertex->position;
    double thetaIJ = std::atan2(u.dot(e2i), u.dot(e1i));
    double thetaJI = std::atan2(u.dot(e2j), u.dot(e1j));
    return -thetaIJ + thetaJI;
}

void NRosyVectorFields::buildIncidence(SparseMatrix<double>& B) const
{
    std::vector<Triplet<double>> trips;
    trips.reserve(nEdges_ * 2);
    for (EdgeCIter e = mesh_.edges.begin(); e != mesh_.edges.end(); ++e) {
        int i = e->he->vertex->index;
        int j = e->he->flip->vertex->index;
        trips.emplace_back(i, e->index,  1.0);
        trips.emplace_back(j, e->index, -1.0);
    }
    B.resize(nVerts_, nEdges_);
    B.setFromTriplets(trips.begin(), trips.end());
}

bool NRosyVectorFields::solveKKT(const SparseMatrix<double>& B,
                                  const VectorXd& rhs,
                                  VectorXd& phi) const
{
    const int Nsys = nEdges_ + nVerts_;

    std::vector<Triplet<double>> trips;
    trips.reserve(nEdges_ + 2 * B.nonZeros() + 1);
    for (int i = 0; i < nEdges_; ++i)
        trips.emplace_back(i, i, 1.0);
    for (int k = 0; k < B.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            trips.emplace_back(it.col(), nEdges_ + it.row(), it.value());
            trips.emplace_back(nEdges_ + it.row(), it.col(), it.value());
        }
    }
    trips.emplace_back(nEdges_, nEdges_, 1e-12); // gauge regularizer

    SparseMatrix<double> KKT(Nsys, Nsys);
    KKT.setFromTriplets(trips.begin(), trips.end());

    VectorXd b = VectorXd::Zero(Nsys);
    b.tail(nVerts_) = rhs;

    SparseLU<SparseMatrix<double>> solver;
    solver.analyzePattern(KKT);
    solver.factorize(KKT);
    if (solver.info() != Success) return false;

    VectorXd sol = solver.solve(b);
    if (solver.info() != Success) return false;

    phi = sol.head(nEdges_);
    return true;
}

// ===================================================================
//  Per-face angle propagation (BFS on dual graph)
// ===================================================================

void NRosyVectorFields::propagateAngles(const VectorXd& phi, double /*period*/)
{
    // Build global→local face index mapping
    std::vector<int> globalToLocal(mesh_.faces.size(), -1);
    int localIdx = 0;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (!f->isBoundary()) {
            globalToLocal[f->index] = localIdx++;
        }
    }

    // Use a temporary vector indexed by global face index for BFS propagation
    std::vector<double> alphaGlobal(mesh_.faces.size(), 0.0);
    std::vector<char> visited(mesh_.faces.size(), 0);

    int root = -1;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (!f->isBoundary()) { root = f->index; break; }
    }
    if (root < 0) return;

    std::queue<int> q;
    visited[root] = 1;
    q.push(root);

    while (!q.empty()) {
        int fi = q.front(); q.pop();
        FaceCIter f = mesh_.faces.begin() + fi;
        HalfEdgeCIter h0 = f->he;
        HalfEdgeCIter h = h0;
        do {
            if (!h->edge->isBoundary()) {
                FaceCIter g = h->flip->face;
                if (!g->isBoundary() && !visited[g->index]) {
                    double tau = transportRotation(h);
                    double sgn = (h == h->edge->he) ? 1.0 : -1.0;
                    double conn = sgn * phi(h->edge->index);
                    alphaGlobal[g->index] = alphaGlobal[fi] + tau - conn;

                    visited[g->index] = 1;
                    q.push(g->index);
                }
            }
            h = h->next;
        } while (h != h0);
    }

    // Map back to interior-only faceTheta_
    for (int gi = 0; gi < static_cast<int>(alphaGlobal.size()); ++gi) {
        int li = globalToLocal[gi];
        if (li >= 0) faceTheta_(li) = alphaGlobal[gi];
    }
}
