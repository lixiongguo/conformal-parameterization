#include "ComplexPolyField.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <unordered_map>
#include <vector>

namespace {

constexpr double kPi  = 3.14159265358979323846;
constexpr double kEps = 1e-20;

uint64_t edgeKey(int a, int b) {
    uint32_t lo = static_cast<uint32_t>(std::min(a, b));
    uint32_t hi = static_cast<uint32_t>(std::max(a, b));
    return (static_cast<uint64_t>(lo) << 32) | hi;
}

} // namespace

// ================================================================
//  Constructor
// ================================================================

ComplexPolyField::ComplexPolyField(Mesh& mesh, int N)
    : mesh_(mesh), N_(N)
{
    // Count interior faces
    nFaces_ = 0;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f)
        if (!f->isBoundary()) ++nFaces_;

    nVerts_ = static_cast<int>(mesh_.vertices.size());
    z_.resize(nFaces_, std::complex<double>(0, 0));

    buildLocalFrames();
    buildAdjacency();
}


// ================================================================
//  buildLocalFrames — per‑face {X, Y, N} from first edge & normal
// ================================================================

void ComplexPolyField::buildLocalFrames() {
    frames_.resize(nFaces_);

    int idx = 0;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;

        HalfEdgeCIter he = f->he;
        const Eigen::Vector3d& p0 = he->vertex->position;
        he = he->next;
        const Eigen::Vector3d& p1 = he->vertex->position;
        he = he->next;
        const Eigen::Vector3d& p2 = he->vertex->position;

        LocalFrame& lf = frames_[idx];

        lf.X = p1 - p0;
        if (lf.X.norm() < kEps) lf.X = Eigen::Vector3d::UnitX();
        lf.X.normalize();

        lf.N = (p1 - p0).cross(p2 - p0);
        if (lf.N.norm() < kEps) lf.N = Eigen::Vector3d::UnitZ();
        lf.N.normalize();

        lf.Y = lf.N.cross(lf.X);
        if (lf.Y.norm() < kEps) lf.Y = Eigen::Vector3d::UnitY();
        lf.Y.normalize();

        ++idx;
    }
}


// ================================================================
//  buildAdjacency — classify each edge as interior or boundary
// ================================================================

void ComplexPolyField::buildAdjacency() {
    interiorEdges_.clear();
    boundaryEdges_.clear();

    // Map: vertex face index → face index in our numbering
    std::unordered_map<int, int> globalToLocal;
    {
        int loc = 0;
        for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
            if (!f->isBoundary()) globalToLocal[f->index] = loc++;
        }
    }

    // For each interior edge, collect the incident interior faces
    std::unordered_map<uint64_t, std::vector<std::pair<int,int>>> edgeFaces;
    //                                                   localFace, {va,vb}

    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        int floc = globalToLocal[f->index];

        HalfEdgeCIter he = f->he;
        for (int e = 0; e < 3; ++e) {
            int va = he->vertex->index;
            int vb = he->next->vertex->index;
            uint64_t k = edgeKey(va, vb);
            edgeFaces[k].emplace_back(floc, va * 100000 + vb);
            he = he->next;
        }
    }

    for (auto& kv : edgeFaces) {
        auto& list = kv.second;
        if (list.size() == 2) {
            int fi = list[0].first, fj = list[1].first;
            int va = list[0].second / 100000;
            int vb = list[0].second % 100000;
            interiorEdges_.push_back({fi, fj, va, vb});
        } else if (list.size() == 1) {
            int fi = list[0].first;
            int va = list[0].second / 100000;
            int vb = list[0].second % 100000;
            boundaryEdges_.push_back({fi, va, vb});
        }
    }
}


// ================================================================
//  edgeComplexConjPowN — conj(edge_in_face) raised to N_th power
// ================================================================

std::complex<double> ComplexPolyField::edgeComplexConjPowN(
    int fi, int va, int vb) const
{
    const auto& posA = mesh_.vertices[va].position;
    const auto& posB = mesh_.vertices[vb].position;
    Eigen::Vector3d e = posA - posB;
    double l = e.norm();
    if (l < kEps) e = frames_[fi].X;
    else e /= l;

    const std::complex<double> c(e.dot(frames_[fi].X),
                                  e.dot(frames_[fi].Y));
    // conj(c) in the face: c̅ = (X·e) - i(Y·e)
    std::complex<double> conjC(c.real(), -c.imag());
    double ang = std::atan2(conjC.imag(), conjC.real()) * static_cast<double>(N_);
    return {std::cos(ang), std::sin(ang)};
}


// ================================================================
//  generateBoundaryConstraints
// ================================================================

void ComplexPolyField::generateBoundaryConstraints(
    std::vector<Eigen::Triplet<std::complex<double>>>& trips,
    Eigen::VectorXcd& b,
    int& row) const
{
    for (const auto& be : boundaryEdges_) {
        // Edge direction e = vb - va
        const auto& pa = mesh_.vertices[be.va].position;
        const auto& pb = mesh_.vertices[be.vb].position;
        Eigen::Vector3d dir = pb - pa;
        if (dir.norm() < kEps) continue;
        dir.normalize();

        // Boundary outward normal: N × edge_direction
        const auto& N = frames_[be.fi].N;
        Eigen::Vector3d outward = N.cross(dir);
        if (outward.norm() < kEps) continue;
        outward.normalize();

        // Express in face local coordinates: conj(dir_perp) ≈  cos + i sin
        std::complex<double> constraint(
            outward.dot(frames_[be.fi].X),
            outward.dot(frames_[be.fi].Y));

        // Raise to N_th power
        double ang = std::atan2(constraint.imag(), constraint.real())
                     * static_cast<double>(N_);
        std::complex<double> cN(std::cos(ang), std::sin(ang));

        trips.emplace_back(row, be.fi, std::complex<double>(1.0, 0.0));
        b[row] = cN;
        ++row;
    }
}


// ================================================================
//  generateInteriorEquations
// ================================================================

void ComplexPolyField::generateInteriorEquations(
    std::vector<Eigen::Triplet<std::complex<double>>>& trips,
    int& row) const
{
    for (const auto& ie : interiorEdges_) {
        auto e_f = edgeComplexConjPowN(ie.fi, ie.va, ie.vb);
        auto e_g = edgeComplexConjPowN(ie.fj, ie.va, ie.vb);

        trips.emplace_back(row, ie.fi,  e_f);
        trips.emplace_back(row, ie.fj, -e_g);
        ++row;
    }
}


// ================================================================
//  compute — main solver
// ================================================================

void ComplexPolyField::compute() {
    if (nFaces_ == 0) return;

    int nEqInterior = static_cast<int>(interiorEdges_.size());
    int nEqBoundary = static_cast<int>(boundaryEdges_.size());
    int totalRows    = nEqInterior + nEqBoundary;

    if (totalRows == 0) {
        // No edges – every face gets the same trivial field
        for (int i = 0; i < nFaces_; ++i)
            z_[i] = std::complex<double>(1.0, 0.0);
        return;
    }

    std::vector<Eigen::Triplet<std::complex<double>>> trips;
    trips.reserve(nEqInterior * 2 + nEqBoundary);

    Eigen::VectorXcd b = Eigen::VectorXcd::Zero(totalRows);
    int row = 0;

    generateInteriorEquations(trips, row);
    generateBoundaryConstraints(trips, b, row);

    // Build sparse matrix
    Eigen::SparseMatrix<std::complex<double>> A(totalRows, nFaces_);
    A.setFromTriplets(trips.begin(), trips.end());

    // Solve normal equations: AᵀA z = Aᵀb
    Eigen::SparseMatrix<std::complex<double>> ATA = A.adjoint() * A;
    Eigen::VectorXcd ATb = A.adjoint() * b;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(ATA);
    if (solver.info() != Eigen::Success) {
        // Fallback: identity field
        for (int i = 0; i < nFaces_; ++i)
            z_[i] = std::complex<double>(1.0, 0.0);
        return;
    }

    Eigen::VectorXcd x = solver.solve(ATb);
    if (solver.info() != Eigen::Success) {
        for (int i = 0; i < nFaces_; ++i)
            z_[i] = std::complex<double>(1.0, 0.0);
        return;
    }

    // Recover per‑face complex representation z = e^{i·N·θ}
    for (int i = 0; i < nFaces_; ++i) {
        double mag = std::abs(x[i]);
        if (mag < kEps) {
            z_[i] = std::complex<double>(1.0, 0.0);
        } else {
            z_[i] = x[i] / mag;   // unit complex
        }
    }
}


// ================================================================
//  theta — per‑face angle
// ================================================================

double ComplexPolyField::theta(int faceIndex) const {
    if (faceIndex < 0 || faceIndex >= nFaces_) return 0.0;
    double ang = std::atan2(z_[faceIndex].imag(), z_[faceIndex].real());
    return ang / static_cast<double>(N_);
}


// ================================================================
//  direction — one of N symmetric 3‑D directions
// ================================================================

Eigen::Vector3d ComplexPolyField::direction(int faceIndex, int k) const {
    if (faceIndex < 0 || faceIndex >= nFaces_) return Eigen::Vector3d::UnitX();
    double th = theta(faceIndex);
    double step = (2.0 * kPi) / static_cast<double>(N_);
    double a = th + k * step;
    const auto& f = frames_[faceIndex];
    return std::cos(a) * f.X + std::sin(a) * f.Y;
}


// ================================================================
//  allDirections
// ================================================================

std::vector<Eigen::Vector3d> ComplexPolyField::allDirections(int faceIndex) const {
    std::vector<Eigen::Vector3d> dirs;
    dirs.reserve(N_);
    for (int k = 0; k < N_; ++k)
        dirs.push_back(direction(faceIndex, k));
    return dirs;
}


// ================================================================
//  windingIndex — singularity detection around a vertex
//  Returns: ±1 = singularity, 0 = regular, -1 = boundary
// ================================================================

int ComplexPolyField::windingIndex(int vertexIndex) const {
    if (vertexIndex < 0 || vertexIndex >= nVerts_) return -1;

    const auto& v = mesh_.vertices[vertexIndex];

    // Check if boundary vertex (no singularity computation)
    for (auto bIt = mesh_.boundaries.begin(); bIt != mesh_.boundaries.end(); ++bIt) {
        HalfEdgeCIter he = *bIt;
        do {
            if (he->vertex->index == vertexIndex) return -1;
            he = he->next;
        } while (he != *bIt);
    }

    // Collect incident interior faces in order around the vertex
    struct FaceEntry { int fi; Eigen::Vector3d center; double angle; };
    std::vector<FaceEntry> incident;

    // Walk around the vertex via half‑edges
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        HalfEdgeCIter he = f->he;
        for (int e = 0; e < 3; ++e) {
            if (he->vertex->index == vertexIndex) {
                // Find the face index in our local numbering
                // Simple linear scan (nFaces_ is small)
                for (int fi = 0; fi < nFaces_; ++fi) {
                    // Check if this is the right face
                    FaceCIter fc = mesh_.faces.begin();
                    int cnt = 0;
                    for (; fc != mesh_.faces.end(); ++fc) {
                        if (!fc->isBoundary()) {
                            if (cnt == fi) break;
                            ++cnt;
                        }
                    }
                    if (fc != mesh_.faces.end() && fc->index == f->index) {
                        // Compute face center
                        HalfEdgeCIter hec = fc->he;
                        const auto& p0 = hec->vertex->position;
                        hec = hec->next;
                        const auto& p1 = hec->vertex->position;
                        hec = hec->next;
                        const auto& p2 = hec->vertex->position;
                        Eigen::Vector3d c = (p0 + p1 + p2) / 3.0;
                        incident.push_back({fi, c, 0.0});
                        break;
                    }
                }
                break;
            }
            he = he->next;
        }
    }

    if (incident.size() < 2) return 0;

    // Sort incident faces by angle around the vertex
    const Eigen::Vector3d& vPos = v.position;

    // Compute average normal for projection
    Eigen::Vector3d avgN = Eigen::Vector3d::Zero();
    for (auto& in : incident) {
        avgN += frames_[in.fi].N;
    }
    if (avgN.norm() < kEps) return 0;
    avgN.normalize();

    // Build tangent basis
    Eigen::Vector3d T1, T2;
    {
        if (std::abs(avgN.x()) < 0.9) T1 = avgN.cross(Eigen::Vector3d::UnitX());
        else T1 = avgN.cross(Eigen::Vector3d::UnitY());
        T1.normalize();
        T2 = avgN.cross(T1);
        T2.normalize();
    }

    for (auto& in : incident) {
        Eigen::Vector3d d = in.center - vPos;
        double u = d.dot(T1);
        double v = d.dot(T2);
        in.angle = std::atan2(v, u);
    }
    std::sort(incident.begin(), incident.end(),
              [](const FaceEntry& a, const FaceEntry& b) { return a.angle < b.angle; });

    // Compute total winding
    double totalAngle = 0.0;
    for (size_t i = 0; i < incident.size(); ++i) {
        size_t j = (i + 1) % incident.size();
        double t0 = theta(incident[i].fi);
        double t1 = theta(incident[j].fi);

        // Find the branch of t1 closest to t0
        double step = (2.0 * kPi) / N_;
        double bestDiff = 1e10, bestT1 = t1;
        for (int m = -N_; m <= N_; ++m) {
            double cand = t1 + m * step;
            double diff = std::abs(cand - t0);
            if (diff < bestDiff) { bestDiff = diff; bestT1 = cand; }
        }
        totalAngle += bestT1 - t0;
    }

    // Winding = totalAngle / (2π) — should be near an integer
    double wind = totalAngle / (2.0 * kPi);
    int w = static_cast<int>(std::round(wind));

    return (std::abs(wind - w) < 0.2) ? w : 0;
}
