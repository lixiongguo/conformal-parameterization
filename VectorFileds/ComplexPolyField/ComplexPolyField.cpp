#include "ComplexPolyField.h"
#include "vector_field_topology.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <vector>

using vf::kEps;
using vf::kPi;

ComplexPolyField::ComplexPolyField(Mesh& mesh, int N)
    : mesh_(mesh), N_(std::max(1, N)) {
    nFaces_ = 0;
    globalToLocalFace_.clear();
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (!f->isBoundary()) {
            globalToLocalFace_[f->index] = nFaces_++;
        }
    }

    nVerts_ = static_cast<int>(mesh_.vertices.size());
    z_.assign(static_cast<size_t>(nFaces_), std::complex<double>(1.0, 0.0));
    zRaw_ = z_;

    buildLocalFrames();
    buildAdjacency();
}

void ComplexPolyField::buildLocalFrames() {
    frames_.resize(static_cast<size_t>(nFaces_));

    int idx = 0;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;

        HalfEdgeCIter he = f->he;
        const Eigen::Vector3d& p0 = he->vertex->position;
        he = he->next;
        const Eigen::Vector3d& p1 = he->vertex->position;
        he = he->next;
        const Eigen::Vector3d& p2 = he->vertex->position;

        const vf::FaceBasis b = vf::buildFaceBasis(p0, p1, p2);
        LocalFrame& lf = frames_[static_cast<size_t>(idx)];
        lf.X = b.x;
        lf.Y = b.y;
        lf.N = b.n;
        ++idx;
    }
}

void ComplexPolyField::buildAdjacency() {
    interiorEdges_.clear();
    boundaryEdges_.clear();

    std::unordered_map<uint64_t, std::vector<std::pair<int, int>>> edgeFaces;

    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        const int floc = globalToLocalFace_.at(f->index);

        HalfEdgeCIter he = f->he;
        for (int e = 0; e < 3; ++e) {
            const int va = he->vertex->index;
            const int vb = he->next->vertex->index;
            const uint64_t k = vf::edgeKey(va, vb);
            edgeFaces[k].emplace_back(floc, va * 100000 + vb);
            he = he->next;
        }
    }

    for (auto& kv : edgeFaces) {
        auto& list = kv.second;
        if (list.size() == 2) {
            const int fi = list[0].first;
            const int fj = list[1].first;
            const int va_f = list[0].second / 100000;
            const int vb_f = list[0].second % 100000;
            const int va_g = list[1].second / 100000;
            const int vb_g = list[1].second % 100000;
            interiorEdges_.push_back({fi, fj, va_f, vb_f, va_g, vb_g});
        } else if (list.size() == 1) {
            const int fi = list[0].first;
            const int va = list[0].second / 100000;
            const int vb = list[0].second % 100000;
            boundaryEdges_.push_back({fi, va, vb});
        }
    }
}

std::complex<double> ComplexPolyField::edgeComplexConjPowN(
    int fi, int va, int vb) const {
    const auto& posA = mesh_.vertices[va].position;
    const auto& posB = mesh_.vertices[vb].position;
    const LocalFrame& lf = frames_[static_cast<size_t>(fi)];
    const vf::FaceBasis b{lf.X, lf.Y, lf.N};
    return vf::edgeConjPowN(posA, posB, b, N_);
}

void ComplexPolyField::generateBoundaryConstraints(
    std::vector<Eigen::Triplet<std::complex<double>>>& trips,
    Eigen::VectorXcd& b,
    int& row) const {
    for (const auto& be : boundaryEdges_) {
        const auto& pa = mesh_.vertices[be.va].position;
        const auto& pb = mesh_.vertices[be.vb].position;
        Eigen::Vector3d edgeTan = pb - pa;
        if (edgeTan.norm() < kEps) continue;
        edgeTan.normalize();

        const LocalFrame& lf = frames_[static_cast<size_t>(be.fi)];
        Eigen::Vector3d inPlane = edgeTan - edgeTan.dot(lf.N) * lf.N;
        if (inPlane.norm() < kEps) continue;
        inPlane.normalize();

        const vf::FaceBasis bfs{lf.X, lf.Y, lf.N};
        const std::complex<double> cN = vf::directionConjPowN(inPlane, bfs, N_);

        trips.emplace_back(row, be.fi, std::complex<double>(1.0, 0.0));
        b[row] = cN;
        ++row;
    }
}

void ComplexPolyField::generateInteriorEquations(
    std::vector<Eigen::Triplet<std::complex<double>>>& trips,
    int& row) const {
    for (const auto& ie : interiorEdges_) {
        const auto e_f = edgeComplexConjPowN(ie.fi, ie.va_f, ie.vb_f);
        const auto e_g = edgeComplexConjPowN(ie.fj, ie.va_g, ie.vb_g);

        trips.emplace_back(row, ie.fi, e_f);
        trips.emplace_back(row, ie.fj, -e_g);
        ++row;
    }
}

void ComplexPolyField::generateGaugeConstraint(
    std::vector<Eigen::Triplet<std::complex<double>>>& trips,
    Eigen::VectorXcd& b,
    int& row) const {
    if (nFaces_ <= 0) return;
    trips.emplace_back(row, 0, std::complex<double>(1.0, 0.0));
    b[row] = std::complex<double>(1.0, 0.0);
    ++row;
}

void ComplexPolyField::compute() {
    if (nFaces_ == 0) return;

    const int nEqInterior = static_cast<int>(interiorEdges_.size());
    const int nEqBoundary = static_cast<int>(boundaryEdges_.size());
    const bool needGauge = boundaryEdges_.empty() && nFaces_ > 0;
    const int totalRows = nEqInterior + nEqBoundary + (needGauge ? 1 : 0);

    if (totalRows == 0) {
        std::fill(z_.begin(), z_.end(), std::complex<double>(1.0, 0.0));
        return;
    }

    std::vector<Eigen::Triplet<std::complex<double>>> trips;
    trips.reserve(static_cast<size_t>(nEqInterior) * 2 + nEqBoundary + 1);

    Eigen::VectorXcd b = Eigen::VectorXcd::Zero(totalRows);
    int row = 0;

    generateInteriorEquations(trips, row);
    generateBoundaryConstraints(trips, b, row);
    if (needGauge) generateGaugeConstraint(trips, b, row);

    Eigen::SparseMatrix<std::complex<double>> A(totalRows, nFaces_);
    A.setFromTriplets(trips.begin(), trips.end());

    Eigen::SparseMatrix<std::complex<double>> ATA = A.adjoint() * A;
    Eigen::VectorXcd ATb = A.adjoint() * b;

    for (int i = 0; i < nFaces_; ++i) {
        ATA.coeffRef(i, i) += std::complex<double>(1e-8, 0.0);
    }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(ATA);
    if (solver.info() != Eigen::Success) {
        std::fill(z_.begin(), z_.end(), std::complex<double>(1.0, 0.0));
        return;
    }

    const Eigen::VectorXcd x = solver.solve(ATb);
    if (solver.info() != Eigen::Success) {
        std::fill(z_.begin(), z_.end(), std::complex<double>(1.0, 0.0));
        return;
    }

    for (int i = 0; i < nFaces_; ++i) {
        zRaw_[static_cast<size_t>(i)] = x[i];
        const double mag = std::abs(x[i]);
        z_[static_cast<size_t>(i)] =
            (mag < kEps) ? std::complex<double>(1.0, 0.0) : x[i] / mag;
    }
}

int ComplexPolyField::localFaceIndex(int globalFaceIndex) const {
    const auto it = globalToLocalFace_.find(globalFaceIndex);
    return (it == globalToLocalFace_.end()) ? -1 : it->second;
}

double ComplexPolyField::theta(int faceIndex) const {
    if (faceIndex < 0 || faceIndex >= nFaces_) return 0.0;
    const double ang = std::atan2(z_[static_cast<size_t>(faceIndex)].imag(),
                                  z_[static_cast<size_t>(faceIndex)].real());
    return ang / static_cast<double>(N_);
}

Eigen::Vector3d ComplexPolyField::direction(int faceIndex, int k) const {
    if (faceIndex < 0 || faceIndex >= nFaces_) return Eigen::Vector3d::UnitX();
    const double th = theta(faceIndex);
    const double step = (2.0 * kPi) / static_cast<double>(N_);
    const double a = th + k * step;
    const auto& f = frames_[static_cast<size_t>(faceIndex)];
    return std::cos(a) * f.X + std::sin(a) * f.Y;
}

std::vector<Eigen::Vector3d> ComplexPolyField::allDirections(int faceIndex) const {
    std::vector<Eigen::Vector3d> dirs;
    dirs.reserve(static_cast<size_t>(N_));
    for (int k = 0; k < N_; ++k) {
        dirs.push_back(direction(faceIndex, k));
    }
    return dirs;
}

int ComplexPolyField::windingIndex(int vertexIndex) const {
    if (vertexIndex < 0 || vertexIndex >= nVerts_) return -1;

    for (auto bIt = mesh_.boundaries.begin(); bIt != mesh_.boundaries.end(); ++bIt) {
        HalfEdgeCIter he = *bIt;
        do {
            if (he->vertex->index == vertexIndex) return -1;
            he = he->next;
        } while (he != *bIt);
    }

    struct FaceHop {
        int fi;
        int va;
        int vb;
    };
    std::vector<FaceHop> cycle;

    HalfEdgeCIter h0 = mesh_.vertices[vertexIndex].he;
    HalfEdgeCIter h = h0;
    do {
        if (!h->onBoundary) {
            const auto it = globalToLocalFace_.find(h->face->index);
            if (it != globalToLocalFace_.end()) {
                cycle.push_back({it->second, h->vertex->index, h->next->vertex->index});
            }
        }
        h = h->flip->next;
    } while (h != h0);

    if (cycle.size() < 2) return 0;

    auto sharedOrientedEdge = [](const FaceHop& a, const FaceHop& b,
                                 int& va, int& vb) -> bool {
        if (a.vb == b.va) {
            va = a.va;
            vb = a.vb;
            return true;
        }
        if (a.vb == b.vb) {
            va = a.va;
            vb = a.vb;
            return true;
        }
        if (a.va == b.va) {
            va = a.vb;
            vb = a.va;
            return true;
        }
        if (a.va == b.vb) {
            va = a.vb;
            vb = a.va;
            return true;
        }
        return false;
    };

    double totalJump = 0.0;

    for (size_t i = 0; i < cycle.size(); ++i) {
        const FaceHop& a = cycle[i];
        const FaceHop& b = cycle[(i + 1) % cycle.size()];

        int ea_va = 0, ea_vb = 0;
        if (!sharedOrientedEdge(a, b, ea_va, ea_vb)) continue;

        const auto e_a = edgeComplexConjPowN(a.fi, ea_va, ea_vb);
        const auto e_b = edgeComplexConjPowN(b.fi, ea_vb, ea_va);
        if (std::abs(e_b) < kEps) continue;

        const std::complex<double> transported =
            zRaw_[static_cast<size_t>(a.fi)] * e_a / e_b;
        const std::complex<double> ratio =
            zRaw_[static_cast<size_t>(b.fi)] / transported;

        const double ang = std::atan2(ratio.imag(), ratio.real());
        double best = ang;
        double bestDist = 1e10;
        for (int m = -2; m <= 2; ++m) {
            const double cand = ang + m * 2.0 * kPi;
            if (std::abs(cand) < bestDist) {
                bestDist = std::abs(cand);
                best = cand;
            }
        }
        totalJump += best;
    }

    const double index = totalJump / (2.0 * kPi);
    if (std::abs(index) < 0.35) return 0;
    const int w = static_cast<int>(std::round(index));
    if (w == 0) return 0;
    return (std::abs(index - w) < 0.12) ? w : 0;
}
