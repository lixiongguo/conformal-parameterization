#include "AbelJacobi.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>

using namespace Eigen;

namespace {

constexpr double kReg = 1e-8;
constexpr double kLatticeTol = 1e-4;

bool isFiniteVector(const VectorXd& v)
{
    for (int i = 0; i < v.size(); ++i) {
        if (!std::isfinite(v(i))) return false;
    }
    return true;
}

} // namespace

AbelJacobi::AbelJacobi(Mesh& mesh)
    : mesh_(mesh)
{
}

void AbelJacobi::extractMesh()
{
    edges_.clear();
    edgeMap_.clear();
    vtxInc_.clear();

    nV_ = static_cast<int>(mesh_.vertices.size());
    nF_ = 0;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (!f->isBoundary()) ++nF_;
    }

    V_.resize(nV_, 3);
    F_.resize(nF_ * 3);
    for (VertexCIter v = mesh_.vertices.begin(); v != mesh_.vertices.end(); ++v) {
        if (v->index < 0 || v->index >= nV_) continue;
        V_.row(v->index) = v->position;
    }

    int fi = 0;
    for (FaceCIter f = mesh_.faces.begin(); f != mesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        const int a = f->he->vertex->index;
        const int b = f->he->next->vertex->index;
        const int c = f->he->next->next->vertex->index;
        if (a < 0 || a >= nV_ || b < 0 || b >= nV_ || c < 0 || c >= nV_) continue;
        F_[fi * 3] = a;
        F_[fi * 3 + 1] = b;
        F_[fi * 3 + 2] = c;
        ++fi;
    }
    nF_ = fi;

    edgeMap_.clear();
    for (int f = 0; f < nF_; ++f) {
        for (int k = 0; k < 3; ++k) {
            const int a = F_[f * 3 + k];
            const int b = F_[f * 3 + (k + 1) % 3];
            if (a < 0 || a >= nV_ || b < 0 || b >= nV_) continue;
            const int v0 = std::min(a, b);
            const int v1 = std::max(a, b);
            const std::pair<int, int> key{v0, v1};
            auto it = edgeMap_.find(key);
            if (it == edgeMap_.end()) {
                EdgeInfo e;
                e.v0 = v0;
                e.v1 = v1;
                e.f0 = f;
                e.f1 = -1;
                e.length = (V_.row(v0) - V_.row(v1)).norm();
                const int idx = static_cast<int>(edges_.size());
                edgeMap_[key] = idx;
                edges_.push_back(e);
            } else {
                edges_[it->second].f1 = f;
            }
        }
    }

    nE_ = static_cast<int>(edges_.size());
    nB_ = static_cast<int>(mesh_.boundaries.size());
    vtxInc_.assign(nV_, {});

    for (int ei = 0; ei < nE_; ++ei) {
        const EdgeInfo& e = edges_[ei];
        if (e.v0 >= 0 && e.v0 < nV_) vtxInc_[e.v0].push_back({ei, +1});
        if (e.v1 >= 0 && e.v1 < nV_) vtxInc_[e.v1].push_back({ei, -1});
    }

    const int chi = nV_ - nE_ + nF_;
    genus_ = std::max(0, (2 - chi - nB_) / 2);
}

int AbelJacobi::findEdge(int a, int b) const
{
    if (a < 0 || b < 0 || a >= nV_ || b >= nV_) return -1;
    const int v0 = std::min(a, b);
    const int v1 = std::max(a, b);
    const auto it = edgeMap_.find({v0, v1});
    return it == edgeMap_.end() ? -1 : it->second;
}

int AbelJacobi::edgeSign(int from, int to) const
{
    const int ei = findEdge(from, to);
    if (ei < 0) return 0;
    const EdgeInfo& e = edges_[ei];
    if (from == e.v0 && to == e.v1) return +1;
    if (from == e.v1 && to == e.v0) return -1;
    return 0;
}

double AbelJacobi::cotan(const Vector3d& a, const Vector3d& b, const Vector3d& c) const
{
    const Vector3d u = a - b;
    const Vector3d v = c - b;
    const double dot = u.dot(v);
    const double cross = u.cross(v).norm();
    if (cross < 1e-12) return 0.0;
    return dot / cross;
}

double AbelJacobi::cotanWeightAtVertex(int vi, int edgeIdx) const
{
    if (vi < 0 || vi >= nV_ || edgeIdx < 0 || edgeIdx >= nE_) return 0.0;
    const EdgeInfo& e = edges_[edgeIdx];
    const int vj = (e.v0 == vi) ? e.v1 : e.v0;
    if (vj < 0 || vj >= nV_) return 0.0;
    double w = 0.0;

    auto accumulate = [&](int fi) {
        if (fi < 0 || fi >= nF_) return;
        int opp = -1;
        for (int k = 0; k < 3; ++k) {
            const int vv = F_[fi * 3 + k];
            if (vv != e.v0 && vv != e.v1) {
                opp = vv;
                break;
            }
        }
        if (opp >= 0 && opp < nV_) {
            w += cotan(V_.row(opp), V_.row(vi), V_.row(vj));
        }
    };

    accumulate(e.f0);
    accumulate(e.f1);
    return -0.5 * w;
}

void AbelJacobi::collectHomologyBasis()
{
    aCycles_.clear();
    bCycles_.clear();
    if (genus_ <= 0 || nB_ != 0) return;

    std::vector<topology::TreeCotreeBasis::Edge> topoEdges;
    topoEdges.reserve(static_cast<size_t>(nE_));
    for (int ei = 0; ei < nE_; ++ei) {
        const EdgeInfo& e = edges_[ei];
        topology::TreeCotreeBasis::Edge te;
        te.v0 = e.v0;
        te.v1 = e.v1;
        te.f0 = e.f0;
        te.f1 = e.f1;
        topoEdges.push_back(te);
    }

    auto basis = topology::TreeCotreeBasis::buildClosedMeshBasis(
        nV_, nF_, topoEdges, genus_,
        [&](int a, int b) { return findEdge(a, b); },
        [&](int from, int to) { return edgeSign(from, to); });

    const int g = genus_;
    for (int i = 0; i < g && i < static_cast<int>(basis.size()); ++i) {
        aCycles_.push_back(basis[static_cast<size_t>(i)]);
    }
    for (int i = g; i < 2 * g && i < static_cast<int>(basis.size()); ++i) {
        bCycles_.push_back(basis[static_cast<size_t>(i)]);
    }
}

bool AbelJacobi::solveHarmonicWithPeriods(
    const std::vector<double>& aPeriodRhs,
    Eigen::VectorXd& omega) const
{
    std::vector<Triplet<double>> trips;
    std::vector<double> rhs;
    int row = 0;

    for (int fi = 0; fi < nF_; ++fi) {
        const int a = F_[fi * 3];
        const int b = F_[fi * 3 + 1];
        const int c = F_[fi * 3 + 2];
        const int eab = findEdge(a, b);
        const int ebc = findEdge(b, c);
        const int eca = findEdge(c, a);
        if (eab < 0 || ebc < 0 || eca < 0) continue;

        trips.emplace_back(row, eab, edgeSign(a, b));
        trips.emplace_back(row, ebc, edgeSign(b, c));
        trips.emplace_back(row, eca, edgeSign(c, a));
        rhs.push_back(0.0);
        ++row;
    }

    for (int vi = 0; vi < nV_; ++vi) {
        for (const auto& inc : vtxInc_[vi]) {
            if (inc.first < 0 || inc.first >= nE_) continue;
            const double w = cotanWeightAtVertex(vi, inc.first);
            trips.emplace_back(row, inc.first, w * static_cast<double>(inc.second));
        }
        rhs.push_back(0.0);
        ++row;
    }

    for (int i = 0; i < static_cast<int>(aCycles_.size()) && i < static_cast<int>(aPeriodRhs.size()); ++i) {
        for (const auto& term : aCycles_[static_cast<size_t>(i)]) {
            if (term.first < 0 || term.first >= nE_) continue;
            trips.emplace_back(row, term.first, static_cast<double>(term.second));
        }
        rhs.push_back(aPeriodRhs[static_cast<size_t>(i)]);
        ++row;
    }

    for (int ei = 0; ei < nE_; ++ei) {
        trips.emplace_back(row, ei, kReg);
        rhs.push_back(0.0);
        ++row;
    }

    if (row == 0 || nE_ == 0) return false;
    if (nE_ > 12000) return false;

    SparseMatrix<double> A(row, nE_);
    A.setFromTriplets(trips.begin(), trips.end());
    VectorXd b = Map<VectorXd>(rhs.data(), static_cast<int>(rhs.size()));

    const SparseMatrix<double> AtA = A.transpose() * A;
    SparseMatrix<double> reg(nE_, nE_);
    reg.setIdentity();
    const SparseMatrix<double> AtAreg = AtA + kReg * reg;
    const VectorXd Atb = A.transpose() * b;

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(AtAreg);
    if (solver.info() != Success) return false;

    omega = solver.solve(Atb);
    return solver.info() == Success && omega.size() == nE_ && isFiniteVector(omega);
}

void AbelJacobi::computeHodgeConjugate(
    const Eigen::VectorXd& omega,
    Eigen::VectorXd& starOmega) const
{
    starOmega = VectorXd::Zero(nE_);
    std::vector<double> weightSum(nE_, 0.0);

    for (int fi = 0; fi < nF_; ++fi) {
        const int a = F_[fi * 3];
        const int b = F_[fi * 3 + 1];
        const int c = F_[fi * 3 + 2];
        if (a < 0 || a >= nV_ || b < 0 || b >= nV_ || c < 0 || c >= nV_) continue;

        const Vector3d p0 = V_.row(a);
        const Vector3d p1 = V_.row(b);
        const Vector3d p2 = V_.row(c);

        const Vector3d e01 = p1 - p0;
        const Vector3d e12 = p2 - p1;
        const Vector3d e20 = p0 - p2;
        const double l0 = e01.norm();
        const double l1 = e12.norm();
        const double l2 = e20.norm();
        const double area = 0.5 * e01.cross(p2 - p0).norm();
        if (area < 1e-14 || l0 < 1e-14 || l1 < 1e-14 || l2 < 1e-14) continue;

        Vector3d wVals = Vector3d::Zero();
        const int eab = findEdge(a, b);
        const int ebc = findEdge(b, c);
        const int eca = findEdge(c, a);
        if (eab >= 0) wVals(0) = edgeSign(a, b) * omega(eab);
        if (ebc >= 0) wVals(1) = edgeSign(b, c) * omega(ebc);
        if (eca >= 0) wVals(2) = edgeSign(c, a) * omega(eca);

        const Vector3d n = e01.cross(p2 - p0);
        const double nLen = n.norm();
        if (nLen < 1e-14) continue;
        const Vector3d nHat = n / nLen;
        Vector3d t1 = e01 / l0;
        Vector3d t2 = nHat.cross(t1);
        const double t2Len = t2.norm();
        if (t2Len < 1e-14) continue;
        t2 /= t2Len;

        const Vector2d q1(l0, 0.0);
        const Vector2d q2(e12.dot(t1), e12.dot(t2));

        Matrix2d matA;
        matA.col(0) = q1 / l0;
        matA.col(1) = q2 / l1;

        const Vector2d bVec(wVals(0) / l0, wVals(1) / l1);
        const Vector2d grad = matA.colPivHouseholderQr().solve(bVec);
        const Vector2d starGrad(-grad.y(), grad.x());

        struct LocalEdge {
            int idx;
            int from;
            int to;
            Vector3d tangent;
            double len;
        };
        const LocalEdge local[3] = {
            {eab, a, b, e01, l0},
            {ebc, b, c, e12, l1},
            {eca, c, a, e20, l2},
        };

        for (int k = 0; k < 3; ++k) {
            if (local[k].idx < 0 || local[k].idx >= nE_) continue;
            const double tLen = local[k].tangent.norm();
            if (tLen < 1e-14) continue;
            const Vector2d t2d(local[k].tangent.dot(t1), local[k].tangent.dot(t2));
            const double t2dLen = t2d.norm();
            if (t2dLen < 1e-14) continue;
            const double starVal = starGrad.dot(t2d / t2dLen) * local[k].len;
            const int s = edgeSign(local[k].from, local[k].to);
            starOmega(local[k].idx) += s * starVal;
            weightSum[local[k].idx] += 1.0;
        }
    }

    for (int ei = 0; ei < nE_; ++ei) {
        if (weightSum[ei] > 0.0) {
            starOmega(ei) /= weightSum[ei];
        }
    }
}

Eigen::Vector2d AbelJacobi::periodOnCycle(
    const HoloForm& form,
    const topology::TreeCotreeBasis::Cycle& cycle) const
{
    double re = 0.0;
    double im = 0.0;
    for (const auto& term : cycle) {
        const int ei = term.first;
        const int s = term.second;
        if (ei < 0 || ei >= nE_) continue;
        re += s * form.omega(ei);
        im += s * form.starOmega(ei);
    }
    return Vector2d(re, im);
}

void AbelJacobi::buildLattice()
{
    const int g = genus_;
    latticeGens_.resize(g, 4 * g);
    latticeGens_.setZero();

    if (static_cast<int>(aCycles_.size()) < g || static_cast<int>(bCycles_.size()) < g) {
        return;
    }

    // Column k = λ_{a_k}, column g+k = λ_{b_k} (each a g-vector in C^g).
    for (int k = 0; k < g; ++k) {
        for (int j = 0; j < g; ++j) {
            const Vector2d pa = periodOnCycle(canonicalBasis_[static_cast<size_t>(j)], aCycles_[static_cast<size_t>(k)]);
            const Vector2d pb = periodOnCycle(canonicalBasis_[static_cast<size_t>(j)], bCycles_[static_cast<size_t>(k)]);
            latticeGens_(j, 2 * k) = pa.x();
            latticeGens_(j, 2 * k + 1) = pa.y();
            latticeGens_(j, 2 * g + 2 * k) = pb.x();
            latticeGens_(j, 2 * g + 2 * k + 1) = pb.y();
        }
    }
}

bool AbelJacobi::build()
{
    ready_ = false;
    extractMesh();

    if (nV_ < 3 || nF_ < 1 || nE_ < 1 || genus_ <= 0 || nB_ != 0) {
        std::cerr << "[AbelJacobi] Requires closed mesh with genus >= 1 (genus="
                  << genus_ << ", boundaries=" << nB_ << ").\n";
        return false;
    }

    collectHomologyBasis();
    if (static_cast<int>(aCycles_.size()) < genus_ || static_cast<int>(bCycles_.size()) < genus_) {
        std::cerr << "[AbelJacobi] Failed to build homology basis.\n";
        return false;
    }

    canonicalBasis_.clear();
    canonicalBasis_.resize(static_cast<size_t>(genus_));

    for (int j = 0; j < genus_; ++j) {
        std::vector<double> rhs(static_cast<size_t>(genus_), 0.0);
        rhs[static_cast<size_t>(j)] = 1.0;

        HoloForm form;
        if (!solveHarmonicWithPeriods(rhs, form.omega)) {
            std::cerr << "[AbelJacobi] Failed to solve canonical form " << j << ".\n";
            return false;
        }
        computeHodgeConjugate(form.omega, form.starOmega);
        if (!isFiniteVector(form.starOmega)) return false;
        canonicalBasis_[static_cast<size_t>(j)] = std::move(form);
    }

    buildLattice();
    ready_ = true;
    return true;
}

bool AbelJacobi::checkPoincareHopf(
    int genus,
    int nBoundary,
    const std::vector<SingularPoint>& singularities)
{
    const int chi = 2 - 2 * genus - nBoundary;
    int sum = 0;
    for (const auto& s : singularities) {
        const int deg = (s.degree != 0) ? s.degree : (4 - s.valence);
        sum += deg;
    }
    return sum == 4 * chi;
}

Eigen::MatrixXd AbelJacobi::mapAlongPath(int baseVertex, int targetVertex) const
{
    MatrixXd mu = MatrixXd::Zero(genus_, 2);
    if (!ready_ || baseVertex < 0 || baseVertex >= nV_ ||
        targetVertex < 0 || targetVertex >= nV_) {
        return mu;
    }
    if (baseVertex == targetVertex) return mu;

    std::vector<int> parent(nV_, -1);
    std::vector<char> visited(nV_, 0);
    std::queue<int> q;
    q.push(baseVertex);
    visited[baseVertex] = 1;

    while (!q.empty()) {
        const int u = q.front();
        q.pop();
        if (u == targetVertex) break;
        for (const auto& inc : vtxInc_[u]) {
            const int ei = inc.first;
            if (ei < 0 || ei >= nE_) continue;
            const EdgeInfo& e = edges_[ei];
            const int v = (e.v0 == u) ? e.v1 : e.v0;
            if (v < 0 || v >= nV_ || visited[v]) continue;
            visited[v] = 1;
            parent[v] = u;
            q.push(v);
        }
    }

    if (!visited[targetVertex]) return mu;

    std::vector<int> path;
    for (int v = targetVertex; v != baseVertex && v >= 0; v = parent[v]) {
        path.push_back(v);
    }
    path.push_back(baseVertex);
    std::reverse(path.begin(), path.end());

    for (size_t i = 0; i + 1 < path.size(); ++i) {
        const int from = path[i];
        const int to = path[i + 1];
        for (int j = 0; j < genus_; ++j) {
            const Vector2d d = integrateEdge(from, to, canonicalBasis_[static_cast<size_t>(j)]);
            mu(j, 0) += d.x();
            mu(j, 1) += d.y();
        }
    }
    return mu;
}

Eigen::Vector2d AbelJacobi::integrateEdge(int from, int to, const HoloForm& form) const
{
    const int ei = findEdge(from, to);
    if (ei < 0 || ei >= nE_) return Vector2d::Zero();
    const int s = edgeSign(from, to);
    return Vector2d(s * form.omega(ei), s * form.starOmega(ei));
}

Eigen::MatrixXd AbelJacobi::abelJacobiMap(int baseVertex, int targetVertex) const
{
    return mapAlongPath(baseVertex, targetVertex);
}

bool AbelJacobi::nearestLatticePoint(
    const Eigen::MatrixXd& mu,
    Eigen::VectorXi& coeffs,
    double& residual) const
{
    const int g = genus_;
    const int nGen = 2 * g;
    if (mu.rows() != g || mu.cols() != 2 || latticeGens_.rows() != g) {
        residual = std::numeric_limits<double>::infinity();
        return false;
    }

    MatrixXd A = MatrixXd::Zero(2 * g, 2 * g);
    VectorXd b(2 * g);
    for (int j = 0; j < g; ++j) {
        b(2 * j) = mu(j, 0);
        b(2 * j + 1) = mu(j, 1);
    }
    for (int k = 0; k < nGen; ++k) {
        const int colRe = (k < g) ? (2 * k) : (2 * g + 2 * (k - g));
        const int colIm = colRe + 1;
        for (int j = 0; j < g; ++j) {
            A(2 * j, k) = latticeGens_(j, colRe);
            A(2 * j + 1, k) = latticeGens_(j, colIm);
        }
    }

    const VectorXd z = A.colPivHouseholderQr().solve(b);
    coeffs.resize(nGen);
    VectorXd zRound(nGen);
    for (int k = 0; k < nGen; ++k) {
        coeffs(k) = static_cast<int>(std::llround(z(k)));
        zRound(k) = static_cast<double>(coeffs(k));
    }
    residual = (A * zRound - b).norm();
    return residual < kLatticeTol;
}

AbelJacobi::DivisorCheck AbelJacobi::checkDivisor(
    int baseVertex,
    const std::vector<SingularPoint>& singularities) const
{
    DivisorCheck out;
    out.poincareHopfOk = checkPoincareHopf(genus_, nB_, singularities);
    out.muD = MatrixXd::Zero(genus_, 2);

    if (!ready_) return out;

    for (const auto& s : singularities) {
        if (s.vertexIndex < 0 || s.vertexIndex >= nV_) continue;
        const int deg = (s.degree != 0) ? s.degree : (4 - s.valence);
        const MatrixXd mp = mapAlongPath(baseVertex, s.vertexIndex);
        out.muD += static_cast<double>(deg) * mp;
    }

    out.abelJacobiOk = nearestLatticePoint(out.muD, out.latticeCoeffs, out.latticeResidual);
    return out;
}

bool AbelJacobi::integrateHolomorphic(
    int baseVertex,
    const Eigen::VectorXd& coeffRe,
    const Eigen::VectorXd& coeffIm,
    std::vector<Eigen::Vector2d>& uv) const
{
    if (!ready_ || coeffRe.size() != genus_ || coeffIm.size() != genus_) return false;

    HoloForm combo;
    combo.omega = VectorXd::Zero(nE_);
    combo.starOmega = VectorXd::Zero(nE_);
    for (int j = 0; j < genus_; ++j) {
        combo.omega += coeffRe(j) * canonicalBasis_[static_cast<size_t>(j)].omega
                     - coeffIm(j) * canonicalBasis_[static_cast<size_t>(j)].starOmega;
        combo.starOmega += coeffRe(j) * canonicalBasis_[static_cast<size_t>(j)].starOmega
                         + coeffIm(j) * canonicalBasis_[static_cast<size_t>(j)].omega;
    }

    uv.assign(static_cast<size_t>(nV_), Vector2d::Zero());
    std::vector<char> visited(nV_, 0);
    std::queue<int> q;
    q.push(baseVertex);
    visited[baseVertex] = 1;

    while (!q.empty()) {
        const int u = q.front();
        q.pop();
        const Vector2d base = uv[static_cast<size_t>(u)];

        for (const auto& inc : vtxInc_[u]) {
            const int ei = inc.first;
            if (ei < 0 || ei >= nE_) continue;
            const EdgeInfo& e = edges_[ei];
            const int v = (e.v0 == u) ? e.v1 : e.v0;
            if (v < 0 || v >= nV_ || visited[v]) continue;

            visited[v] = 1;
            const Vector2d d = integrateEdge(u, v, combo);
            uv[static_cast<size_t>(v)] = base + d;
            q.push(v);
        }
    }
    return true;
}

bool AbelJacobi::quantizeAndParameterize(int baseVertex)
{
    if (!ready_) return false;

    VectorXd cRe = VectorXd::Ones(genus_);
    VectorXd cIm = VectorXd::Zero(genus_);

    if (genus_ == 1 && !bCycles_.empty()) {
        const Vector2d pb = periodOnCycle(canonicalBasis_[0], bCycles_[0]);
        const double target = std::max(std::abs(pb.x()), std::abs(pb.y()));
        if (target > 1e-12) {
            const double scale = std::round(target) / target;
            cRe(0) = scale;
        }
    } else {
        for (int j = 0; j < genus_; ++j) {
            for (int k = 0; k < genus_ && k < static_cast<int>(bCycles_.size()); ++k) {
                const Vector2d pb = periodOnCycle(canonicalBasis_[static_cast<size_t>(j)], bCycles_[static_cast<size_t>(k)]);
                if (std::abs(pb.x()) > 1e-12) {
                    cRe(j) = std::round(pb.x()) / pb.x();
                }
                if (std::abs(pb.y()) > 1e-12) {
                    cIm(j) = std::round(pb.y()) / pb.y();
                }
            }
        }
    }

    std::vector<Vector2d> uv;
    if (!integrateHolomorphic(baseVertex, cRe, cIm, uv)) return false;

    for (VertexIter v = mesh_.vertices.begin(); v != mesh_.vertices.end(); ++v) {
        const int i = v->index;
        if (i >= 0 && i < nV_) {
            v->uv = uv[static_cast<size_t>(i)];
        }
    }
    return true;
}
