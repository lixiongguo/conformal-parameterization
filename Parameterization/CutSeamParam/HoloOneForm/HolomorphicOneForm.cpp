#include "HolomorphicOneForm.h"
#include "Lscm.h"
#include "TreeCotreeBasis.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <queue>
#include <set>
#include <utility>

using namespace Eigen;

namespace {

constexpr double kReg = 1e-8;

bool isFiniteVector(const VectorXd& v)
{
    for (int i = 0; i < v.size(); ++i) {
        if (!std::isfinite(v(i))) return false;
    }
    return true;
}

} // namespace

HolomorphicOneForm::HolomorphicOneForm(Mesh& mesh0)
    : CutSeamParameterization(mesh0)
{
}

double HolomorphicOneForm::cotan(const Vector3d& a, const Vector3d& b, const Vector3d& c) const
{
    const Vector3d u = a - b;
    const Vector3d v = c - b;
    const double dot = u.dot(v);
    const double cross = u.cross(v).norm();
    if (cross < 1e-12) return 0.0;
    return dot / cross;
}

void HolomorphicOneForm::extractMesh()
{
    edges.clear();
    edgeMap.clear();
    vtxInc.clear();

    nV = static_cast<int>(cutMesh_.vertices.size());
    nF = 0;
    for (FaceCIter f = cutMesh_.faces.begin(); f != cutMesh_.faces.end(); ++f) {
        if (!f->isBoundary()) ++nF;
    }

    V.resize(nV, 3);
    F.resize(nF * 3);
    for (VertexCIter v = cutMesh_.vertices.begin(); v != cutMesh_.vertices.end(); ++v) {
        if (v->index < 0 || v->index >= nV) continue;
        V.row(v->index) = v->position;
    }

    int fi = 0;
    for (FaceCIter f = cutMesh_.faces.begin(); f != cutMesh_.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        const int a = f->he->vertex->index;
        const int b = f->he->next->vertex->index;
        const int c = f->he->next->next->vertex->index;
        if (a < 0 || a >= nV || b < 0 || b >= nV || c < 0 || c >= nV) continue;
        F[fi * 3] = a;
        F[fi * 3 + 1] = b;
        F[fi * 3 + 2] = c;
        ++fi;
    }
    nF = fi;

    for (int f = 0; f < nF; ++f) {
        for (int k = 0; k < 3; ++k) {
            int a = F[f * 3 + k];
            int b = F[f * 3 + (k + 1) % 3];
            if (a < 0 || a >= nV || b < 0 || b >= nV) continue;
            int v0 = std::min(a, b);
            int v1 = std::max(a, b);
            const std::pair<int, int> key{v0, v1};
            auto it = edgeMap.find(key);
            if (it == edgeMap.end()) {
                EdgeInfo e;
                e.v0 = v0;
                e.v1 = v1;
                e.f0 = f;
                e.f1 = -1;
                e.length = (V.row(v0) - V.row(v1)).norm();
                const int idx = static_cast<int>(edges.size());
                edgeMap[key] = idx;
                edges.push_back(e);
            } else {
                edges[it->second].f1 = f;
            }
        }
    }

    nE = static_cast<int>(edges.size());
    nB = static_cast<int>(cutMesh_.boundaries.size());
    vtxInc.assign(nV, {});

    for (int ei = 0; ei < nE; ++ei) {
        const EdgeInfo& e = edges[ei];
        if (e.v0 >= 0 && e.v0 < nV) vtxInc[e.v0].push_back({ei, +1});
        if (e.v1 >= 0 && e.v1 < nV) vtxInc[e.v1].push_back({ei, -1});
    }

    const int chi = nV - nE + nF;
    genus = std::max(0, (2 - chi - nB) / 2);

    normalizeEdge = 0;
    double bestLen = 0.0;
    for (int ei = 0; ei < nE; ++ei) {
        if (edges[ei].length > bestLen) {
            bestLen = edges[ei].length;
            normalizeEdge = ei;
        }
    }
}

int HolomorphicOneForm::findEdge(int a, int b) const
{
    if (a < 0 || b < 0 || a >= nV || b >= nV) return -1;
    const int v0 = std::min(a, b);
    const int v1 = std::max(a, b);
    const auto it = edgeMap.find({v0, v1});
    return it == edgeMap.end() ? -1 : it->second;
}

int HolomorphicOneForm::edgeSign(int from, int to) const
{
    const int ei = findEdge(from, to);
    if (ei < 0) return 0;
    const EdgeInfo& e = edges[ei];
    if (from == e.v0 && to == e.v1) return +1;
    if (from == e.v1 && to == e.v0) return -1;
    return 0;
}

double HolomorphicOneForm::cotanWeightAtVertex(int vi, int edgeIdx) const
{
    if (vi < 0 || vi >= nV || edgeIdx < 0 || edgeIdx >= nE) return 0.0;
    const EdgeInfo& e = edges[edgeIdx];
    const int vj = (e.v0 == vi) ? e.v1 : e.v0;
    if (vj < 0 || vj >= nV) return 0.0;
    double w = 0.0;

    auto accumulate = [&](int fi) {
        if (fi < 0 || fi >= nF) return;
        int opp = -1;
        for (int k = 0; k < 3; ++k) {
            const int vv = F[fi * 3 + k];
            if (vv != e.v0 && vv != e.v1) {
                opp = vv;
                break;
            }
        }
        if (opp >= 0 && opp < nV) {
            w += cotan(V.row(opp), V.row(vi), V.row(vj));
        }
    };

    accumulate(e.f0);
    accumulate(e.f1);
    return -0.5 * w;
}

void HolomorphicOneForm::collectBoundaryCycles(
    std::vector<std::vector<std::pair<int, int>>>& cycles) const
{
    cycles.clear();
    for (HalfEdgeIter start : cutMesh_.boundaries) {
        std::vector<std::pair<int, int>> cycle;
        HalfEdgeCIter he = start;
        do {
            const int from = he->vertex->index;
            const int to = he->next->vertex->index;
            const int ei = findEdge(from, to);
            if (ei >= 0) {
                cycle.push_back({ei, edgeSign(from, to)});
            }
            he = he->next;
        } while (he != start);
        if (!cycle.empty()) cycles.push_back(cycle);
    }
}

void HolomorphicOneForm::collectTreeCycles(
    std::vector<std::vector<std::pair<int, int>>>& cycles) const
{
    cycles.clear();
    if (nV == 0 || nF == 0 || nE == 0 || genus <= 0) return;
    if (nB != 0) return; // This helper is for closed meshes only.

    std::vector<topology::TreeCotreeBasis::Edge> topoEdges;
    topoEdges.reserve(static_cast<size_t>(nE));
    for (int ei = 0; ei < nE; ++ei) {
        const EdgeInfo& e = edges[ei];
        topology::TreeCotreeBasis::Edge te;
        te.v0 = e.v0;
        te.v1 = e.v1;
        te.f0 = e.f0;
        te.f1 = e.f1;
        topoEdges.push_back(te);
    }

    auto basis = topology::TreeCotreeBasis::buildClosedMeshBasis(
        nV, nF, topoEdges, genus,
        [&](int a, int b) { return findEdge(a, b); },
        [&](int from, int to) { return edgeSign(from, to); });

    if ((int)basis.size() != 2 * genus) {
        std::cerr << "[HolomorphicOneForm] tree-cotree produced " << basis.size()
                  << " cycles, expected " << (2 * genus) << " for genus=" << genus << ".\n";
    }

    cycles.reserve(basis.size());
    for (const auto& cyc : basis) {
        std::vector<std::pair<int, int>> out;
        out.reserve(cyc.size());
        for (const auto& se : cyc) out.push_back(se);
        if (!out.empty()) cycles.push_back(std::move(out));
    }
}

bool HolomorphicOneForm::solveHarmonic1Form(Eigen::VectorXd& omega)
{
    std::vector<Triplet<double>> trips;
    std::vector<double> rhs;
    int row = 0;

    for (int fi = 0; fi < nF; ++fi) {
        const int a = F[fi * 3];
        const int b = F[fi * 3 + 1];
        const int c = F[fi * 3 + 2];
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

    for (int vi = 0; vi < nV; ++vi) {
        for (const auto& inc : vtxInc[vi]) {
            if (inc.first < 0 || inc.first >= nE) continue;
            const double w = cotanWeightAtVertex(vi, inc.first);
            trips.emplace_back(row, inc.first, w * static_cast<double>(inc.second));
        }
        rhs.push_back(0.0);
        ++row;
    }

    std::vector<std::vector<std::pair<int, int>>> cycles;
    if (genus > 0) {
        collectBoundaryCycles(cycles);
        if (cycles.empty()) {
            collectTreeCycles(cycles);
        }
    }

    const int targetPeriodRows = std::max(0, 2 * genus);
    for (int i = 0; i < targetPeriodRows && i < static_cast<int>(cycles.size()); ++i) {
        for (const auto& term : cycles[i]) {
            if (term.first < 0 || term.first >= nE) continue;
            trips.emplace_back(row, term.first, static_cast<double>(term.second));
        }
        rhs.push_back(0.0);
        ++row;
    }

    if (normalizeEdge >= 0 && normalizeEdge < nE) {
        trips.emplace_back(row, normalizeEdge, 1.0);
        rhs.push_back(1.0);
        ++row;
    }

    for (int ei = 0; ei < nE; ++ei) {
        trips.emplace_back(row, ei, kReg);
        rhs.push_back(0.0);
        ++row;
    }

    if (row == 0 || nE == 0) return false;

    // Avoid forming A^T A for very large meshes (WASM memory limit).
    if (nE > 8000) return false;

    SparseMatrix<double> A(row, nE);
    A.setFromTriplets(trips.begin(), trips.end());
    VectorXd b = Map<VectorXd>(rhs.data(), static_cast<int>(rhs.size()));

    const SparseMatrix<double> AtA = A.transpose() * A;
    SparseMatrix<double> reg(nE, nE);
    reg.setIdentity();
    const SparseMatrix<double> AtAreg = AtA + kReg * reg;
    const VectorXd Atb = A.transpose() * b;

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(AtAreg);
    if (solver.info() != Success) return false;

    omega = solver.solve(Atb);
    return solver.info() == Success && omega.size() == nE && isFiniteVector(omega);
}

double HolomorphicOneForm::triangleArea(int fi) const
{
    if (fi < 0 || fi >= nF) return 0.0;
    const Vector3d p0 = V.row(F[fi * 3]);
    const Vector3d p1 = V.row(F[fi * 3 + 1]);
    const Vector3d p2 = V.row(F[fi * 3 + 2]);
    return 0.5 * (p1 - p0).cross(p2 - p0).norm();
}

void HolomorphicOneForm::faceEdgeValues(int fi, const Eigen::VectorXd& omega, Eigen::Vector3d& vals) const
{
    const int a = F[fi * 3];
    const int b = F[fi * 3 + 1];
    const int c = F[fi * 3 + 2];
    const int eab = findEdge(a, b);
    const int ebc = findEdge(b, c);
    const int eca = findEdge(c, a);
    vals = Vector3d::Zero();
    if (eab >= 0 && eab < omega.size()) vals(0) = edgeSign(a, b) * omega(eab);
    if (ebc >= 0 && ebc < omega.size()) vals(1) = edgeSign(b, c) * omega(ebc);
    if (eca >= 0 && eca < omega.size()) vals(2) = edgeSign(c, a) * omega(eca);
}

void HolomorphicOneForm::computeHodgeConjugate(
    const Eigen::VectorXd& omega,
    Eigen::VectorXd& starOmega) const
{
    starOmega = VectorXd::Zero(nE);
    std::vector<double> weightSum(nE, 0.0);

    for (int fi = 0; fi < nF; ++fi) {
        const int a = F[fi * 3];
        const int b = F[fi * 3 + 1];
        const int c = F[fi * 3 + 2];
        if (a < 0 || a >= nV || b < 0 || b >= nV || c < 0 || c >= nV) continue;

        const Vector3d p0 = V.row(a);
        const Vector3d p1 = V.row(b);
        const Vector3d p2 = V.row(c);

        const Vector3d e01 = p1 - p0;
        const Vector3d e12 = p2 - p1;
        const Vector3d e20 = p0 - p2;
        const double l0 = e01.norm();
        const double l1 = e12.norm();
        const double l2 = e20.norm();
        const double area = triangleArea(fi);
        if (area < 1e-14 || l0 < 1e-14 || l1 < 1e-14 || l2 < 1e-14) continue;

        Vector3d wVals;
        faceEdgeValues(fi, omega, wVals);

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

        const int eab = findEdge(a, b);
        const int ebc = findEdge(b, c);
        const int eca = findEdge(c, a);

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
            if (local[k].idx < 0 || local[k].idx >= nE) continue;
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

    for (int ei = 0; ei < nE; ++ei) {
        if (weightSum[ei] > 0.0) {
            starOmega(ei) /= weightSum[ei];
        }
    }
}

void HolomorphicOneForm::integrateTree(
    const Eigen::VectorXd& omega,
    const Eigen::VectorXd& starOmega)
{
    std::vector<Vector2d> uv(static_cast<size_t>(nV), Vector2d::Zero());
    std::vector<char> visited(nV, 0);
    std::queue<int> q;

    q.push(0);
    visited[0] = 1;

    while (!q.empty()) {
        const int u = q.front();
        q.pop();
        if (u < 0 || u >= nV) continue;
        const Vector2d base = uv[static_cast<size_t>(u)];

        for (const auto& inc : vtxInc[u]) {
            const int ei = inc.first;
            if (ei < 0 || ei >= nE || ei >= omega.size() || ei >= starOmega.size()) continue;
            const EdgeInfo& e = edges[ei];
            const int v = (e.v0 == u) ? e.v1 : e.v0;
            if (v < 0 || v >= nV || visited[v]) continue;

            visited[v] = 1;
            const int s = (e.v0 == u) ? +1 : -1;
            uv[static_cast<size_t>(v)] = base + Vector2d(s * omega(ei), s * starOmega(ei));
            q.push(v);
        }
    }

    for (VertexIter v = cutMesh_.vertices.begin(); v != cutMesh_.vertices.end(); ++v) {
        const int i = v->index;
        if (i >= 0 && i < nV) {
            v->uv = uv[static_cast<size_t>(i)];
        }
    }
}

void HolomorphicOneForm::parameterize()
{
    if (!prepareCutMesh()) {
        return;
    }

    extractMesh();
    if (nV < 3 || nF < 1 || nE < 1) {
        Lscm fallback(cutMesh_);
        fallback.parameterize();
        copyCutMeshUvsToOriginal();
        return;
    }

    VectorXd omega;
    if (!solveHarmonic1Form(omega)) {
        Lscm fallback(cutMesh_);
        fallback.parameterize();
        copyCutMeshUvsToOriginal();
        return;
    }

    VectorXd starOmega;
    computeHodgeConjugate(omega, starOmega);
    if (!isFiniteVector(starOmega)) {
        Lscm fallback(cutMesh_);
        fallback.parameterize();
        copyCutMeshUvsToOriginal();
        return;
    }

    integrateTree(omega, starOmega);

    double spread = 0.0;
    for (VertexCIter v = cutMesh_.vertices.begin(); v != cutMesh_.vertices.end(); ++v) {
        spread = std::max(spread, std::max(std::abs(v->uv.x()), std::abs(v->uv.y())));
    }
    if (spread < 1e-12) {
        Lscm fallback(cutMesh_);
        fallback.parameterize();
        copyCutMeshUvsToOriginal();
        return;
    }

    copyCutMeshUvsToOriginal();
    normalize();
}
