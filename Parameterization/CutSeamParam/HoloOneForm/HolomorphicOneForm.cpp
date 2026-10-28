#include "HolomorphicOneForm.h"
#include "TreeCotreeBasis.h"
#include <Eigen/Dense>
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
            // Cotangent weight uses the angle opposite the edge, not the angle at vi.
            w += cotan(V.row(vi), V.row(opp), V.row(vj));
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
    std::vector<std::vector<std::pair<int, int>>>& cycles)
{
    cycles.clear();
    treeCotree_ = TreeCotree();
    treeCotree_.vertexCount = nV;
    treeCotree_.faceCount = nF;
    treeCotree_.genus = genus;
    if (nV == 0 || nF == 0 || nE == 0 || genus <= 0) return;
    if (nB != 0) return;

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

    topology::TreeCotreeBasis::Parts parts;
    auto basis = topology::TreeCotreeBasis::buildClosedMeshBasis(
        nV, nF, topoEdges, genus,
        [&](int a, int b) { return findEdge(a, b); },
        [&](int from, int to) { return edgeSign(from, to); },
        &parts);

    auto endpoints = [&](const std::vector<int>& edgeIndices, std::vector<std::pair<int, int>>& dest) {
        dest.reserve(edgeIndices.size());
        for (int edgeIndex : edgeIndices) {
            if (edgeIndex < 0 || edgeIndex >= nE) continue;
            dest.emplace_back(edges[edgeIndex].v0, edges[edgeIndex].v1);
        }
    };
    endpoints(parts.primalTree, treeCotree_.primalTree);
    endpoints(parts.dualTree, treeCotree_.dualTree);
    endpoints(parts.homologyEdges, treeCotree_.homologyEdges);
    treeCotree_.visitedVertices = parts.visitedVertices;
    treeCotree_.visitedFaces = parts.visitedFaces;

    if ((int)basis.size() != 2 * genus) {
        std::cerr << "[HolomorphicOneForm] tree-cotree produced " << basis.size()
                  << " cycles, expected " << (2 * genus) << " for genus=" << genus << ".\n";
    }

    cycles.reserve(basis.size());
    treeCotree_.cycles.reserve(basis.size());
    for (const auto& cyc : basis) {
        std::vector<std::pair<int, int>> out;
        std::vector<std::pair<int, int>> loop;
        out.reserve(cyc.size());
        loop.reserve(cyc.size());
        for (const auto& se : cyc) {
            out.push_back(se);
            if (se.first >= 0 && se.first < nE) loop.emplace_back(edges[se.first].v0, edges[se.first].v1);
        }
        if (!out.empty()) cycles.push_back(std::move(out));
        if (!loop.empty()) treeCotree_.cycles.push_back(std::move(loop));
    }
}

bool HolomorphicOneForm::solveHarmonic1Form(
    const std::vector<std::vector<std::pair<int, int>>>& cycles,
    const std::vector<double>& periodRhs,
    Eigen::VectorXd& omega)
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

    const int nPeriods = std::min(static_cast<int>(cycles.size()), static_cast<int>(periodRhs.size()));
    for (int i = 0; i < nPeriods; ++i) {
        bool wrote = false;
        for (const auto& term : cycles[static_cast<size_t>(i)]) {
            if (term.first < 0 || term.first >= nE) continue;
            trips.emplace_back(row, term.first, static_cast<double>(term.second));
            wrote = true;
        }
        if (!wrote) continue;
        rhs.push_back(periodRhs[static_cast<size_t>(i)]);
        ++row;
    }

    for (int ei = 0; ei < nE; ++ei) {
        trips.emplace_back(row, ei, kReg);
        rhs.push_back(0.0);
        ++row;
    }

    if (row == 0 || nE == 0) return false;

    // Avoid forming A^T A for very large meshes (WASM memory limit).
    // 与 AbelJacobi 的上限对齐(12000)：允许 torus/double-torus/3holes 等中等网格。
    if (nE > 12000) return false;

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

double HolomorphicOneForm::faceWedge(const Vector3d& left, const Vector3d& right) const
{
    return (left(0) * (right(1) - right(2)) + left(1) * (right(2) - right(0)) + left(2) * (right(0) - right(1))) / 6.0;
}

double HolomorphicOneForm::faceStarWedge(int face, const Vector3d& left, const Vector3d& right) const
{
    if (face < 0 || face >= nF) return 0.0;
    const int a = F[face * 3];
    const int b = F[face * 3 + 1];
    const int c = F[face * 3 + 2];
    if (a < 0 || a >= nV || b < 0 || b >= nV || c < 0 || c >= nV) return 0.0;
    const double area = triangleArea(face);
    if (area < 1e-14) return 0.0;
    const double l0 = (V.row(b) - V.row(a)).squaredNorm();
    const double l1 = (V.row(c) - V.row(b)).squaredNorm();
    const double l2 = (V.row(a) - V.row(c)).squaredNorm();
    Matrix3d metric;
    metric << -4.0 * l0, l0 + l1 - l2, l0 + l2 - l1,
              l1 + l0 - l2, -4.0 * l1, l1 + l2 - l0,
              l2 + l0 - l1, l2 + l1 - l0, -4.0 * l2;
    metric /= 24.0 * area;
    return left.transpose() * metric * right;
}

bool HolomorphicOneForm::computeHodgeConjugate(
    const std::vector<VectorXd>& basis,
    const VectorXd& omega,
    VectorXd& starOmega) const
{
    const int nBasis = static_cast<int>(basis.size());
    starOmega = VectorXd::Zero(nE);
    if (nBasis == 0) return false;

    MatrixXd wedge(nBasis, nBasis);
    VectorXd starAgainst(nBasis);
    wedge.setZero();
    starAgainst.setZero();

    for (int face = 0; face < nF; ++face) {
        std::vector<Vector3d> values(static_cast<size_t>(nBasis));
        for (int i = 0; i < nBasis; ++i) faceEdgeValues(face, basis[static_cast<size_t>(i)], values[static_cast<size_t>(i)]);
        Vector3d target;
        faceEdgeValues(face, omega, target);
        for (int i = 0; i < nBasis; ++i) {
            starAgainst(i) += faceStarWedge(face, values[static_cast<size_t>(i)], target);
            for (int j = 0; j < nBasis; ++j) {
                wedge(i, j) += faceWedge(values[static_cast<size_t>(i)], values[static_cast<size_t>(j)]);
            }
        }
    }

    const VectorXd lambda = wedge.fullPivLu().solve(starAgainst);
    if (!isFiniteVector(lambda)) return false;
    for (int i = 0; i < nBasis; ++i) starOmega += lambda(i) * basis[static_cast<size_t>(i)];
    return isFiniteVector(starOmega);
}

void HolomorphicOneForm::integrateTree(
    const Eigen::VectorXd& omega,
    const Eigen::VectorXd& starOmega)
{
    int nCut = 0;
    for (VertexCIter vertex = cutMesh_.vertices.begin(); vertex != cutMesh_.vertices.end(); ++vertex) {
        nCut = std::max(nCut, vertex->index + 1);
    }
    if (nCut <= 0) return;

    std::vector<std::vector<std::pair<int, Vector2d>>> neighbors(static_cast<size_t>(nCut));
    for (HalfEdgeCIter he = cutMesh_.halfEdges.begin(); he != cutMesh_.halfEdges.end(); ++he) {
        if (cutMesh_.isCutBoundaryEdge(he->edge)) continue;
        const int from = he->vertex->index;
        const int to = he->next->vertex->index;
        if (from < 0 || to < 0 || from >= nCut || to >= nCut) continue;
        const int originalFrom = cutMesh_.originalVertexIndex(from);
        const int originalTo = cutMesh_.originalVertexIndex(to);
        const int edge = findEdge(originalFrom, originalTo);
        if (edge < 0 || edge >= omega.size() || edge >= starOmega.size()) continue;
        const int sign = edgeSign(originalFrom, originalTo);
        neighbors[static_cast<size_t>(from)].push_back(
            {to, Vector2d(sign * omega(edge), sign * starOmega(edge))});
    }

    std::vector<Vector2d> uv(static_cast<size_t>(nCut), Vector2d::Zero());
    std::vector<char> visited(static_cast<size_t>(nCut), 0);
    std::queue<int> q;
    for (int start = 0; start < nCut; ++start) {
        if (visited[static_cast<size_t>(start)]) continue;
        visited[static_cast<size_t>(start)] = 1;
        q.push(start);
        while (!q.empty()) {
            const int u = q.front();
            q.pop();
            const Vector2d base = uv[static_cast<size_t>(u)];
            for (const auto& step : neighbors[static_cast<size_t>(u)]) {
                if (visited[static_cast<size_t>(step.first)]) continue;
                visited[static_cast<size_t>(step.first)] = 1;
                uv[static_cast<size_t>(step.first)] = base + step.second;
                q.push(step.first);
            }
        }
    }

    double signedArea = 0.0;
    for (FaceCIter face = cutMesh_.faces.begin(); face != cutMesh_.faces.end(); ++face) {
        if (face->isBoundary()) continue;
        const int a = face->he->vertex->index;
        const int b = face->he->next->vertex->index;
        const int c = face->he->next->next->vertex->index;
        if (a < 0 || b < 0 || c < 0 || a >= nCut || b >= nCut || c >= nCut) continue;
        const Vector2d u = uv[static_cast<size_t>(b)] - uv[static_cast<size_t>(a)];
        const Vector2d v = uv[static_cast<size_t>(c)] - uv[static_cast<size_t>(a)];
        signedArea += u.x() * v.y() - u.y() * v.x();
    }
    if (signedArea < 0.0) {
        for (Vector2d& value : uv) value.y() = -value.y();
    }

    for (VertexIter vertex = cutMesh_.vertices.begin(); vertex != cutMesh_.vertices.end(); ++vertex) {
        if (vertex->index >= 0 && vertex->index < nCut) {
            vertex->uv = uv[static_cast<size_t>(vertex->index)];
        }
    }
}

void HolomorphicOneForm::parameterize()
{
    treeCotree_ = TreeCotree();
    extractMesh();
    if (nV < 3 || nF < 1 || nE < 1 || genus <= 0) {
        return;
    }

    std::vector<std::vector<std::pair<int, int>>> cycles;
    collectTreeCycles(cycles);
    const int nBasis = 2 * genus;
    if (static_cast<int>(cycles.size()) != nBasis) {
        return;
    }

    std::vector<VectorXd> basis(static_cast<size_t>(nBasis));
    for (int j = 0; j < nBasis; ++j) {
        std::vector<double> periods(static_cast<size_t>(nBasis), 0.0);
        periods[static_cast<size_t>(j)] = 1.0;
        if (!solveHarmonic1Form(cycles, periods, basis[static_cast<size_t>(j)])) {
            return;
        }
    }

    VectorXd starOmega;
    if (!computeHodgeConjugate(basis, basis[0], starOmega)) {
        return;
    }

    std::vector<std::pair<int, int>> seams;
    for (const auto& cycle : cycles) {
        for (const auto& term : cycle) {
            if (term.first < 0 || term.first >= nE) continue;
            seams.emplace_back(edges[term.first].v0, edges[term.first].v1);
        }
    }
    setSeamVertexPairs(seams);
    if (!prepareCutMesh()) {
        return;
    }

    integrateTree(basis[0], starOmega);

    double spread = 0.0;
    for (VertexCIter v = cutMesh_.vertices.begin(); v != cutMesh_.vertices.end(); ++v) {
        spread = std::max(spread, std::max(std::abs(v->uv.x()), std::abs(v->uv.y())));
    }
    if (spread < 1e-12) {
        return;
    }

    copyCutMeshUvsToOriginal();
    normalize();
}
