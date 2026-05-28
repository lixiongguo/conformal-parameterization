#include "HolomorphicOneForm.h"
#include <algorithm>
#include <cmath>
#include <queue>
#include <set>
#include <utility>

using namespace Eigen;

HolomorphicOneForm::HolomorphicOneForm(Mesh& mesh0)
    : Parameterization(mesh0)
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

    nV = static_cast<int>(mesh.vertices.size());
    nF = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (!f->isBoundary()) ++nF;
    }

    V.resize(nV, 3);
    F.resize(nF * 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        V.row(v->index) = v->position;
    }

    int fi = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        F[fi * 3] = f->he->vertex->index;
        F[fi * 3 + 1] = f->he->next->vertex->index;
        F[fi * 3 + 2] = f->he->next->next->vertex->index;
        ++fi;
    }

    for (int f = 0; f < nF; ++f) {
        for (int k = 0; k < 3; ++k) {
            int a = F[f * 3 + k];
            int b = F[f * 3 + (k + 1) % 3];
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
    nB = static_cast<int>(mesh.boundaries.size());
    vtxInc.assign(nV, {});

    for (int ei = 0; ei < nE; ++ei) {
        const EdgeInfo& e = edges[ei];
        vtxInc[e.v0].push_back({ei, +1});
        vtxInc[e.v1].push_back({ei, -1});
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
    const EdgeInfo& e = edges[edgeIdx];
    const int vj = (e.v0 == vi) ? e.v1 : e.v0;
    double w = 0.0;

    auto accumulate = [&](int fi) {
        if (fi < 0) return;
        int opp = -1;
        for (int k = 0; k < 3; ++k) {
            const int vv = F[fi * 3 + k];
            if (vv != e.v0 && vv != e.v1) {
                opp = vv;
                break;
            }
        }
        if (opp >= 0) {
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
    for (HalfEdgeIter start : mesh.boundaries) {
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
    if (nV == 0 || nE == 0) return;

    std::vector<int> parent(nV, -1);
    std::vector<int> parentEdge(nV, -1);
    std::vector<int> parentSign(nV, 0);
    std::vector<char> visited(nV, 0);
    std::queue<int> q;

    q.push(0);
    visited[0] = 1;

    while (!q.empty()) {
        const int u = q.front();
        q.pop();
        for (const auto& inc : vtxInc[u]) {
            const int ei = inc.first;
            const int s = inc.second;
            const EdgeInfo& e = edges[ei];
            const int v = (e.v0 == u) ? e.v1 : e.v0;
            if (visited[v]) continue;
            visited[v] = 1;
            parent[v] = u;
            parentEdge[v] = ei;
            parentSign[v] = (e.v0 == u) ? +1 : -1;
            q.push(v);
        }
    }

    for (int ei = 0; ei < nE; ++ei) {
        const EdgeInfo& e = edges[ei];
        const int a = e.v0;
        const int b = e.v1;
        if (parent[a] == b || parent[b] == a) continue;

        std::vector<std::pair<int, int>> cycle;
        cycle.push_back({ei, +1});

        int x = a;
        int y = b;
        std::set<int> seenA;
        while (x != 0 && seenA.insert(x).second) {
            if (parentEdge[x] < 0) break;
            cycle.push_back({parentEdge[x], parentSign[x]});
            x = parent[x];
        }
        std::set<int> seenB;
        while (y != 0 && seenB.insert(y).second) {
            if (parentEdge[y] < 0) break;
            cycle.push_back({parentEdge[y], -parentSign[y]});
            y = parent[y];
        }
        if (!cycle.empty()) cycles.push_back(cycle);
        if (static_cast<int>(cycles.size()) >= std::max(1, 2 * genus)) break;
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
            const double w = cotanWeightAtVertex(vi, inc.first);
            trips.emplace_back(row, inc.first, w * static_cast<double>(inc.second));
        }
        rhs.push_back(0.0);
        ++row;
    }

    std::vector<std::vector<std::pair<int, int>>> cycles;
    collectBoundaryCycles(cycles);
    if (cycles.empty() && genus > 0) {
        collectTreeCycles(cycles);
    }

    const int targetPeriodRows = std::max(0, 2 * genus);
    for (int i = 0; i < targetPeriodRows && i < static_cast<int>(cycles.size()); ++i) {
        for (const auto& term : cycles[i]) {
            trips.emplace_back(row, term.first, static_cast<double>(term.second));
        }
        rhs.push_back(0.0);
        ++row;
    }

    if (genus == 0) {
        trips.emplace_back(row, normalizeEdge, 1.0);
        rhs.push_back(1.0);
        ++row;
    } else if (row < nE) {
        trips.emplace_back(row, normalizeEdge, 1.0);
        rhs.push_back(1.0);
        ++row;
    }

    if (row == 0 || nE == 0) return false;

    SparseMatrix<double> A(row, nE);
    A.setFromTriplets(trips.begin(), trips.end());
    VectorXd b = Map<VectorXd>(rhs.data(), static_cast<int>(rhs.size()));

    const SparseMatrix<double> AtA = A.transpose() * A;
    const VectorXd Atb = A.transpose() * b;

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(AtA);
    if (solver.info() != Success) return false;

    omega = solver.solve(Atb);
    return solver.info() == Success;
}

double HolomorphicOneForm::triangleArea(int fi) const
{
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
    if (eab >= 0) vals(0) = edgeSign(a, b) * omega(eab);
    if (ebc >= 0) vals(1) = edgeSign(b, c) * omega(ebc);
    if (eca >= 0) vals(2) = edgeSign(c, a) * omega(eca);
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

        const Vector3d n = e01.cross(p2 - p0).normalized();
        Vector3d t1 = e01.normalized();
        Vector3d t2 = n.cross(t1).normalized();

        const Vector2d q1(l0, 0.0);
        const Vector2d q2(e12.dot(t1), e12.dot(t2));
        const Vector2d q0(0.0, 0.0);

        Matrix2d A;
        A.col(0) = (q1 - q0) / l0;
        A.col(1) = (q2 - q0) / l1;

        const Vector2d bVec(
            wVals(0) / l0,
            wVals(1) / l1);

        Vector2d grad = A.colPivHouseholderQr().solve(bVec);
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
            if (local[k].idx < 0) continue;
            const Vector2d t2d(local[k].tangent.dot(t1), local[k].tangent.dot(t2));
            const double starVal = starGrad.dot(t2d.normalized()) * local[k].len;
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
    std::vector<char> visited(nV, 0);
    std::queue<int> q;
    q.push(0);
    visited[0] = 1;
    mesh.vertices[0].uv = Vector2d::Zero();

    while (!q.empty()) {
        const int u = q.front();
        q.pop();
        const Vector2d base = mesh.vertices[u].uv;

        for (const auto& inc : vtxInc[u]) {
            const int ei = inc.first;
            const EdgeInfo& e = edges[ei];
            const int v = (e.v0 == u) ? e.v1 : e.v0;
            if (visited[v]) continue;

            visited[v] = 1;
            const int s = (e.v0 == u) ? +1 : -1;
            mesh.vertices[v].uv = base + Vector2d(s * omega(ei), s * starOmega(ei));
            q.push(v);
        }
    }

    for (int i = 0; i < nV; ++i) {
        if (!visited[i]) {
            mesh.vertices[i].uv = Vector2d::Zero();
        }
    }
}

void HolomorphicOneForm::parameterize()
{
    extractMesh();
    if (nV < 3 || nF < 1 || nE < 1) return;

    VectorXd omega;
    if (!solveHarmonic1Form(omega)) return;

    VectorXd starOmega;
    computeHodgeConjugate(omega, starOmega);

    integrateTree(omega, starOmega);
    normalize();
}
