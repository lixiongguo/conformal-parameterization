#include "IncrementalFlattening.h"
#include "GaussianCurvature.h"
#include "TreeCotreeBasis.h"
#include "Tutte.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <set>
#include <vector>

using namespace Eigen;

namespace {

constexpr double kReg = 1e-10;
constexpr double kAngleEps = 1e-12;

bool isFiniteVec(const VectorXd& v)
{
    for (int i = 0; i < v.size(); ++i) {
        if (!std::isfinite(v(i))) return false;
    }
    return true;
}

} // namespace

IncrementalFlattening::IncrementalFlattening(Mesh& mesh0)
    : Parameterization(mesh0)
{
}

int IncrementalFlattening::seamEdgeCount() const
{
    int c = 0;
    for (char b : m_isSeamEdge) {
        if (b) ++c;
    }
    return c;
}

int IncrementalFlattening::cutVertexCount() const
{
    return m_cutMesh.nCutVerts;
}

bool IncrementalFlattening::isMeshBoundaryVertex(int vi) const
{
    if (vi < 0 || vi >= nV) return false;
    return mesh.vertices[vi].isBoundary();
}

Eigen::Matrix2d IncrementalFlattening::quarterTurnMatrix(int k)
{
  k = ((k % 4) + 4) % 4;
  static const Matrix2d R90 = (Matrix2d() << 0.0, -1.0, 1.0, 0.0).finished();
  Matrix2d R = Matrix2d::Identity();
  for (int i = 0; i < k; ++i) {
    R = R90 * R;
  }
  return R;
}

int IncrementalFlattening::findEdge(int a, int b) const
{
    if (a < 0 || b < 0 || a >= nV || b >= nV) return -1;
    const int v0 = std::min(a, b);
    const int v1 = std::max(a, b);
    const auto it = edgeMap.find({v0, v1});
    return it == edgeMap.end() ? -1 : it->second;
}

int IncrementalFlattening::edgeSign(int from, int to) const
{
    const int ei = findEdge(from, to);
    if (ei < 0) return 0;
    const EdgeInfo& e = edges[ei];
    if (from == e.v0 && to == e.v1) return +1;
    if (from == e.v1 && to == e.v0) return -1;
    return 0;
}

double IncrementalFlattening::cotanWeight(int vi, int vj, int opp) const
{
    if (vi < 0 || vi >= nV || vj < 0 || vj >= nV || opp < 0 || opp >= nV) return 0.0;
    const Vector3d a = V.row(vi) - V.row(opp);
    const Vector3d b = V.row(vj) - V.row(opp);
    const double cross = a.cross(b).norm();
    if (cross < kAngleEps) return 0.0;
    return a.dot(b) / cross;
}

bool IncrementalFlattening::extractMeshData()
{
    edges.clear();
    edgeMap.clear();
    vtxInc.clear();

    nV = static_cast<int>(mesh.vertices.size());
    nF = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (!f->isBoundary()) ++nF;
    }
    nB = static_cast<int>(mesh.boundaries.size());

    if (nV < 3 || nF < 1) return false;

    V.resize(nV, 3);
    F.resize(nF * 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        if (v->index >= 0 && v->index < nV) {
            V.row(v->index) = v->position;
        }
    }

    int fi = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        const int a = f->he->vertex->index;
        const int b = f->he->next->vertex->index;
        const int c = f->he->next->next->vertex->index;
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
            const int v0 = std::min(a, b);
            const int v1 = std::max(a, b);
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
    vtxInc.assign(nV, {});
    for (int ei = 0; ei < nE; ++ei) {
        const EdgeInfo& e = edges[ei];
        vtxInc[e.v0].push_back({ei, +1});
        vtxInc[e.v1].push_back({ei, -1});
    }

    const int chi = nV - nE + nF;
    genus = std::max(0, (2 - chi - nB) / 2);

    buildCotanLaplacian(m_laplace);
    buildVertexAreas(m_vertexAreas);
    m_K0 = geometry::gaussianCurvatureAngleDeficit(mesh);
    m_K = m_K0;
    m_phi = VectorXd::Zero(nV);
    return true;
}

void IncrementalFlattening::buildVertexAreas(Eigen::VectorXd& areas) const
{
    areas = VectorXd::Zero(nV);
    for (int f = 0; f < nF; ++f) {
        const int a = F[f * 3];
        const int b = F[f * 3 + 1];
        const int c = F[f * 3 + 2];
        const double area = 0.5 * (V.row(b) - V.row(a)).cross(V.row(c) - V.row(a)).norm();
        if (!std::isfinite(area) || area <= 0.0) continue;
        const double third = area / 3.0;
        areas(a) += third;
        areas(b) += third;
        areas(c) += third;
    }
    for (int i = 0; i < nV; ++i) {
        areas(i) = std::max(areas(i), 1e-16);
    }
}

void IncrementalFlattening::buildCotanLaplacian(Eigen::SparseMatrix<double>& L) const
{
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nV);

    for (int ei = 0; ei < nE; ++ei) {
        const EdgeInfo& e = edges[ei];
        const int vi = e.v0;
        const int vj = e.v1;

        auto opposite = [&](int f) {
            if (f < 0) return -1;
            for (int k = 0; k < 3; ++k) {
                const int vv = F[f * 3 + k];
                if (vv != vi && vv != vj) return vv;
            }
            return -1;
        };

        const int o1 = opposite(e.f0);
        const int o2 = opposite(e.f1);
        double w = 0.0;
        if (o1 >= 0) w += 0.5 * cotanWeight(vi, vj, o1);
        if (o2 >= 0) w += 0.5 * cotanWeight(vi, vj, o2);
        if (w <= 0.0) w = kReg;

        diag(vi) += w;
        diag(vj) += w;
        trips.emplace_back(vi, vj, -w);
        trips.emplace_back(vj, vi, -w);
    }

    for (int i = 0; i < nV; ++i) {
        trips.emplace_back(i, i, diag(i));
    }

    L.resize(nV, nV);
    L.setFromTriplets(trips.begin(), trips.end());
}

bool IncrementalFlattening::solveConstrainedPhi(const std::vector<int>& constraintVerts,
                                                const Eigen::VectorXd& rhs,
                                                Eigen::VectorXd& phi) const
{
    const int m = static_cast<int>(constraintVerts.size());
    if (m == 0) return false;

    std::vector<Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(m) * 8);
    for (int r = 0; r < m; ++r) {
        const int vi = constraintVerts[r];
        if (vi < 0 || vi >= nV) continue;
        for (SparseMatrix<double>::InnerIterator it(m_laplace, vi); it; ++it) {
            trips.emplace_back(r, it.col(), it.value());
        }
    }

    SparseMatrix<double> Lc(m, nV);
    Lc.setFromTriplets(trips.begin(), trips.end());

    VectorXd invA = VectorXd::Ones(nV);
    for (int i = 0; i < nV; ++i) {
        invA(i) = 1.0 / (2.0 * m_vertexAreas(i));
    }

    const SparseMatrix<double> LcInvA = Lc * invA.asDiagonal();
    const SparseMatrix<double> S = LcInvA * Lc.transpose();

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(S);
    if (solver.info() != Success) return false;

    VectorXd lambda = solver.solve(rhs);
    if (!isFiniteVec(lambda)) return false;

    phi = invA.asDiagonal() * (Lc.transpose() * lambda);
    return isFiniteVec(phi);
}

double IncrementalFlattening::roundToHalfPi(double k) const
{
    return 0.5 * M_PI * std::round(2.0 * k / M_PI);
}

bool IncrementalFlattening::identifyIsolatedCones(const std::vector<bool>& inOmega,
                                                  std::vector<int>& cones) const
{
    cones.clear();
    for (int vi = 0; vi < nV; ++vi) {
        if (inOmega[vi]) continue;
        if (mesh.vertices[vi].isBoundary()) continue;

        bool isolated = true;
        HalfEdgeCIter h = mesh.vertices[vi].he;
        do {
            const int nb = h->flip->vertex->index;
            if (!inOmega[nb] && nb != vi) {
                isolated = false;
                break;
            }
            h = h->flip->next;
        } while (h != mesh.vertices[vi].he);

        if (isolated) cones.push_back(vi);
    }
    return !cones.empty();
}

void IncrementalFlattening::flatteningPhase()
{
    double maxAbsK = 0.0;
    for (int i = 0; i < nV; ++i) {
        maxAbsK = std::max(maxAbsK, std::abs(m_K0(i)));
    }
    double epsilon = (m_epsilon0 > 0.0) ? m_epsilon0 : std::max(1e-4, 0.02 * maxAbsK);

    m_phi.setZero();
    m_K = m_K0;
    m_coneVerts.clear();

    for (int iter = 0; iter < m_maxFlattenIters; ++iter) {
        std::vector<bool> inOmega(nV, false);
        std::vector<int> omegaVerts;
        omegaVerts.reserve(nV);

        for (int i = 0; i < nV; ++i) {
            if (mesh.vertices[i].isBoundary()) continue;
            if (std::abs(m_K(i)) < epsilon) {
                inOmega[i] = true;
                omegaVerts.push_back(i);
            }
        }

        if (omegaVerts.empty()) break;

        VectorXd rhs(omegaVerts.size());
        for (size_t r = 0; r < omegaVerts.size(); ++r) {
            rhs(static_cast<int>(r)) = -m_K(omegaVerts[r]);
        }

        VectorXd phiStep;
        if (!solveConstrainedPhi(omegaVerts, rhs, phiStep)) break;
        m_phi += phiStep;

        const VectorXd lapPhi = m_laplace * m_phi;
        m_K = m_K0 - lapPhi;

        std::vector<int> cones;
        if (identifyIsolatedCones(inOmega, cones)) {
            m_coneVerts = cones;
            break;
        }

        epsilon *= m_epsilonGrowth;
    }

    if (m_coneVerts.empty()) {
        int best = -1;
        double bestVal = -1.0;
        for (int i = 0; i < nV; ++i) {
            if (mesh.vertices[i].isBoundary()) continue;
            const double val = std::abs(m_K(i));
            if (val > bestVal) {
                bestVal = val;
                best = i;
            }
        }
        if (best >= 0) m_coneVerts.push_back(best);
    }
}

void IncrementalFlattening::roundingPhase()
{
    if (m_coneVerts.empty()) return;

    std::set<int> fixed;
    std::map<int, double> targetK;
    for (int vi : m_coneVerts) {
        targetK[vi] = m_K(vi);
    }

    auto buildConstraints = [&](std::vector<int>& constraints, VectorXd& rhs) {
        constraints.clear();
        for (int vi = 0; vi < nV; ++vi) {
            if (mesh.vertices[vi].isBoundary()) continue;
            const bool isCone = targetK.count(vi) > 0;
            if (isCone && fixed.count(vi)) {
                constraints.push_back(vi);
                rhs.conservativeResize(rhs.size() + 1);
                rhs(rhs.size() - 1) = targetK[vi] - m_K0(vi);
            } else if (!isCone) {
                constraints.push_back(vi);
                rhs.conservativeResize(rhs.size() + 1);
                rhs(rhs.size() - 1) = -m_K0(vi);
            }
        }
    };

    const int maxRound = static_cast<int>(m_coneVerts.size()) + 5;
    for (int roundIter = 0; roundIter < maxRound; ++roundIter) {
        int pick = -1;
        double bestErr = std::numeric_limits<double>::max();
        for (const auto& kv : targetK) {
            if (fixed.count(kv.first)) continue;
            const double rounded = roundToHalfPi(kv.second);
            const double err = std::abs(kv.second - rounded);
            if (err < bestErr) {
                bestErr = err;
                pick = kv.first;
            }
        }
        if (pick < 0) break;

        fixed.insert(pick);
        targetK[pick] = roundToHalfPi(targetK[pick]);

        std::vector<int> constraints;
        VectorXd rhs(0);
        buildConstraints(constraints, rhs);

        VectorXd phiStep;
        if (!solveConstrainedPhi(constraints, rhs, phiStep)) break;
        m_phi += phiStep;
        m_K = m_K0 - m_laplace * m_phi;

        for (const auto& kv : targetK) {
            if (fixed.count(kv.first)) {
                targetK[kv.first] = m_K(kv.first);
            }
        }
    }

    m_coneK.resize(m_coneVerts.size());
    for (size_t i = 0; i < m_coneVerts.size(); ++i) {
        const int vi = m_coneVerts[i];
        m_coneK(static_cast<int>(i)) = roundToHalfPi(m_K(vi));
        m_K(vi) = m_coneK(static_cast<int>(i));
        targetK[vi] = m_coneK(static_cast<int>(i));
    }

    std::vector<std::vector<std::pair<int, int>>> homology;
    collectHomologyCycles(homology);

    VectorXd vertexHolonomy = VectorXd::Zero(nV);
    for (size_t i = 0; i < m_coneVerts.size(); ++i) {
        vertexHolonomy(m_coneVerts[i]) = m_coneK(static_cast<int>(i));
    }

    VectorXd omegaProbe;
    solveEdgeOmegas(vertexHolonomy, {}, omegaProbe, false);

    std::vector<std::pair<std::vector<std::pair<int, int>>, double>> extra;
    for (const auto& cyc : homology) {
        double h = 0.0;
        for (const auto& se : cyc) {
            if (se.first >= 0 && se.first < omegaProbe.size()) {
                h += static_cast<double>(se.second) * omegaProbe(se.first);
            }
        }
        const double hRound = roundToHalfPi(h);
        if (std::abs(h - hRound) > 1e-3) {
            extra.push_back({cyc, hRound});
        }
    }

    if (!extra.empty()) {
        solveEdgeOmegas(vertexHolonomy, extra, m_omega, false);
    }
}

void IncrementalFlattening::collectHomologyCycles(
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
            if (ei >= 0) cycle.push_back({ei, edgeSign(from, to)});
            he = he->next;
        } while (he != start);
        if (!cycle.empty()) cycles.push_back(std::move(cycle));
    }

    if (genus <= 0 || nB != 0) return;

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

    for (const auto& cyc : basis) {
        std::vector<std::pair<int, int>> out;
        out.reserve(cyc.size());
        for (const auto& se : cyc) out.push_back(se);
        if (!out.empty()) cycles.push_back(std::move(out));
    }
}

bool IncrementalFlattening::dijkstraBetween(int source,
                                            int target,
                                            std::vector<int>& pathEdges,
                                            std::vector<char>& pathForward) const
{
    pathEdges.clear();
    pathForward.clear();
    if (source < 0 || source >= nV || target < 0 || target >= nV) return false;

    std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
    std::vector<int> parent(nV, -1);
    std::vector<int> parentEdge(nV, -1);
    using Node = std::pair<double, int>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

    dist[source] = 0.0;
    pq.push({0.0, source});

    while (!pq.empty()) {
        const Node top = pq.top();
        pq.pop();
        const double d = top.first;
        const int u = top.second;
        if (d > dist[u]) continue;
        if (u == target) break;

        for (const auto& inc : vtxInc[u]) {
            const int ei = inc.first;
            if (ei < 0 || ei >= nE) continue;
            const EdgeInfo& e = edges[ei];
            const int v = (e.v0 == u) ? e.v1 : e.v0;
            if (v < 0 || v >= nV) continue;
            const double w = e.length;
            if (w <= 0.0) continue;
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                parent[v] = u;
                parentEdge[v] = ei;
                pq.push({dist[v], v});
            }
        }
    }

    if (parent[target] < 0 && source != target) return false;

    int cur = target;
    while (cur != source && parent[cur] >= 0) {
        const int ei = parentEdge[cur];
        if (ei < 0) break;
        const int from = parent[cur];
        pathEdges.push_back(ei);
        pathForward.push_back(static_cast<char>(edgeSign(from, cur) > 0 ? 1 : 0));
        cur = from;
    }
    if (cur != source) return false;
    std::reverse(pathEdges.begin(), pathEdges.end());
    std::reverse(pathForward.begin(), pathForward.end());
    return !pathEdges.empty();
}

bool IncrementalFlattening::dijkstraToBoundary(int source,
                                               std::vector<int>& pathEdges,
                                               std::vector<char>& pathForward) const
{
    int target = -1;
    double best = std::numeric_limits<double>::infinity();
    std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
    using Node = std::pair<double, int>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
    dist[source] = 0.0;
    pq.push({0.0, source});
    std::vector<int> parent(nV, -1);
    while (!pq.empty()) {
        const Node top = pq.top();
        pq.pop();
        const double d = top.first;
        const int u = top.second;
        if (d > dist[u]) continue;
        for (const auto& inc : vtxInc[u]) {
            const int ei = inc.first;
            if (ei < 0 || ei >= nE) continue;
            const EdgeInfo& e = edges[ei];
            const int v = (e.v0 == u) ? e.v1 : e.v0;
            if (v < 0 || v >= nV) continue;
            if (dist[u] + e.length < dist[v]) {
                dist[v] = dist[u] + e.length;
                parent[v] = u;
                pq.push({dist[v], v});
            }
        }
    }
    for (int vi = 0; vi < nV; ++vi) {
        if (!isMeshBoundaryVertex(vi)) continue;
        if (dist[vi] < best) {
            best = dist[vi];
            target = vi;
        }
    }
    if (target < 0) return false;
    return dijkstraBetween(source, target, pathEdges, pathForward);
}

void IncrementalFlattening::buildTreeCotreeSeamPaths()
{
    if (genus <= 0 || nB != 0) return;

    std::vector<std::vector<std::pair<int, int>>> cycles;
    collectHomologyCycles(cycles);
    for (const auto& cyc : cycles) {
        if (cyc.empty()) continue;
        SeamPath path;
        for (const auto& se : cyc) {
            path.edges.push_back(se.first);
            path.forward.push_back(static_cast<char>(se.second > 0 ? 1 : 0));
            if (se.first >= 0 && se.first < nE) {
                m_isSeamEdge[se.first] = 1;
            }
        }
        m_cutPaths.push_back(std::move(path));
    }
}

void IncrementalFlattening::buildCutSeams()
{
    m_cutPaths.clear();
    m_isSeamEdge.assign(nE, 0);

    std::vector<int> pathE;
    std::vector<char> pathF;

    if (nB > 0) {
        for (int cone : m_coneVerts) {
            if (!dijkstraToBoundary(cone, pathE, pathF)) continue;
            SeamPath sp;
            sp.edges = pathE;
            sp.forward = pathF;
            m_cutPaths.push_back(sp);
            for (int ei : pathE) {
                if (ei >= 0 && ei < nE) m_isSeamEdge[ei] = 1;
            }
        }
    } else {
        buildTreeCotreeSeamPaths();
        if (m_coneVerts.empty()) return;
        int farthest = -1;
        double best = -1.0;
        std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
        std::vector<char> vis(nV, 0);
        std::queue<int> q;
        const int src = m_coneVerts[0];
        dist[src] = 0.0;
        q.push(src);
        vis[src] = 1;
        while (!q.empty()) {
            const int u = q.front();
            q.pop();
            for (const auto& inc : vtxInc[u]) {
                const int ei = inc.first;
                if (ei < 0 || ei >= nE) continue;
                const EdgeInfo& e = edges[ei];
                const int v = (e.v0 == u) ? e.v1 : e.v0;
                if (v < 0 || v >= nV || vis[v]) continue;
                vis[v] = 1;
                dist[v] = dist[u] + e.length;
                q.push(v);
            }
        }
        for (int vi = 0; vi < nV; ++vi) {
            if (dist[vi] > best) {
                best = dist[vi];
                farthest = vi;
            }
        }
        if (farthest >= 0 && farthest != src) {
            if (dijkstraBetween(src, farthest, pathE, pathF)) {
                SeamPath sp;
                sp.edges = pathE;
                sp.forward = pathF;
                m_cutPaths.push_back(sp);
                for (int ei : pathE) {
                    if (ei >= 0 && ei < nE) m_isSeamEdge[ei] = 1;
                }
            }
        }
    }

    if (std::none_of(m_isSeamEdge.begin(), m_isSeamEdge.end(), [](char b) { return b != 0; })) {
        for (int ei = 0; ei < nE; ++ei) {
            if (!edges[ei].isInterior()) continue;
            m_isSeamEdge[ei] = 1;
            SeamPath sp;
            sp.edges.push_back(ei);
            sp.forward.push_back(1);
            m_cutPaths.push_back(sp);
            break;
        }
    }
}

bool IncrementalFlattening::buildCutMesh()
{
    CutMesh& cm = m_cutMesh;
    cm.cutPos.resize(0, 3);
    cm.cutFaces.resize(0);
    cm.origVert.clear();
    cm.seamMate.clear();
    cm.faceOrig.clear();
    cm.sheetOfFace.clear();
    cm.cutUV.clear();
    cm.isCutBoundary.clear();
    cm.nCutVerts = 0;
    cm.nCutFaces = nF;

    if (nF < 1) return false;

    std::vector<char> vtxOnSeam(nV, 0);
    for (int ei = 0; ei < nE; ++ei) {
        if (!m_isSeamEdge[ei]) continue;
        vtxOnSeam[edges[ei].v0] = 1;
        vtxOnSeam[edges[ei].v1] = 1;
    }

    cm.sheetOfFace.assign(nF, -1);
    std::queue<std::pair<int, int>> fq;
    cm.sheetOfFace[0] = 0;
    fq.push({0, 0});

    std::vector<std::vector<std::pair<int, int>>> faceAdj(nF);
    for (int ei = 0; ei < nE; ++ei) {
        const EdgeInfo& e = edges[ei];
        if (!e.isInterior()) continue;
        faceAdj[e.f0].push_back({e.f1, ei});
        faceAdj[e.f1].push_back({e.f0, ei});
    }

    while (!fq.empty()) {
        const int f = fq.front().first;
        const int s = fq.front().second;
        fq.pop();
        for (const auto& nb : faceAdj[f]) {
            const int g = nb.first;
            const int ei = nb.second;
            if (g < 0 || g >= nF) continue;
            if (cm.sheetOfFace[g] >= 0) continue;
            const int ns = m_isSeamEdge[ei] ? (1 - s) : s;
            cm.sheetOfFace[g] = ns;
            fq.push({g, ns});
        }
    }
    for (int f = 0; f < nF; ++f) {
        if (cm.sheetOfFace[f] < 0) cm.sheetOfFace[f] = 0;
    }

    std::map<std::pair<int, int>, int> vtxSheetToCut;
    auto getCutVert = [&](int origV, int sheet) -> int {
        if (origV < 0 || origV >= nV) return -1;
        const bool seamInterior = vtxOnSeam[origV] && !isMeshBoundaryVertex(origV);
        const int s = seamInterior ? sheet : 0;
        const std::pair<int, int> key{origV, s};
        auto it = vtxSheetToCut.find(key);
        if (it != vtxSheetToCut.end()) return it->second;

        const int idx = static_cast<int>(cm.origVert.size());
        cm.origVert.push_back(origV);
        cm.seamMate.push_back(-1);
        cm.cutUV.push_back(mesh.vertices[origV].uv);
        cm.isCutBoundary.push_back(isMeshBoundaryVertex(origV));
        cm.cutPos.conservativeResize(idx + 1, 3);
        cm.cutPos.row(idx) = V.row(origV);
        vtxSheetToCut[key] = idx;

        if (seamInterior) {
            const int otherSheet = 1 - s;
            const std::pair<int, int> key2{origV, otherSheet};
            if (vtxSheetToCut.count(key2) == 0) {
                const int idx2 = static_cast<int>(cm.origVert.size());
                cm.origVert.push_back(origV);
                cm.seamMate.push_back(idx);
                cm.cutUV.push_back(mesh.vertices[origV].uv);
                cm.isCutBoundary.push_back(false);
                cm.cutPos.conservativeResize(idx2 + 1, 3);
                cm.cutPos.row(idx2) = V.row(origV);
                vtxSheetToCut[key2] = idx2;
                cm.seamMate[idx] = idx2;
            }
        }
        return idx;
    };

    cm.cutFaces.resize(nF * 3);
    cm.faceOrig.resize(nF);
    for (int f = 0; f < nF; ++f) {
        cm.faceOrig[f] = f;
        const int sheet = cm.sheetOfFace[f];
        for (int k = 0; k < 3; ++k) {
            const int ov = F[f * 3 + k];
            cm.cutFaces[f * 3 + k] = getCutVert(ov, sheet);
        }
    }

    cm.nCutVerts = static_cast<int>(cm.origVert.size());
    cm.nCutFaces = nF;
    return cm.nCutVerts >= 3;
}

void IncrementalFlattening::computeSeamRotations(const Eigen::VectorXd& omega)
{
    m_seamRotation.assign(nE, Matrix2d::Identity());
    for (int ei = 0; ei < nE; ++ei) {
        if (!m_isSeamEdge[ei] || !edges[ei].isInterior()) continue;
        const int k = static_cast<int>(std::round(2.0 * omega(ei) / M_PI));
        m_seamRotation[ei] = quarterTurnMatrix(k);
    }
}

bool IncrementalFlattening::solveEdgeOmegas(
    const Eigen::VectorXd& vertexHolonomy,
    const std::vector<std::pair<std::vector<std::pair<int, int>>, double>>& extraCycles,
    Eigen::VectorXd& omega,
    bool seamOnly) const
{
    std::vector<int> interiorEdges;
    for (int ei = 0; ei < nE; ++ei) {
        if (!edges[ei].isInterior()) continue;
        if (seamOnly && !m_isSeamEdge[ei]) continue;
        interiorEdges.push_back(ei);
    }
    const int nOmega = static_cast<int>(interiorEdges.size());
    if (nOmega == 0) return false;

    std::map<int, int> edgeToOmega;
    for (int i = 0; i < nOmega; ++i) {
        edgeToOmega[interiorEdges[i]] = i;
    }

    std::vector<Triplet<double>> trips;
    std::vector<double> rhs;
    int row = 0;

    for (int vi = 0; vi < nV; ++vi) {
        if (std::abs(vertexHolonomy(vi)) < 1e-14) continue;
        for (const auto& inc : vtxInc[vi]) {
            const auto it = edgeToOmega.find(inc.first);
            if (it == edgeToOmega.end()) continue;
            trips.emplace_back(row, it->second, static_cast<double>(inc.second));
        }
        rhs.push_back(vertexHolonomy(vi));
        ++row;
    }

    for (const auto& extra : extraCycles) {
        for (const auto& se : extra.first) {
            const auto it = edgeToOmega.find(se.first);
            if (it == edgeToOmega.end()) continue;
            trips.emplace_back(row, it->second, static_cast<double>(se.second));
        }
        rhs.push_back(extra.second);
        ++row;
    }

    if (row == 0) return false;

    SparseMatrix<double> A(row, nOmega);
    A.setFromTriplets(trips.begin(), trips.end());
    VectorXd b = Map<const VectorXd>(rhs.data(), static_cast<int>(rhs.size()));

    const SparseMatrix<double> AtA = A.transpose() * A + kReg * SparseMatrix<double>(Identity(nOmega, nOmega));
    const VectorXd Atb = A.transpose() * b;

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(AtA);
    if (solver.info() != Success) return false;

    omega = VectorXd::Zero(nE);
    const VectorXd sol = solver.solve(Atb);
    if (!isFiniteVec(sol)) return false;

    for (int i = 0; i < nOmega; ++i) {
        omega(interiorEdges[i]) = sol(i);
    }
    return true;
}

double IncrementalFlattening::transportAngle(int fFrom, int fTo, int edgeIdx) const
{
    if (fFrom < 0 || fTo < 0 || edgeIdx < 0) return 0.0;
    const EdgeInfo& e = edges[edgeIdx];

    auto edgeDirInFace = [&](int f, int& vTail, int& vHead) {
        for (int k = 0; k < 3; ++k) {
            const int a = F[f * 3 + k];
            const int b = F[f * 3 + (k + 1) % 3];
            if ((a == e.v0 && b == e.v1) || (a == e.v1 && b == e.v0)) {
                vTail = a;
                vHead = b;
                return;
            }
        }
        vTail = vHead = -1;
    };

    int t0 = -1, h0 = -1, t1 = -1, h1 = -1;
    edgeDirInFace(fFrom, t0, h0);
    edgeDirInFace(fTo, t1, h1);
    if (t0 < 0 || t1 < 0) return 0.0;

    auto tangentAt = [&](int f, int from, int to) {
        const Vector3d p = V.row(from);
        const Vector3d q = V.row(to);
        Vector3d n = (V.row(F[f * 3 + 1]) - V.row(F[f * 3]))
                         .cross(V.row(F[f * 3 + 2]) - V.row(F[f * 3]));
        n.normalize();
        Vector3d t = q - p;
        t -= t.dot(n) * n;
        t.normalize();
        return t;
    };

    const Vector3d u0 = tangentAt(fFrom, t0, h0);
    const Vector3d u1 = tangentAt(fTo, t1, h1);
    Vector3d n0 = (V.row(F[fFrom * 3 + 1]) - V.row(F[fFrom * 3]))
                      .cross(V.row(F[fFrom * 3 + 2]) - V.row(F[fFrom * 3]));
    n0.normalize();
    const double x = u0.dot(u1);
    const double y = n0.dot(u0.cross(u1));
    return std::atan2(y, x);
}

void IncrementalFlattening::solveCrossField(const Eigen::VectorXd& omega,
                                            Eigen::VectorXd& theta) const
{
    std::vector<std::vector<std::pair<int, int>>> faceAdj(nF);
    for (int ei = 0; ei < nE; ++ei) {
        const EdgeInfo& e = edges[ei];
        if (!e.isInterior()) continue;
        faceAdj[e.f0].push_back({e.f1, ei});
        faceAdj[e.f1].push_back({e.f0, ei});
    }

    std::vector<Triplet<double>> trips;
    std::vector<double> rhs;
    int row = 0;

    for (int f = 0; f < nF; ++f) {
        for (const auto& nb : faceAdj[f]) {
            const int g = nb.first;
            const int ei = nb.second;
            if (f > g) continue;

            const double kappa = transportAngle(f, g, ei);
            const double w = 1.0;
            trips.emplace_back(row, f, w);
            trips.emplace_back(row, g, -w);
            rhs.push_back(-kappa - omega(ei));
            ++row;
        }
    }

    if (row == 0) {
        theta = VectorXd::Zero(nF);
        return;
    }

    trips.emplace_back(row, 0, 1.0);
    rhs.push_back(0.0);
    ++row;

    SparseMatrix<double> A(row, nF);
    A.setFromTriplets(trips.begin(), trips.end());
    VectorXd b = Map<const VectorXd>(rhs.data(), static_cast<int>(rhs.size()));

    const SparseMatrix<double> AtA = A.transpose() * A + kReg * SparseMatrix<double>(Identity(nF, nF));
    const VectorXd Atb = A.transpose() * b;

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(AtA);
    theta = solver.solve(Atb);
}

void IncrementalFlattening::initCutMeshUV()
{
    Tutte tutte(mesh, TutteBoundary::CIRCLE);
    tutte.parameterize();

    CutMesh& cm = m_cutMesh;
    for (int cv = 0; cv < cm.nCutVerts; ++cv) {
        const int ov = cm.origVert[cv];
        if (ov >= 0 && ov < nV) {
            cm.cutUV[cv] = mesh.vertices[ov].uv;
        }
    }
}

void IncrementalFlattening::cutArapLocalStep()
{
    const CutMesh& cm = m_cutMesh;
    m_cutRotations.resize(cm.nCutFaces);

    for (int cf = 0; cf < cm.nCutFaces; ++cf) {
        const int f = cm.faceOrig[cf];
        const int ca = cm.cutFaces[cf * 3];
        const int cb = cm.cutFaces[cf * 3 + 1];
        const int cc = cm.cutFaces[cf * 3 + 2];

        const Vector3d p0 = cm.cutPos.row(ca);
        const Vector3d p1 = cm.cutPos.row(cb);
        const Vector3d p2 = cm.cutPos.row(cc);
        Vector3d e1 = p1 - p0;
        e1.normalize();
        Vector3d n = e1.cross(p2 - p0);
        n.normalize();
        Vector3d e2 = n.cross(e1);
        e2.normalize();

        const Vector2d x0(0, 0);
        const Vector2d x1(e1.dot(p1 - p0), e2.dot(p1 - p0));
        const Vector2d x2(e1.dot(p2 - p0), e2.dot(p2 - p0));

        const Vector2d u0 = cm.cutUV[ca];
        const Vector2d u1 = cm.cutUV[cb];
        const Vector2d u2 = cm.cutUV[cc];

        Matrix2d S = Matrix2d::Zero();
        const int edgeIdx[3][2] = {{1, 2}, {2, 0}, {0, 1}};
        const int opp[3] = {0, 1, 2};
        const Vector2d x[3] = {x0, x1, x2};
        const Vector2d u[3] = {u0, u1, u2};

        for (int e = 0; e < 3; ++e) {
            const int j = edgeIdx[e][0];
            const int k = edgeIdx[e][1];
            const int i = opp[e];
            const Vector2d dx = x[j] - x[k];
            const Vector2d du = u[j] - u[k];
            const Vector2d ea = x[j] - x[i];
            const Vector2d eb = x[k] - x[i];
            const double cotVal = ea.dot(eb) / std::max(std::abs(ea.x() * eb.y() - ea.y() * eb.x()), kAngleEps);
            const double w = std::max(cotVal, kReg);
            S += w * du * dx.transpose();
        }

        const double a00 = S(0, 0), a01 = S(0, 1), a10 = S(1, 0), a11 = S(1, 1);
        const double denom = std::sqrt((a00 + a11) * (a00 + a11) + (a10 - a01) * (a10 - a01));
        if (denom < kAngleEps) {
            m_cutRotations[cf] = Matrix2d::Identity();
        } else {
            const double c0 = (a00 + a11) / denom;
            const double s0 = (a10 - a01) / denom;
            m_cutRotations[cf] << c0, -s0, s0, c0;
        }
        (void)f;
    }
}

void IncrementalFlattening::enforceSeamRotations()
{
    const CutMesh& cm = m_cutMesh;
    for (int ei = 0; ei < nE; ++ei) {
        if (!m_isSeamEdge[ei] || !edges[ei].isInterior()) continue;
        const EdgeInfo& e = edges[ei];
        const int f0 = e.f0;
        const int f1 = e.f1;
        if (f0 < 0 || f1 < 0) continue;

        const int s0 = cm.sheetOfFace[f0];
        const int s1 = cm.sheetOfFace[f1];
        const int cf0 = f0;
        const int cf1 = f1;
        const Matrix2d& r = m_seamRotation[ei];

        if (s0 == 0 && s1 == 1) {
            m_cutRotations[cf1] = r * m_cutRotations[cf0];
        } else if (s1 == 0 && s0 == 1) {
            m_cutRotations[cf0] = r * m_cutRotations[cf1];
        }
    }
}

void IncrementalFlattening::cutArapGlobalStep()
{
    CutMesh& cm = m_cutMesh;
    const int nCV = cm.nCutVerts;
    if (nCV < 3) return;

    std::vector<int> interiorIdx(nCV, -1);
    int nInterior = 0;
    for (int i = 0; i < nCV; ++i) {
        if (!cm.isCutBoundary[i]) interiorIdx[i] = nInterior++;
    }
    if (nInterior == 0) return;

    std::vector<Triplet<double>> triplets;
    VectorXd bx = VectorXd::Zero(nInterior);
    VectorXd by = VectorXd::Zero(nInterior);

    for (int cf = 0; cf < cm.nCutFaces; ++cf) {
        const Matrix2d& R = m_cutRotations[cf];
        const int vIdx[3] = {
            cm.cutFaces[cf * 3],
            cm.cutFaces[cf * 3 + 1],
            cm.cutFaces[cf * 3 + 2]
        };

        const Vector3d p0 = cm.cutPos.row(vIdx[0]);
        const Vector3d p1 = cm.cutPos.row(vIdx[1]);
        const Vector3d p2 = cm.cutPos.row(vIdx[2]);
        Vector3d e1 = p1 - p0;
        e1.normalize();
        Vector3d n = e1.cross(p2 - p0);
        n.normalize();
        Vector3d e2 = n.cross(e1);
        e2.normalize();

        const Vector2d x[3] = {
            Vector2d(0, 0),
            Vector2d(e1.dot(p1 - p0), e2.dot(p1 - p0)),
            Vector2d(e1.dot(p2 - p0), e2.dot(p2 - p0))
        };

        const int edgeIdx[3][2] = {{1, 2}, {2, 0}, {0, 1}};
        const int opp[3] = {0, 1, 2};

        for (int e = 0; e < 3; ++e) {
            const int a = edgeIdx[e][0];
            const int b = edgeIdx[e][1];
            const int o = opp[e];
            const Vector2d ea = x[a] - x[o];
            const Vector2d eb = x[b] - x[o];
            const double cotVal = ea.dot(eb) / std::max(std::abs(ea.x() * eb.y() - ea.y() * eb.x()), kAngleEps);
            const double w = std::max(cotVal, kReg);
            const Vector2d rotDx = R * (x[a] - x[b]);

            const int vi = vIdx[a];
            const int vj = vIdx[b];

            if (!cm.isCutBoundary[vi]) {
                const int ri = interiorIdx[vi];
                triplets.emplace_back(ri, ri, w);
                bx(ri) += w * rotDx.x();
                by(ri) += w * rotDx.y();
                if (cm.isCutBoundary[vj]) {
                    bx(ri) += w * cm.cutUV[vj].x();
                    by(ri) += w * cm.cutUV[vj].y();
                } else {
                    triplets.emplace_back(ri, interiorIdx[vj], -w);
                }
            }
            if (!cm.isCutBoundary[vj]) {
                const int rj = interiorIdx[vj];
                triplets.emplace_back(rj, rj, w);
                bx(rj) -= w * rotDx.x();
                by(rj) -= w * rotDx.y();
                if (cm.isCutBoundary[vi]) {
                    bx(rj) += w * cm.cutUV[vi].x();
                    by(rj) += w * cm.cutUV[vi].y();
                } else {
                    triplets.emplace_back(rj, interiorIdx[vi], -w);
                }
            }
        }
    }

    SparseMatrix<double> L(nInterior, nInterior);
    L.setFromTriplets(triplets.begin(), triplets.end());

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Success) return;

    const VectorXd ux = solver.solve(bx);
    const VectorXd uy = solver.solve(by);
    if (!isFiniteVec(ux) || !isFiniteVec(uy)) return;

    for (int i = 0; i < nCV; ++i) {
        if (!cm.isCutBoundary[i]) {
            cm.cutUV[i] = Vector2d(ux(interiorIdx[i]), uy(interiorIdx[i]));
        }
    }
}

void IncrementalFlattening::projectCutUVToMesh()
{
    const CutMesh& cm = m_cutMesh;
    std::vector<char> assigned(nV, 0);
    for (int f = 0; f < cm.nCutFaces; ++f) {
        if (cm.sheetOfFace[f] != 0) continue;
        for (int k = 0; k < 3; ++k) {
            const int cv = cm.cutFaces[f * 3 + k];
            const int ov = cm.origVert[cv];
            if (ov < 0 || ov >= nV) continue;
            mesh.vertices[ov].uv = cm.cutUV[cv];
            assigned[ov] = 1;
        }
    }
    for (int cv = 0; cv < cm.nCutVerts; ++cv) {
        const int ov = cm.origVert[cv];
        if (ov < 0 || ov >= nV || assigned[ov]) continue;
        mesh.vertices[ov].uv = cm.cutUV[cv];
        assigned[ov] = 1;
    }
}

void IncrementalFlattening::globalArapPhase()
{
    buildCutSeams();
    if (!buildCutMesh()) {
        std::cerr << "[IncrementalFlattening] cut mesh build failed, fallback Tutte\n";
        Tutte tutte(mesh, TutteBoundary::CIRCLE);
        tutte.parameterize();
        return;
    }

    VectorXd vertexHolonomy = VectorXd::Zero(nV);
    for (size_t i = 0; i < m_coneVerts.size(); ++i) {
        vertexHolonomy(m_coneVerts[i]) = m_coneK(static_cast<int>(i));
    }

    std::vector<std::vector<std::pair<int, int>>> homology;
    collectHomologyCycles(homology);
    std::vector<std::pair<std::vector<std::pair<int, int>>, double>> extra;
    VectorXd omegaProbe;
    if (solveEdgeOmegas(vertexHolonomy, {}, omegaProbe, true)) {
        for (const auto& cyc : homology) {
            double h = 0.0;
            for (const auto& se : cyc) {
                if (se.first >= 0 && se.first < omegaProbe.size()) {
                    h += static_cast<double>(se.second) * omegaProbe(se.first);
                }
            }
            const double hRound = roundToHalfPi(h);
            if (std::abs(h - hRound) > 1e-3) {
                extra.push_back({cyc, hRound});
            }
        }
    }

    m_omega = VectorXd::Zero(nE);
    if (!solveEdgeOmegas(vertexHolonomy, extra, m_omega, true)) {
        m_omega.setZero();
    }

    computeSeamRotations(m_omega);

    VectorXd theta;
    solveCrossField(m_omega, theta);

    m_cutRotations.resize(m_cutMesh.nCutFaces);
    for (int f = 0; f < nF; ++f) {
        const double c = std::cos(theta(f));
        const double s = std::sin(theta(f));
        m_cutRotations[f] << c, -s, s, c;
    }

    initCutMeshUV();

    for (int iter = 0; iter < m_maxArapIters; ++iter) {
        cutArapLocalStep();
        enforceSeamRotations();
        cutArapGlobalStep();
    }

    projectCutUVToMesh();
}

void IncrementalFlattening::parameterize()
{
    if (!extractMeshData()) {
        std::cerr << "[IncrementalFlattening] mesh extraction failed\n";
        return;
    }

    flatteningPhase();
    roundingPhase();
    globalArapPhase();
    normalize();

    std::cout << "[IncrementalFlattening] cones: ";
    for (int vi : m_coneVerts) std::cout << vi << " ";
    std::cout << " | seam edges: " << seamEdgeCount()
              << " | cut verts: " << cutVertexCount() << "\n";
}
