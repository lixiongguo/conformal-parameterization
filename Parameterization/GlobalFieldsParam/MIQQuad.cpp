#include "MIQQuad.h"
#include <cmath>
#include <map>
#include <queue>
#include <limits>

using namespace Eigen;

namespace {
const double kQuarterTurn = std::acos(-1.0) / 2.0;
const double kPi = std::acos(-1.0);

int mod4(int x) { int r = x % 4; return r < 0 ? r + 4 : r; }
} // namespace

MIQQuad::MIQQuad(Mesh& mesh0)
    : GlobalFieldsParameterization(mesh0)
{
}

void MIQQuad::buildLocalFrame(const Vector3d& n, Vector3d& t1, Vector3d& t2) const
{
    Vector3d ref(1, 0, 0);
    if (std::fabs(n.dot(ref)) > 0.9) ref = Vector3d(0, 1, 0);
    t1 = (ref - ref.dot(n) * n).normalized();
    t2 = n.cross(t1).normalized();
}

double MIQQuad::cotan(const Vector3d& a, const Vector3d& b, const Vector3d& c) const
{
    const Vector3d u = a - b;
    const Vector3d v = c - b;
    const double d = u.dot(v);
    const double cr = u.cross(v).norm();
    if (cr < 1e-12) return 0.0;
    return d / cr;
}

int MIQQuad::findSharedEdge(int f1, int f2) const
{
    for (int ei = 0; ei < nE; ++ei) {
        const Edge& e = edgeList[ei];
        if ((e.f1 == f1 && e.f2 == f2) || (e.f1 == f2 && e.f2 == f1))
            return ei;
    }
    return -1;
}

// ============================== Phase 1 =====================================

bool MIQQuad::initMeshData()
{
    nV = static_cast<int>(mesh.vertices.size());
    nF = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (!f->isBoundary()) ++nF;
    }
    if (nV < 3 || nF < 1) return false;

    vertPos.resize(nV, 3);
    faces.resize(nF * 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        vertPos.row(v->index) = v->position;
    }

    int fi = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        faces[fi * 3] = f->he->vertex->index;
        faces[fi * 3 + 1] = f->he->next->vertex->index;
        faces[fi * 3 + 2] = f->he->next->next->vertex->index;
        ++fi;
    }

    edgeList.clear();
    std::map<std::pair<int, int>, int> em;
    for (int f = 0; f < nF; ++f) {
        for (int k = 0; k < 3; ++k) {
            int v1 = faces[f * 3 + k];
            int v2 = faces[f * 3 + (k + 1) % 3];
            if (v1 > v2) std::swap(v1, v2);
            const auto key = std::make_pair(v1, v2);
            auto it = em.find(key);
            if (it == em.end()) {
                Edge e;
                e.v1 = v1;
                e.v2 = v2;
                e.f1 = f;
                e.f2 = -1;
                e.idx = static_cast<int>(edgeList.size());
                em[key] = e.idx;
                edgeList.push_back(e);
            } else {
                edgeList[it->second].f2 = f;
            }
        }
    }
    nE = static_cast<int>(edgeList.size());
    jump.resize(nE);
    jump.setZero();

    mixedIntegerProgram_.reset();
    return true;
}

void MIQQuad::initCrossField()
{
    faceN.resize(nF, 3);
    faceT1.resize(nF, 3);
    faceT2.resize(nF, 3);
    theta.resize(nF);
    theta.setZero();

    for (int fi = 0; fi < nF; ++fi) {
        const int v0 = faces[fi * 3];
        const int v1 = faces[fi * 3 + 1];
        const int v2 = faces[fi * 3 + 2];
        const Vector3d p0 = vertPos.row(v0);
        const Vector3d p1 = vertPos.row(v1);
        const Vector3d p2 = vertPos.row(v2);
        const Vector3d fn = (p1 - p0).cross(p2 - p0).normalized();
        faceN.row(fi) = fn;

        Vector3d t1, t2;
        buildLocalFrame(fn, t1, t2);
        faceT1.row(fi) = t1;
        faceT2.row(fi) = t2;

        const Vector3d e0 = p1 - p0;
        const Vector3d e1 = p2 - p1;
        const Vector3d e2 = p0 - p2;
        const double l0 = e0.squaredNorm();
        const double l1 = e1.squaredNorm();
        const double l2 = e2.squaredNorm();
        Vector3d d;
        if (l0 > l1 && l0 > l2) d = e0;
        else if (l1 > l0 && l1 > l2) d = e1;
        else d = e2;
        const double dx = d.dot(t1);
        const double dy = d.dot(t2);
        theta(fi) = std::atan2(dy, dx);
    }
}

std::vector<MixedIntegerProgram::Constraint> MIQQuad::buildFaceAdjacencyConstraints() const
{
    std::vector<MixedIntegerProgram::Constraint> constraints;
    constraints.reserve(edgeList.size());
    for (const auto& e : edgeList) {
        if (e.f1 < 0 || e.f2 < 0) continue;

        const Vector3d edgeDir = (vertPos.row(e.v2) - vertPos.row(e.v1)).normalized();

        const Vector3d t1i = faceT1.row(e.f1);
        const Vector3d t2i = faceT2.row(e.f1);
        const double alpha_i = std::atan2(edgeDir.dot(t2i), edgeDir.dot(t1i));

        const Vector3d t1j = faceT1.row(e.f2);
        const Vector3d t2j = faceT2.row(e.f2);
        const double alpha_j = std::atan2(edgeDir.dot(t2j), edgeDir.dot(t1j));

        const double kappa = alpha_j - alpha_i;

        MixedIntegerProgram::Constraint c;
        c.i = e.f1;
        c.j = e.f2;
        c.idx = e.idx;
        c.kappa = kappa;
        constraints.push_back(c);
    }
    return constraints;
}

double MIQQuad::crossFieldEnergy() const
{
    if (!mixedIntegerProgram_) return 0.0;
    return mixedIntegerProgram_->energy(theta, jump);
}

bool MIQQuad::solveCrossFieldIP()
{
    auto constraints = buildFaceAdjacencyConstraints();
    mixedIntegerProgram_ = std::make_unique<MixedIntegerProgram>(
        nF,
        std::move(constraints),
        kQuarterTurn);
    mixedIntegerProgram_->setIntegerBounds(jumpLo_, jumpHi_);
    mixedIntegerProgram_->setAlternatingIterations(crossIters_);
    mixedIntegerProgram_->setRefinePasses(jumpRefinePasses_);
    if (!mixedIntegerProgram_->solve(theta, jump)) {
        return false;
    }

    wrapCrossFieldAngles();
    return true;
}

void MIQQuad::wrapCrossFieldAngles()
{
    for (int fi = 0; fi < theta.size(); ++fi) {
        theta(fi) = std::fmod(theta(fi), kQuarterTurn);
        if (theta(fi) < 0.0) {
            theta(fi) += kQuarterTurn;
        }
    }
}

void MIQQuad::buildTargetDirs()
{
    faceD1.resize(nF, 3);
    faceD2.resize(nF, 3);

    for (int fi = 0; fi < nF; ++fi) {
        const Vector3d n = faceN.row(fi);
        Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);
        const double c = std::cos(theta(fi));
        const double s = std::sin(theta(fi));
        faceD1.row(fi) = (c * t1 + s * t2).normalized();
        faceD2.row(fi) = (-s * t1 + c * t2).normalized();
    }
}

void MIQQuad::buildCotLaplacian(SparseMatrix<double>& L)
{
    L.resize(nV, nV);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nV);

    for (const auto& e : edgeList) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        const int vi = e.v1;
        const int vj = e.v2;

        int o1 = -1, o2 = -1;
        for (int k = 0; k < 3; ++k) {
            const int v = faces[e.f1 * 3 + k];
            if (v != vi && v != vj) { o1 = v; break; }
        }
        for (int k = 0; k < 3; ++k) {
            const int v = faces[e.f2 * 3 + k];
            if (v != vi && v != vj) { o2 = v; break; }
        }

        double w = 0.5 * (cotan(vertPos.row(o1), vertPos.row(vi), vertPos.row(vj)) +
                          cotan(vertPos.row(o2), vertPos.row(vi), vertPos.row(vj)));
        if (w < 0) w = 0.0;

        diag(vi) += w;
        diag(vj) += w;
        trips.emplace_back(vi, vj, -w);
        trips.emplace_back(vj, vi, -w);
    }
    for (int i = 0; i < nV; ++i) {
        trips.emplace_back(i, i, diag(i) + 1e-8);
    }
    L.setFromTriplets(trips.begin(), trips.end());
}

void MIQQuad::solvePoisson()
{
    VectorXd div1 = VectorXd::Zero(nV);
    VectorXd div2 = VectorXd::Zero(nV);

    for (int fi = 0; fi < nF; ++fi) {
        const int v0 = faces[fi * 3];
        const int v1 = faces[fi * 3 + 1];
        const int v2 = faces[fi * 3 + 2];
        const Vector3d p0 = vertPos.row(v0);
        const Vector3d p1 = vertPos.row(v1);
        const Vector3d p2 = vertPos.row(v2);

        const Vector3d n = (p1 - p0).cross(p2 - p0);
        const Vector3d e12 = p2 - p1;
        const Vector3d e20 = p0 - p2;
        const Vector3d e01 = p1 - p0;

        const Vector3d fd1 = faceD1.row(fi);
        const Vector3d fd2 = faceD2.row(fi);

        div1(v0) += 0.5 * fd1.dot(n.cross(e12));
        div1(v1) += 0.5 * fd1.dot(n.cross(e20));
        div1(v2) += 0.5 * fd1.dot(n.cross(e01));

        div2(v0) += 0.5 * fd2.dot(n.cross(e12));
        div2(v1) += 0.5 * fd2.dot(n.cross(e20));
        div2(v2) += 0.5 * fd2.dot(n.cross(e01));
    }

    SparseMatrix<double> L;
    buildCotLaplacian(L);

    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        if (v->isBoundary()) {
            L.coeffRef(v->index, v->index) += 1e6;
            div1(v->index) = 0.0;
            div2(v->index) = 0.0;
        }
    }
    bool hasBnd = false;
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        if (v->isBoundary()) { hasBnd = true; break; }
    }
    if (!hasBnd) {
        L.coeffRef(0, 0) += 1e6;
        div1(0) = 0.0;
        div2(0) = 0.0;
    }

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Success) { UV = MatrixXd::Zero(nV, 2); return; }
    const VectorXd u = solver.solve(div1);
    const VectorXd v = solver.solve(div2);

    const double uMin = u.minCoeff(), uMax = u.maxCoeff();
    const double vMin = v.minCoeff(), vMax = v.maxCoeff();
    double uR = uMax - uMin, vR = vMax - vMin;
    if (uR < 1e-10) uR = 1.0;
    if (vR < 1e-10) vR = 1.0;

    UV.resize(nV, 2);
    for (int i = 0; i < nV; ++i) {
        UV(i, 0) = (u(i) - uMin) / uR;
        UV(i, 1) = (v(i) - vMin) / vR;
    }
}

// ============================== Phase 2 =====================================

int MIQQuad::countSingularities()
{
    singularVerts_.clear();
    for (int v = 0; v < nV; ++v) {
        // collect incident faces
        std::vector<int> vf;
        for (int fi = 0; fi < nF; ++fi)
            for (int k = 0; k < 3; ++k)
                if (faces[fi * 3 + k] == v) { vf.push_back(fi); break; }
        if (vf.size() < 3) continue;

        // angle defect
        double sumAng = 0.0;
        for (int fi : vf) {
            int k0 = -1;
            for (int k = 0; k < 3; ++k)
                if (faces[fi * 3 + k] == v) { k0 = k; break; }
            int k1 = (k0 + 1) % 3, k2 = (k0 + 2) % 3;
            Vector3d e1 = vertPos.row(faces[fi * 3 + k1]) - vertPos.row(v);
            Vector3d e2 = vertPos.row(faces[fi * 3 + k2]) - vertPos.row(v);
            double denom = e1.norm() * e2.norm();
            if (denom < 1e-12) continue;
            double c = std::max(-1.0, std::min(1.0, e1.dot(e2) / denom));
            sumAng += std::acos(c);
        }
        double I0 = (2.0 * kPi - sumAng) / (2.0 * kPi);

        // accumulated jumps (signed) around v from Phase 1
        double jSum = 0.0;
        // order faces around v via angle in tangent plane
        Vector3d n = faceN.row(vf[0]);
        if (n.norm() < 1e-12) continue;
        Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);
        const Vector3d pv = vertPos.row(v);
        std::vector<std::pair<double, int>> ord;
        for (int fi : vf) {
            Vector3d c(0, 0, 0);
            for (int k = 0; k < 3; ++k) c += vertPos.row(faces[fi * 3 + k]);
            c = c / 3.0 - pv;
            c -= c.dot(n) * n;
            if (c.norm() < 1e-12) continue;
            ord.emplace_back(std::atan2(c.dot(t2), c.dot(t1)), fi);
        }
        if (ord.size() < 3) continue;
        std::sort(ord.begin(), ord.end());

        for (size_t i = 0; i < ord.size(); ++i) {
            int f0 = ord[i].second;
            int f1 = ord[(i + 1) % ord.size()].second;
            int ei = findSharedEdge(f0, f1);
            if (ei < 0) continue;
            const Edge& e = edgeList[ei];
            // signed jump: +jump if f0→f1 matches e.f1→e.f2, else −jump
            int sign = (e.f1 == f0 && e.f2 == f1) ? 1 : -1;
            if (ei < jump.size()) jSum += sign * jump(ei);
        }

        double I = I0 + 0.25 * jSum;
        if (std::fabs(I) > 1e-6)
            singularVerts_.push_back(v);
    }
    return static_cast<int>(singularVerts_.size());
}

void MIQQuad::buildCutGraph()
{
    cutEdgeMap_.assign(nE, 0);
    cutEdges_.clear();

    // face adjacency for dual spanning tree
    std::vector<std::vector<std::pair<int, int>>> faceAdj(nF);
    for (int ei = 0; ei < nE; ++ei) {
        const Edge& e = edgeList[ei];
        if (e.f2 < 0) continue;
        faceAdj[e.f1].emplace_back(e.f2, ei);
        faceAdj[e.f2].emplace_back(e.f1, ei);
    }

    // BFS spanning tree on dual faces
    std::vector<char> inTree(nE, 0);
    std::vector<char> visited(nF, 0);
    std::queue<int> q;
    q.push(0);
    visited[0] = 1;
    while (!q.empty()) {
        int f = q.front(); q.pop();
        for (const auto& nb : faceAdj[f]) {
            if (visited[nb.first]) continue;
            visited[nb.first] = 1;
            inTree[nb.second] = 1;         // this dual edge is in the tree
            q.push(nb.first);
        }
    }

    // Cotree edges → cut
    for (int ei = 0; ei < nE; ++ei) {
        if (edgeList[ei].f2 < 0) continue;                // boundary
        if (inTree[ei]) continue;                          // tree edge
        cutEdgeMap_[ei] = 1;
        cutEdges_.push_back(ei);
    }

    // --- remove leaf cut-edges (non-separating) ---
    // Build valence of each vertex within the cut-graph subgraph.
    std::vector<int> cutDeg(nV, 0);
    for (int ei : cutEdges_) {
        cutDeg[edgeList[ei].v1]++;
        cutDeg[edgeList[ei].v2]++;
    }

    int nSings = countSingularities();
    std::vector<char> isSing(nV, 0);
    for (int sv : singularVerts_) isSing[sv] = 1;

    bool removed = true;
    while (removed && !cutEdges_.empty()) {
        removed = false;
        for (size_t k = 0; k < cutEdges_.size(); ) {
            int ei = cutEdges_[k];
            const Edge& e = edgeList[ei];
            // can remove if at least one endpoint is a leaf (deg==1) and NOT a singularity
            bool leafA = (cutDeg[e.v1] == 1 && !isSing[e.v1]);
            bool leafB = (cutDeg[e.v2] == 1 && !isSing[e.v2]);
            if (leafA || leafB) {
                if (leafA) cutDeg[e.v1]--;
                if (leafB) cutDeg[e.v2]--;
                cutEdgeMap_[ei] = 0;
                cutEdges_.erase(cutEdges_.begin() + k);
                removed = true;
            } else {
                ++k;
            }
        }
    }

    // --- connect remaining singularities to cut graph (Dijkstra) ---
    if (nSings > 0 && !cutEdges_.empty()) {
        // mark vertices on current cut
        std::vector<char> onCut(nV, 0);
        for (int ei : cutEdges_) {
            onCut[edgeList[ei].v1] = 1;
            onCut[edgeList[ei].v2] = 1;
        }

        // vertex adjacency
        std::vector<std::vector<int>> vAdj(nV);
        for (int ei = 0; ei < nE; ++ei) {
            const Edge& e = edgeList[ei];
            vAdj[e.v1].push_back(e.v2);
            vAdj[e.v2].push_back(e.v1);
        }

        for (int sv : singularVerts_) {
            if (onCut[sv]) continue;

            // Dijkstra to nearest cut vertex
            std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
            std::vector<int> parent(nV, -1);
            std::vector<int> parentEdge(nV, -1);
            using P = std::pair<double, int>;
            std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
            dist[sv] = 0;
            pq.emplace(0, sv);

            int target = -1;
            while (!pq.empty()) {
                auto [d, v] = pq.top(); pq.pop();
                if (d != dist[v]) continue;
                if (onCut[v]) { target = v; break; }
                for (int nb : vAdj[v]) {
                    double nd = d + (vertPos.row(nb) - vertPos.row(v)).norm();
                    if (nd < dist[nb]) {
                        dist[nb] = nd;
                        parent[nb] = v;
                        // find edge between v and nb
                        for (int ei = 0; ei < nE; ++ei) {
                            const Edge& e = edgeList[ei];
                            if ((e.v1 == v && e.v2 == nb) || (e.v1 == nb && e.v2 == v)) {
                                parentEdge[nb] = ei; break;
                            }
                        }
                        pq.emplace(nd, nb);
                    }
                }
            }

            if (target < 0) continue;

            // walk back and add edges to cut
            int cur = target;
            while (cur != sv && parent[cur] >= 0) {
                int ei = parentEdge[cur];
                if (ei >= 0 && !cutEdgeMap_[ei]) {
                    cutEdgeMap_[ei] = 1;
                    cutEdges_.push_back(ei);
                    onCut[edgeList[ei].v1] = 1;
                    onCut[edgeList[ei].v2] = 1;
                }
                cur = parent[cur];
            }
        }
    }
}

bool MIQQuad::buildCutMesh()
{
    // Union-Find on (face, localCorner) to duplicate vertices along cut edges.
    const int nCorners = nF * 3;
    std::vector<int> parent(nCorners);
    for (int i = 0; i < nCorners; ++i) parent[i] = i;

    auto find = [&](int x) {
        int root = x;
        while (parent[root] != root) root = parent[root];
        while (parent[x] != root) { int nxt = parent[x]; parent[x] = root; x = nxt; }
        return root;
    };
    auto unite = [&](int a, int b) {
        a = find(a); b = find(b);
        if (a != b) parent[b] = a;
    };

    // Unite corners across TREE edges (non-cut interior edges).
    for (int ei = 0; ei < nE; ++ei) {
        if (edgeList[ei].f2 < 0) continue;  // boundary
        if (cutEdgeMap_[ei]) continue;       // cut edge → do NOT unite
        const Edge& e = edgeList[ei];
        // find local corner indices for v1, v2 in each face
        auto loc = [&](int fi, int vert) {
            for (int k = 0; k < 3; ++k)
                if (faces[fi * 3 + k] == vert) return k;
            return -1;
        };
        int c1_f1 = loc(e.f1, e.v1), c2_f1 = loc(e.f1, e.v2);
        int c1_f2 = loc(e.f2, e.v1), c2_f2 = loc(e.f2, e.v2);
        if (c1_f1 < 0 || c2_f1 < 0 || c1_f2 < 0 || c2_f2 < 0) return false;
        unite(e.f1 * 3 + c1_f1, e.f2 * 3 + c1_f2);
        unite(e.f1 * 3 + c2_f1, e.f2 * 3 + c2_f2);
    }

    // map DSU root → new vertex index
    std::map<int, int> root2idx;
    std::vector<Vector3d> cvpos;
    cutParent.clear();
    cutFaces.resize(nF * 3);
    cutFaces.setConstant(-1);

    for (int fi = 0; fi < nF; ++fi) {
        for (int k = 0; k < 3; ++k) {
            int r = find(fi * 3 + k);
            auto it = root2idx.find(r);
            int ci;
            if (it == root2idx.end()) {
                ci = static_cast<int>(cvpos.size());
                root2idx[r] = ci;
                int bv = faces[fi * 3 + k];
                cvpos.push_back(vertPos.row(bv));
                cutParent.push_back(bv);
            } else {
                ci = it->second;
            }
            cutFaces[fi * 3 + k] = ci;
        }
    }

    nCV = static_cast<int>(cvpos.size());
    cutVertPos.resize(nCV, 3);
    for (int i = 0; i < nCV; ++i)
        cutVertPos.row(i) = cvpos[i];

    cutRot_.assign(cutEdges_.size(), 0);
    cutTx.resize(static_cast<int>(cutEdges_.size())); cutTx.setZero();
    cutTy.resize(static_cast<int>(cutEdges_.size())); cutTy.setZero();

    return nCV >= 3;
}

void MIQQuad::determineCutRotations()
{
    // Propagate cross-field orientation from face 0 via tree edges.
    // The rotation on a cut edge is determined by comparing the per-face
    // orientation propagated to both incident faces.

    const int nCuts = static_cast<int>(cutEdges_.size());
    if (nCuts < 1) return;

    // BFS on face adjacency (tree edges only)
    std::vector<int> faceRot(nF, -1);   // orientation π/2 rotation from seed
    std::queue<int> q;
    faceRot[0] = 0;
    q.push(0);

    // build face adjacency via non-cut interior edges
    std::vector<std::vector<std::pair<int, int>>> treeAdj(nF);
    for (int ei = 0; ei < nE; ++ei) {
        const Edge& e = edgeList[ei];
        if (e.f2 < 0) continue;
        if (cutEdgeMap_[ei]) continue;
        treeAdj[e.f1].emplace_back(e.f2, ei);
        treeAdj[e.f2].emplace_back(e.f1, ei);
    }

    while (!q.empty()) {
        int f = q.front(); q.pop();
        for (const auto& nb : treeAdj[f]) {
            if (faceRot[nb.first] >= 0) continue;
            // compute matching across this tree edge
            const Edge& e = edgeList[nb.second];
            int matchingVal = 0;
            // compare faceDirs across the edge (same logic as computeMatching in QuadCover)
            Vector3d n1 = faceN.row(e.f1);
            Vector3d n2 = faceN.row(e.f2);
            Vector3d t11, t12, t21, t22;
            buildLocalFrame(n1, t11, t12);
            buildLocalFrame(n2, t21, t22);
            Vector2d d1(faceD1.row(e.f1).dot(t11), faceD1.row(e.f1).dot(t12));
            Vector2d d2(faceD1.row(e.f2).dot(t21), faceD1.row(e.f2).dot(t22));
            if (d1.norm() > 1e-12) d1.normalize();
            if (d2.norm() > 1e-12) d2.normalize();
            double best = 1e10;
            for (int k = 0; k < 4; ++k) {
                double cs = std::cos(k * kQuarterTurn), sn = std::sin(k * kQuarterTurn);
                Vector2d rot(cs * d2.x() - sn * d2.y(), sn * d2.x() + cs * d2.y());
                double dot = std::max(-1.0, std::min(1.0, d1.dot(rot)));
                double ang = std::acos(std::fabs(dot));
                if (ang < best) { best = ang; matchingVal = k; }
            }
            // The propagated rotation = parent rotation + matching
            // Since faceRot[f] is in units of π/2, matchingVal is also
            int rot = mod4(faceRot[f] + matchingVal);
            faceRot[nb.first] = rot;
            q.push(nb.first);
        }
    }

    // For each cut edge, the rotation between the two sides.
    // cutRot = (orientation_f2 − orientation_f1) mod 4.
    for (int ci = 0; ci < nCuts; ++ci) {
        const Edge& e = edgeList[cutEdges_[ci]];
        if (faceRot[e.f1] >= 0 && faceRot[e.f2] >= 0) {
            cutRot_[ci] = mod4(faceRot[e.f2] - faceRot[e.f1]);
        } else {
            cutRot_[ci] = 0;
        }
    }
}

void MIQQuad::buildCotLaplacianCut(SparseMatrix<double>& L) const
{
    L.resize(nCV, nCV);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nCV);

    for (int fi = 0; fi < nF; ++fi) {
        int i0 = cutFaces[fi * 3], i1 = cutFaces[fi * 3 + 1], i2 = cutFaces[fi * 3 + 2];
        if (i0 < 0 || i1 < 0 || i2 < 0) continue;
        if (i0 == i1 || i1 == i2 || i2 == i0) continue;

        const Vector3d p0 = cutVertPos.row(i0);
        const Vector3d p1 = cutVertPos.row(i1);
        const Vector3d p2 = cutVertPos.row(i2);

        auto addW = [&](int a, int b, double w) {
            if (!std::isfinite(w) || w <= 0.0) return;
            diag(a) += w; diag(b) += w;
            trips.emplace_back(a, b, -w);
            trips.emplace_back(b, a, -w);
        };
        addW(i1, i2, 0.5 * cotan(p1, p0, p2));
        addW(i2, i0, 0.5 * cotan(p2, p1, p0));
        addW(i0, i1, 0.5 * cotan(p0, p2, p1));
    }
    for (int i = 0; i < nCV; ++i)
        trips.emplace_back(i, i, diag(i) + 1e-8);

    L.setFromTriplets(trips.begin(), trips.end());
}

bool MIQQuad::solveSeamlessParam()
{
    buildCutGraph();
    if (!buildCutMesh()) return false;

    const int nCuts = static_cast<int>(cutEdges_.size());

    // If no cuts, the mesh is a topological disk → plain Poisson is enough.
    if (nCuts == 0) {
        // plain Poisson on the cut mesh (identical to original)
        SparseMatrix<double> L;
        buildCotLaplacianCut(L);
        VectorXd div1 = VectorXd::Zero(nCV);
        VectorXd div2 = VectorXd::Zero(nCV);
        for (int fi = 0; fi < nF; ++fi) {
            int i0 = cutFaces[fi * 3], i1 = cutFaces[fi * 3 + 1], i2 = cutFaces[fi * 3 + 2];
            if (i0 < 0 || i1 < 0 || i2 < 0) continue;
            const Vector3d p0 = cutVertPos.row(i0), p1 = cutVertPos.row(i1), p2 = cutVertPos.row(i2);
            const Vector3d n = (p1 - p0).cross(p2 - p0);
            const Vector3d fd1 = faceD1.row(fi), fd2 = faceD2.row(fi);
            div1(i0) += 0.5 * fd1.dot(n.cross(p2 - p1));
            div1(i1) += 0.5 * fd1.dot(n.cross(p0 - p2));
            div1(i2) += 0.5 * fd1.dot(n.cross(p1 - p0));
            div2(i0) += 0.5 * fd2.dot(n.cross(p2 - p1));
            div2(i1) += 0.5 * fd2.dot(n.cross(p0 - p2));
            div2(i2) += 0.5 * fd2.dot(n.cross(p1 - p0));
        }
        L.coeffRef(0, 0) += 1e6; div1(0) = 0.0; div2(0) = 0.0;

        SimplicialLDLT<SparseMatrix<double>> solver;
        solver.compute(L);
        if (solver.info() != Success) return false;
        VectorXd cu = solver.solve(div1);
        VectorXd cv = solver.solve(div2);
        if (solver.info() != Success) return false;

        // project back to base mesh
        UV.resize(nV, 2);
        UV.setZero();
        std::vector<int> cnt(nV, 0);
        for (int i = 0; i < nCV; ++i) {
            int bv = cutParent[i];
            UV(bv, 0) += cu(i); UV(bv, 1) += cv(i);
            cnt[bv]++;
        }
        for (int v = 0; v < nV; ++v)
            if (cnt[v] > 0) UV.row(v) /= static_cast<double>(cnt[v]);

        double uMin = UV.col(0).minCoeff(), uMax = UV.col(0).maxCoeff();
        double vMin = UV.col(1).minCoeff(), vMax = UV.col(1).maxCoeff();
        double s = std::max(uMax - uMin, vMax - vMin);
        if (s < 1e-10) s = 1.0;
        UV.col(0) = (UV.col(0).array() - uMin) / s;
        UV.col(1) = (UV.col(1).array() - vMin) / s;
        return true;
    }

    determineCutRotations();

    // --- stage 1: Poisson on cut mesh (continuous solve) ---
    SparseMatrix<double> L;
    buildCotLaplacianCut(L);

    VectorXd div1 = VectorXd::Zero(nCV);
    VectorXd div2 = VectorXd::Zero(nCV);
    for (int fi = 0; fi < nF; ++fi) {
        int i0 = cutFaces[fi * 3], i1 = cutFaces[fi * 3 + 1], i2 = cutFaces[fi * 3 + 2];
        if (i0 < 0 || i1 < 0 || i2 < 0 || i0 == i1 || i1 == i2 || i2 == i0) continue;
        const Vector3d p0 = cutVertPos.row(i0), p1 = cutVertPos.row(i1), p2 = cutVertPos.row(i2);
        const Vector3d n = (p1 - p0).cross(p2 - p0);
        const Vector3d fd1 = faceD1.row(fi), fd2 = faceD2.row(fi);
        div1(i0) += 0.5 * fd1.dot(n.cross(p2 - p1));
        div1(i1) += 0.5 * fd1.dot(n.cross(p0 - p2));
        div1(i2) += 0.5 * fd1.dot(n.cross(p1 - p0));
        div2(i0) += 0.5 * fd2.dot(n.cross(p2 - p1));
        div2(i1) += 0.5 * fd2.dot(n.cross(p0 - p2));
        div2(i2) += 0.5 * fd2.dot(n.cross(p1 - p0));
    }

    // Pin one vertex to fix nullspace
    L.coeffRef(0, 0) += 1e6; div1(0) = 0.0; div2(0) = 0.0;

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Success) return false;
    VectorXd cuTilde = solver.solve(div1);
    VectorXd cvTilde = solver.solve(div2);
    if (solver.info() != Success) return false;

    // --- stage 2: period constraints (joint U-V KKT) ---
    // For each cut edge ci with rotation r = cutRot_[ci]:
    //   (u_B, v_B) = Rot(r·π/2)(u_A, v_A) + (tx, ty)
    // r=0: u_B − u_A = tx ,  v_B − v_A = ty
    // r=1: u_B + v_A = tx ,  v_B − u_A = ty   ← U,V coupled
    // r=2: u_B + u_A = tx ,  v_B + v_A = ty
    // r=3: u_B − v_A = tx ,  v_B + u_A = ty   ← U,V coupled

    // Build cut-corner map
    struct CutCorners { int ca1, ca2, cb1, cb2; };
    std::vector<CutCorners> cc(nCuts);
    for (int ci = 0; ci < nCuts; ++ci) {
        const Edge& e = edgeList[cutEdges_[ci]];
        auto loc = [&](int fi, int vert) {
            for (int k = 0; k < 3; ++k) if (faces[fi * 3 + k] == vert) return fi * 3 + k;
            return -1;
        };
        cc[ci].ca1 = cutFaces[loc(e.f1, e.v1)];
        cc[ci].ca2 = cutFaces[loc(e.f1, e.v2)];
        cc[ci].cb1 = cutFaces[loc(e.f2, e.v1)];
        cc[ci].cb2 = cutFaces[loc(e.f2, e.v2)];
        if (cc[ci].ca1 < 0 || cc[ci].cb1 < 0) return false;
    }

    // Compute continuous-period mismatch mu
    VectorXd muU(nCuts), muV(nCuts);
    for (int ci = 0; ci < nCuts; ++ci) {
        int r = cutRot_[ci];
        double ca_u = cuTilde(cc[ci].ca1), ca_v = cvTilde(cc[ci].ca1);
        double cb_u = cuTilde(cc[ci].cb1), cb_v = cvTilde(cc[ci].cb1);
        double ru, rv;
        if      (r == 0) { ru =  ca_u;      rv =  ca_v; }
        else if (r == 1) { ru = -ca_v;      rv =  ca_u; }
        else if (r == 2) { ru = -ca_u;      rv = -ca_v; }
        else             { ru =  ca_v;      rv = -ca_u; }
        muU(ci) = cb_u - ru;
        muV(ci) = cb_v - rv;
    }

    // Round to integer
    for (int ci = 0; ci < nCuts; ++ci) {
        cutTx(ci) = static_cast<int>(std::round(muU(ci)));
        cutTy(ci) = static_cast<int>(std::round(muV(ci)));
    }

    // --- stage 3: joint U-V KKT harmonic correction ---
    // Variable order: [u(0..nCV-1), v(0..nCV-1), λ_u(0..nCuts-1), λ_v(0..nCuts-1)]
    const int M = nCV;
    const int J = 2 * M + 2 * nCuts;
    std::vector<Triplet<double>> kkT;

    // Copy L into top-left 2×2 blocks
    for (int k = 0; k < L.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
            kkT.emplace_back(it.row(),       it.col(),       it.value());  // L_uu
            kkT.emplace_back(M + it.row(),   M + it.col(),   it.value());  // L_vv
        }
    kkT.emplace_back(0,     0,     1e6);  // pin  ψ_u(0)=0
    kkT.emplace_back(M,     M,     1e6);  // pin  ψ_v(0)=0

    // Period constraint rows (2·nCuts constraints, rows 2M .. 2M+2nCuts-1)
    for (int ci = 0; ci < nCuts; ++ci) {
        int r = cutRot_[ci];
        int ca = cc[ci].ca1, cb = cc[ci].cb1;
        int ru = 2 * M + ci;          // λ_u row
        int rv = 2 * M + nCuts + ci;  // λ_v row

        // U constraint:  u(cb) + c_uu·u(ca) + c_uv·v(ca) = tx
        // V constraint:  v(cb) + c_vu·u(ca) + c_vv·v(ca) = ty
        double cuu=0, cuv=0, cvu=0, cvv=0;
        //  (u_B,v_B) − Rot(r)(u_A,v_A) = (tx,ty)
        //  u_B − [Rot·(u_A,v_A)]_u = tx
        //  v_B − [Rot·(u_A,v_A)]_v = ty
        if      (r == 0) { cuu = -1;               cvv = -1;              }
        else if (r == 1) {             cuv =  1;   cvu = -1;              }
        else if (r == 2) { cuu =  1;               cvv =  1;              }
        else             {             cuv = -1;   cvu =  1;              }

        if (std::fabs(cuu) > 1e-12) {
            kkT.emplace_back(ru, ca,  cuu);  kkT.emplace_back(ca, ru,  cuu);
        }
        if (std::fabs(cuv) > 1e-12) {
            kkT.emplace_back(ru, M+ca, cuv); kkT.emplace_back(M+ca, ru, cuv);
        }
        if (std::fabs(cvu) > 1e-12) {
            kkT.emplace_back(rv, ca,  cvu);  kkT.emplace_back(ca, rv,  cvu);
        }
        if (std::fabs(cvv) > 1e-12) {
            kkT.emplace_back(rv, M+ca, cvv); kkT.emplace_back(M+ca, rv, cvv);
        }
        // cb column: always +1
        kkT.emplace_back(ru, cb,  1.0);  kkT.emplace_back(cb,  ru, 1.0);
        kkT.emplace_back(rv, M+cb, 1.0); kkT.emplace_back(M+cb, rv, 1.0);
    }

    SparseMatrix<double> KKTmat(J, J);
    KKTmat.setFromTriplets(kkT.begin(), kkT.end());

    VectorXd rhs = VectorXd::Zero(J);
    for (int ci = 0; ci < nCuts; ++ci) {
        rhs(2 * M + ci)            = cutTx(ci) - muU(ci);
        rhs(2 * M + nCuts + ci)    = cutTy(ci) - muV(ci);
    }

    SimplicialLDLT<SparseMatrix<double>> kktSolver;
    kktSolver.compute(KKTmat);
    if (kktSolver.info() != Success) return false;
    VectorXd sol = kktSolver.solve(rhs);
    if (kktSolver.info() != Success) return false;

    VectorXd cu = cuTilde + sol.head(M);
    VectorXd cv = cvTilde + sol.segment(M, M);

    // --- project back to base mesh ---
    UV.resize(nV, 2);
    UV.setZero();
    std::vector<int> cnt(nV, 0);
    for (int i = 0; i < M; ++i) {
        int bv = cutParent[i];
        UV(bv, 0) += cu(i); UV(bv, 1) += cv(i);
        cnt[bv]++;
    }
    for (int v = 0; v < nV; ++v)
        if (cnt[v] > 0) UV.row(v) /= static_cast<double>(cnt[v]);

    // Normalize
    double uMin = UV.col(0).minCoeff(), uMax = UV.col(0).maxCoeff();
    double vMin = UV.col(1).minCoeff(), vMax = UV.col(1).maxCoeff();
    double s = std::max(uMax - uMin, vMax - vMin);
    if (s < 1e-10) s = 1.0;
    UV.col(0) = (UV.col(0).array() - uMin) / s;
    UV.col(1) = (UV.col(1).array() - vMin) / s;
    return true;
}

void MIQQuad::parameterize()
{
    if (!initMeshData()) return;

    initCrossField();
    if (!solveCrossFieldIP()) return;
    buildTargetDirs();

    // Phase 2: seamless parameterization with cut graph
    if (solveSeamlessParam()) {
        for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
            v->uv = Vector2d(UV(v->index, 0), UV(v->index, 1));
        }
        normalize();
    } else {
        // Fallback: plain Poisson (existing behaviour)
        solvePoisson();
        for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
            v->uv = Vector2d(UV(v->index, 0), UV(v->index, 1));
        }
        normalize();
    }
}
