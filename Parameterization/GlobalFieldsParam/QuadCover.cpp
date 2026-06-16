#include "QuadCover.h"
#include "PrincipalCurvatureField.h"
#include <map>
#include <algorithm>
#include <queue>
#include <cmath>
#include <Eigen/IterativeLinearSolvers>

using namespace Eigen;

namespace {
int mod4(int x) {
    int r = x % 4;
    return r < 0 ? r + 4 : r;
}

struct DisjointSet {
    std::vector<int> parent, rank;
    explicit DisjointSet(int n = 0) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    int find(int x) {
        if (parent[x] == x) return x;
        parent[x] = find(parent[x]);
        return parent[x];
    }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) rank[a]++;
    }
};

double normalizeAngleHalfPi(double a) {
    const double k = 0.5 * M_PI;
    a = std::fmod(a, k);
    if (a < 0.0) a += k;
    return a;
}
}

QuadCover::QuadCover(Mesh& mesh0)
    : GlobalFieldsParameterization(mesh0),
      smoothIters(3),
      requirePureQuads(true) {}

void QuadCover::setCrossField(const MatrixXd& dirs) {
    GlobalFieldsParameterization::setCrossField(dirs);
    externalFaceDirs = dirs;  // preserve original input for principal curvature estimation
}

void QuadCover::setSmoothIterations(int iters) {
    smoothIters = std::max(0, iters);
}

void QuadCover::setRequirePureQuads(bool pureQuads) {
    requirePureQuads = pureQuads;
}

void QuadCover::buildLocalFrame(const Vector3d& n, Vector3d& t1, Vector3d& t2) {
    Vector3d ref(1, 0, 0);
    if (std::fabs(n.dot(ref)) > 0.9) ref = Vector3d(0, 1, 0);
    t1 = (ref - ref.dot(n) * n).normalized();
    t2 = n.cross(t1).normalized();
}

double QuadCover::cotan(const Vector3d& a, const Vector3d& b, const Vector3d& c) const {
    Vector3d u = a - b, v = c - b;
    double d = u.dot(v), cr = u.cross(v).norm();
    if (cr < 1e-12) return 0.0;
    return d / cr;
}

void QuadCover::buildCotLaplacian(SparseMatrix<double>& L) {
    L.resize(nVerts, nVerts);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nVerts);

    for (const Edge& e : edgeList) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        int vi = e.v1, vj = e.v2;

        int o1 = -1, o2 = -1;
        int f1v0 = faces[e.f1 * 3], f1v1 = faces[e.f1 * 3 + 1], f1v2 = faces[e.f1 * 3 + 2];
        int f2v0 = faces[e.f2 * 3], f2v1 = faces[e.f2 * 3 + 1], f2v2 = faces[e.f2 * 3 + 2];

        if (f1v0 != vi && f1v0 != vj) o1 = f1v0;
        else if (f1v1 != vi && f1v1 != vj) o1 = f1v1;
        else o1 = f1v2;

        if (f2v0 != vi && f2v0 != vj) o2 = f2v0;
        else if (f2v1 != vi && f2v1 != vj) o2 = f2v1;
        else o2 = f2v2;

        double w = 0.5 * (cotan(vertPos.row(o1), vertPos.row(vi), vertPos.row(vj)) +
                          cotan(vertPos.row(o2), vertPos.row(vi), vertPos.row(vj)));
        if (w < 0) w = 0;

        diag(vi) += w; diag(vj) += w;
        trips.emplace_back(vi, vj, -w);
        trips.emplace_back(vj, vi, -w);
    }
    for (int i = 0; i < nVerts; ++i)
        trips.emplace_back(i, i, diag(i) + 1e-8);
    L.setFromTriplets(trips.begin(), trips.end());
}

bool QuadCover::initMeshData() {
    nVerts = (int)mesh.vertices.size();
    nFaces = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f)
        if (!f->isBoundary()) nFaces++;

    if (nVerts < 3 || nFaces < 1) return false;

    vertPos.resize(nVerts, 3);
    faces.resize(nFaces * 3);

    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v)
        vertPos.row(v->index) = v->position;

    int fi = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (!f->isBoundary()) {
            faces[fi * 3] = f->he->vertex->index;
            faces[fi * 3 + 1] = f->he->next->vertex->index;
            faces[fi * 3 + 2] = f->he->next->next->vertex->index;
            fi++;
        }
    }

    std::map<std::pair<int, int>, int> edgeMap;
    edgeList.clear();
    for (int fi2 = 0; fi2 < nFaces; ++fi2) {
        for (int k = 0; k < 3; ++k) {
            int v1 = faces[fi2 * 3 + k], v2 = faces[fi2 * 3 + (k + 1) % 3];
            if (v1 > v2) std::swap(v1, v2);
            auto key = std::make_pair(v1, v2);
            auto it = edgeMap.find(key);
            if (it == edgeMap.end()) {
                edgeMap[key] = (int)edgeList.size();
                Edge e; e.v1 = v1; e.v2 = v2; e.f1 = fi2; e.f2 = -1;
                edgeList.push_back(e);
            } else {
                edgeList[it->second].f2 = fi2;
            }
        }
    }
    nEdges = (int)edgeList.size();
    return true;
}

void QuadCover::computeVertexNormals() {
    vertexNormals.resize(nVerts, 3);
    vertexNormals.setZero();

    for (int fi = 0; fi < nFaces; ++fi) {
        int v0 = faces[fi * 3], v1 = faces[fi * 3 + 1], v2 = faces[fi * 3 + 2];
        Vector3d p0 = vertPos.row(v0), p1 = vertPos.row(v1), p2 = vertPos.row(v2);
        Vector3d fn = (p1 - p0).cross(p2 - p0);
        double area = 0.5 * fn.norm();
        fn.normalize();
        for (int k = 0; k < 3; ++k)
            vertexNormals.row(faces[fi * 3 + k]) += area * fn;
    }
    for (int i = 0; i < nVerts; ++i) {
        double len = vertexNormals.row(i).norm();
        if (len > 1e-12) vertexNormals.row(i) /= len;
    }
}

void QuadCover::estimatePrincipalCurvature() {
    faceNormals.resize(nFaces, 3);
    faceDirs_.resize(nFaces, 3);

    for (int fi = 0; fi < nFaces; ++fi) {
        int v0 = faces[fi * 3], v1 = faces[fi * 3 + 1], v2 = faces[fi * 3 + 2];
        Vector3d p0 = vertPos.row(v0), p1 = vertPos.row(v1), p2 = vertPos.row(v2);

        Vector3d fn = (p1 - p0).cross(p2 - p0);
        double area = 0.5 * fn.norm();
        fn /= std::max(2.0 * area, 1e-12);
        faceNormals.row(fi) = fn;

        if (hasExternalField_ && externalFaceDirs.rows() == nFaces) {
            Vector3d dir = externalFaceDirs.row(fi);
            dir -= dir.dot(fn) * fn;
            double dlen = dir.norm();
            if (dlen > 1e-12) faceDirs_.row(fi) = dir / dlen;
            else {
                Vector3d t1, t2;
                buildLocalFrame(fn, t1, t2);
                faceDirs_.row(fi) = t1;
            }
            continue;
        }
    }

    // If the user didn't provide an external cross field, estimate per-face principal directions.
    if (!(hasExternalField_ && externalFaceDirs.rows() == nFaces)) {
        PrincipalCurvatureField::estimateFromWeingarten(
            vertPos, faces, vertexNormals, faceNormals, faceDirs_);
    }
}

void QuadCover::computeFaceTheta() {
    faceTheta.resize(nFaces);
    for (int fi = 0; fi < nFaces; ++fi) {
        Vector3d n = faceNormals.row(fi).normalized();
        Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);
        Vector3d d = faceDirs_.row(fi);
        d -= d.dot(n) * n;
        if (d.norm() < 1e-12) d = t1;
        else d.normalize();
        faceTheta(fi) = normalizeAngleHalfPi(std::atan2(d.dot(t2), d.dot(t1)));
    }
}

void QuadCover::syncFaceDirsFromTheta() {
    for (int fi = 0; fi < nFaces; ++fi) {
        Vector3d n = faceNormals.row(fi).normalized();
        Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);
        double c = std::cos(faceTheta(fi)), s = std::sin(faceTheta(fi));
        faceDirs_.row(fi) = (c * t1 + s * t2).normalized();
    }
}

void QuadCover::smoothCrossField(int iters) {
    if (iters <= 0 || nFaces < 1) return;
    const double kQuarter = 0.5 * M_PI;

    for (int iter = 0; iter < iters; ++iter) {
        VectorXd newTheta = faceTheta;
        for (int eIdx = 0; eIdx < nEdges; ++eIdx) {
            const Edge& e = edgeList[eIdx];
            if (e.f1 < 0 || e.f2 < 0) continue;

            double target = faceTheta(e.f2) + kQuarter * matching[eIdx];
            target = normalizeAngleHalfPi(target);

            const double alpha = 0.5;
            newTheta(e.f1) = normalizeAngleHalfPi(alpha * faceTheta(e.f1) + (1.0 - alpha) * target);
        }
        faceTheta = newTheta;
    }
    syncFaceDirsFromTheta();
}

int QuadCover::signedMatching(int fromFace, int toFace) const {
    for (int ei = 0; ei < nEdges; ++ei) {
        const Edge& e = edgeList[ei];
        if (e.f1 == fromFace && e.f2 == toFace) return matching[ei];
        if (e.f2 == fromFace && e.f1 == toFace) return -matching[ei];
    }
    return 0;
}

void QuadCover::computeMatching() {
    matching.resize(nEdges);
    matching.setZero();

    for (int eIdx = 0; eIdx < nEdges; ++eIdx) {
        const Edge& e = edgeList[eIdx];
        if (e.f1 < 0 || e.f2 < 0) continue;

        Vector3d n1 = faceNormals.row(e.f1);
        Vector3d n2 = faceNormals.row(e.f2);
        Vector3d t11, t12, t21, t22;
        buildLocalFrame(n1, t11, t12);
        buildLocalFrame(n2, t21, t22);

        Vector2d d1_2d(faceDirs_.row(e.f1).dot(t11), faceDirs_.row(e.f1).dot(t12));
        Vector2d d2_2d(faceDirs_.row(e.f2).dot(t21), faceDirs_.row(e.f2).dot(t12));
        if (d1_2d.norm() > 1e-12) d1_2d.normalize();
        if (d2_2d.norm() > 1e-12) d2_2d.normalize();

        double best = 1e10;
        int bestK = 0;
        for (int k = 0; k < 4; ++k) {
            double cs = std::cos(k * M_PI_2), sn = std::sin(k * M_PI_2);
            Vector2d rot(cs * d2_2d.x() - sn * d2_2d.y(), sn * d2_2d.x() + cs * d2_2d.y());
            double dot = d1_2d.dot(rot);
            dot = std::max(-1.0, std::min(1.0, dot));
            double ang = std::acos(std::fabs(dot));
            if (ang < best) { best = ang; bestK = k; }
        }
        matching[eIdx] = bestK;
    }
}

void QuadCover::computeLayerShift() {
    layerShift = VectorXd::Zero(nVerts);

    std::vector<std::vector<int>> incidentFaces(nVerts);
    std::vector<char> boundaryVertex(nVerts, 0);
    for (int fi = 0; fi < nFaces; ++fi)
        for (int k = 0; k < 3; ++k)
            incidentFaces[faces[fi * 3 + k]].push_back(fi);

    for (int ei = 0; ei < nEdges; ++ei) {
        const Edge& e = edgeList[ei];
        if (e.f1 < 0 || e.f2 < 0) {
            boundaryVertex[e.v1] = 1;
            boundaryVertex[e.v2] = 1;
        }
    }

    for (int v = 0; v < nVerts; ++v) {
        if (boundaryVertex[v] || incidentFaces[v].size() < 3) continue;

        Vector3d n = vertexNormals.row(v);
        if (n.norm() < 1e-12) continue;
        n.normalize();
        Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);
        const Vector3d pv = vertPos.row(v);

        std::vector<std::pair<double, int>> ordered;
        for (int fi : incidentFaces[v]) {
            Vector3d c(0, 0, 0);
            for (int k = 0; k < 3; ++k) c += vertPos.row(faces[fi * 3 + k]);
            c = c / 3.0 - pv;
            c -= c.dot(n) * n;
            if (c.norm() < 1e-12) continue;
            ordered.emplace_back(std::atan2(c.dot(t2), c.dot(t1)), fi);
        }
        if (ordered.size() < 3) continue;
        std::sort(ordered.begin(), ordered.end());

        int holonomy = 0;
        for (size_t i = 0; i < ordered.size(); ++i) {
            int f0 = ordered[i].second;
            int f1 = ordered[(i + 1) % ordered.size()].second;
            holonomy += signedMatching(f0, f1);
        }
        int h = mod4(holonomy);
        if (h > 2) h -= 4;
        layerShift(v) = (double)h / 4.0;
    }
}

bool QuadCover::buildBranchCover(CoverData& cover) {
    const int sheets = 4;
    const int rawCoverCorners = nFaces * sheets * 3;
    DisjointSet dsu(rawCoverCorners);

    auto cornerId = [&](int fi, int local, int s) {
        return (fi * sheets + mod4(s)) * 3 + local;
    };
    auto localCorner = [&](int fi, int vertex) {
        for (int k = 0; k < 3; ++k)
            if (faces[fi * 3 + k] == vertex) return k;
        return -1;
    };

    for (int ei = 0; ei < nEdges; ++ei) {
        const Edge& e = edgeList[ei];
        if (e.f1 < 0 || e.f2 < 0) continue;
        const int r = matching[ei];
        const int f1v1 = localCorner(e.f1, e.v1);
        const int f1v2 = localCorner(e.f1, e.v2);
        const int f2v1 = localCorner(e.f2, e.v1);
        const int f2v2 = localCorner(e.f2, e.v2);
        if (f1v1 < 0 || f1v2 < 0 || f2v1 < 0 || f2v2 < 0) return false;

        for (int s = 0; s < sheets; ++s) {
            const int t = mod4(s + r);
            dsu.unite(cornerId(e.f1, f1v1, s), cornerId(e.f2, f2v1, t));
            dsu.unite(cornerId(e.f1, f1v2, s), cornerId(e.f2, f2v2, t));
        }
    }

    std::map<int, int> repToIndex;
    std::vector<int> coverIndex(rawCoverCorners, -1);
    std::vector<Vector3d> cpos;
    cover.baseToCover.assign(nVerts, {});
    cover.baseSheet0.assign(nVerts, -1);

    for (int fi = 0; fi < nFaces; ++fi) {
        for (int s = 0; s < sheets; ++s) {
            for (int k = 0; k < 3; ++k) {
                const int raw = cornerId(fi, k, s);
                const int v = faces[fi * 3 + k];
                const int rep = dsu.find(raw);
                auto it = repToIndex.find(rep);
                int idx;
                if (it == repToIndex.end()) {
                    idx = (int)cpos.size();
                    repToIndex[rep] = idx;
                    coverIndex[raw] = idx;
                    cpos.push_back(vertPos.row(v));
                } else {
                    idx = it->second;
                    coverIndex[raw] = idx;
                }
                auto& copies = cover.baseToCover[v];
                if (std::find(copies.begin(), copies.end(), idx) == copies.end())
                    copies.push_back(idx);
                if (s == 0 && cover.baseSheet0[v] < 0)
                    cover.baseSheet0[v] = idx;
            }
        }
    }

    cover.nCoverVerts = (int)cpos.size();
    cover.nCoverFaces = nFaces * sheets;
    cover.coverPos.resize(cover.nCoverVerts, 3);
    for (int i = 0; i < cover.nCoverVerts; ++i)
        cover.coverPos.row(i) = cpos[i];

    cover.coverFaces.resize(cover.nCoverFaces, 3);
    cover.coverD1.resize(cover.nCoverFaces, 3);
    cover.coverD2.resize(cover.nCoverFaces, 3);

    for (int fi = 0; fi < nFaces; ++fi) {
        Vector3d n = faceNormals.row(fi).normalized();
        Vector3d base = faceDirs_.row(fi);
        base -= base.dot(n) * n;
        if (base.norm() < 1e-12) {
            Vector3d t1, t2;
            buildLocalFrame(n, t1, t2);
            base = t1;
        } else base.normalize();
        Vector3d ortho = n.cross(base).normalized();

        for (int s = 0; s < sheets; ++s) {
            const int row = fi * sheets + s;
            cover.coverFaces(row, 0) = coverIndex[cornerId(fi, 0, s)];
            cover.coverFaces(row, 1) = coverIndex[cornerId(fi, 1, s)];
            cover.coverFaces(row, 2) = coverIndex[cornerId(fi, 2, s)];

            const double a = s * 0.5 * M_PI;
            cover.coverD1.row(row) = (std::cos(a) * base + std::sin(a) * ortho).normalized();
            cover.coverD2.row(row) = (-std::sin(a) * base + std::cos(a) * ortho).normalized();
        }
    }
    return cover.nCoverVerts >= 3;
}

bool QuadCover::solveCoverPoisson(const MatrixXd& coverPos,
                                  const MatrixXi& coverFaces,
                                  const MatrixXd& coverD1,
                                  const MatrixXd& coverD2,
                                  MatrixXd& coverUV) {
    const int nV = (int)coverPos.rows();
    const int nF = (int)coverFaces.rows();
    if (nV < 3 || nF < 1) return false;

    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nV);
    VectorXd rhsU = VectorXd::Zero(nV);
    VectorXd rhsV = VectorXd::Zero(nV);
    std::vector<std::vector<int>> adj(nV);

    auto addWeight = [&](int i, int j, double w) {
        if (!std::isfinite(w) || w <= 0.0) return;
        diag(i) += w; diag(j) += w;
        trips.emplace_back(i, j, -w);
        trips.emplace_back(j, i, -w);
        adj[i].push_back(j);
        adj[j].push_back(i);
    };

    for (int fi = 0; fi < nF; ++fi) {
        const int i0 = coverFaces(fi, 0), i1 = coverFaces(fi, 1), i2 = coverFaces(fi, 2);
        if (i0 < 0 || i0 >= nV || i1 < 0 || i1 >= nV || i2 < 0 || i2 >= nV) return false;
        if (i0 == i1 || i1 == i2 || i2 == i0) continue;

        const Vector3d p0 = coverPos.row(i0);
        const Vector3d p1 = coverPos.row(i1);
        const Vector3d p2 = coverPos.row(i2);
        Vector3d n = (p1 - p0).cross(p2 - p0);
        const double dblArea = n.norm();
        if (dblArea < 1e-14) continue;
        n /= dblArea;
        const double area = 0.5 * dblArea;

        addWeight(i1, i2, 0.5 * cotan(p1, p0, p2));
        addWeight(i2, i0, 0.5 * cotan(p2, p1, p0));
        addWeight(i0, i1, 0.5 * cotan(p0, p2, p1));

        const Vector3d grad0 = n.cross(p2 - p1) / dblArea;
        const Vector3d grad1 = n.cross(p0 - p2) / dblArea;
        const Vector3d grad2 = n.cross(p1 - p0) / dblArea;
        const Vector3d xu = coverD1.row(fi);
        const Vector3d xv = coverD2.row(fi);

        rhsU(i0) += area * xu.dot(grad0);
        rhsU(i1) += area * xu.dot(grad1);
        rhsU(i2) += area * xu.dot(grad2);
        rhsV(i0) += area * xv.dot(grad0);
        rhsV(i1) += area * xv.dot(grad1);
        rhsV(i2) += area * xv.dot(grad2);
    }

    for (int i = 0; i < nV; ++i)
        trips.emplace_back(i, i, diag(i) + 1e-8);

    SparseMatrix<double> L(nV, nV);
    L.setFromTriplets(trips.begin(), trips.end());

    std::vector<char> seen(nV, 0);
    for (int seed = 0; seed < nV; ++seed) {
        if (seen[seed]) continue;
        std::queue<int> q;
        q.push(seed);
        seen[seed] = 1;
        L.coeffRef(seed, seed) += 1e8;
        rhsU(seed) = 0.0;
        rhsV(seed) = 0.0;
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (int nb : adj[v]) {
                if (!seen[nb]) { seen[nb] = 1; q.push(nb); }
            }
        }
    }
    L.makeCompressed();

    ConjugateGradient<SparseMatrix<double>, Lower | Upper, DiagonalPreconditioner<double>> solver;
    solver.setMaxIterations(std::max(200, nV * 2));
    solver.setTolerance(1e-8);
    solver.compute(L);
    if (solver.info() != Success) return false;

    VectorXd u = solver.solve(rhsU);
    VectorXd v = solver.solve(rhsV);
    if (solver.info() != Success) return false;

    coverUV.resize(nV, 2);
    coverUV.col(0) = u;
    coverUV.col(1) = v;
    return true;
}

bool QuadCover::integrateOnCover(CoverData& cover) {
    return solveCoverPoisson(cover.coverPos, cover.coverFaces,
                             cover.coverD1, cover.coverD2, cover.coverUV);
}

void QuadCover::buildCoverLaplacian(const CoverData& cover,
                                    SparseMatrix<double>& L,
                                    std::vector<std::vector<int>>& adj) const {
    const int nV = cover.nCoverVerts;
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nV);
    adj.assign(nV, {});

    auto addWeight = [&](int i, int j, double w) {
        if (!std::isfinite(w) || w <= 0.0) return;
        diag(i) += w; diag(j) += w;
        trips.emplace_back(i, j, -w);
        trips.emplace_back(j, i, -w);
        adj[i].push_back(j);
        adj[j].push_back(i);
    };

    for (int fi = 0; fi < cover.nCoverFaces; ++fi) {
        const int i0 = cover.coverFaces(fi, 0);
        const int i1 = cover.coverFaces(fi, 1);
        const int i2 = cover.coverFaces(fi, 2);
        const Vector3d p0 = cover.coverPos.row(i0);
        const Vector3d p1 = cover.coverPos.row(i1);
        const Vector3d p2 = cover.coverPos.row(i2);
        addWeight(i1, i2, 0.5 * cotan(p1, p0, p2));
        addWeight(i2, i0, 0.5 * cotan(p2, p1, p0));
        addWeight(i0, i1, 0.5 * cotan(p0, p2, p1));
    }
    for (int i = 0; i < nV; ++i)
        trips.emplace_back(i, i, diag(i) + 1e-8);

    L.resize(nV, nV);
    L.setFromTriplets(trips.begin(), trips.end());
}

int QuadCover::coverVertexFor(const CoverData& cover, int baseVert) const {
    if (baseVert < 0 || baseVert >= nVerts) return -1;
    if (cover.baseSheet0[baseVert] >= 0) return cover.baseSheet0[baseVert];
    if (!cover.baseToCover[baseVert].empty()) return cover.baseToCover[baseVert][0];
    return -1;
}

int QuadCover::findSharedEdge(int f1, int f2) const {
    for (int ei = 0; ei < nEdges; ++ei) {
        const Edge& e = edgeList[ei];
        if ((e.f1 == f1 && e.f2 == f2) || (e.f1 == f2 && e.f2 == f1))
            return ei;
    }
    return -1;
}

int QuadCover::countBoundaryComponents() const {
    std::vector<char> isBnd(nVerts, 0);
    for (int ei = 0; ei < nEdges; ++ei) {
        const Edge& e = edgeList[ei];
        if (e.f2 < 0) { isBnd[e.v1] = 1; isBnd[e.v2] = 1; }
    }

    std::vector<std::vector<int>> adj(nVerts);
    for (int ei = 0; ei < nEdges; ++ei) {
        const Edge& e = edgeList[ei];
        if (e.f2 < 0) {
            adj[e.v1].push_back(e.v2);
            adj[e.v2].push_back(e.v1);
        }
    }

    std::vector<char> seen(nVerts, 0);
    int components = 0;
    for (int v = 0; v < nVerts; ++v) {
        if (!isBnd[v] || seen[v]) continue;
        components++;
        std::queue<int> q;
        q.push(v);
        seen[v] = 1;
        while (!q.empty()) {
            int cur = q.front(); q.pop();
            for (int nb : adj[cur]) {
                if (!seen[nb]) { seen[nb] = 1; q.push(nb); }
            }
        }
    }
    return components;
}

void QuadCover::buildDualCyclePath(int cotreeEdgeIdx,
                                   const std::vector<int>& parentFace,
                                   CutPath& path) const {
    path.edges.clear();
    path.forward.clear();

    const Edge& cotree = edgeList[cotreeEdgeIdx];
    int fA = cotree.f1, fB = cotree.f2;

    auto addOrientedEdge = [&](int ei, int fromFace, int toFace) {
        const Edge& e = edgeList[ei];
        int va = -1, vb = -1;
        if (e.f1 == fromFace && e.f2 == toFace) { va = e.v1; vb = e.v2; }
        else if (e.f2 == fromFace && e.f1 == toFace) { va = e.v2; vb = e.v1; }
        else {
            for (int k = 0; k < 3; ++k) {
                int v0 = faces[fromFace * 3 + k];
                int v1 = faces[fromFace * 3 + (k + 1) % 3];
                if ((v0 == e.v1 && v1 == e.v2) || (v0 == e.v2 && v1 == e.v1)) {
                    va = v0; vb = v1; break;
                }
            }
        }
        if (va < 0) return;
        path.edges.push_back(ei);
        path.forward.push_back((e.v1 == va && e.v2 == vb) ? 1 : 0);
    };

    std::vector<int> treeFaces;
    int cur = fB;
    while (cur != fA && parentFace[cur] != cur) {
        treeFaces.push_back(cur);
        cur = parentFace[cur];
    }
    std::reverse(treeFaces.begin(), treeFaces.end());

    int prevF = fA;
    for (int nextF : treeFaces) {
        int ei = findSharedEdge(prevF, nextF);
        if (ei >= 0) addOrientedEdge(ei, prevF, nextF);
        prevF = nextF;
    }
    addOrientedEdge(cotreeEdgeIdx, fB, fA);
}

void QuadCover::buildCutGraphAndPaths() {
    cutPaths.clear();
    if (nFaces < 1) return;

    const int nb = countBoundaryComponents();
    const int chi = nVerts - nEdges + nFaces;
    int genus = (2 - chi - nb) / 2;
    if (genus < 0) genus = 0;

    std::vector<std::vector<std::pair<int, int>>> faceAdj(nFaces);
    for (int ei = 0; ei < nEdges; ++ei) {
        const Edge& e = edgeList[ei];
        if (e.f1 < 0 || e.f2 < 0) continue;
        faceAdj[e.f1].emplace_back(e.f2, ei);
        faceAdj[e.f2].emplace_back(e.f1, ei);
    }

    std::vector<int> parentFace(nFaces, -1);
    std::vector<char> inTree(nEdges, 0);
    std::queue<int> q;
    parentFace[0] = 0;
    q.push(0);
    while (!q.empty()) {
        int f = q.front(); q.pop();
        for (const auto& nb : faceAdj[f]) {
            if (parentFace[nb.first] >= 0) continue;
            parentFace[nb.first] = f;
            inTree[nb.second] = 1;
            q.push(nb.first);
        }
    }

    for (int ei = 0; ei < nEdges; ++ei) {
        if (inTree[ei] || edgeList[ei].f2 < 0) continue;
        CutPath path;
        buildDualCyclePath(ei, parentFace, path);
        if (!path.edges.empty()) cutPaths.push_back(path);
    }

    if (nb > 1) {
        std::vector<char> isBnd(nVerts, 0);
        for (int ei = 0; ei < nEdges; ++ei) {
            const Edge& e = edgeList[ei];
            if (e.f2 < 0) { isBnd[e.v1] = 1; isBnd[e.v2] = 1; }
        }
        std::vector<int> bndVerts;
        for (int v = 0; v < nVerts; ++v)
            if (isBnd[v]) bndVerts.push_back(v);

        std::vector<std::vector<int>> vAdj(nVerts);
        for (int ei = 0; ei < nEdges; ++ei) {
            const Edge& e = edgeList[ei];
            vAdj[e.v1].push_back(e.v2);
            vAdj[e.v2].push_back(e.v1);
        }

        std::vector<char> seen(nVerts, 0);
        int start = bndVerts.empty() ? 0 : bndVerts[0];
        for (int target : bndVerts) {
            if (target == start || seen[target]) continue;
            std::vector<int> parent(nVerts, -1);
            std::vector<int> parentEdge(nVerts, -1);
            std::queue<int> bfs;
            bfs.push(start);
            seen[start] = 1;
            parent[start] = start;
            bool found = false;
            while (!bfs.empty() && !found) {
                int v = bfs.front(); bfs.pop();
                for (int nb : vAdj[v]) {
                    if (seen[nb]) continue;
                    seen[nb] = 1;
                    parent[nb] = v;
                    for (int ei = 0; ei < nEdges; ++ei) {
                        const Edge& e = edgeList[ei];
                        if ((e.v1 == v && e.v2 == nb) || (e.v1 == nb && e.v2 == v)) {
                            parentEdge[nb] = ei; break;
                        }
                    }
                    if (nb == target) { found = true; break; }
                    bfs.push(nb);
                }
            }
            if (!found) continue;

            CutPath path;
            int cur = target;
            while (cur != start && parent[cur] != cur) {
                int ei = parentEdge[cur];
                if (ei >= 0) {
                    const Edge& e = edgeList[ei];
                    path.edges.push_back(ei);
                    path.forward.push_back(e.v1 == parent[cur] ? 1 : 0);
                }
                cur = parent[cur];
            }
            std::reverse(path.edges.begin(), path.edges.end());
            std::reverse(path.forward.begin(), path.forward.end());
            if (!path.edges.empty()) cutPaths.push_back(path);
            start = target;
        }
    }

    if (cutPaths.empty()) {
        CutPath dummy;
        for (int ei = 0; ei < nEdges; ++ei) {
            if (edgeList[ei].f2 < 0) continue;
            dummy.edges.push_back(ei);
            dummy.forward.push_back(1);
            break;
        }
        if (!dummy.edges.empty()) cutPaths.push_back(dummy);
    }
}

double QuadCover::computePathPeriod(const CoverData& cover,
                                    const CutPath& path,
                                    int component) const {
    double period = 0.0;
    for (size_t i = 0; i < path.edges.size(); ++i) {
        const Edge& e = edgeList[path.edges[i]];
        int va = path.forward[i] ? e.v1 : e.v2;
        int vb = path.forward[i] ? e.v2 : e.v1;
        int cva = coverVertexFor(cover, va);
        int cvb = coverVertexFor(cover, vb);
        if (cva < 0 || cvb < 0) continue;
        period += cover.coverUV(cvb, component) - cover.coverUV(cva, component);
    }
    return period;
}

void QuadCover::computeMuCoefficients(const CoverData& cover,
                                      VectorXd& muU,
                                      VectorXd& muV) const {
    const int np = (int)cutPaths.size();
    muU.resize(np);
    muV.resize(np);
    for (int j = 0; j < np; ++j) {
        muU(j) = computePathPeriod(cover, cutPaths[j], 0);
        muV(j) = computePathPeriod(cover, cutPaths[j], 1);
    }
}

double QuadCover::roundMu(double value) const {
    const double step = requirePureQuads ? 2.0 : 1.0;
    return step * std::round(value / step);
}

bool QuadCover::buildHarmonicBasis(const CoverData& cover,
                                   const SparseMatrix<double>& L,
                                   std::vector<VectorXd>& basisU,
                                   std::vector<VectorXd>& basisV) const {
    const int nV = cover.nCoverVerts;
    const int np = (int)cutPaths.size();
    if (nV < 1 || np < 1) return false;

    MatrixXd P = MatrixXd::Zero(np, nV);
    for (int j = 0; j < np; ++j) {
        for (size_t k = 0; k < cutPaths[j].edges.size(); ++k) {
            const Edge& e = edgeList[cutPaths[j].edges[k]];
            int va = cutPaths[j].forward[k] ? e.v1 : e.v2;
            int vb = cutPaths[j].forward[k] ? e.v2 : e.v1;
            int cva = coverVertexFor(cover, va);
            int cvb = coverVertexFor(cover, vb);
            if (cva < 0 || cvb < 0) continue;
            P(j, cvb) += 1.0;
            P(j, cva) -= 1.0;
        }
    }

    std::vector<std::vector<int>> adj;
    SparseMatrix<double> Lcopy = L;
    buildCoverLaplacian(cover, Lcopy, adj);

    std::vector<char> seen(nV, 0);
    std::vector<int> anchors;
    for (int seed = 0; seed < nV; ++seed) {
        if (seen[seed]) continue;
        anchors.push_back(seed);
        std::queue<int> bfs;
        bfs.push(seed);
        seen[seed] = 1;
        while (!bfs.empty()) {
            int v = bfs.front(); bfs.pop();
            for (int nb : adj[v]) {
                if (!seen[nb]) { seen[nb] = 1; bfs.push(nb); }
            }
        }
    }

    basisU.assign(np, VectorXd::Zero(nV));
    basisV.assign(np, VectorXd::Zero(nV));

    const int nA = (int)anchors.size();
    const int sysSize = nV + np + nA;
    for (int j = 0; j < np; ++j) {
        std::vector<Triplet<double>> trips;
        trips.reserve(L.nonZeros() + np + nA);

        for (int k = 0; k < L.outerSize(); ++k) {
            for (SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
                trips.emplace_back(it.row(), it.col(), it.value());
        }
        for (int r = 0; r < np; ++r) {
            for (int c = 0; c < nV; ++c) {
                if (std::fabs(P(r, c)) > 1e-12) {
                    trips.emplace_back(nV + r, c, P(r, c));
                    trips.emplace_back(c, nV + r, P(r, c));
                }
            }
        }
        for (int a = 0; a < nA; ++a)
            trips.emplace_back(nV + np + a, anchors[a], 1.0);

        SparseMatrix<double> KKT(sysSize, sysSize);
        KKT.setFromTriplets(trips.begin(), trips.end());

        VectorXd rhs = VectorXd::Zero(sysSize);
        rhs(nV + j) = 1.0;

        SimplicialLDLT<SparseMatrix<double>> solver;
        solver.compute(KKT);
        if (solver.info() != Success) return false;

        VectorXd sol = solver.solve(rhs);
        if (solver.info() != Success) return false;

        basisU[j] = sol.head(nV);
        basisV[j] = sol.head(nV);
    }
    return true;
}

bool QuadCover::enforceGlobalContinuity(CoverData& cover) {
    buildCutGraphAndPaths();
    if (cutPaths.empty()) return true;

    computeMuCoefficients(cover, muU, muV);

    SparseMatrix<double> L;
    std::vector<std::vector<int>> adj;
    buildCoverLaplacian(cover, L, adj);

    std::vector<VectorXd> basisU, basisV;
    if (!buildHarmonicBasis(cover, L, basisU, basisV)) return false;

    const int np = (int)cutPaths.size();
    VectorXd deltaMuU(np), deltaMuV(np);
    for (int j = 0; j < np; ++j) {
        deltaMuU(j) = roundMu(muU(j)) - muU(j);
        deltaMuV(j) = roundMu(muV(j)) - muV(j);
    }

    VectorXd psiU = VectorXd::Zero(cover.nCoverVerts);
    VectorXd psiV = VectorXd::Zero(cover.nCoverVerts);
    for (int j = 0; j < np; ++j) {
        psiU += deltaMuU(j) * basisU[j];
        psiV += deltaMuV(j) * basisV[j];
    }

    cover.coverUV.col(0) += psiU;
    cover.coverUV.col(1) += psiV;

    computeMuCoefficients(cover, muU, muV);
    return true;
}

void QuadCover::projectCoverUVToMesh(const CoverData& cover) {
    UV.resize(nVerts, 2);
    for (int v = 0; v < nVerts; ++v) {
        const bool branchVertex = std::abs(layerShift(v)) > 1e-6;
        if (branchVertex && !cover.baseToCover[v].empty()) {
            RowVector2d avg(0.0, 0.0);
            for (int ci : cover.baseToCover[v]) avg += cover.coverUV.row(ci);
            UV.row(v) = avg / (double)cover.baseToCover[v].size();
        } else {
            int ci = cover.baseSheet0[v];
            if (ci < 0 && !cover.baseToCover[v].empty()) ci = cover.baseToCover[v][0];
            if (ci < 0) { UV.row(v).setZero(); continue; }
            UV.row(v) = cover.coverUV.row(ci);
        }
    }

    double uMin = UV.col(0).minCoeff(), uMax = UV.col(0).maxCoeff();
    double vMin = UV.col(1).minCoeff(), vMax = UV.col(1).maxCoeff();
    const double scale = std::max(uMax - uMin, vMax - vMin);
    if (scale > 1e-12) {
        UV.col(0) = (UV.col(0).array() - uMin) / scale;
        UV.col(1) = (UV.col(1).array() - vMin) / scale;
    }
}

void QuadCover::writeUVToMesh() {
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v)
        v->uv = Vector2d(UV(v->index, 0), UV(v->index, 1));
}

bool QuadCover::runFullPipeline() {
    // Input: triangle mesh (constructor) + cross field (setCrossField or principal curvature).
    computeVertexNormals();
    estimatePrincipalCurvature();  // uses externalFaceDirs when setCrossField() was called
    computeFaceTheta();

    computeMatching();
    if (smoothIters > 0) smoothCrossField(smoothIters);
    syncFaceDirsFromTheta();
    computeMatching();

    computeLayerShift();

    CoverData cover;
    if (!buildBranchCover(cover)) {
        UV = MatrixXd::Zero(nVerts, 2);
        return false;
    }
    if (!integrateOnCover(cover)) {
        UV = MatrixXd::Zero(nVerts, 2);
        return false;
    }
    if (!enforceGlobalContinuity(cover)) {
        UV = MatrixXd::Zero(nVerts, 2);
        return false;
    }

    projectCoverUVToMesh(cover);
    return true;
}

void QuadCover::parameterize() {
    if (!initMeshData()) return;
    if (!runFullPipeline()) {
        UV = MatrixXd::Zero(nVerts, 2);
    }
    writeUVToMesh();
}

bool QuadCover::parameterizeFull() {
    if (!initMeshData()) return false;
    const bool ok = runFullPipeline();
    writeUVToMesh();
    return ok;
}

void QuadCover::computeTransportMatching() {
    computeMatching();
}

void QuadCover::buildCoveringField() {
    faceD1.resize(nFaces, 3);
    faceD2.resize(nFaces, 3);
    for (int fi = 0; fi < nFaces; ++fi) {
        Vector3d n = faceNormals.row(fi);
        Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);
        Vector3d d = faceDirs_.row(fi);
        double dx = d.dot(t1), dy = d.dot(t2);
        faceD1.row(fi) = (dx * t1 + dy * t2).normalized();
        faceD2.row(fi) = (-dy * t1 + dx * t2).normalized();
    }
}

void QuadCover::solvePoisson() {
    buildCoveringField();
    MatrixXd Vd1(nVerts, 3), Vd2(nVerts, 3);
    Vd1.setZero(); Vd2.setZero();
    VectorXd vertArea = VectorXd::Zero(nVerts);

    for (int fi = 0; fi < nFaces; ++fi) {
        int v0 = faces[fi * 3], v1 = faces[fi * 3 + 1], v2 = faces[fi * 3 + 2];
        Vector3d p0 = vertPos.row(v0), p1 = vertPos.row(v1), p2 = vertPos.row(v2);
        double area = 0.5 * (p1 - p0).cross(p2 - p0).norm();
        for (int k = 0; k < 3; ++k) {
            int v = faces[fi * 3 + k];
            vertArea(v) += area;
            Vd1.row(v) += area * faceD1.row(fi);
            Vd2.row(v) += area * faceD2.row(fi);
        }
    }
    for (int i = 0; i < nVerts; ++i) {
        if (vertArea(i) > 1e-12) {
            Vd1.row(i) /= vertArea(i);
            Vd2.row(i) /= vertArea(i);
        }
    }

    VectorXd div1 = VectorXd::Zero(nVerts);
    VectorXd div2 = VectorXd::Zero(nVerts);
    for (int fi = 0; fi < nFaces; ++fi) {
        int v0 = faces[fi * 3], v1 = faces[fi * 3 + 1], v2 = faces[fi * 3 + 2];
        Vector3d p0 = vertPos.row(v0), p1 = vertPos.row(v1), p2 = vertPos.row(v2);
        Vector3d fn = (p1 - p0).cross(p2 - p0);
        Vector3d e01 = p1 - p0, e12 = p2 - p1, e20 = p0 - p2;
        Vector3d en01 = fn.cross(e01).normalized();
        Vector3d en12 = fn.cross(e12).normalized();
        Vector3d en20 = fn.cross(e20).normalized();
        double cot0 = cotan(p2, p0, p1);
        double cot1 = cotan(p0, p1, p2);
        double cot2 = cotan(p1, p2, p0);
        Vector3d fd1 = faceD1.row(fi), fd2 = faceD2.row(fi);
        double d1 = 0.5 * (fd1.dot(en01) * cot2 + fd1.dot(en12) * cot0 + fd1.dot(en20) * cot1);
        double d2 = 0.5 * (fd2.dot(en01) * cot2 + fd2.dot(en12) * cot0 + fd2.dot(en20) * cot1);
        div1(v0) += d1; div1(v1) += d1; div1(v2) += d1;
        div2(v0) += d2; div2(v1) += d2; div2(v2) += d2;
    }

    SparseMatrix<double> L;
    buildCotLaplacian(L);

    std::vector<int> bndVerts;
    for (int i = 0; i < nVerts; ++i) {
        for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
            if ((int)v->index == i && v->isBoundary()) {
                bndVerts.push_back(i);
                break;
            }
        }
    }
    if (bndVerts.empty()) bndVerts.push_back(0);

    SparseMatrix<double> Lm = L;
    VectorXd b1(div1), b2(div2);
    for (int bi : bndVerts) {
        Lm.coeffRef(bi, bi) += 1e6;
        b1(bi) = 0;
        b2(bi) = 0;
    }

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(Lm);
    if (solver.info() != Success) { UV = MatrixXd::Zero(nVerts, 2); return; }

    VectorXd u = solver.solve(b1);
    VectorXd v = solver.solve(b2);

    double uMin = u.minCoeff(), uMax = u.maxCoeff();
    double vMin = v.minCoeff(), vMax = v.maxCoeff();
    double uRng = uMax - uMin, vRng = vMax - vMin;
    if (uRng < 1e-10) uRng = 1.0;
    if (vRng < 1e-10) vRng = 1.0;

    UV.resize(nVerts, 2);
    for (int i = 0; i < nVerts; ++i) {
        UV(i, 0) = (u(i) - uMin) / uRng;
        UV(i, 1) = (v(i) - vMin) / vRng;
    }
}

