#include "QuadCover.h"
#include <map>
#include <algorithm>
#include <queue>
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
}

QuadCover::QuadCover(Mesh& mesh0): Parameterization(mesh0), hasExternalFaceDirs(false) {}

void QuadCover::setFaceDirections(const Eigen::MatrixXd& dirs) {
    externalFaceDirs = dirs;
    hasExternalFaceDirs = (dirs.cols() == 3 && dirs.rows() > 0);
}

// ============================================================
// Helpers
// ============================================================

void QuadCover::buildLocalFrame(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2) {
    Eigen::Vector3d ref(1,0,0);
    if (fabs(n.dot(ref)) > 0.9) ref = Eigen::Vector3d(0,1,0);
    t1 = (ref - ref.dot(n)*n).normalized();
    t2 = n.cross(t1).normalized();
}

double QuadCover::cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c) {
    // cot(angle at b) = ((a-b)·(c-b)) / |(a-b)×(c-b)|
    Eigen::Vector3d u = a - b, v = c - b;
    double d = u.dot(v), cr = u.cross(v).norm();
    if (cr < 1e-12) return 0.0;
    return d / cr;
}

void QuadCover::buildCotLaplacian(Eigen::SparseMatrix<double>& L) {
    L.resize(nVerts, nVerts);
    std::vector<Eigen::Triplet<double>> trips;
    Eigen::VectorXd diag = Eigen::VectorXd::Zero(nVerts);
    
    for (auto& e : edgeList) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        int vi = e.v1, vj = e.v2;
        
        // Find the two vertices opposite to the edge
        int o1 = -1, o2 = -1;
        int f1v0 = faces[e.f1*3], f1v1 = faces[e.f1*3+1], f1v2 = faces[e.f1*3+2];
        int f2v0 = faces[e.f2*3], f2v1 = faces[e.f2*3+1], f2v2 = faces[e.f2*3+2];
        
        if (f1v0 != vi && f1v0 != vj) o1 = f1v0;
        else if (f1v1 != vi && f1v1 != vj) o1 = f1v1;
        else o1 = f1v2;
        
        if (f2v0 != vi && f2v0 != vj) o2 = f2v0;
        else if (f2v1 != vi && f2v1 != vj) o2 = f2v1;
        else o2 = f2v2;
        
        double w = 0.5 * (cotan(vertPos.row(o1), vertPos.row(vi), vertPos.row(vj)) +
                          cotan(vertPos.row(o2), vertPos.row(vi), vertPos.row(vj)));
        if (w < 0) w = 0;  // ensure non-negative for Delaunay meshes
        
        diag(vi) += w; diag(vj) += w;
        trips.push_back(Eigen::Triplet<double>(vi, vj, -w));
        trips.push_back(Eigen::Triplet<double>(vj, vi, -w));
    }
    for (int i = 0; i < nVerts; i++)
        trips.push_back(Eigen::Triplet<double>(i, i, diag(i) + 1e-8));
    L.setFromTriplets(trips.begin(), trips.end());
}

// ============================================================
// Phase 0: Vertex normals (area-weighted face normals)
// ============================================================

void QuadCover::computeVertexNormals() {
    vertexNormals.resize(nVerts, 3);
    vertexNormals.setZero();
    
    for (int fi = 0; fi < nFaces; fi++) {
        int v0 = faces[fi*3], v1 = faces[fi*3+1], v2 = faces[fi*3+2];
        Eigen::Vector3d p0 = vertPos.row(v0), p1 = vertPos.row(v1), p2 = vertPos.row(v2);
        
        Eigen::Vector3d fn = (p1 - p0).cross(p2 - p0);
        double area = 0.5 * fn.norm();
        fn.normalize();
        
        // Accumulate area-weighted face normal to each vertex
        for (int k = 0; k < 3; k++) {
            int v = faces[fi*3+k];
            vertexNormals.row(v) += area * fn;
        }
    }
    
    // Normalize
    for (int i = 0; i < nVerts; i++) {
        double len = vertexNormals.row(i).norm();
        if (len > 1e-12) vertexNormals.row(i) /= len;
    }
}

// ============================================================
// Phase 0: Weingarten map W = I⁻¹·II → principal curvature direction
// ============================================================

void QuadCover::estimatePrincipalCurvature() {
    faceNormals.resize(nFaces, 3);
    faceDirs.resize(nFaces, 3);
    
    for (int fi = 0; fi < nFaces; fi++) {
        int v0 = faces[fi*3], v1 = faces[fi*3+1], v2 = faces[fi*3+2];
        Eigen::Vector3d p0 = vertPos.row(v0), p1 = vertPos.row(v1), p2 = vertPos.row(v2);
        
        // 1. Face normal
        Eigen::Vector3d fn = (p1 - p0).cross(p2 - p0);
        double area = 0.5 * fn.norm();
        fn /= std::max(2.0 * area, 1e-12);
        faceNormals.row(fi) = fn;

        if (hasExternalFaceDirs && externalFaceDirs.rows() == nFaces) {
            Eigen::Vector3d dir = externalFaceDirs.row(fi);
            dir -= dir.dot(fn) * fn;
            double dlen = dir.norm();
            if (dlen > 1e-12) {
                faceDirs.row(fi) = dir / dlen;
            } else {
                Eigen::Vector3d t1, t2;
                buildLocalFrame(fn, t1, t2);
                faceDirs.row(fi) = t1;
            }
            continue;
        }
        
        // 2. Tangent plane local frame {t1, t2}
        Eigen::Vector3d t1, t2;
        buildLocalFrame(fn, t1, t2);
        
        // 3. Vertex normals (area-weighted, precomputed)
        Eigen::Vector3d n0 = vertexNormals.row(v0);
        Eigen::Vector3d n1 = vertexNormals.row(v1);
        Eigen::Vector3d n2 = vertexNormals.row(v2);
        
        // 4. Local 2D coordinates of edges in tangent plane
        Eigen::Vector3d e01 = p1 - p0, e02 = p2 - p0;
        double u1 = e01.dot(t1), v1_coord = e01.dot(t2);
        double u2 = e02.dot(t1), v2_coord = e02.dot(t2);
        
        // 5. First fundamental form I (constant for the face geometry)
        Matrix2d I_mat;
        I_mat << u1*u1 + v1_coord*v1_coord,  u1*u2 + v1_coord*v2_coord,
                 u1*u2 + v1_coord*v2_coord,  u2*u2 + v2_coord*v2_coord;
        
        // 6. Second fundamental form II via normal variation
        // Δn·t₁ = (n₁-n₀)·t₁,  Δn·t₂ = (n₂-n₀)·t₂  projected onto tangent plane
        Eigen::Vector3d dn1 = n1 - n0;
        Eigen::Vector3d dn2 = n2 - n0;
        
        // II = [e f; f g]  where:
        //   e = -∂²p/∂u² · n  approximated by finite differences
        //   f = -∂²p/∂u∂v · n
        //   g = -∂²p/∂v² · n
        //
        // Using: dn = -II · du  (neglecting higher-order terms)
        // Solve 2×2 system:  [u1 v1; u2 v2]ᵀ [e;f] ≈ [-dn1·t1; -dn2·t1]
        // but for symmetric II we estimate e,f,g via least squares
        
        // Build RHS: normal variation dotted with tangent basis
        Eigen::Vector2d rhs1, rhs2;
        rhs1 << -dn1.dot(t1), -dn2.dot(t1);  // for e, f (first row of II)
        rhs2 << -dn1.dot(t2), -dn2.dot(t2);  // for f, g (second row of II)
        
        // Solve least squares: 2 equations for 2 unknowns each row
        // M = [u1 v1; u2 v2],  solve Mᵀ [e;f] = rhs1,  Mᵀ [f;g] = rhs2
        Matrix2d M;
        M << u1, v1_coord,
             u2, v2_coord;
        
        Matrix2d Mt = M.transpose();
        
        // Solve Mt * x = rhs → normal equations: M*Mt*x = M*rhs
        Eigen::Vector2d ef = (Mt * M).ldlt().solve(Mt * rhs1);
        Eigen::Vector2d fg = (Mt * M).ldlt().solve(Mt * rhs2);
        
        double e_coef = ef(0), f_coef = ef(1);
        double fg_f = fg(0), g_coef = fg(1);
        
        // Symmetrize f = (f_coef + fg_f) / 2
        double f_sym = 0.5 * (f_coef + fg_f);
        
        // 7. Shape operator (Weingarten map) W = I⁻¹ · II
        Matrix2d II_mat;
        II_mat << e_coef, f_sym,
                  f_sym,  g_coef;
        
        Matrix2d I_inv = I_mat.inverse();
        Matrix2d W = I_inv * II_mat;
        
        // Symmetrize W = (W + Wᵀ)/2 for real eigenvalues
        Matrix2d W_sym = 0.5 * (W + W.transpose());
        
        // 8. Eigen-decomposition → principal curvature directions
        Eigen::SelfAdjointEigenSolver<Matrix2d> es(W_sym);
        
        // Use direction of max |eigenvalue| (maximal curvature magnitude)
        Eigen::Vector2d dir2d;
        if (fabs(es.eigenvalues()(0)) >= fabs(es.eigenvalues()(1)))
            dir2d = es.eigenvectors().col(0);
        else
            dir2d = es.eigenvectors().col(1);
        
        // 9. Map back to 3D
        Eigen::Vector3d dir3d = dir2d(0) * t1 + dir2d(1) * t2;
        double dlen = dir3d.norm();
        if (dlen > 1e-12) dir3d /= dlen;
        else dir3d = t1;  // fallback for degenerate faces
        
        faceDirs.row(fi) = dir3d;
    }
}

// ============================================================
// Phase 1: Matching
// ============================================================

void QuadCover::computeMatching() {
    matching.resize(nEdges);
    matching.setZero();
    
    for (int eIdx = 0; eIdx < nEdges; eIdx++) {
        auto& e = edgeList[eIdx];
        if (e.f1 < 0 || e.f2 < 0) { matching[eIdx] = 0; continue; }
        
        Eigen::Vector3d n1 = faceNormals.row(e.f1);
        Eigen::Vector3d n2 = faceNormals.row(e.f2);
        
        Eigen::Vector3d t11, t12, t21, t22;
        buildLocalFrame(n1, t11, t12);
        buildLocalFrame(n2, t21, t22);
        
        Eigen::Vector3d d1 = faceDirs.row(e.f1);
        Eigen::Vector3d d2 = faceDirs.row(e.f2);
        
        // Convert to 2D local coordinates
        Eigen::Vector2d d1_2d(d1.dot(t11), d1.dot(t12));
        Eigen::Vector2d d2_2d(d2.dot(t21), d2.dot(t22));
        d1_2d.normalize(); d2_2d.normalize();
        
        // Try 4 rotations, pick min angle
        double best = 1e10;
        int bestK = 0;
        for (int k = 0; k < 4; k++) {
            double cs = cos(k * M_PI_2), sn = sin(k * M_PI_2);
            Eigen::Vector2d rot(cs*d2_2d.x() - sn*d2_2d.y(),
                                sn*d2_2d.x() + cs*d2_2d.y());
            double dot = d1_2d.dot(rot);
            dot = std::max(-1.0, std::min(1.0, dot));
            double ang = acos(fabs(dot));  // 4-RoSy symmetry uses abs
            if (ang < best) { best = ang; bestK = k; }
        }
        matching[eIdx] = bestK;
    }
}

// ============================================================
// Phase 2: Layer shift
// ============================================================

void QuadCover::computeLayerShift() {
    layerShift = Eigen::VectorXd::Zero(nVerts);

    std::vector<std::vector<int>> incidentFaces(nVerts);
    std::vector<char> boundaryVertex(nVerts, 0);
    for (int fi = 0; fi < nFaces; ++fi) {
        for (int k = 0; k < 3; ++k) incidentFaces[faces[fi * 3 + k]].push_back(fi);
    }
    for (int ei = 0; ei < nEdges; ++ei) {
        const Edge& e = edgeList[ei];
        if (e.f1 < 0 || e.f2 < 0) {
            boundaryVertex[e.v1] = 1;
            boundaryVertex[e.v2] = 1;
        }
    }

    auto signedMatch = [&](int fromFace, int toFace) {
        for (int ei = 0; ei < nEdges; ++ei) {
            const Edge& e = edgeList[ei];
            if (e.f1 == fromFace && e.f2 == toFace) return matching[ei];
            if (e.f2 == fromFace && e.f1 == toFace) return -matching[ei];
        }
        return 0;
    };

    for (int v = 0; v < nVerts; ++v) {
        if (boundaryVertex[v] || incidentFaces[v].size() < 3) continue;

        Eigen::Vector3d n = vertexNormals.row(v);
        if (n.norm() < 1e-12) continue;
        n.normalize();
        Eigen::Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);
        const Eigen::Vector3d pv = vertPos.row(v);

        std::vector<std::pair<double, int>> ordered;
        for (int fi : incidentFaces[v]) {
            Eigen::Vector3d c(0, 0, 0);
            for (int k = 0; k < 3; ++k) c += vertPos.row(faces[fi * 3 + k]);
            c = c / 3.0 - pv;
            c -= c.dot(n) * n;
            if (c.norm() < 1e-12) continue;
            ordered.push_back(std::make_pair(std::atan2(c.dot(t2), c.dot(t1)), fi));
        }
        if (ordered.size() < 3) continue;
        std::sort(ordered.begin(), ordered.end());

        int holonomy = 0;
        for (size_t i = 0; i < ordered.size(); ++i) {
            int f0 = ordered[i].second;
            int f1 = ordered[(i + 1) % ordered.size()].second;
            holonomy += signedMatch(f0, f1);
        }
        int h = mod4(holonomy);
        if (h > 2) h -= 4;
        layerShift(v) = (double)h / 4.0;
    }
}

// ============================================================
// Phase 3: Covering vector field
// ============================================================

void QuadCover::buildCoveringField() {
    faceD1.resize(nFaces, 3);
    faceD2.resize(nFaces, 3);
    
    for (int fi = 0; fi < nFaces; fi++) {
        Eigen::Vector3d n = faceNormals.row(fi);
        Eigen::Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);
        
        Eigen::Vector3d d = faceDirs.row(fi);
        double dx = d.dot(t1), dy = d.dot(t2);
        
        // Two orthogonal directions
        faceD1.row(fi) = ( dx * t1 + dy * t2).normalized();
        faceD2.row(fi) = (-dy * t1 + dx * t2).normalized();
    }
}

// ============================================================
// Phase 4: Poisson solve (Hodge decomposition)
// ============================================================

void QuadCover::solvePoisson() {
    // 1. Interpolate face directions to vertices (area-weighted)
    Eigen::MatrixXd Vd1(nVerts, 3), Vd2(nVerts, 3);
    Vd1.setZero(); Vd2.setZero();
    Eigen::VectorXd vertArea = Eigen::VectorXd::Zero(nVerts);
    
    for (int fi = 0; fi < nFaces; fi++) {
        int v0 = faces[fi*3], v1 = faces[fi*3+1], v2 = faces[fi*3+2];
        Eigen::Vector3d p0 = vertPos.row(v0), p1 = vertPos.row(v1), p2 = vertPos.row(v2);
        double area = 0.5 * (p1-p0).cross(p2-p0).norm();
        
        for (int k = 0; k < 3; k++) {
            int v = faces[fi*3+k];
            vertArea(v) += area;
            Vd1.row(v) += area * faceD1.row(fi);
            Vd2.row(v) += area * faceD2.row(fi);
        }
    }
    for (int i = 0; i < nVerts; i++) {
        if (vertArea(i) > 1e-12) {
            Vd1.row(i) /= vertArea(i);
            Vd2.row(i) /= vertArea(i);
        }
    }
    
    // 2. Compute divergence: div = G^T * M_faces * d (integrated discrete divergence)
    Eigen::VectorXd div1 = Eigen::VectorXd::Zero(nVerts);
    Eigen::VectorXd div2 = Eigen::VectorXd::Zero(nVerts);
    
    for (int fi = 0; fi < nFaces; fi++) {
        int v0 = faces[fi*3], v1 = faces[fi*3+1], v2 = faces[fi*3+2];
        Eigen::Vector3d p0 = vertPos.row(v0), p1 = vertPos.row(v1), p2 = vertPos.row(v2);
        Eigen::Vector3d fn = (p1-p0).cross(p2-p0);
        
        // Edge vectors
        Eigen::Vector3d e01 = p1-p0, e12 = p2-p1, e20 = p0-p2;
        
        // Edge normals in face plane (90° rotated inward)
        Eigen::Vector3d en01 = fn.cross(e01).normalized();
        Eigen::Vector3d en12 = fn.cross(e12).normalized();
        Eigen::Vector3d en20 = fn.cross(e20).normalized();
        
        // Cotan weights
        double cot0 = cotan(p2, p0, p1);  // opposite v0, at edge v1-v2
        double cot1 = cotan(p0, p1, p2);  // opposite v1, at edge v2-v0
        double cot2 = cotan(p1, p2, p0);  // opposite v2, at edge v0-v1
        
        // Face average direction
        Eigen::Vector3d fd1 = faceD1.row(fi), fd2 = faceD2.row(fi);
        
        // Discrete divergence: sum over edges of (d·edge_normal)*cot/2
        double d1 = 0.5 * (fd1.dot(en01)*cot2 + fd1.dot(en12)*cot0 + fd1.dot(en20)*cot1);
        double d2 = 0.5 * (fd2.dot(en01)*cot2 + fd2.dot(en12)*cot0 + fd2.dot(en20)*cot1);
        
        div1(v0) += d1; div1(v1) += d1; div1(v2) += d1;
        div2(v0) += d2; div2(v1) += d2; div2(v2) += d2;
    }
    
    // 3. Build cotangent Laplacian
    Eigen::SparseMatrix<double> L;
    buildCotLaplacian(L);
    
    // 4. Fix boundary: pick first boundary vertex as anchor
    // Find boundary vertices
    std::vector<int> bndVerts;
    for (int i = 0; i < nVerts; i++) {
        for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
            if ((int)v->index == i && v->isBoundary()) {
                bndVerts.push_back(i);
                break;
            }
        }
    }
    
    // If closed mesh, fix vertex 0
    bool closed = bndVerts.empty();
    if (closed) bndVerts.push_back(0);
    
    // Fix first boundary vertex to (0,0)
    int anchor = bndVerts[0];
    
    // 5. Solve constrained linear system: eliminate anchor row/col
    // We modify the RHS: for each anchored vertex, set RHS = 0 and
    // add penalty to diagonal
    Eigen::SparseMatrix<double> Lmod(L);
    Eigen::VectorXd b1(div1), b2(div2);
    
    // Zero out row/col for anchor and set to identity
    for (int k = 0; k < Lmod.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Lmod, k); it; ++it) {
            if (it.row() == anchor || it.col() == anchor) {
                // Will be handled by the reset below
            }
        }
    }
    
    // Better approach: use penalty method
    Eigen::SparseMatrix<double> Lm = L;
    for (int bi : bndVerts) {
        Lm.coeffRef(bi, bi) += 1e6;
        b1(bi) = 0;
        b2(bi) = 0;
    }
    
    // Solve: L*u = b
    // Remove zero eigenvalues by shifting slightly
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(Lm);
    
    if (solver.info() != Eigen::Success) {
        UV = Eigen::MatrixXd::Zero(nVerts, 2);
        return;
    }
    
    Eigen::VectorXd u = solver.solve(b1);
    Eigen::VectorXd v = solver.solve(b2);
    
    // 6. Normalize
    double uMin = u.minCoeff(), uMax = u.maxCoeff();
    double vMin = v.minCoeff(), vMax = v.maxCoeff();
    double uRng = uMax - uMin, vRng = vMax - vMin;
    if (uRng < 1e-10) uRng = 1.0;
    if (vRng < 1e-10) vRng = 1.0;
    
    UV.resize(nVerts, 2);
    for (int i = 0; i < nVerts; i++) {
        UV(i, 0) = (u(i) - uMin) / uRng;
        UV(i, 1) = (v(i) - vMin) / vRng;
    }
}

// ============================================================
// Full QuadCover: matching via shared-edge parallel transport
// ============================================================

void QuadCover::computeTransportMatching() {
    matching.resize(nEdges);
    matching.setZero();

    for (int eIdx = 0; eIdx < nEdges; eIdx++) {
        Edge& e = edgeList[eIdx];
        if (e.f1 < 0 || e.f2 < 0) {
            matching[eIdx] = 0;
            continue;
        }

        Eigen::Vector3d edgeDir = vertPos.row(e.v2) - vertPos.row(e.v1);
        double elen = edgeDir.norm();
        if (elen < 1e-12) {
            matching[eIdx] = 0;
            continue;
        }
        edgeDir /= elen;

        auto faceAngle = [&](int fi) {
            Eigen::Vector3d n = faceNormals.row(fi).normalized();
            Eigen::Vector3d x = edgeDir - edgeDir.dot(n) * n;
            double xlen = x.norm();
            if (xlen < 1e-12) {
                Eigen::Vector3d yFallback;
                buildLocalFrame(n, x, yFallback);
            } else {
                x /= xlen;
            }
            Eigen::Vector3d y = n.cross(x).normalized();

            Eigen::Vector3d d = faceDirs.row(fi);
            d -= d.dot(n) * n;
            double dlen = d.norm();
            if (dlen < 1e-12) d = x;
            else d /= dlen;
            return std::atan2(d.dot(y), d.dot(x));
        };

        const double a1 = faceAngle(e.f1);
        const double a2 = faceAngle(e.f2);
        int r = (int)std::floor((a1 - a2) / (0.5 * M_PI) + 0.5);
        matching[eIdx] = mod4(r);
    }
}

bool QuadCover::solveCoverPoisson(const Eigen::MatrixXd& coverPos,
                                  const Eigen::MatrixXi& coverFaces,
                                  const Eigen::MatrixXd& coverD1,
                                  const Eigen::MatrixXd& coverD2,
                                  Eigen::MatrixXd& coverUV) {
    const int nV = (int)coverPos.rows();
    const int nF = (int)coverFaces.rows();
    if (nV < 3 || nF < 1) return false;

    std::vector<Eigen::Triplet<double>> trips;
    Eigen::VectorXd diag = Eigen::VectorXd::Zero(nV);
    Eigen::VectorXd rhsU = Eigen::VectorXd::Zero(nV);
    Eigen::VectorXd rhsV = Eigen::VectorXd::Zero(nV);
    std::vector<std::vector<int>> adj(nV);

    auto addWeight = [&](int i, int j, double w) {
        if (!std::isfinite(w) || w <= 0.0) return;
        diag(i) += w;
        diag(j) += w;
        trips.push_back(Eigen::Triplet<double>(i, j, -w));
        trips.push_back(Eigen::Triplet<double>(j, i, -w));
        adj[i].push_back(j);
        adj[j].push_back(i);
    };

    for (int fi = 0; fi < nF; ++fi) {
        const int i0 = coverFaces(fi, 0), i1 = coverFaces(fi, 1), i2 = coverFaces(fi, 2);
        if (i0 < 0 || i0 >= nV || i1 < 0 || i1 >= nV || i2 < 0 || i2 >= nV) return false;
        if (i0 == i1 || i1 == i2 || i2 == i0) continue;
        const Eigen::Vector3d p0 = coverPos.row(i0);
        const Eigen::Vector3d p1 = coverPos.row(i1);
        const Eigen::Vector3d p2 = coverPos.row(i2);
        Eigen::Vector3d n = (p1 - p0).cross(p2 - p0);
        const double dblArea = n.norm();
        if (dblArea < 1e-14) continue;
        n /= dblArea;
        const double area = 0.5 * dblArea;

        const double c0 = cotan(p1, p0, p2);
        const double c1 = cotan(p2, p1, p0);
        const double c2 = cotan(p0, p2, p1);
        addWeight(i1, i2, 0.5 * c0);
        addWeight(i2, i0, 0.5 * c1);
        addWeight(i0, i1, 0.5 * c2);

        const Eigen::Vector3d grad0 = n.cross(p2 - p1) / dblArea;
        const Eigen::Vector3d grad1 = n.cross(p0 - p2) / dblArea;
        const Eigen::Vector3d grad2 = n.cross(p1 - p0) / dblArea;
        const Eigen::Vector3d xu = coverD1.row(fi);
        const Eigen::Vector3d xv = coverD2.row(fi);

        rhsU(i0) += area * xu.dot(grad0);
        rhsU(i1) += area * xu.dot(grad1);
        rhsU(i2) += area * xu.dot(grad2);
        rhsV(i0) += area * xv.dot(grad0);
        rhsV(i1) += area * xv.dot(grad1);
        rhsV(i2) += area * xv.dot(grad2);
    }

    for (int i = 0; i < nV; ++i) {
        trips.push_back(Eigen::Triplet<double>(i, i, diag(i) + 1e-8));
    }
    Eigen::SparseMatrix<double> L(nV, nV);
    L.setFromTriplets(trips.begin(), trips.end());

    // Anchor one vertex per connected component to remove the constant nullspace.
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
            int v = q.front();
            q.pop();
            for (int nb : adj[v]) {
                if (!seen[nb]) {
                    seen[nb] = 1;
                    q.push(nb);
                }
            }
        }
    }
    L.makeCompressed();

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
                             Eigen::DiagonalPreconditioner<double>> solver;
    solver.setMaxIterations(std::max(200, nV * 2));
    solver.setTolerance(1e-8);
    solver.compute(L);
    if (solver.info() != Eigen::Success) return false;

    Eigen::VectorXd u = solver.solve(rhsU);
    Eigen::VectorXd v = solver.solve(rhsV);
    if (solver.info() != Eigen::Success) return false;

    coverUV.resize(nV, 2);
    coverUV.col(0) = u;
    coverUV.col(1) = v;
    return true;
}

bool QuadCover::buildBranchCoverAndIntegrate() {
    const int sheets = 4;
    const int rawCoverCorners = nFaces * sheets * 3;
    DisjointSet dsu(rawCoverCorners);

    auto cornerId = [&](int fi, int local, int s) {
        return (fi * sheets + mod4(s)) * 3 + local;
    };
    auto localCorner = [&](int fi, int vertex) {
        for (int k = 0; k < 3; ++k) {
            if (faces[fi * 3 + k] == vertex) return k;
        }
        return -1;
    };

    // Glue face-corner copies according to matching. Around a vertex, the
    // accumulated transition is exactly layerShift/holonomy: ordinary vertices
    // keep four independent preimages, while non-zero layerShift vertices
    // become branch points where the sheet cycle does not close on itself.
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
    std::vector<Eigen::Vector3d> cpos;
    std::vector<std::vector<int>> baseToCover(nVerts);
    std::vector<int> baseSheet0(nVerts, -1);

    for (int fi = 0; fi < nFaces; ++fi) {
        for (int s = 0; s < sheets; ++s) {
            for (int k = 0; k < 3; ++k) {
                const int raw = cornerId(fi, k, s);
                const int v = faces[fi * 3 + k];
                const int rep = dsu.find(raw);
                auto it = repToIndex.find(rep);
                if (it == repToIndex.end()) {
                    const int idx = (int)cpos.size();
                    repToIndex[rep] = idx;
                    coverIndex[raw] = idx;
                    cpos.push_back(vertPos.row(v));
                } else {
                    coverIndex[raw] = it->second;
                }
                if (std::find(baseToCover[v].begin(), baseToCover[v].end(), coverIndex[raw]) == baseToCover[v].end()) {
                    baseToCover[v].push_back(coverIndex[raw]);
                }
                if (s == 0 && baseSheet0[v] < 0) baseSheet0[v] = coverIndex[raw];
            }
        }
    }

    Eigen::MatrixXd coverPos((int)cpos.size(), 3);
    for (int i = 0; i < (int)cpos.size(); ++i) coverPos.row(i) = cpos[i];

    Eigen::MatrixXi coverFaces(nFaces * sheets, 3);
    Eigen::MatrixXd coverD1(nFaces * sheets, 3), coverD2(nFaces * sheets, 3);
    for (int fi = 0; fi < nFaces; ++fi) {
        Eigen::Vector3d n = faceNormals.row(fi).normalized();
        Eigen::Vector3d base = faceDirs.row(fi);
        base -= base.dot(n) * n;
        double blen = base.norm();
        if (blen < 1e-12) {
            Eigen::Vector3d t1, t2;
            buildLocalFrame(n, t1, t2);
            base = t1;
        } else {
            base /= blen;
        }
        Eigen::Vector3d ortho = n.cross(base).normalized();

        for (int s = 0; s < sheets; ++s) {
            const int row = fi * sheets + s;
            coverFaces(row, 0) = coverIndex[cornerId(fi, 0, s)];
            coverFaces(row, 1) = coverIndex[cornerId(fi, 1, s)];
            coverFaces(row, 2) = coverIndex[cornerId(fi, 2, s)];

            const double a = s * 0.5 * M_PI;
            Eigen::Vector3d d1 = std::cos(a) * base + std::sin(a) * ortho;
            Eigen::Vector3d d2 = -std::sin(a) * base + std::cos(a) * ortho;
            coverD1.row(row) = d1.normalized();
            coverD2.row(row) = d2.normalized();
        }
    }

    Eigen::MatrixXd coverUV;
    if (!solveCoverPoisson(coverPos, coverFaces, coverD1, coverD2, coverUV)) return false;

    UV.resize(nVerts, 2);
    for (int v = 0; v < nVerts; ++v) {
        const bool branchVertex = std::abs(layerShift(v)) > 1e-6;
        if (branchVertex && !baseToCover[v].empty()) {
            Eigen::RowVector2d avg(0.0, 0.0);
            for (int ci : baseToCover[v]) avg += coverUV.row(ci);
            UV.row(v) = avg / (double)baseToCover[v].size();
        } else {
            int ci = baseSheet0[v];
            if (ci < 0 && !baseToCover[v].empty()) ci = baseToCover[v][0];
            if (ci < 0) return false;
            UV.row(v) = coverUV.row(ci);
        }
    }

    double uMin = UV.col(0).minCoeff(), uMax = UV.col(0).maxCoeff();
    double vMin = UV.col(1).minCoeff(), vMax = UV.col(1).maxCoeff();
    const double scale = std::max(uMax - uMin, vMax - vMin);
    if (scale > 1e-12) {
        UV.col(0) = (UV.col(0).array() - uMin) / scale;
        UV.col(1) = (UV.col(1).array() - vMin) / scale;
    }
    return true;
}

// ============================================================
// Main entry point
// ============================================================

void QuadCover::parameterize() {
    // Extract vertex positions and face indices from mesh
    nVerts = (int)mesh.vertices.size();
    nFaces = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++)
        if (!f->isBoundary()) nFaces++;
    
    vertPos.resize(nVerts, 3);
    faces.resize(nFaces * 3);
    
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        vertPos.row(v->index) = v->position;
    }
    
    int fi = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            faces[fi*3]   = f->he->vertex->index;
            faces[fi*3+1] = f->he->next->vertex->index;
            faces[fi*3+2] = f->he->next->next->vertex->index;
            fi++;
        }
    }
    
    // Build edge list
    std::map<std::pair<int,int>, int> edgeMap;
    edgeList.clear();
    for (int fi2 = 0; fi2 < nFaces; fi2++) {
        for (int k = 0; k < 3; k++) {
            int v1 = faces[fi2*3+k], v2 = faces[fi2*3+(k+1)%3];
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
    
    // Phase 0: Compute vertex normals, then Weingarten map → principal curvature
    computeVertexNormals();
    estimatePrincipalCurvature();
    
    // Phase 1: Matching
    computeMatching();
    
    // Phase 2: Layer shift
    computeLayerShift();
    
    // Phase 3: Covering field
    buildCoveringField();
    
    // Phase 4: Poisson solve
    solvePoisson();
    
    // Copy UV back to mesh vertices
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->uv = Eigen::Vector2d(UV(v->index, 0), UV(v->index, 1));
    }
}

bool QuadCover::parameterizeFull() {
    nVerts = (int)mesh.vertices.size();
    nFaces = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++)
        if (!f->isBoundary()) nFaces++;

    vertPos.resize(nVerts, 3);
    faces.resize(nFaces * 3);

    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        vertPos.row(v->index) = v->position;
    }

    int fi = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            faces[fi * 3] = f->he->vertex->index;
            faces[fi * 3 + 1] = f->he->next->vertex->index;
            faces[fi * 3 + 2] = f->he->next->next->vertex->index;
            fi++;
        }
    }

    std::map<std::pair<int, int>, int> edgeMap;
    edgeList.clear();
    for (int fi2 = 0; fi2 < nFaces; fi2++) {
        for (int k = 0; k < 3; k++) {
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

    computeVertexNormals();
    estimatePrincipalCurvature();
    computeTransportMatching();
    computeLayerShift();

    if (!buildBranchCoverAndIntegrate()) {
        UV = Eigen::MatrixXd::Zero(nVerts, 2);
        return false;
    }

    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->uv = Eigen::Vector2d(UV(v->index, 0), UV(v->index, 1));
    }
    return true;
}
