#include "QuadCover.h"
#include <map>
#include <algorithm>

using namespace Eigen;

QuadCover::QuadCover(Mesh& mesh0): Parameterization(mesh0) {}

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
    for (int ei = 0; ei < nEdges; ei++) {
        auto& e = edgeList[ei];
        if (e.f1 < 0 || e.f2 < 0) continue;
        layerShift(e.v1) += matching[ei];
        layerShift(e.v2) += matching[ei];
    }
    layerShift /= 4.0;
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
