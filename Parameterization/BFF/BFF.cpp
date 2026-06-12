#include "BFF.h"
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Eigen;

BFF::BFF(Mesh& mesh0) : Parameterization(mesh0), mode(CURVATURE), scale(1.0) {}

// ============================================================
// Helpers
// ============================================================

void BFF::buildLocalFrame(const Vector3d& n, Vector3d& t1, Vector3d& t2) {
    Vector3d ref(1, 0, 0);
    if (fabs(n.dot(ref)) > 0.9) ref = Vector3d(0, 1, 0);
    t1 = (ref - ref.dot(n) * n).normalized();
    t2 = n.cross(t1).normalized();
}

double BFF::cotan(const Vector3d& a, const Vector3d& b, const Vector3d& c) {
    Vector3d u = a - b, v = c - b;
    double d = u.dot(v), cr = u.cross(v).norm();
    if (cr < 1e-12) return 0.0;
    return d / cr;
}

// ============================================================
// Step 1: Extract mesh data
// ============================================================

void BFF::extractMeshData() {
    // Count vertices and faces
    nV = (int)mesh.vertices.size();
    nF = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++)
        if (!f->isBoundary()) nF++;

    V.resize(nV, 3);
    F.resize(nF * 3);

    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
        V.row(v->index) = v->position;

    int fi = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            F[fi*3]   = f->he->vertex->index;
            F[fi*3+1] = f->he->next->vertex->index;
            F[fi*3+2] = f->he->next->next->vertex->index;
            fi++;
        }
    }
}

// ============================================================
// Step 1b: Build cotan-Laplacian + edge list + mass matrix
// ============================================================

void BFF::computeCotLaplacian() {
    // Build edge list
    std::map<std::pair<int,int>, int> edgeMap;
    edges.clear();
    for (int fi = 0; fi < nF; fi++) {
        for (int k = 0; k < 3; k++) {
            int v1 = F[fi*3+k], v2 = F[fi*3+(k+1)%3];
            if (v1 > v2) std::swap(v1, v2);
            auto key = std::make_pair(v1, v2);
            auto it = edgeMap.find(key);
            if (it == edgeMap.end()) {
                edgeMap[key] = (int)edges.size();
                Edge e; e.v1 = v1; e.v2 = v2; e.f1 = fi; e.f2 = -1;
                edges.push_back(e);
            } else {
                edges[it->second].f2 = fi;
            }
        }
    }
    nEdges = (int)edges.size();

    // Build cotan-Laplacian L
    L.resize(nV, nV);
    massVec = VectorXd::Zero(nV);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nV);

    for (auto& e : edges) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        int vi = e.v1, vj = e.v2;

        // Find opposite vertices for cotan computation
        int o1 = -1, o2 = -1;
        for (int k = 0; k < 3; k++) {
            int v = F[e.f1*3+k];
            if (v != vi && v != vj) { o1 = v; break; }
        }
        for (int k = 0; k < 3; k++) {
            int v = F[e.f2*3+k];
            if (v != vi && v != vj) { o2 = v; break; }
        }

        double w = 0.5 * (cotan(V.row(o1), V.row(vi), V.row(vj)) +
                          cotan(V.row(o2), V.row(vi), V.row(vj)));
        if (w < 0) w = 0;

        diag(vi) += w; diag(vj) += w;
        trips.push_back(Triplet<double>(vi, vj, -w));
        trips.push_back(Triplet<double>(vj, vi, -w));
    }

    // Diagonal + small regularization
    for (int i = 0; i < nV; i++)
        trips.push_back(Triplet<double>(i, i, diag(i) + 1e-8));

    L.setFromTriplets(trips.begin(), trips.end());

    // Build mass matrix (lumped: vertex area = 1/3 sum of adjacent face areas)
    for (int fi = 0; fi < nF; fi++) {
        int v0 = F[fi*3], v1 = F[fi*3+1], v2 = F[fi*3+2];
        Vector3d p0 = V.row(v0), p1 = V.row(v1), p2 = V.row(v2);
        double area = 0.5 * (p1 - p0).cross(p2 - p0).norm();
        double third = area / 3.0;
        massVec(v0) += third; massVec(v1) += third; massVec(v2) += third;
    }
}

// ============================================================
// Step 2: Compute discrete boundary curvature
// ============================================================

void BFF::computeBoundaryCurvature() {
    boundaryVertices.clear();
    boundaryCurvature.clear();

    // Find boundary half-edges and order them
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->isBoundary()) {
            boundaryVertices.push_back(v->index);
        }
    }

    // Order boundary vertices by walking along the boundary
    // Find a boundary edge to start
    EdgeCIter startEdge = mesh.edges.end();
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (e->isBoundary()) { startEdge = e; break; }
    }
    if (startEdge == mesh.edges.end()) return;

    // Walk boundary in order
    std::vector<int> orderedBnd;
    HalfEdgeCIter he = startEdge->he;
    if (!he->onBoundary) he = he->flip;
    
    VertexCIter startV = he->vertex;
    orderedBnd.push_back(startV->index);
    
    HalfEdgeCIter cur = he->next;
    while (cur->vertex != startV) {
        if (cur->onBoundary) {
            orderedBnd.push_back(cur->vertex->index);
            cur = cur->next;
        } else {
            cur = cur->flip->next;
        }
        if (orderedBnd.size() > nV) break; // safety
    }

    if (orderedBnd.size() < 3) {
        boundaryVertices = orderedBnd;
        for (size_t i = 0; i < orderedBnd.size(); i++)
            boundaryCurvature.push_back(0.0);
        return;
    }

    boundaryVertices = orderedBnd;
    int B = (int)orderedBnd.size();

    // Compute discrete curvature: π - angle at boundary vertex
    // (angle formed by two incident boundary edges)
    boundaryCurvature.resize(B, 0.0);
    
    for (int i = 0; i < B; i++) {
        int prev = orderedBnd[(i - 1 + B) % B];
        int curr = orderedBnd[i];
        int next = orderedBnd[(i + 1) % B];
        
        Vector3d pPrev = V.row(prev), pCurr = V.row(curr), pNext = V.row(next);
        
        Vector3d e1 = (pPrev - pCurr).normalized();
        Vector3d e2 = (pNext - pCurr).normalized();
        
        double dot = std::max(-1.0, std::min(1.0, e1.dot(e2)));
        double angle = acos(dot);
        
        // Curvature = π - interior angle at boundary vertex
        // For boundary vertex, the angle is the external angle defect
        boundaryCurvature[i] = M_PI - angle;
    }
}

// ============================================================
// Step 3: Compute Neumann data from target curvature (CURVATURE mode)
// ============================================================

void BFF::computeNeumannData() {
    int B = (int)boundaryVertices.size();
    if (B < 3 || targetCurvature.size() != (size_t)B) return;

    // Compute Neumann data: h_i = κ_i - κ̃_i  (difference between current and target curvature)
    // h is the Neumann boundary data for the Poisson equation
    // h vector will be stored and used in solvePoissonNeumann
}

// ============================================================
// Step 4: Solve Poisson with Neumann BC → harmonic a (CURVATURE mode)
// ============================================================

void BFF::solvePoissonNeumann() {
    int B = (int)boundaryVertices.size();
    if (B < 3) return;

    UV.resize(nV, 2);
    
    // Target curvature: default to uniform (circle boundary)
    // Total boundary curvature = 2π for simply-connected domain
    // Uniform distribution: 2π/B per boundary vertex
    double uniformCurv = 2.0 * M_PI / B;
    
    // Neumann data h = current curvature - target curvature
    // Solve L*a = 0 with Neumann BC ∂a/∂n = h
    VectorXd h(B);
    for (int i = 0; i < B; i++) {
        double target = (targetCurvature.empty()) ? uniformCurv : targetCurvature[i];
        h(i) = boundaryCurvature[i] - target;
    }

    // Incorporate Neumann BC into the Laplacian system:
    // The boundary condition ∂a/∂n = h modifies the RHS for boundary vertices
    // For each boundary edge, the Neumann integral contributes to boundary vertices
    
    // Build RHS: b = 0 interior, b = -mass*(curvature correction) on boundary
    VectorXd rhs = VectorXd::Zero(nV);
    
    // Boundary curvature correction via Neumann condition
    // For boundary vertex i, the Neumann boundary integral over incident edges
    // contributes to the RHS. The discrete Neumann BC is encoded as:
    //   (L*a)_i = -∫_{∂Ω_i} h ds  for boundary vertex i
    // where the integral is over the boundary edges incident to vertex i.
    for (int i = 0; i < B; i++) {
        int v = boundaryVertices[i];
        int vPrev = boundaryVertices[(i - 1 + B) % B];
        int vNext = boundaryVertices[(i + 1) % B];
        
        // Edge lengths for boundary integration
        double len1 = (V.row(v) - V.row(vPrev)).norm();
        double len2 = (V.row(vNext) - V.row(v)).norm();
        double totalLen = len1 + len2;
        
        if (totalLen > 1e-12) {
            // Neumann integral: ∫ h ds over half-edge on each side
            double integral = 0.5 * len1 * h((i - 1 + B) % B) + 0.5 * len2 * h(i);
            rhs(v) = -integral * scale;
        }
    }

    // Fix one vertex to remove nullspace (adds constraint for unique solution)
    int anchor = boundaryVertices[0];
    
    // Penalty method: add large weight to anchor
    SparseMatrix<double> Lm = L;
    for (int k = 0; k < Lm.outerSize(); k++) {
        for (SparseMatrix<double>::InnerIterator it(Lm, k); it; ++it) {
            if (it.row() == anchor || it.col() == anchor) {
                // handled below
            }
        }
    }
    Lm.coeffRef(anchor, anchor) += 1e8;
    rhs(anchor) = 0.0;

    // Solve
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(Lm);
    if (solver.info() != Success) return;

    VectorXd a = solver.solve(rhs);

    // Store as u-coordinate
    for (int i = 0; i < nV; i++)
        UV(i, 0) = a(i);
}

// ============================================================
// Step 5: Hilbert transform → conjugate harmonic b
// ============================================================

void BFF::hilbertTransform() {
    if (UV.rows() == 0) return;

    // The Hilbert transform: given harmonic a, find conjugate harmonic b
    // such that a + i*b is holomorphic (i.e., conformal)
    // 
    // Method: Compute gradient of a per face, rotate by 90° to get gradient of b,
    // then solve Poisson equation for b from the divergence of the rotated gradient.
    //
    // ∇b = (∇a)^⊥  where ⊥ is 90° counterclockwise rotation in the tangent plane

    // Step 1: Per-face gradient of a
    MatrixXd gradA(nF, 3);
    gradA.setZero();

    for (int fi = 0; fi < nF; fi++) {
        int v0 = F[fi*3], v1 = F[fi*3+1], v2 = F[fi*3+2];
        
        Vector3d p0 = V.row(v0), p1 = V.row(v1), p2 = V.row(v2);
        double a0 = UV(v0, 0), a1 = UV(v1, 0), a2 = UV(v2, 0);
        
        // Face normal and area
        Vector3d fn = (p1 - p0).cross(p2 - p0);
        double area = 0.5 * fn.norm();
        if (area < 1e-12) continue;
        
        // Gradients in 3D: ∇a = (1/2A) * Σ_i a_i * (n × e_i)
        // where e_i are edges opposite to vertex i
        Vector3d e01 = p1 - p0, e12 = p2 - p1, e20 = p0 - p2;
        
        Vector3d ga = (1.0 / (2.0 * area)) * (
            a0 * fn.cross(e12) + a1 * fn.cross(e20) + a2 * fn.cross(e01)
        );
        
        // Rotate by 90° CCW in the face: gb = fn × ga
        Vector3d gb = fn.normalized().cross(ga);
        
        gradA.row(fi) = gb; // store rotated gradient (= gradient of b)
    }

    // Step 2: Interpolate face gradients to vertices (area-weighted)
    VectorXd vertArea = VectorXd::Zero(nV);
    MatrixXd Vgrad(nV, 3);
    Vgrad.setZero();
    
    for (int fi = 0; fi < nF; fi++) {
        Vector3d p0 = V.row(F[fi*3]), p1 = V.row(F[fi*3+1]), p2 = V.row(F[fi*3+2]);
        double area = 0.5 * (p1 - p0).cross(p2 - p0).norm();
        
        for (int k = 0; k < 3; k++) {
            int v = F[fi*3+k];
            vertArea(v) += area;
            Vgrad.row(v) += area * gradA.row(fi);
        }
    }
    for (int i = 0; i < nV; i++) {
        if (vertArea(i) > 1e-12)
            Vgrad.row(i) /= vertArea(i);
    }

    // Step 3: Compute divergence of the vertex-interpolated gradient field
    VectorXd div = VectorXd::Zero(nV);
    
    for (int fi = 0; fi < nF; fi++) {
        int v0 = F[fi*3], v1 = F[fi*3+1], v2 = F[fi*3+2];
        Vector3d p0 = V.row(v0), p1 = V.row(v1), p2 = V.row(v2);
        Vector3d fn = (p1 - p0).cross(p2 - p0);
        
        Vector3d e01 = p1 - p0, e12 = p2 - p1, e20 = p0 - p2;
        Vector3d en01 = fn.cross(e01).normalized();
        Vector3d en12 = fn.cross(e12).normalized();
        Vector3d en20 = fn.cross(e20).normalized();
        
        double cot0 = cotan(p2, p0, p1);
        double cot1 = cotan(p0, p1, p2);
        double cot2 = cotan(p1, p2, p0);
        
        Vector3d g = gradA.row(fi);
        
        double d = 0.5 * (g.dot(en01) * cot2 + g.dot(en12) * cot0 + g.dot(en20) * cot1);
        
        div(v0) += d; div(v1) += d; div(v2) += d;
    }

    // Step 4: Solve Poisson equation L*b = div
    // Fix one vertex
    int anchor = boundaryVertices.empty() ? 0 : boundaryVertices[0];
    
    SparseMatrix<double> Lm = L;
    Lm.coeffRef(anchor, anchor) += 1e8;
    VectorXd bRhs = div;
    bRhs(anchor) = 0.0;

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(Lm);
    if (solver.info() != Success) return;

    VectorXd b = solver.solve(bRhs);

    // Store as v-coordinate
    for (int i = 0; i < nV; i++)
        UV(i, 1) = b(i);
}

// ============================================================
// Step 6: Solve Laplace with Dirichlet BC (POSITION mode)
// ============================================================

void BFF::solveDirichlet() {
    int B = (int)boundaryVertices.size();
    if (B < 3) return;

    UV.resize(nV, 2);
    
    // Find interior vertices
    std::vector<bool> isBoundary(nV, false);
    for (int idx : boundaryVertices) isBoundary[idx] = true;

    // Build interior-only Laplacian: L_II
    std::vector<int> interiorToGlobal;
    std::map<int, int> globalToInterior;
    for (int i = 0; i < nV; i++) {
        if (!isBoundary[i]) {
            globalToInterior[i] = (int)interiorToGlobal.size();
            interiorToGlobal.push_back(i);
        }
    }
    int nI = (int)interiorToGlobal.size();
    if (nI == 0) return;

    SparseMatrix<double> L_II(nI, nI);
    VectorXd bu = VectorXd::Zero(nI), bv = VectorXd::Zero(nI);
    std::vector<Triplet<double>> trips;

    for (auto& e : edges) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        int vi = e.v1, vj = e.v2;

        int o1 = -1, o2 = -1;
        for (int k = 0; k < 3; k++) { int v = F[e.f1*3+k]; if (v != vi && v != vj) { o1 = v; break; } }
        for (int k = 0; k < 3; k++) { int v = F[e.f2*3+k]; if (v != vi && v != vj) { o2 = v; break; } }

        double w = 0.5 * (cotan(V.row(o1), V.row(vi), V.row(vj)) +
                          cotan(V.row(o2), V.row(vi), V.row(vj)));
        if (w < 0) w = 0;

        auto itI_i = globalToInterior.find(vi);
        auto itI_j = globalToInterior.find(vj);

        if (itI_i != globalToInterior.end() && itI_j != globalToInterior.end()) {
            // Both interior: L_II contribution
            int ii = itI_i->second, ij = itI_j->second;
            trips.push_back(Triplet<double>(ii, ij, -w));
            trips.push_back(Triplet<double>(ij, ii, -w));
            trips.push_back(Triplet<double>(ii, ii, w));
            trips.push_back(Triplet<double>(ij, ij, w));
        } else if (itI_i != globalToInterior.end()) {
            // vi interior, vj boundary: -L_IB * g contribution to RHS
            int ii = itI_i->second;
            trips.push_back(Triplet<double>(ii, ii, w));
            if (isBoundary[vi]) continue; // shouldn't happen
            // Move contribution: w * g(vj) to RHS
            if ((size_t)vj < targetBoundaryU.size()) {
                bu(ii) += w * targetBoundaryU[vj];
                bv(ii) += w * targetBoundaryV[vj];
            }
        } else if (itI_j != globalToInterior.end()) {
            int ij = itI_j->second;
            trips.push_back(Triplet<double>(ij, ij, w));
            if ((size_t)vi < targetBoundaryU.size()) {
                bu(ij) += w * targetBoundaryU[vi];
                bv(ij) += w * targetBoundaryV[vi];
            }
        }
    }
    
    // Regularization
    for (int i = 0; i < nI; i++)
        trips.push_back(Triplet<double>(i, i, 1e-8));

    L_II.setFromTriplets(trips.begin(), trips.end());

    // Solve
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(L_II);
    if (solver.info() != Success) return;

    VectorXd uInterior = solver.solve(bu);
    VectorXd vInterior = solver.solve(bv);

    // Assemble full UV
    for (int i = 0; i < nV; i++) {
        if (isBoundary[i]) {
            UV(i, 0) = targetBoundaryU[i];
            UV(i, 1) = targetBoundaryV[i];
        } else {
            int ii = globalToInterior[i];
            UV(i, 0) = uInterior(ii);
            UV(i, 1) = vInterior(ii);
        }
    }
}

// ============================================================
// Public setters
// ============================================================

void BFF::setTargetBoundaryCurvature(const std::vector<double>& curvatures) {
    targetCurvature = curvatures;
}

void BFF::setTargetBoundaryPositions(const std::vector<double>& u, const std::vector<double>& v) {
    targetBoundaryU = u;
    targetBoundaryV = v;
}

// ============================================================
// Main entry point
// ============================================================

void BFF::parameterize() {
    // Extract data
    extractMeshData();
    computeCotLaplacian();
    computeBoundaryCurvature();

    if (mode == POSITION) {
        // Dirichlet mode: user specifies target boundary shape
        solveDirichlet();
    } else {
        // Curvature mode: specify target curvature → Neumann solve + Hilbert
        solvePoissonNeumann();
        hilbertTransform();
    }

    // Copy UV back to mesh
    // Normalize to [0,1] × [0,1]
    double uMin = UV.col(0).minCoeff(), uMax = UV.col(0).maxCoeff();
    double vMin = UV.col(1).minCoeff(), vMax = UV.col(1).maxCoeff();
    double uRng = uMax - uMin, vRng = vMax - vMin;
    if (uRng < 1e-10) uRng = 1.0;
    if (vRng < 1e-10) vRng = 1.0;

    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int idx = v->index;
        v->uv = Vector2d((UV(idx, 0) - uMin) / uRng, (UV(idx, 1) - vMin) / vRng);
    }
}
