#include "Tutte.h"
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <unordered_set>

Tutte::Tutte(Mesh& mesh0, TutteBoundary boundaryShape)
: Parameterization(mesh0), m_shape(boundaryShape) {}

void Tutte::findBoundaryLoop(std::vector<int>& boundaryVerts) const
{
    if (mesh.boundaries.empty()) return;
    
    // Walk along the first boundary loop
    HalfEdgeCIter heStart = mesh.boundaries[0];
    HalfEdgeCIter he = heStart;
    do {
        boundaryVerts.push_back(he->vertex->index);
        he = he->next;
    } while (he != heStart);
}

void Tutte::pinBoundary()
{
    std::vector<int> boundaryVerts;
    findBoundaryLoop(boundaryVerts);
    int n = (int)boundaryVerts.size();
    if (n < 3) return;
    
    for (int i = 0; i < n; i++) {
        double t = (double)i / n;  // [0, 1)
        double angle = t * 2.0 * M_PI;
        
        double u, v;
        if (m_shape == TutteBoundary::CIRCLE) {
            u = 0.5 + 0.5 * cos(angle);
            v = 0.5 + 0.5 * sin(angle);
        } else {  // SQUARE
            double s = t * 4.0; // [0, 4)
            if (s < 1.0)       { u = s;        v = 0; }
            else if (s < 2.0)  { u = 1;        v = s - 1; }
            else if (s < 3.0)  { u = 1 - (s-2); v = 1; }
            else               { u = 0;        v = 1 - (s-3); }
        }
        
        mesh.vertices[boundaryVerts[i]].uv = Eigen::Vector2d(u, v);
    }
}

void Tutte::solveLaplacian(const std::vector<int>& boundaryVerts)
{
    int N = (int)mesh.vertices.size();
    int nInterior = N - (int)boundaryVerts.size();
    if (nInterior <= 0) return;
    
    // Map: global -> interior idx (-1 for boundary)
    std::vector<int> iIdx(N, -1);
    std::unordered_set<int> bSet(boundaryVerts.begin(), boundaryVerts.end());
    int cnt = 0;
    for (int i = 0; i < N; i++) if (!bSet.count(i)) iIdx[i] = cnt++;
    
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd bx = Eigen::VectorXd::Zero(nInterior);
    Eigen::VectorXd by = Eigen::VectorXd::Zero(nInterior);
    
    // Iterate over faces to build cot-Laplacian (more robust than vertex iteration)
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        
        HalfEdgeCIter he = f->he;
        int vi[3]; 
        for (int k = 0; k < 3; k++) { vi[k] = he->vertex->index; he = he->next; }
        
        for (int k = 0; k < 3; k++) {
            int a = vi[k], b = vi[(k+1)%3], c = vi[(k+2)%3]; // edge (b,c) opposite to a
            
            // Cotangent of angle at vertex a
            Eigen::Vector3d e1 = mesh.vertices[b].position - mesh.vertices[a].position;
            Eigen::Vector3d e2 = mesh.vertices[c].position - mesh.vertices[a].position;
            double cotVal = e1.dot(e2) / std::max(e1.cross(e2).norm(), 1e-12);
            double w = std::max(cotVal, 1e-8);
            
            // Contribution to vertices b and c (the edge endpoints)
            for (int side = 0; side < 2; side++) {
                int p = (side==0) ? b : c;
                int q = (side==0) ? c : b;
                
                if (!bSet.count(p)) {
                    int rp = iIdx[p];
                    triplets.push_back(Eigen::Triplet<double>(rp, rp, w));
                    
                    if (bSet.count(q)) {
                        bx(rp) += w * mesh.vertices[q].uv.x();
                        by(rp) += w * mesh.vertices[q].uv.y();
                    } else {
                        int rq = iIdx[q];
                        triplets.push_back(Eigen::Triplet<double>(rp, rq, -w));
                    }
                }
            }
        }
    }
    
    Eigen::SparseMatrix<double> L(nInterior, nInterior);
    L.setFromTriplets(triplets.begin(), triplets.end());
    
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(L);
    if (solver.info() == Eigen::Success) {
        Eigen::VectorXd ux = solver.solve(bx);
        Eigen::VectorXd uy = solver.solve(by);
        for (int i = 0; i < N; i++) {
            if (!bSet.count(i)) {
                mesh.vertices[i].uv = Eigen::Vector2d(ux(iIdx[i]), uy(iIdx[i]));
            }
        }
    }
}

void Tutte::parameterize()
{
    pinBoundary();
    std::vector<int> boundaryVerts;
    findBoundaryLoop(boundaryVerts);
    solveLaplacian(boundaryVerts);
    normalize();
}
