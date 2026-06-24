#include "Tutte.h"
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <unordered_set>

Tutte::Tutte(Mesh& mesh0, TutteBoundary boundaryShape, TutteWeight weight)
: Parameterization(mesh0), m_shape(boundaryShape), m_weight(weight) {}

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

void Tutte::pinTwoVertices(std::vector<int>& pinVerts)
{
    pinVerts.clear();
    int bestA = 0, bestB = 1;
    double maxDist2 = 0.0;
    
    // Farthest pair on boundary loops (same heuristic as LSCM)
    for (HalfEdgeCIter heStart : mesh.boundaries) {
        HalfEdgeCIter he1 = heStart;
        do {
            const int vIdx1 = he1->vertex->index;
            const Eigen::Vector3d& p1 = he1->vertex->position;
            
            HalfEdgeCIter he2 = heStart;
            do {
                const int vIdx2 = he2->vertex->index;
                if (vIdx1 >= vIdx2) { he2 = he2->next; continue; }
                
                const Eigen::Vector3d& p2 = he2->vertex->position;
                const double d2 = (p2 - p1).squaredNorm();
                if (d2 > maxDist2) {
                    maxDist2 = d2;
                    bestA = vIdx1;
                    bestB = vIdx2;
                }
                he2 = he2->next;
            } while (he2 != heStart);
            
            he1 = he1->next;
        } while (he1 != heStart);
    }
    
    // Fallback: no boundary (should not happen for disk meshes)
    if (maxDist2 < 1e-24) {
        const int N = (int)mesh.vertices.size();
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                const double d2 = (mesh.vertices[j].position - mesh.vertices[i].position).squaredNorm();
                if (d2 > maxDist2) {
                    maxDist2 = d2;
                    bestA = i;
                    bestB = j;
                }
            }
        }
    }
    
    const double d = std::sqrt(std::max(maxDist2, 1e-24));
    mesh.vertices[bestA].uv = Eigen::Vector2d(0.0, 0.0);
    mesh.vertices[bestB].uv = Eigen::Vector2d(d, 0.0);
    pinVerts = {bestA, bestB};
}

void Tutte::solveLaplacian(const std::vector<int>& fixedVerts)
{
    int N = (int)mesh.vertices.size();
    int nFree = N - (int)fixedVerts.size();
    if (nFree <= 0) return;
    
    // Map: global -> free idx (-1 for fixed)
    std::vector<int> iIdx(N, -1);
    std::unordered_set<int> fixedSet(fixedVerts.begin(), fixedVerts.end());
    int cnt = 0;
    for (int i = 0; i < N; i++) if (!fixedSet.count(i)) iIdx[i] = cnt++;
    
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd bx = Eigen::VectorXd::Zero(nFree);
    Eigen::VectorXd by = Eigen::VectorXd::Zero(nFree);
    
    // Iterate over faces to build Laplacian (cotan or uniform weights)
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        
        HalfEdgeCIter he = f->he;
        int vi[3]; 
        for (int k = 0; k < 3; k++) { vi[k] = he->vertex->index; he = he->next; }
        
        for (int k = 0; k < 3; k++) {
            int a = vi[k], b = vi[(k+1)%3], c = vi[(k+2)%3]; // edge (b,c) opposite to a
            
            double w;
            if (m_weight == TutteWeight::UNIFORM) {
                w = 1.0;
            } else {
                Eigen::Vector3d e1 = mesh.vertices[b].position - mesh.vertices[a].position;
                Eigen::Vector3d e2 = mesh.vertices[c].position - mesh.vertices[a].position;
                double cotVal = e1.dot(e2) / std::max(e1.cross(e2).norm(), 1e-12);
                w = std::max(cotVal, 1e-8);
            }
            
            for (int side = 0; side < 2; side++) {
                int p = (side==0) ? b : c;
                int q = (side==0) ? c : b;
                
                if (!fixedSet.count(p)) {
                    int rp = iIdx[p];
                    triplets.push_back(Eigen::Triplet<double>(rp, rp, w));
                    
                    if (fixedSet.count(q)) {
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
    
    Eigen::SparseMatrix<double> L(nFree, nFree);
    L.setFromTriplets(triplets.begin(), triplets.end());
    
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(L);
    if (solver.info() == Eigen::Success) {
        Eigen::VectorXd ux = solver.solve(bx);
        Eigen::VectorXd uy = solver.solve(by);
        for (int i = 0; i < N; i++) {
            if (!fixedSet.count(i)) {
                mesh.vertices[i].uv = Eigen::Vector2d(ux(iIdx[i]), uy(iIdx[i]));
            }
        }
    }
}

void Tutte::parameterize()
{
    std::vector<int> fixedVerts;
    if (m_shape == TutteBoundary::FREE) {
        pinTwoVertices(fixedVerts);
    } else {
        pinBoundary();
        findBoundaryLoop(fixedVerts);
    }
    solveLaplacian(fixedVerts);
    normalize();
}
