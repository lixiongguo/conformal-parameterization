#include "ARAP.h"
#include "Tutte.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SVD>
#include <cmath>

ARAP::ARAP(Mesh& mesh0, int maxIter)
: Parameterization(mesh0), m_maxIter(maxIter) {}

void ARAP::computeLocalFrames()
{
    m_localRefs.clear();
    m_localRefs.resize(mesh.faces.size());
    m_rotations.resize(mesh.faces.size());
    
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        int fi = f->index;
        
        // Get three vertices of the triangle
        HalfEdgeCIter he = f->he;
        std::vector<Eigen::Vector3d> pos3D(3);
        for (int i = 0; i < 3; i++) {
            pos3D[i] = he->vertex->position;
            he = he->next;
        }
        
        // Compute local 2D orthonormal basis (tangent plane)
        Eigen::Vector3d e1 = pos3D[1] - pos3D[0];
        Eigen::Vector3d e2 = pos3D[2] - pos3D[0];
        
        Eigen::Vector3d xAxis = e1.normalized();
        Eigen::Vector3d zAxis = xAxis.cross(e2);
        zAxis.normalize();
        Eigen::Vector3d yAxis = zAxis.cross(xAxis);
        yAxis.normalize();
        
        // Project 3D triangle to 2D local frame
        m_localRefs[fi].resize(3);
        m_localRefs[fi][0] = Eigen::Vector2d(0, 0);
        m_localRefs[fi][1] = Eigen::Vector2d(e1.dot(xAxis), e1.dot(yAxis));
        m_localRefs[fi][2] = Eigen::Vector2d(e2.dot(xAxis), e2.dot(yAxis));
        
        // Initialize rotations to identity
        m_rotations[fi] = Eigen::Matrix2d::Identity();
    }
}

void ARAP::localStep()
{
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        int fi = f->index;
        
        // Get UV coordinates of this triangle
        HalfEdgeCIter he = f->he;
        std::vector<Eigen::Vector2d> uvs(3);
        for (int i = 0; i < 3; i++) {
            uvs[i] = he->vertex->uv;
            he = he->next;
        }
        
        // Compute cotangent weights for edges
        // For edge between vertex j and k, the weight is cot(angle at i)
        const auto& x = m_localRefs[fi];
        
        Eigen::Matrix2d S = Eigen::Matrix2d::Zero();
        
        const int edges[3][2] = {{1,2},{2,0},{0,1}};
        const int opposite[3]  = {0, 1, 2};
        
        for (int e = 0; e < 3; e++) {
            int j = edges[e][0];
            int k = edges[e][1];
            int i = opposite[e];
            
            Eigen::Vector2d dx = x[j] - x[k];
            Eigen::Vector2d du = uvs[j] - uvs[k];
            
            // Cotangent of angle at vertex i
            Eigen::Vector2d a = x[j] - x[i];
            Eigen::Vector2d b = x[k] - x[i];
            double cotVal = (a.x()*b.x() + a.y()*b.y()) / 
                            std::max(std::abs(a.x()*b.y() - a.y()*b.x()), 1e-12);
            double w = std::max(cotVal, 1e-8);  // ensure positive for convex combination
            
            // Cross-covariance contribution
            S += w * du * dx.transpose();
        }
        
        // Extract optimal 2D rotation from S
        double a = S(0,0), b = S(0,1), c = S(1,0), d = S(1,1);
        double denom = std::sqrt((a + d)*(a + d) + (c - b)*(c - b));
        if (denom < 1e-12) {
            m_rotations[fi] = Eigen::Matrix2d::Identity();
        } else {
            double cosTheta = (a + d) / denom;
            double sinTheta = (c - b) / denom;
            m_rotations[fi] << cosTheta, -sinTheta,
                               sinTheta,  cosTheta;
        }
    }
}

void ARAP::globalStep()
{
    int N = (int)mesh.vertices.size();
    
    // Find boundary vertices (to fix them)
    std::vector<bool> isBoundary(N, false);
    for (auto& he : mesh.boundaries) {
        HalfEdgeCIter h = he;
        do {
            isBoundary[h->vertex->index] = true;
            h = h->next;
        } while (h != he);
    }
    
    // Map vertex -> interior index
    std::vector<int> interiorIdx(N, -1);
    int nInterior = 0;
    for (int i = 0; i < N; i++) {
        if (!isBoundary[i]) interiorIdx[i] = nInterior++;
    }
    
    if (nInterior == 0) return;
    
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd bx = Eigen::VectorXd::Zero(nInterior);
    Eigen::VectorXd by = Eigen::VectorXd::Zero(nInterior);
    
    // Accumulate contributions from each face
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        int fi = f->index;
        const auto& R = m_rotations[fi];
        const auto& x = m_localRefs[fi];
        
        HalfEdgeCIter he = f->he;
        std::vector<int> vIdx(3);
        for (int a = 0; a < 3; a++) {
            vIdx[a] = he->vertex->index;
            he = he->next;
        }
        
        // For each edge in this triangle
        const int edges[3][2] = {{1,2},{2,0},{0,1}};
        const int opposite[3]  = {0, 1, 2};
        
        for (int e = 0; e < 3; e++) {
            int a = edges[e][0];
            int b_idx = edges[e][1];
            int opp = opposite[e];
            
            Eigen::Vector2d dx = x[a] - x[b_idx];
            
            // Cotangent weight
            Eigen::Vector2d ea = x[a] - x[opp];
            Eigen::Vector2d eb = x[b_idx] - x[opp];
            double cotVal = (ea.x()*eb.x() + ea.y()*eb.y()) /
                            std::max(std::abs(ea.x()*eb.y() - ea.y()*eb.x()), 1e-12);
            double w = std::max(cotVal, 1e-8);
            
            // Rotated reference edge
            Eigen::Vector2d rotDx = R * (x[a] - x[b_idx]);  // R * (x_i - x_j)
            
            int vi = vIdx[a];
            int vj = vIdx[b_idx];
            
            // Normal equations for ||(u_i - u_j) - R(x_i - x_j)||^2.
            // Both endpoints of the edge must contribute with opposite RHS signs.
            if (!isBoundary[vi]) {
                int ri = interiorIdx[vi];
                triplets.push_back(Eigen::Triplet<double>(ri, ri, w));
                bx(ri) += w * rotDx.x();
                by(ri) += w * rotDx.y();
                
                if (isBoundary[vj]) {
                    bx(ri) += w * mesh.vertices[vj].uv.x();
                    by(ri) += w * mesh.vertices[vj].uv.y();
                } else {
                    int cj = interiorIdx[vj];
                    triplets.push_back(Eigen::Triplet<double>(ri, cj, -w));
                }
            }

            if (!isBoundary[vj]) {
                int rj = interiorIdx[vj];
                triplets.push_back(Eigen::Triplet<double>(rj, rj, w));
                bx(rj) -= w * rotDx.x();
                by(rj) -= w * rotDx.y();
                
                if (isBoundary[vi]) {
                    bx(rj) += w * mesh.vertices[vi].uv.x();
                    by(rj) += w * mesh.vertices[vi].uv.y();
                } else {
                    int ci = interiorIdx[vi];
                    triplets.push_back(Eigen::Triplet<double>(rj, ci, -w));
                }
            }
        }
    }
    
    // Solve
    Eigen::SparseMatrix<double> L(nInterior, nInterior);
    L.setFromTriplets(triplets.begin(), triplets.end());
    
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L);
    
    if (solver.info() == Eigen::Success) {
        Eigen::VectorXd ux = solver.solve(bx);
        Eigen::VectorXd uy = solver.solve(by);
        
        for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
            int vi = v->index;
            if (!isBoundary[vi]) {
                int ri = interiorIdx[vi];
                v->uv = Eigen::Vector2d(ux(ri), uy(ri));
            }
        }
    }
}

void ARAP::initTutte()
{
    Tutte tutte(mesh, TutteBoundary::CIRCLE);
    tutte.parameterize();
}

void ARAP::parameterize()
{
    computeLocalFrames();
    initTutte();
    
    for (int iter = 0; iter < m_maxIter; iter++) {
        localStep();
        globalStep();
    }
    
    normalize();
}
