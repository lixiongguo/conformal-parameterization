#ifndef QUADCOVER_H
#define QUADCOVER_H

#include "Parameterization.h"
#include <Eigen/SparseCholesky>

/**
 * QuadCover: Global seamless quad parameterization via branch covering.
 * 
 * Algorithm phases (Kälberer et al. 2007):
 *   0. Compute vertex normals, then Weingarten map W = I⁻¹·II per face
 *      → eigen-decomposition yields principal curvature direction
 *   1. Compute matching r_ij ∈ {0,1,2,3} per edge
 *   2. Compute layer shift ls(v) = Σ r / 4, identify singularities
 *   3. Build covering vector field (two orthogonal directions per face)
 *   4. Hodge decomposition → Poisson solve: L·u = div(X), L·v = div(Y)
 *   5. Global continuity & normalization
 */
class QuadCover : public Parameterization {
public:
    QuadCover(Mesh& mesh0);
    void parameterize() override;
    void setFaceDirections(const Eigen::MatrixXd& dirs);
    
protected:
    // Phase 0: Compute vertex normals (area-weighted face normals)
    void computeVertexNormals();
    
    // Phase 0: Weingarten map W = I⁻¹·II → eigen-decomposition for principal directions
    void estimatePrincipalCurvature();
    
    // Phase 1: Compute matching per edge
    void computeMatching();
    
    // Phase 2: Compute layer shift per vertex
    void computeLayerShift();
    
    // Phase 3: Build two orthogonal direction fields per face
    void buildCoveringField();
    
    // Phase 4: Poisson solve for UV (Hodge decomposition)
    void solvePoisson();
    
    // helpers
    void buildLocalFrame(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2);
    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);
    void buildCotLaplacian(Eigen::SparseMatrix<double>& L);
    
    // member variables
    Eigen::MatrixXd vertPos;          // V rows × 3
    Eigen::VectorXi faces;            // F rows × 3, flattened
    int nVerts, nFaces;
    
    Eigen::MatrixXd vertexNormals;    // per-vertex normal (nVerts × 3), area-weighted
    Eigen::MatrixXd faceNormals;      // per-face normal (nFaces × 3)
    Eigen::MatrixXd faceDirs;         // per-face reference direction (nFaces × 3)
    Eigen::MatrixXd externalFaceDirs; // optional input frame field directions
    Eigen::MatrixXd faceD1, faceD2;   // covering field: two orthogonal dirs per face
    bool hasExternalFaceDirs;
    
    Eigen::VectorXi matching;         // r_ij per edge, resized to nEdges
    Eigen::VectorXd layerShift;       // ls(v) per vertex
    Eigen::MatrixXd UV;               // output UV (nVerts × 2)
    
    // edge list
    struct Edge { int v1, v2, f1, f2; };
    std::vector<Edge> edgeList;
    int nEdges;
};

#endif
