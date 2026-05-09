#ifndef BFF_H
#define BFF_H

#include "Parameterization.h"
#include <Eigen/SparseCholesky>

/**
 * BFF: Boundary First Flattening (Sawhney & Crane 2017)
 * 
 * Conformal parameterization by separating boundary and interior problems:
 *   1. Determine target boundary curvature (or shape)
 *   2. Solve Poisson equation for interior via Poincaré-Steklov operator
 *   3. Hilbert transform for conjugate harmonic function
 * 
 * Supports two modes:
 *   - CURVATURE:  Specify target boundary curvature → get conformal flattening
 *   - POSITION:   Specify target boundary positions → solve directly
 */
class BFF : public Parameterization {
public:
    enum Mode { CURVATURE = 0, POSITION = 1 };
    
    BFF(Mesh& mesh0);
    void parameterize() override;
    
    // Set target mode and parameters (call before parameterize())
    void setMode(Mode m) { mode = m; }
    
    // For CURVATURE mode: set target boundary curvature (uniform = 2π/B for circle)
    // angles[i] = target angle sum for boundary vertex i
    void setTargetBoundaryCurvature(const std::vector<double>& curvatures);
    
    // For POSITION mode: set target boundary UV positions
    void setTargetBoundaryPositions(const std::vector<double>& u, const std::vector<double>& v);
    
    // Set scale factor for boundary (controls boundary length in parameter domain)
    void setScale(double s) { scale = s; }

protected:
    // Step 1: Extract mesh data and compute cotan-Laplacian
    void extractMeshData();
    void computeCotLaplacian();
    
    // Step 2: Compute discrete boundary curvature
    void computeBoundaryCurvature();
    
    // Step 3: Compute Neumann boundary data from target curvature
    void computeNeumannData();
    
    // Step 4: Solve Poisson equation with Neumann BC → harmonic function a
    void solvePoissonNeumann();
    
    // Step 5: Hilbert transform → conjugate harmonic function b
    void hilbertTransform();
    
    // Step 6 (POSITION mode): Solve Laplace with Dirichlet BC
    void solveDirichlet();
    
    // Helper: cotan of angle at vertex c in triangle (a,b,c)
    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);
    
    // Helper: build local 2D frame
    void buildLocalFrame(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2);
    
    // Member variables
    Mode mode;
    double scale;
    
    // Mesh data
    Eigen::MatrixXd V;      // vertices (nV × 3)
    Eigen::VectorXi F;      // faces (nF × 3, flattened)
    int nV, nF;
    
    // Boundary data
    std::vector<int> boundaryVertices;      // ordered boundary vertex indices
    std::vector<double> boundaryCurvature;   // current discrete curvature per boundary vertex
    std::vector<double> targetCurvature;     // target curvature (CURVATURE mode)
    std::vector<double> targetBoundaryU;     // target boundary positions (POSITION mode)
    std::vector<double> targetBoundaryV;
    
    // Laplacian data
    Eigen::SparseMatrix<double> L;           // cotan-Laplacian
    Eigen::VectorXd massVec;                 // vertex areas (lumped mass)
    
    // Edge list for cotan weights
    struct Edge { int v1, v2, f1, f2; };
    std::vector<Edge> edges;
    int nEdges;
    
    // Output
    Eigen::MatrixXd UV;     // (nV × 2)
};

#endif
