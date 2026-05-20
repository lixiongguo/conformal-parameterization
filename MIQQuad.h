#ifndef MIQQUAD_H
#define MIQQUAD_H

#include "Parameterization.h"
#include <Eigen/SparseCholesky>

/**
 * MIQ (Mixed-Integer Quadrangulation) global parameterization.
 * 
 * Two-phase algorithm (Bommes et al. 2009):
 *   Phase 1: Cross field optimization via alternating MIQP
 *     - Alternating between continuous θ (angle per face) and
 *       integer p_ij (quarter-turn jumps per edge)
 *     - min Σ ||θ_i - θ_j + (π/2)·p_ij||²
 *   Phase 2: Poisson-based UV parameterization
 *     - Given optimized θ, compute target directions per face
 *     - Solve L·u = div_x, L·v = div_y (Hodge decomposition)
 */
class MIQQuad : public Parameterization {
public:
    MIQQuad(Mesh& mesh0);
    void parameterize() override;

protected:
    // Phase 1: Cross field optimization
    void initCrossField();                  // estimate initial θ from geometry
    void buildFaceLaplacian(Eigen::SparseMatrix<double>& L);
    bool optimizeCrossField(int maxIter=10);// alternating MIQP solver
    void buildTargetDirs();                 // target u,v per face from θ

    // Phase 2: Poisson solve for UV
    void buildCotLaplacian(Eigen::SparseMatrix<double>& L);
    void solvePoisson();                    // L·u = div_x, L·v = div_y

    // Helpers
    void buildLocalFrame(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2);
    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);
    int roundToInt(double x);

    // Data
    Eigen::MatrixXd vertPos;
    Eigen::VectorXi faces;       // flattened face indices
    int nV, nF, nE;
    struct Edge { int v1, v2, f1, f2; };
    std::vector<Edge> edgeList;

    Eigen::MatrixXd faceN;       // face normals
    Eigen::VectorXd theta;       // optimized angle per face (mod π/2)
    Eigen::VectorXi jump;        // integer p_ij per edge
    Eigen::MatrixXd faceD1, faceD2; // two orthogonal dirs per face
    Eigen::MatrixXd UV;
};

#endif
