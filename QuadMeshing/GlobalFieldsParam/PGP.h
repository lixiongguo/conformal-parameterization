#ifndef PGP_H
#define PGP_H

#include "GlobalFieldsParameterization.h"
#include <Eigen/SparseCholesky>
#include <vector>

/**
 * PGP: Periodic Global Parameterization (Ray et al. 2006).
 *
 * Uses alternative variables U_i=(cosθ_i, sinθ_i), V_i=(cosφ_i, sinφ_i)
 * on vertices to bake in 2π-periodicity.  π/2 rotational ambiguity is
 * handled by locally reorienting control vector fields per triangle.
 *
 * Algorithm:
 *   1. Mesh extraction, local frames, cross field
 *   2. Per-vertex rotation indices r_i from vector field alignment
 *   3. Build quadratic energy in U,V (edge-based with triangle integration)
 *   4. Linear solve → Newton iteration with norm-enforcing penalty
 *   5. Reconstruct (θ,φ) per triangle (Algorithm 1 from paper)
 *   6. BFS-propagate consistent UV and write to mesh
 */
class PGP : public GlobalFieldsParameterization {
public:
    PGP(Mesh& mesh0);

    void parameterize() override;

    /** Chart size parameter ω (default 10× average edge length). */
    void setOmega(double w) { omega = w; }

    /** Penalty weight for unit-norm constraint (default 1e-3). */
    void setPenaltyWeight(double eps) { penaltyEps = eps; }

    /** Maximum Newton outer iterations (default 5). */
    void setMaxNewtonIters(int iters) { maxNewtonIters = iters; }

protected:
    struct InternalEdge {
        int v1, v2;
        int f1, f2;
        int idx;
    };

    // --- Mesh extraction ---
    bool initMeshData();
    void buildLocalFrames();
    void buildLocalFrame(const Eigen::Vector3d& n,
                         Eigen::Vector3d& t1,
                         Eigen::Vector3d& t2);
    void computeCrossField();

    // --- Vertex rotations ---
    void computeVertexRotations();

    // --- Quadratic energy in (U,V) space ---
    void assembleSystemUV(
        Eigen::SparseMatrix<double>&  A,
        Eigen::VectorXd&              b,
        const std::vector<Eigen::Vector2d>& Ucur,
        const std::vector<Eigen::Vector2d>& Vcur);

    // --- Linear solve for initial guess ---
    void solveLinearInitial();

    // --- Newton iteration with penalty ---
    double computePenalty();
    void   addPenaltyGradientHessian(
        Eigen::VectorXd& grad,
        std::vector<Eigen::Triplet<double>>& hessTrips);
    void   newtonSolve();

    // --- Triangle integration weights (Appendix B, Eq 24) ---
    Eigen::Vector3d computeLambda(int f) const;

    // --- Per-triangle reconstruction (Algorithm 1) ---
    void reconstructPerTriangle();
    int  determineEdgeRotation(int f, int vA, int vB) const;

    // --- Consistent UV propagation ---
    void propagateUV();

    // --- Output ---
    void writeUV();

    // --- Helpers ---
    static Eigen::Vector2d rotateK(int k, const Eigen::Vector2d& v);
    static Eigen::Matrix2d rotMat2(double angle);
    Eigen::Vector2d  vecFieldAtVertexInFace(int f, int localVi, bool perp) const;

    // --- Mesh data ---
    int nV, nF;
    Eigen::MatrixXd vertPos;
    Eigen::MatrixXi faceIndices;
    std::vector<InternalEdge> edges;

    // --- Cross field ---
    std::vector<Eigen::Vector2d> faceD1;   // per-face d1 in 2D local frame
    std::vector<Eigen::Vector2d> faceD2;   // per-face d2 = rot90(d1)

    // --- Per-face geometry ---
    Eigen::MatrixXd faceNormals;
    std::vector<Eigen::Matrix2d> faceGrad;  // G_T ∈ R^{2×3}
    std::vector<double>  faceAreas;
    // edge vectors in 2D local frame (nF × 3 edges)
    std::vector<Eigen::Vector2d> faceEdges[3];

    // --- PGP alternative variables ---
    std::vector<Eigen::Vector2d> U;       // U_i = (cos θ_i, sin θ_i)
    std::vector<Eigen::Vector2d> V;       // V_i = (cos φ_i, sin φ_i)

    // --- Per-vertex rotation index ---
    std::vector<int> vertexR;             // r_i ∈ {0,1,2,3}

    // --- Reconstructed per-triangle coordinates ---
    std::vector<Eigen::Vector3d> triTheta;  // per-triangle (θ_1,θ_2,θ_3)
    std::vector<Eigen::Vector3d> triPhi;    // per-triangle (φ_1,φ_2,φ_3)
    std::vector<int>            triRot;     // per-triangle r_T

    // --- Sparse solver ---
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    bool solverReady;

    // --- Options ---
    double omega;
    double penaltyEps;
    int    maxNewtonIters;
};

#endif
