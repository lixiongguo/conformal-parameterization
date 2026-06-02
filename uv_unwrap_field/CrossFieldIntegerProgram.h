#ifndef CROSSFIELDINTEGERPROGRAM_H
#define CROSSFIELDINTEGERPROGRAM_H

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <vector>

/**
 * Discrete MIQP surrogate for 4-RoSy cross fields (MIQ phase 1, Bommes et al. 2009).
 *
 *   min_{θ,p} Σ_e (θ_{f1} - θ_{f2} + (π/2)·p_e)²,   p_e ∈ ℤ
 *
 * Solver:
 *   1. Alternating optimization — fix p, linear solve for θ; fix θ, round p.
 *   2. Coordinate descent on p with θ re-solved per candidate.
 *
 * Mesh-agnostic: supply internal face–face edges only (f1, f2 >= 0).
 */
class CrossFieldIntegerProgram {
public:
    struct Edge {
        int f1 = -1;
        int f2 = -1;
        int idx = -1;
    };

    CrossFieldIntegerProgram(int numFaces, std::vector<Edge> internalEdges);

    void setJumpBounds(int lo, int hi);
    void setAlternatingIterations(int iters) { alternatingIters_ = iters; }
    void setRefinePasses(int passes) { refinePasses_ = passes; }

    /** Run alternating + refinement; updates theta and jump in place. */
    bool solve(Eigen::VectorXd& theta, Eigen::VectorXi& jump);

    double energy(const Eigen::VectorXd& theta, const Eigen::VectorXi& jump) const;

    int numFaces() const { return nF_; }
    int numEdges() const { return nE_; }

private:
    void buildFaceLaplacian(Eigen::SparseMatrix<double>& L) const;
    bool solveThetaGivenJumps(
        const Eigen::VectorXi& jump,
        Eigen::VectorXd& outTheta) const;
    bool optimizeAlternating(Eigen::VectorXd& theta, Eigen::VectorXi& jump);
    void refineJumpsCoordinateDescent(Eigen::VectorXd& theta, Eigen::VectorXi& jump);
    void wrapThetaToQuarterTurn(Eigen::VectorXd& theta) const;

    static int roundToInt(double x);

    int nF_ = 0;
    int nE_ = 0;
    std::vector<Edge> edges_;

    int jumpLo_ = -4;
    int jumpHi_ = 4;
    int alternatingIters_ = 12;
    int refinePasses_ = 2;

    mutable Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> faceSolver_;
    mutable bool faceSolverReady_ = false;
    mutable Eigen::SparseMatrix<double> faceLaplacian_;
};

#endif
