#ifndef BOUNDED_LSCM_H
#define BOUNDED_LSCM_H

#include "Parameterization.h"
#include <Eigen/Sparse>
#include <vector>

// ================================================================
//  BoundedLscm — LSCM with Lipman's bounded conformal distortion
//
//  Reference:
//    Lipman, Y. "Bounded Distortion Mapping Spaces for Triangular
//    Meshes", ACM Trans. Graph. (SIGGRAPH), 2012.
//
//  Approach:
//    Start from standard LSCM, then iteratively apply per-face
//    linear constraints derived from Lipman's convexification:
//      1. Half-plane constraint  |Im(α_j)| ≤ tan(θ_max)·Re(α_j)
//         where θ_max = cos⁻¹(4/5) ≈ 36.87°
//      2. Bounded-distortion constraint  |β_j| ≤ (C-1)/(C+1)·|α_j|
//         (linearised around previous solution)
//      3. Local frame ⊕_j is realigned to α_j after each iteration.
//
//  The per-iteration problem is a linearly-constrained QP solved
//  via KKT system with LDLT decomposition.
// ================================================================

class BoundedLscm : public Parameterization {
public:
    /// @param mesh0       mesh to parameterise (modified in-place)
    /// @param C            conformal distortion upper bound (C > 1)
    /// @param maxIter      maximum outer-iteration count
    BoundedLscm(Mesh& mesh0, double C = 5.0, int maxIter = 12);

    void parameterize() override;

private:
    double m_C;           // distortion bound K ≤ C
    int    m_maxIter;     // outer iterations (frame update)
    double m_kappa;       // = (C-1)/(C+1)

    // ------ per-face angle data ------
    int    nFaces_;
    std::vector<int> faceIndices_;    // interior-face indices

    // ------ local frame realignment ------
    void buildLocalFrames(const Eigen::VectorXd& uv,
                          std::vector<std::complex<double>>& frameDir) const;

    // ------ QP build per iteration ------
    // Build objective  min Σ |β_j|²·area_j  →  xᵀP x
    // plus linear constraints  Aineq·x ≤ bineq
    bool buildQP(const Eigen::VectorXd& uvPrev,
                 const std::vector<std::complex<double>>& frameDir,
                 Eigen::SparseMatrix<double>& P,
                 Eigen::SparseMatrix<double>& Aineq,
                 Eigen::VectorXd& bineq,
                 int& nVars);

    // Solve linearly-constrained QP via KKT (active-set heuristic)
    bool solveLscmQP(const Eigen::SparseMatrix<double>& P,
                     const Eigen::SparseMatrix<double>& A,
                     const Eigen::VectorXd& b,
                     Eigen::VectorXd& uv);

    // ------ helpers ------
    void extractUV(const Eigen::VectorXd& uvAll);
    void computeAlphaBeta(const Eigen::VectorXd& uv,
                          std::vector<std::complex<double>>& alpha,
                          std::vector<std::complex<double>>& beta) const;
    double maxConformalDistortion(const Eigen::VectorXd& uv) const;
};

#endif
