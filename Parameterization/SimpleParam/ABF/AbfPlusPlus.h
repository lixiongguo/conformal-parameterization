#ifndef ABF_PLUS_PLUS_H
#define ABF_PLUS_PLUS_H

#include "Parameterization.h"
#include "AugmentedLagrangian.h"
#include <Eigen/Sparse>
#include <array>
#include <vector>
#include <memory>

class AbfPlusPlus : public Parameterization {
public:
    explicit AbfPlusPlus(Mesh& mesh0, int maxIterations = 200);
    void parameterize() override;

private:
    int m_maxIterations;

    // --- angle variable data ---
    // faceAngles_[fi][k] is the global angle index for corner k.
    std::vector<std::array<int, 3>> faceAngles_;
    std::vector<int> angleToVertex_;          // [angleIdx] → vertex index
    Eigen::VectorXd beta_;                    // target angles (3D geometry)
    int nAngles_       = 0;
    int nInteriorFaces_ = 0;

    // --- constraint data (serves AugmentedLagrangian callbacks) ---
    Eigen::SparseMatrix<double> constraintJacobian_;   // M × nAngles  (constant, linear)
    Eigen::VectorXd              constraintRhs_;       // M
    int nConstraints_ = 0;

    // --- augmented Lagrangian solver ---
    std::unique_ptr<AugmentedLagrangian> al_;

    void buildAngleMapping();
    void buildConstraints();

    // Augmented Lagrangian callbacks (ABF++ has linear constraints: Jα = rhs)
    Eigen::VectorXd evalConstraint(const Eigen::VectorXd& alpha) const;
    Eigen::SparseMatrix<double> evalJacobian(const Eigen::VectorXd& alpha) const;
    Eigen::SparseMatrix<double> evalConstraintHessian(
        const Eigen::VectorXd& alpha,
        const Eigen::VectorXd& lambda,
        double rho) const;

    void angleToUv(const Eigen::VectorXd& alpha);

    // helpers
    double boundaryCurvature(int v0, int v1, int v2) const;
};

#endif
