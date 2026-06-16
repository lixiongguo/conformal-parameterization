#ifndef MIXEDINTEGERPROGRAM_H
#define MIXEDINTEGERPROGRAM_H

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <vector>

/**
 * Small mixed-integer least-squares solver for pairwise difference constraints.
 *
 *   min_{x,z} sum_c (x_i - x_j + kappa_c + scale * z_c)^2,   z_c in Z
 *
 * kappa_c is the base rotation between the local coordinate frames of face i
 * and face j.  It is a fixed constant per constraint — callers must compute
 * it from geometry (e.g. the angle of the shared edge expressed in each face's
 * tangent frame).
 *
 * This class deliberately knows nothing about application-specific geometry.
 * Callers are responsible for mapping their problem to constraints.
 */
class MixedIntegerProgram {
public:
    struct Constraint {
        int i = -1;
        int j = -1;
        int idx = -1;
        double kappa = 0.0;   ///< base rotation  i→j  (θ_i + κ_ij − θ_j + s·z ≈ 0)
    };

    MixedIntegerProgram(int numVariables,
                        std::vector<Constraint> constraints,
                        double integerScale);

    void setIntegerBounds(int lo, int hi);
    void setAlternatingIterations(int iters) { alternatingIters_ = iters; }
    void setRefinePasses(int passes) { refinePasses_ = passes; }

    /** Run alternating + refinement; updates variables and integers in place. */
    bool solve(Eigen::VectorXd& variables, Eigen::VectorXi& integers);

    double energy(const Eigen::VectorXd& variables, const Eigen::VectorXi& integers) const;

    int numVariables() const { return nVars_; }
    int numConstraints() const { return nConstraints_; }

private:
    void buildNormalMatrix(Eigen::SparseMatrix<double>& A) const;
    bool solveVariablesGivenIntegers(
        const Eigen::VectorXi& integers,
        Eigen::VectorXd& outVariables) const;
    bool optimizeAlternating(Eigen::VectorXd& variables, Eigen::VectorXi& integers);
    void refineIntegersCoordinateDescent(Eigen::VectorXd& variables, Eigen::VectorXi& integers);

    static int roundToInt(double x);

    int nVars_ = 0;
    int nConstraints_ = 0;
    double integerScale_ = 1.0;
    std::vector<Constraint> constraints_;

    int integerLo_ = -4;
    int integerHi_ = 4;
    int alternatingIters_ = 12;
    int refinePasses_ = 2;

    mutable Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> variableSolver_;
    mutable bool variableSolverReady_ = false;
    mutable Eigen::SparseMatrix<double> normalMatrix_;
};

#endif
