#ifndef CROSSFIELDINTEGERPROGRAM_H
#define CROSSFIELDINTEGERPROGRAM_H

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <vector>

/**
 * Small mixed-integer least-squares solver for pairwise difference constraints.
 *
 *   min_{x,z} Σ_c (x_i - x_j + scale * z_c)^2,   z_c ∈ Z
 *
 * This class deliberately knows nothing about application-specific geometry.
 * Callers are responsible for mapping their problem to constraints.
 */
class CrossFieldIntegerProgram {
public:
    struct Constraint {
        int i = -1;
        int j = -1;
        int idx = -1;
    };

    CrossFieldIntegerProgram(int numVariables,
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
    void buildDifferenceLaplacian(Eigen::SparseMatrix<double>& L) const;
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
    mutable Eigen::SparseMatrix<double> variableLaplacian_;
};

#endif
