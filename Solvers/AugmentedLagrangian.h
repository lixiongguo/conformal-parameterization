#ifndef AUGMENTED_LAGRANGIAN_H
#define AUGMENTED_LAGRANGIAN_H

#include <Eigen/Sparse>
#include <functional>

// ============================================================
//  Augmented Lagrangian constrained optimisation solver
//
//  Problem:
//      min  f(x)
//      s.t. c(x) = 0        (M nonlinear constraints)
//           x_i > 0         (optional box constraint)
//
//  f(x) is separable quadratic:  f(x) = Σ (x_i − β_i)² / β_i
//  – gradient      ∇f  = D·(x − β)  with  D_ii = 2/β_i
//  – Hessian       ∇²f = D           (constant, positive definite)
//
//  c(x), ∇c(x), and the combined constraint‑Hessian term
//      Σ_j (λ_j + ρ·c_j) · ∇²c_j
//  are supplied by the caller via callbacks.
// ============================================================

// Moved outside the class to avoid clang's nested-struct
// default-member-initializer parsing issue with Emscripten.
struct AugmentedLagrangianConfig {
    double rhoInit   = 1e-2;          // initial penalty parameter
    double rhoMax    = 1e8;           // upper bound on ρ
    double rhoScale  = 3.16227766;    // = √10, penalty increase factor
    double tolOuter  = 1e-6;          // constraint violation tolerance
    double tolInner  = 1e-8;          // gradient norm tolerance (Newton)
    int    maxOuter  = 40;            // max AL outer iterations
    int    maxInner  = 25;            // max Newton inner iterations
    bool   positiveConstraint = true; // x_i ≥ 0 box constraint
};

class AugmentedLagrangian {
public:
    using Config = AugmentedLagrangianConfig;

    // ---------- constraint callbacks ----------

    /// Evaluate all M constraint values at x.
    using ConstraintFunc = std::function<Eigen::VectorXd(const Eigen::VectorXd& x)>;

    /// Evaluate constraint Jacobian ∇c(x)  (M × N sparse).
    using JacobianFunc = std::function<Eigen::SparseMatrix<double>(const Eigen::VectorXd& x)>;

    /// Evaluate the combined constraint‑Hessian contribution:
    ///   Σ_j (λ_j + ρ·c_j) · ∇²c_j(x)
    /// This is added to the Lagrangian Hessian ∇²L.
    /// For linear constraints (∇²c_j ≡ 0) simply return an empty matrix.
    using ConstraintHessianFunc = std::function<Eigen::SparseMatrix<double>(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& lambda,
        double rho)>;

    // ---------- constructor ----------

    /// @param N          number of variables
    /// @param M          number of constraints
    /// @param conFunc    c(x) evaluator
    /// @param jacFunc    ∇c(x) evaluator
    /// @param hessFunc   Σ (λ_j + ρ·c_j) ∇²c_j  evaluator
    /// @param cfg        solver configuration
    AugmentedLagrangian(int N,
                        int M,
                        ConstraintFunc  conFunc,
                        JacobianFunc    jacFunc,
                        ConstraintHessianFunc hessFunc,
                        const Config&   cfg = Config{});

    // ---------- solve ----------

    /// Run augmented Lagrangian optimisation.
    /// @param  x      input: initial guess; output: optimal solution
    /// @param  beta   target values  β_i > 0  for the separable objective
    /// @return true on success
    bool solve(Eigen::VectorXd& x, const Eigen::VectorXd& beta);

private:
    // dimensions
    int N_;   // variable count
    int M_;   // constraint count

    // callbacks
    ConstraintFunc         conFunc_;
    JacobianFunc           jacFunc_;
    ConstraintHessianFunc  hessFunc_;

    // configuration
    Config cfg_;

    // helpers during inner Newton loop
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_;

    // ---- internal helpers ----
    double evaluateLagrangian(const Eigen::VectorXd& x,
                              const Eigen::VectorXd& beta,
                              const Eigen::VectorXd& lambda,
                              double rho) const;

    void   evaluateGradient(const Eigen::VectorXd& x,
                            const Eigen::VectorXd& beta,
                            const Eigen::VectorXd& diagGrad,
                            const Eigen::VectorXd& lambda,
                            double rho,
                            const Eigen::VectorXd& constraintVals,
                            const Eigen::SparseMatrix<double>& jac,
                            Eigen::VectorXd& grad) const;

    bool   assembleHessian(const Eigen::VectorXd& x,
                           const Eigen::VectorXd& beta,
                           const Eigen::VectorXd& diagGrad,
                           const Eigen::VectorXd& lambda,
                           double rho,
                           const Eigen::VectorXd& constraintVals,
                           const Eigen::SparseMatrix<double>& jac);
};

#endif
