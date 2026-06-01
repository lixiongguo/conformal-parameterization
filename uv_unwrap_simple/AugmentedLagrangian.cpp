#include "AugmentedLagrangian.h"

// ============================================================
//  Constructor
// ============================================================

AugmentedLagrangian::AugmentedLagrangian(int N,
                                         int M,
                                         ConstraintFunc  conFunc,
                                         JacobianFunc    jacFunc,
                                         ConstraintHessianFunc hessFunc,
                                         const Config&   cfg)
    : N_(N)
    , M_(M)
    , conFunc_(std::move(conFunc))
    , jacFunc_(std::move(jacFunc))
    , hessFunc_(std::move(hessFunc))
    , cfg_(cfg)
{}


// ============================================================
//  Evaluate augmented Lagrangian:  L(x) = f(x) + λᵀc + ρ/2‖c‖²
// ============================================================

double AugmentedLagrangian::evaluateLagrangian(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& lambda,
    double rho) const
{
    // objective  f = Σ (x_i − β_i)² / β_i
    double L = 0.0;
    for (int i = 0; i < N_; i++) {
        double diff = x(i) - beta(i);
        L += diff * diff / beta(i);
    }

    // constraints
    Eigen::VectorXd c = conFunc_(x);
    L += lambda.dot(c) + 0.5 * rho * c.squaredNorm();

    return L;
}


// ============================================================
//  Evaluate gradient  ∇L = ∇f + ∇cᵀ (λ + ρ·c)
// ============================================================

void AugmentedLagrangian::evaluateGradient(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& diagGrad,          // D_ii = 2/β_i
    const Eigen::VectorXd& lambda,
    double rho,
    const Eigen::VectorXd& constraintVals,
    const Eigen::SparseMatrix<double>& jac,
    Eigen::VectorXd& grad) const
{
    // ∇f = D·(x − β)
    grad.resize(N_);
    for (int i = 0; i < N_; i++)
        grad(i) = diagGrad(i) * (x(i) - beta(i));

    // ∇cᵀ · (λ + ρ·c)
    Eigen::VectorXd rhs = lambda + rho * constraintVals;
    grad += jac.transpose() * rhs;
}


// ============================================================
//  Assemble full Hessian  ∇²L = D + ρ·∇cᵀ·∇c + Σ(λ_j+ρc_j)∇²c_j
//  and factorise it with SimplicialLDLT.
// ============================================================

bool AugmentedLagrangian::assembleHessian(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& diagGrad,          // D_ii = 2/β_i
    const Eigen::VectorXd& lambda,
    double rho,
    const Eigen::VectorXd& constraintVals,
    const Eigen::SparseMatrix<double>& jac)
{
    const int N = N_;

    // ----  ∇²f = D  (constant diagonal)  ----
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(N + jac.nonZeros() * 2);

    for (int i = 0; i < N; i++)
        triplets.emplace_back(i, i, diagGrad(i));

    // ----  ρ · ∇cᵀ·∇c  (Gauss‑Newton term)  ----
    Eigen::SparseMatrix<double> JtJ = jac.transpose() * jac;
    for (int k = 0; k < JtJ.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(JtJ, k); it; ++it)
            triplets.emplace_back(it.row(), it.col(), rho * it.value());
    }

    // ----  Σ (λ_j + ρ·c_j) ∇²c_j  (nonlinear Hessian contribution)  ----
    // This is the critical term: handles genuinely nonlinear constraints.
    // For linear constraints ∇²c_j = 0, hessFunc_ returns an empty matrix.
    Eigen::SparseMatrix<double> conHess = hessFunc_(x, lambda, rho);
    for (int k = 0; k < conHess.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(conHess, k); it; ++it)
            triplets.emplace_back(it.row(), it.col(), it.value());
    }

    // Assemble and factorise
    Eigen::SparseMatrix<double> H(N, N);
    H.setFromTriplets(triplets.begin(), triplets.end());

    solver_.compute(H);
    return (solver_.info() == Eigen::Success);
}


// ============================================================
//  solve() – main augmented Lagrangian loop
// ============================================================

bool AugmentedLagrangian::solve(Eigen::VectorXd& x, const Eigen::VectorXd& beta)
{
    const int N = N_;
    const int M = M_;

    // ---- precompute constant D = diag(2/β_i) for gradient & Hessian ----
    Eigen::VectorXd diagGrad(N);
    for (int i = 0; i < N; i++)
        diagGrad(i) = 2.0 / beta(i);

    // ---- AL state ----
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(M);
    double rho = cfg_.rhoInit;

    // ---- outer loop (multiplier / penalty update) ----
    for (int outer = 0; outer < cfg_.maxOuter; outer++) {

        // evaluate constraints & Jacobian at current x
        Eigen::VectorXd cVal           = conFunc_(x);
        Eigen::SparseMatrix<double> JC = jacFunc_(x);

        // assemble & factorise Hessian for this outer iteration
        if (!assembleHessian(x, beta, diagGrad, lambda, rho, cVal, JC))
            return false;

        // ---- inner loop (Newton on L₍ρ₎( · ; λ)) ----
        for (int inner = 0; inner < cfg_.maxInner; inner++) {

            // re‑evaluate constraints at current x
            cVal = conFunc_(x);
            JC   = jacFunc_(x);

            // gradient
            Eigen::VectorXd grad;
            evaluateGradient(x, beta, diagGrad, lambda, rho, cVal, JC, grad);

            if (grad.lpNorm<Eigen::Infinity>() < cfg_.tolInner)
                break;

            // Newton direction
            Eigen::VectorXd delta = solver_.solve(-grad);
            if (solver_.info() != Eigen::Success)
                return false;

            // ---- box constraint: x_i > 0 → clip step length ----
            double maxT = 1.0;
            if (cfg_.positiveConstraint) {
                for (int i = 0; i < N; i++) {
                    if (delta(i) < 0.0) {
                        double ti = -0.95 * x(i) / delta(i);
                        if (ti < maxT) maxT = ti;
                    }
                }
            }
            double t = std::min(1.0, maxT);

            // ---- Armijo backtracking line search ----
            double curL     = evaluateLagrangian(x, beta, lambda, rho);
            double gradDotD = grad.dot(delta);
            const double armijoC = 1e-3;

            bool accepted = false;
            while (t > 1e-10) {
                Eigen::VectorXd xTry = x + t * delta;
                double tryL = evaluateLagrangian(xTry, beta, lambda, rho);

                if (tryL <= curL + armijoC * t * gradDotD) {
                    x        = xTry;
                    accepted = true;
                    break;
                }
                t *= 0.5;
            }

            if (!accepted || t <= 1e-10)
                break;       // Newton converged (step too small)
        }

        // ---- constraint convergence check ----
        cVal = conFunc_(x);
        double maxViolation = cVal.lpNorm<Eigen::Infinity>();
        if (maxViolation < cfg_.tolOuter)
            break;

        // ---- multiplier update ----
        lambda += rho * cVal;

        // ---- penalty increase ----
        if (rho < cfg_.rhoMax)
            rho = std::min(rho * cfg_.rhoScale, cfg_.rhoMax);

        // After multiplier update, re‑factorise Hessian with new ρ
        JC = jacFunc_(x);
        if (!assembleHessian(x, beta, diagGrad, lambda, rho, cVal, JC))
            return false;
    }

    return true;
}
