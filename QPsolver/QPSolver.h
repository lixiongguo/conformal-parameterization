#ifndef QPSOLVER_QPSOLVER_H
#define QPSOLVER_QPSOLVER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace QPsolver {

inline double clampd(double v, double lo, double hi) {
    return std::min(std::max(v, lo), hi);
}

inline Eigen::VectorXd projectBox(const Eigen::VectorXd& v,
                                    const Eigen::VectorXd& lo,
                                    const Eigen::VectorXd& hi) {
    Eigen::VectorXd out = v;
    out = out.cwiseMax(lo).cwiseMin(hi);
    return out;
}

// Problem:
//   min 0.5 * x^T * I * x
//   s.t. lowerY <= A * x <= upperY
//        lowerX <= x <= upperX
struct MinNormBoxQPProblem {
    int n = 0; // variables
    int m = 0; // constraints on y = A*x

    // A is sparse m x n
    Eigen::SparseMatrix<double> A;

    Eigen::VectorXd lowerX;
    Eigen::VectorXd upperX;

    Eigen::VectorXd lowerY;
    Eigen::VectorXd upperY;
};

struct MinNormBoxQPOptions {
    int admmMaxIters = 200;
    int cgMaxIters = 80; // per ADMM iteration

    // Augmented Lagrangian penalties
    double rhoY = 1.0;
    double rhoX = 1.0;

    // CG stopping
    double cgRelTol = 1e-6;

    // Simple feasibility tolerance for early exit
    double feasTol = 1e-5;
};

class QPSolver {
public:
    // Returns min-norm solution x, clamped within variable bounds.
    static bool solveMinNormBoxQPIdentityQ(const MinNormBoxQPProblem& prob,
                                            Eigen::VectorXd& outX,
                                            const MinNormBoxQPOptions& opts = MinNormBoxQPOptions()) {
        if (prob.n <= 0 || prob.m < 0) return false;
        if (prob.A.rows() != prob.m || prob.A.cols() != prob.n) return false;

        if (prob.lowerX.size() != prob.n || prob.upperX.size() != prob.n) return false;
        if (prob.lowerY.size() != prob.m || prob.upperY.size() != prob.m) return false;

        // Basic bound sanity.
        for (int i = 0; i < prob.n; ++i) {
            if (!(prob.lowerX[i] <= prob.upperX[i])) return false;
        }
        for (int j = 0; j < prob.m; ++j) {
            if (!(prob.lowerY[j] <= prob.upperY[j])) return false;
        }

        const int n = prob.n;
        const int m = prob.m;

        // Splitting variables:
        //   s = x, enforced by projection onto [lowerX, upperX]
        //   y = A*x, enforced by projection onto [lowerY, upperY]
        Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd s = projectBox(x, prob.lowerX, prob.upperX);

        Eigen::VectorXd Ax = prob.A * x;
        Eigen::VectorXd y = projectBox(Ax, prob.lowerY, prob.upperY);

        // Scaled dual variables (ADMM)
        Eigen::VectorXd uX = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd uY = Eigen::VectorXd::Zero(m);

        // Precompute diagonal for Jacobi preconditioning of M = (1+rhoX)I + rhoY*A^T*A
        Eigen::VectorXd diagA2 = Eigen::VectorXd::Zero(n); // diag(A^T A)
        for (int col = 0; col < prob.A.outerSize(); ++col) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(prob.A, col); it; ++it) {
                const double v = it.value();
                diagA2[it.col()] += v * v;
            }
        }
        Eigen::VectorXd diagM = (1.0 + opts.rhoX) * Eigen::VectorXd::Ones(n) + opts.rhoY * diagA2;
        // Prevent division by zero.
        for (int i = 0; i < n; ++i) {
            if (diagM[i] <= 0.0) diagM[i] = 1e-12;
        }

        auto matvecM = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
            Eigen::VectorXd Av = prob.A * v;               // m
            Eigen::VectorXd AtAv = prob.A.transpose() * Av; // n
            return (1.0 + opts.rhoX) * v + opts.rhoY * AtAv;
        };

        auto cgSolve = [&](const Eigen::VectorXd& rhs,
                            const Eigen::VectorXd& x0) -> Eigen::VectorXd {
            Eigen::VectorXd xk = x0;
            Eigen::VectorXd r = rhs - matvecM(xk);
            Eigen::VectorXd z = r.cwiseQuotient(diagM);
            Eigen::VectorXd p = z;
            double rz = r.dot(z);

            const double rhsNorm = rhs.norm();
            const double tol = std::max(opts.cgRelTol * rhsNorm, 1e-14);

            for (int it = 0; it < opts.cgMaxIters; ++it) {
                Eigen::VectorXd Ap = matvecM(p);
                const double pAp = p.dot(Ap);
                if (!(pAp > 0.0)) break; // SPD expected; stop if numerical issues.

                const double alpha = rz / pAp;
                xk = xk + alpha * p;
                r = r - alpha * Ap;

                if (r.norm() <= tol) break;

                z = r.cwiseQuotient(diagM);
                const double rzNew = r.dot(z);
                const double beta = rzNew / rz;
                p = z + beta * p;
                rz = rzNew;
            }
            return xk;
        };

        for (int iter = 0; iter < opts.admmMaxIters; ++iter) {
            // x-update: solve Mx = rhoY * A^T (y - uY) + rhoX * (s - uX)
            Eigen::VectorXd yMinusU = y - uY;
            Eigen::VectorXd rhs = opts.rhoY * (prob.A.transpose() * yMinusU) + opts.rhoX * (s - uX);

            x = cgSolve(rhs, x);

            // s-update (box projection)
            s = projectBox(x + uX, prob.lowerX, prob.upperX);

            // y-update (box projection)
            Eigen::VectorXd Ax2 = prob.A * x;
            y = projectBox(Ax2 + uY, prob.lowerY, prob.upperY);

            // Dual updates
            uX = uX + (x - s);
            uY = uY + (Ax2 - y);

            // Early feasibility check using current s and y.
            // (This is a heuristic; ADMM may still be improving.)
            const Eigen::VectorXd As = prob.A * s;
            double maxVx = 0.0;
            for (int i = 0; i < n; ++i) {
                const double below = prob.lowerX[i] - s[i];
                const double above = s[i] - prob.upperX[i];
                maxVx = std::max(maxVx, std::max(below, above));
            }
            double maxVy = 0.0;
            for (int j = 0; j < m; ++j) {
                const double below = prob.lowerY[j] - As[j];
                const double above = As[j] - prob.upperY[j];
                maxVy = std::max(maxVy, std::max(below, above));
            }

            if (maxVx <= opts.feasTol && maxVy <= opts.feasTol) {
                outX = s;
                return true;
            }
        }

        outX = projectBox(s, prob.lowerX, prob.upperX);

        // Final feasibility check (strict enough to avoid silent failures).
        const Eigen::VectorXd As = prob.A * outX;
        double maxVy = 0.0;
        for (int j = 0; j < m; ++j) {
            const double below = prob.lowerY[j] - As[j];
            const double above = As[j] - prob.upperY[j];
            maxVy = std::max(maxVy, std::max(below, above));
        }

        return maxVy <= (10.0 * opts.feasTol);
    }
};

} // namespace QPsolver

#endif // QPSOLVER_QPSOLVER_H

