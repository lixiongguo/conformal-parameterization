#define _USE_MATH_DEFINES
#include "FastHGPNumerics.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace FastHGPNumerics {

namespace {

constexpr double kTol = 1e-10;
constexpr double kEps = 1e-14;

double safeAbs(double v)
{
    return std::max(std::abs(v), kEps);
}

Eigen::VectorXcd safeDiv(const Eigen::VectorXcd& num, const Eigen::VectorXd& den)
{
    Eigen::VectorXcd out(num.size());
    for (int i = 0; i < num.size(); ++i) {
        out(i) = num(i) / safeAbs(den(i));
    }
    return out;
}

// Active-set QP: min 0.5||p - target||^2 s.t. A p >= b (2D).
Eigen::Vector2d solveHalfPlaneQP(const Eigen::Vector2d& target,
                                 const Eigen::MatrixXd& A,
                                 const Eigen::VectorXd& b)
{
    const int m = static_cast<int>(A.rows());
    if (m == 0) {
        return target;
    }

    Eigen::Vector2d p = target;
    std::vector<int> active;
    const int maxIter = m + 4;

    for (int iter = 0; iter < maxIter; ++iter) {
        bool feasible = true;
        for (int i = 0; i < m; ++i) {
            if (A.row(i).dot(p) < b(i) - 1e-12) {
                feasible = false;
                break;
            }
        }
        if (feasible) {
            return p;
        }

        int worst = -1;
        double worstViol = 0.0;
        for (int i = 0; i < m; ++i) {
            const double viol = b(i) - A.row(i).dot(p);
            if (viol > worstViol) {
                worstViol = viol;
                worst = i;
            }
        }
        if (worst < 0) {
            return p;
        }

        if (std::find(active.begin(), active.end(), worst) == active.end()) {
            active.push_back(worst);
        }

        if (active.size() == 1) {
            const Eigen::Vector2d n = A.row(active[0]).transpose();
            const double nn = n.squaredNorm();
            if (nn < kEps) {
                return target;
            }
            const double t = (b(active[0]) - n.dot(target)) / nn;
            p = target + t * n;
            continue;
        }

        if (active.size() >= 2) {
            const int i0 = active[active.size() - 2];
            const int i1 = active[active.size() - 1];
            Eigen::Matrix2d M;
            M.row(0) = A.row(i0);
            M.row(1) = A.row(i1);
            Eigen::Vector2d rhs(b(i0), b(i1));
            if (std::abs(M.determinant()) < kEps) {
                const Eigen::Vector2d n = A.row(worst).transpose();
                const double nn = n.squaredNorm();
                if (nn > kEps) {
                    const double t = (b(worst) - n.dot(target)) / nn;
                    p = target + t * n;
                }
            } else {
                p = M.fullPivLu().solve(rhs);
            }
        }
    }

    bool allFeasible = true;
    for (int i = 0; i < m; ++i) {
        if (A.row(i).dot(p) < b(i) - 1e-8) {
            allFeasible = false;
            break;
        }
    }
    return allFeasible ? p : target;
}

} // namespace

void computeLocalBasis(const Eigen::MatrixXd& v1,
                       const Eigen::MatrixXd& v2,
                       Eigen::MatrixXd& e1,
                       Eigen::MatrixXd& e2)
{
    const int n = static_cast<int>(v1.rows());
    e1.resize(n, 3);
    e2.resize(n, 3);
    for (int i = 0; i < n; ++i) {
        const double len1 = v1.row(i).norm();
        e1.row(i) = (len1 > kEps) ? (v1.row(i) / len1) : Eigen::RowVector3d::Zero();
        const double proj = e1.row(i).dot(v2.row(i));
        Eigen::RowVector3d e2row = v2.row(i) - proj * e1.row(i);
        const double len2 = e2row.norm();
        e2.row(i) = (len2 > kEps) ? (e2row / len2) : Eigen::RowVector3d::Zero();
    }
}

PerpsResult computePerps(const Eigen::MatrixXd& V,
                         const Eigen::MatrixXi& F,
                         int xlb)
{
    const int numF = static_cast<int>(F.rows());
    PerpsResult result;
    result.tc.resize(numF, 3);
    result.dblAc.resize(numF);

    Eigen::MatrixXd v1(numF, 3), v2(numF, 3), v3(numF, 3);
    for (int f = 0; f < numF; ++f) {
        const int i1 = F(f, 0) - 1;
        const int i2 = F(f, 1) - 1;
        const int i3 = F(f, 2) - 1;
        v1.row(f) = V.row(i3) - V.row(i2);
        v2.row(f) = V.row(i1) - V.row(i3);
        v3.row(f) = V.row(i2) - V.row(i1);
    }

    Eigen::MatrixXd vlb1 = v1, vlb2 = v2;
    if (xlb == 2) {
        vlb1 = v2;
        vlb2 = v3;
    } else if (xlb == 3) {
        vlb1 = v3;
        vlb2 = v1;
    }

    Eigen::MatrixXd e1, e2;
    computeLocalBasis(vlb1, vlb2, e1, e2);

    for (int f = 0; f < numF; ++f) {
        auto dotC = [&](const Eigen::RowVector3d& a, const Eigen::RowVector3d& b) {
            return Complex(a.dot(b), 0.0);
        };
        auto dotImag = [&](const Eigen::RowVector3d& e, const Eigen::RowVector3d& v) {
            return Complex(0.0, e.dot(v));
        };

        Complex cv1 = dotC(e1.row(f), v1.row(f)) + dotImag(e2.row(f), v1.row(f));
        Complex cv2 = dotC(e1.row(f), v2.row(f)) + dotImag(e2.row(f), v2.row(f));
        Complex cv3 = -cv1 - cv2;

        result.dblAc(f) = std::abs((std::conj(cv1) * cv2).imag());

        Complex t1r = Complex(0.0, -1.0) * cv1;
        Complex t2r = Complex(0.0, -1.0) * cv2;
        Complex t3r = Complex(0.0, -1.0) * cv3;
        result.tc(f, 0) = t1r;
        result.tc(f, 1) = t2r;
        result.tc(f, 2) = t3r;
    }

    return result;
}

void createJmatrix(const Eigen::SparseMatrix<Complex>& harmonicBasis,
                   const Eigen::MatrixXd& V,
                   const Eigen::MatrixXi& F_he,
                   Eigen::MatrixXcd& J_fz,
                   Eigen::MatrixXcd& J_fbz,
                   Eigen::VectorXd& Area)
{
    const int numF = static_cast<int>(F_he.rows());
    const int numDOF = static_cast<int>(harmonicBasis.cols());

    PerpsResult perps = computePerps(V, F_he, 3);

    J_fz.resize(numF, numDOF);
    J_fbz.resize(numF, numDOF);

    auto harmonicRow = [&harmonicBasis](int row) {
        Eigen::RowVectorXcd hb(harmonicBasis.cols());
        hb.setZero();
        for (Eigen::SparseMatrix<Complex>::InnerIterator it(harmonicBasis, row); it; ++it) {
            hb(it.col()) = it.value();
        }
        return hb;
    };

    for (int f = 0; f < numF; ++f) {
        const double cf = -0.5 / std::max(perps.dblAc(f), kEps);
        Eigen::RowVectorXcd hb1 = harmonicRow(F_he(f, 0) - 1);
        Eigen::RowVectorXcd hb2 = harmonicRow(F_he(f, 1) - 1);
        Eigen::RowVectorXcd hb3 = harmonicRow(F_he(f, 2) - 1);

        for (int j = 0; j < numDOF; ++j) {
            const Complex jfz = cf * (std::conj(perps.tc(f, 0)) * hb1(j)
                                    + std::conj(perps.tc(f, 1)) * hb2(j)
                                    + std::conj(perps.tc(f, 2)) * hb3(j));
            const Complex jfbz = cf * (perps.tc(f, 0) * hb1(j)
                                     + perps.tc(f, 1) * hb2(j)
                                     + perps.tc(f, 2) * hb3(j));
            J_fz(f, j) = jfz;
            J_fbz(f, j) = jfbz;
        }
    }

    Area = perps.dblAc / 2.0;
    const double sumA = Area.sum();
    if (sumA > kEps) {
        Area /= sumA;
    }
}

void localStep(const Eigen::VectorXcd& frames,
               double sigma2Eps,
               double k,
               const Eigen::VectorXcd& fz,
               const Eigen::VectorXcd& fbz,
               Eigen::VectorXcd& fz_local,
               Eigen::VectorXcd& fbz_local)
{
    const int n = static_cast<int>(fz.size());
    Eigen::VectorXd y(n), x(n);
    for (int i = 0; i < n; ++i) {
        y(i) = (fz(i) * frames(i)).real();
        x(i) = std::abs(fbz(i));
    }

    Eigen::VectorXd y1 = y;
    Eigen::VectorXd x1 = x;

    const double eps12 = sigma2Eps * ((1.0 + k * k) / (1.0 - k));
    const double eps34 = sigma2Eps * (1.0 + k) / (1.0 - k);

    for (int i = 0; i < n; ++i) {
        const bool I2 = (y(i) < x(i) / k) && (y(i) > -k * x(i) + eps12);
        const bool I3 = (y(i) <= -k * x(i) + eps12) && (y(i) >= -x(i) + eps34);
        const bool I4 = (y(i) < -x(i) + eps34) && (y(i) > -x(i) + sigma2Eps) && (y(i) < sigma2Eps + x(i));
        const bool I5 = (y(i) <= -x(i) + sigma2Eps);

        if (I2) {
            y1(i) = (y(i) + k * x(i)) / (1.0 + k * k);
            x1(i) = (k * y(i) + k * k * x(i)) / (1.0 + k * k);
        } else if (I3) {
            y1(i) = sigma2Eps / (1.0 - k);
            x1(i) = k * sigma2Eps / (1.0 - k);
        } else if (I4) {
            y1(i) = (y(i) + x(i) + sigma2Eps) / 2.0;
            x1(i) = (y(i) + x(i) - sigma2Eps) / 2.0;
        } else if (I5) {
            y1(i) = sigma2Eps;
            x1(i) = 0.0;
        }
    }

    fz_local.resize(n);
    fbz_local.resize(n);
    for (int i = 0; i < n; ++i) {
        if (std::abs(fbz(i)) > kEps) {
            fbz_local(i) = x1(i) * fbz(i) / std::abs(fbz(i));
        } else {
            fbz_local(i) = Complex(0.0, 0.0);
        }
        const Complex fzFrame = fz(i) * frames(i);
        fz_local(i) = Complex(y1(i), fzFrame.imag()) / frames(i);
    }
}

void globalStepMAP(const Eigen::VectorXcd& a_l_in,
                   const Eigen::MatrixXcd& MInvJtrans,
                   const Eigen::MatrixXcd& JWithOnes,
                   Eigen::VectorXcd& c,
                   Eigen::VectorXcd& a_l_out)
{
    c = MInvJtrans * a_l_in;
    a_l_out = JWithOnes * c;
    a_l_out.conservativeResize(a_l_out.size() - 1);
}

void globalStepATP(const Eigen::MatrixXcd& MInvJtrans,
                   const Eigen::MatrixXcd& JWithOnes,
                   const Eigen::VectorXcd& c_i,
                   const Eigen::VectorXcd& ni,
                   double norm_ni_sqr,
                   const Eigen::VectorXd& W,
                   Eigen::VectorXcd& c_i_p_1,
                   Eigen::VectorXcd& a_l)
{
    const Eigen::VectorXcd M_J_ni = MInvJtrans * ni;
    Eigen::RowVectorXcd row = ni.adjoint().cwiseProduct(W.transpose().cast<double>());
    const Complex denom = row * JWithOnes * M_J_ni;
    const double denomReal = std::max(std::abs(denom), kEps);
    const double scalar = norm_ni_sqr / denomReal;
    c_i_p_1 = c_i - scalar * M_J_ni;
    a_l = JWithOnes * c_i_p_1;
    a_l.conservativeResize(a_l.size() - 1);
}

double symDirEnergyByfzfbz(const Eigen::VectorXcd& fz,
                           const Eigen::VectorXcd& fbz,
                           const Eigen::VectorXd& Area)
{
    double E = 0.0;
    for (int i = 0; i < fz.size(); ++i) {
        const double fz2 = std::norm(fz(i));
        const double fbz2 = std::norm(fbz(i));
        const double g1 = fz2 + fbz2;
        const double g2 = fz2 - fbz2;
        const double g2sq = std::max(g2 * g2, kEps);
        E += Area(i) * (g1 + g1 / g2sq);
    }
    return E;
}

Eigen::VectorXcd expandReducedX(const Eigen::VectorXd& xSmall,
                              int n,
                              const std::vector<int>& fixedIndices,
                              const Eigen::VectorXd& fixedValues)
{
    Eigen::VectorXd full(2 * n);
    full.setZero();
    std::vector<bool> isFixed(2 * n, false);
    for (int idx : fixedIndices) {
        if (idx >= 0 && idx < 2 * n) {
            isFixed[idx] = true;
        }
    }
    int smallIdx = 0;
    for (int i = 0; i < 2 * n; ++i) {
        if (!isFixed[i]) {
            full(i) = xSmall(smallIdx++);
        }
    }
    for (size_t k = 0; k < fixedIndices.size(); ++k) {
        const int fi = fixedIndices[k];
        if (fi >= 0 && fi < 2 * n) {
            full(fi) = fixedValues(k);
        }
    }

    Eigen::VectorXcd xComplex(n);
    for (int i = 0; i < n; ++i) {
        xComplex(i) = Complex(full(i), full(i + n));
    }
    return xComplex;
}

Eigen::VectorXd packReducedX(const Eigen::VectorXd& fullX,
                             const std::vector<int>& fixedIndices)
{
    std::vector<bool> isFixed(static_cast<size_t>(fullX.size()), false);
    for (int idx : fixedIndices) {
        if (idx >= 0 && idx < fullX.size()) {
            isFixed[static_cast<size_t>(idx)] = true;
        }
    }
    int count = 0;
    for (bool f : isFixed) {
        if (!f) ++count;
    }
    Eigen::VectorXd small(count);
    int j = 0;
    for (int i = 0; i < fullX.size(); ++i) {
        if (!isFixed[static_cast<size_t>(i)]) {
            small(j++) = fullX(i);
        }
    }
    return small;
}

SymDirEnergyResult symDirEnergyGradHess(const Eigen::VectorXd& xSmall,
                                        const Eigen::MatrixXcd& J_fz,
                                        const Eigen::MatrixXcd& J_fbz,
                                        const Eigen::VectorXd& Area,
                                        const std::vector<int>& fixedIndices,
                                        const Eigen::VectorXd& fixedValues,
                                        bool spdHessian)
{
    const int n = static_cast<int>(J_fz.cols());
    SymDirEnergyResult out;
    const Eigen::VectorXcd xComplex = expandReducedX(xSmall, n, fixedIndices, fixedValues);

    out.fz = J_fz * xComplex;
    out.fbz = J_fbz * xComplex;

    const int numF = static_cast<int>(out.fz.size());
    Eigen::VectorXd abs_fz(numF), abs_fbz(numF), SIGMA(numF), sigma(numF);
    for (int i = 0; i < numF; ++i) {
        abs_fz(i) = std::abs(out.fz(i));
        abs_fbz(i) = std::abs(out.fbz(i));
        SIGMA(i) = abs_fz(i) + abs_fbz(i);
        sigma(i) = abs_fz(i) - abs_fbz(i);
    }

    out.energy = 0.0;
    for (int i = 0; i < numF; ++i) {
        const double s = std::max(SIGMA(i), kEps);
        const double sig = std::max(std::abs(sigma(i)), kEps);
        out.energy += 0.5 * Area(i) * (s * s + 1.0 / (s * s) + sig * sig + 1.0 / (sig * sig));
    }

    Eigen::VectorXd grad_h_SIGMA(numF), grad_h_sigma(numF);
    for (int i = 0; i < numF; ++i) {
        const double s = std::max(SIGMA(i), kEps);
        const double sig = std::max(std::abs(sigma(i)), kEps);
        grad_h_SIGMA(i) = Area(i) * ((std::pow(s, 4) - 1.0) / std::pow(s, 3));
        grad_h_sigma(i) = Area(i) * ((std::pow(sig, 4) - 1.0) / std::pow(sig, 3));
    }

    Eigen::MatrixXcd grad_abs_fz(numF, n);
    Eigen::MatrixXcd grad_abs_fbz(numF, n);
    for (int i = 0; i < numF; ++i) {
        const Complex scaleFz = out.fz(i) / safeAbs(abs_fz(i));
        const Complex scaleFbz = out.fbz(i) / safeAbs(abs_fbz(i));
        for (int j = 0; j < n; ++j) {
            grad_abs_fz(i, j) = std::conj(J_fz(i, j)) * scaleFz;
            grad_abs_fbz(i, j) = std::conj(J_fbz(i, j)) * scaleFbz;
        }
    }

    Eigen::MatrixXcd grad_SIGMA_x = grad_abs_fz + grad_abs_fbz;
    Eigen::MatrixXcd grad_sigma_x = grad_abs_fz - grad_abs_fbz;

    Eigen::VectorXcd G_complex = grad_SIGMA_x.adjoint() * grad_h_SIGMA
                               + grad_sigma_x.adjoint() * grad_h_sigma;

    Eigen::VectorXd Gfull(2 * n);
    for (int i = 0; i < n; ++i) {
        Gfull(i) = G_complex(i).real();
        Gfull(i + n) = G_complex(i).imag();
    }
    out.gradient = packReducedX(Gfull, fixedIndices);

    const int nFree = static_cast<int>(out.gradient.size());
    out.hessian = Eigen::MatrixXd::Zero(nFree, nFree);

    Eigen::MatrixXd grad_SIGMA_x_real(2 * n, numF);
    Eigen::MatrixXd grad_sigma_x_real(2 * n, numF);
    Eigen::MatrixXd grad_abs_fz_real(2 * n, numF);
    Eigen::MatrixXd grad_abs_fbz_real(2 * n, numF);

    for (int i = 0; i < numF; ++i) {
        for (int j = 0; j < n; ++j) {
            grad_SIGMA_x_real(j, i) = grad_SIGMA_x(i, j).real();
            grad_SIGMA_x_real(j + n, i) = grad_SIGMA_x(i, j).imag();
            grad_sigma_x_real(j, i) = grad_sigma_x(i, j).real();
            grad_sigma_x_real(j + n, i) = grad_sigma_x(i, j).imag();
            grad_abs_fz_real(j, i) = grad_abs_fz(i, j).real();
            grad_abs_fz_real(j + n, i) = grad_abs_fz(i, j).imag();
            grad_abs_fbz_real(j, i) = grad_abs_fbz(i, j).real();
            grad_abs_fbz_real(j + n, i) = grad_abs_fbz(i, j).imag();
        }
    }

    Eigen::VectorXd hess_h11 = Area.array() * (1.0 + 3.0 * SIGMA.array().pow(-4));
    Eigen::VectorXd hess_h22 = Area.array() * (1.0 + 3.0 * sigma.array().pow(-4));

    Eigen::MatrixXd H1 = Eigen::MatrixXd::Zero(2 * n, 2 * n);
    H1 += grad_SIGMA_x_real * hess_h11.asDiagonal() * grad_SIGMA_x_real.transpose();
    H1 += grad_sigma_x_real * hess_h22.asDiagonal() * grad_sigma_x_real.transpose();

    Eigen::MatrixXd Ar1(2 * n, numF), Ar2(2 * n, numF);
    Eigen::MatrixXd Br1(2 * n, numF), Br2(2 * n, numF);
    for (int i = 0; i < numF; ++i) {
        for (int j = 0; j < n; ++j) {
            const Complex a = J_fz(i, j);
            const Complex b = J_fbz(i, j);
            Ar1(j, i) = a.real();
            Ar1(j + n, i) = a.imag();
            Ar2(j, i) = -a.imag();
            Ar2(j + n, i) = a.real();
            Br1(j, i) = b.real();
            Br1(j + n, i) = b.imag();
            Br2(j, i) = -b.imag();
            Br2(j + n, i) = b.real();
        }
    }

    Eigen::VectorXd alpha = grad_h_SIGMA + grad_h_sigma;
    Eigen::VectorXd beta = grad_h_SIGMA - grad_h_sigma;
    if (spdHessian) {
        for (int i = 0; i < alpha.size(); ++i) {
            alpha(i) = std::max(alpha(i), 0.0);
        }
    }

    Eigen::MatrixXd H2 = Eigen::MatrixXd::Zero(2 * n, 2 * n);
    for (int i = 0; i < numF; ++i) {
        const double sA = alpha(i) / safeAbs(abs_fz(i));
        const double sB = beta(i) / safeAbs(abs_fbz(i));
        H2 += sA * (Ar1.col(i) * Ar1.col(i).transpose() + Ar2.col(i) * Ar2.col(i).transpose());
        H2 -= sA * (grad_abs_fz_real.col(i) * grad_abs_fz_real.col(i).transpose());
        H2 += sB * (Br1.col(i) * Br1.col(i).transpose() + Br2.col(i) * Br2.col(i).transpose());
        H2 -= sB * (grad_abs_fbz_real.col(i) * grad_abs_fbz_real.col(i).transpose());
    }

    const Eigen::MatrixXd Hfull = H1 + H2;

    std::vector<int> freeIdx;
    std::vector<bool> isFixed(2 * n, false);
    for (int idx : fixedIndices) {
        if (idx >= 0 && idx < 2 * n) isFixed[idx] = true;
    }
    for (int i = 0; i < 2 * n; ++i) {
        if (!isFixed[i]) freeIdx.push_back(i);
    }
    for (int i = 0; i < nFree; ++i) {
        for (int j = 0; j < nFree; ++j) {
            out.hessian(i, j) = Hfull(freeIdx[i], freeIdx[j]);
        }
    }

    return out;
}

double lineSearchLocalInjectivity(const Eigen::VectorXcd& fz,
                                  const Eigen::VectorXcd& fbz,
                                  const Eigen::VectorXcd& d_fz,
                                  const Eigen::VectorXcd& d_fbz)
{
    const int n = static_cast<int>(fz.size());
    double tGlobal = std::numeric_limits<double>::infinity();

    for (int i = 0; i < n; ++i) {
        const double A = std::norm(d_fz(i)) - std::norm(d_fbz(i));
        const double B = 2.0 * ((d_fz(i).real() * fz(i).real() + d_fz(i).imag() * fz(i).imag())
                              - (d_fbz(i).real() * fbz(i).real() + d_fbz(i).imag() * fbz(i).imag()));
        const double C = std::norm(fz(i)) - std::norm(fbz(i));
        const double Delta = B * B - 4.0 * A * C;

        double t1 = std::numeric_limits<double>::infinity();
        double t2 = std::numeric_limits<double>::infinity();

        if (std::abs(A) <= kTol) {
            if (std::abs(B) > kTol) {
                t1 = -C / B;
                t2 = t1;
            }
        } else if (Delta > -kTol) {
            const double sqrtD = (Delta > 0.0) ? std::sqrt(Delta) : 0.0;
            if (B > 0.0) {
                const double Q = -0.5 * (B + sqrtD);
                t1 = Q / A;
                t2 = C / Q;
            } else {
                const double Q = -0.5 * (B - sqrtD);
                t1 = Q / A;
                t2 = C / Q;
            }
        }

        if (t1 < 0.0) t1 = std::numeric_limits<double>::infinity();
        if (t2 < 0.0) t2 = std::numeric_limits<double>::infinity();
        tGlobal = std::min(tGlobal, std::min(t1, t2));
    }

    if (!std::isfinite(tGlobal)) {
        tGlobal = 1.0;
    }
    return std::min(1.0, 0.9 * tGlobal);
}

double lineSearchDecreasingEnergy(const Eigen::VectorXcd& fz,
                                  const Eigen::VectorXcd& fbz,
                                  const Eigen::VectorXcd& d_fz,
                                  const Eigen::VectorXcd& d_fbz,
                                  const Eigen::VectorXd& Area,
                                  double E,
                                  const Eigen::VectorXd& G,
                                  const Eigen::VectorXd& d,
                                  double tMax)
{
    double t = tMax;
    const double alpha = 1e-4;
    const double beta = 0.5;
    const double alphaGd = alpha * G.dot(d);

    for (int ii = 0; ii < 30; ++ii) {
        const Eigen::VectorXcd fz_t = fz + t * d_fz;
        const Eigen::VectorXcd fbz_t = fbz + t * d_fbz;
        const double E_t = symDirEnergyByfzfbz(fz_t, fbz_t, Area);
        if (E_t < E + t * alphaGd) {
            return t;
        }
        t *= beta;
    }
    return t;
}

double optimizeSymDirEnergyByGlobalScaling(const Eigen::VectorXcd& fz,
                                           const Eigen::VectorXcd& fbz,
                                           const Eigen::VectorXd& Area)
{
    double x = 0.0, y = 0.0;
    for (int i = 0; i < fz.size(); ++i) {
        const double a = std::norm(fz(i));
        const double b = std::norm(fbz(i));
        const double g2 = a - b;
        const double g2sq = std::max(g2 * g2, kEps);
        x += Area(i) * ((a + b) / g2sq);
        y += Area(i) * (a + b);
    }
    if (y < kEps) {
        return 1.0;
    }
    const double s = x / y;
    return std::pow(s, 0.25);
}

void fixFirstCone(const Eigen::VectorXd& initialValue,
                  Eigen::VectorXd& x,
                  std::vector<int>& fixedIndices,
                  Eigen::VectorXd& fixedValues)
{
    const int n = static_cast<int>(initialValue.size()) / 2;
    fixedIndices = {0, n};
    fixedValues.resize(2);
    fixedValues(0) = initialValue(0);
    fixedValues(1) = initialValue(n);

    std::vector<bool> isFixed(2 * n, false);
    isFixed[0] = true;
    isFixed[n] = true;
    int count = 0;
    for (bool f : isFixed) if (!f) ++count;
    x.resize(count);
    int j = 0;
    for (int i = 0; i < 2 * n; ++i) {
        if (!isFixed[i]) {
            x(j++) = initialValue(i);
        }
    }
}

Eigen::Vector2d putVertexInKernel(const Eigen::Vector2d& UV,
                                  const Eigen::Matrix2Xd& oneRing)
{
    const int m = static_cast<int>(oneRing.cols());
    if (m < 2) {
        return UV;
    }

    Eigen::Matrix2Xd edges(2, m);
    for (int i = 0; i < m; ++i) {
        const int next = (i + 1) % m;
        edges.col(i) = oneRing.col(next) - oneRing.col(i);
    }

    Eigen::MatrixXd A(m, 2);
    Eigen::VectorXd b(m);
    double minCross = std::numeric_limits<double>::infinity();
    for (int i = 0; i < m; ++i) {
        const Eigen::Vector2d e = edges.col(i);
        const Eigen::Vector2d v = UV - oneRing.col(i);
        const double cross = e.x() * v.y() - e.y() * v.x();
        minCross = std::min(minCross, std::abs(cross));
        A(i, 0) = e.y();
        A(i, 1) = -e.x();
        b(i) = cross + (minCross / 100.0);
    }

    const double crossEps = minCross / 100.0;
    for (int i = 0; i < m; ++i) {
        const Eigen::Vector2d e = edges.col(i);
        const double rhs = crossEps + e.y() * oneRing(0, i) - e.x() * oneRing(1, i);
        A(i, 0) = e.y();
        A(i, 1) = -e.x();
        b(i) = rhs;
    }

    Eigen::Vector2d newUV = solveHalfPlaneQP(UV, A, b);

    bool valid = true;
    for (int i = 0; i < m; ++i) {
        const Eigen::Vector2d e = edges.col(i);
        const Eigen::Vector2d v = newUV - oneRing.col(i);
        const double cross = e.x() * v.y() - e.y() * v.x();
        if (cross < 0.0) {
            valid = false;
            break;
        }
    }
    return valid ? newUV : UV;
}

ATPResult ATPForInitialValue(const Eigen::MatrixXcd& J_fz,
                             const Eigen::MatrixXcd& J_fbz,
                             const Eigen::VectorXd& Area,
                             const Eigen::VectorXcd& frames)
{
    ATPResult result;
    const int numCones = static_cast<int>(J_fz.cols());
    const int numTriangles = static_cast<int>(J_fz.rows());

    Eigen::MatrixXcd J(2 * numTriangles, numCones);
    J.topRows(numTriangles) = J_fz;
    J.bottomRows(numTriangles) = J_fbz;

    Eigen::RowVectorXcd onesRow = Eigen::RowVectorXcd::Ones(numCones);
    Eigen::MatrixXcd JWithOnes(2 * numTriangles + 1, numCones);
    JWithOnes.topRows(2 * numTriangles) = J;
    JWithOnes.bottomRows(1) = onesRow;

    Eigen::VectorXd W_small(2 * numTriangles);
    for (int i = 0; i < numTriangles; ++i) {
        W_small(i) = Area(i);
        W_small(i + numTriangles) = Area(i);
    }
    Eigen::VectorXd W(W_small.size() + 1);
    W.head(W_small.size()) = W_small;
    W(W.size() - 1) = Area.mean();

    Eigen::MatrixXcd JW = JWithOnes.adjoint();
    for (int j = 0; j < JW.cols(); ++j) {
        JW.col(j) *= W(j);
    }
    const Eigen::MatrixXcd MInv = (JW * JWithOnes).inverse();
    const Eigen::MatrixXcd MInvJtrans = MInv * JW;

    Eigen::VectorXcd start_fz(numTriangles), start_fbz(numTriangles);
    for (int i = 0; i < numTriangles; ++i) {
        start_fz(i) = Complex(1.0, 0.0) / frames(i);
        start_fbz(i) = Complex(0.0, 0.0);
    }

    Eigen::VectorXcd a_start(2 * numTriangles + 1);
    a_start.head(numTriangles) = start_fz;
    a_start.segment(numTriangles, numTriangles) = start_fbz;
    a_start(a_start.size() - 1) = Complex(0.0, 0.0);

    Eigen::VectorXcd c, a_l;
    globalStepMAP(a_start, MInvJtrans, JWithOnes, c, a_l);

    Eigen::VectorXcd fz_global = a_l.head(numTriangles);
    Eigen::VectorXcd fbz_global = a_l.tail(numTriangles);

    const int maxIter = 500;
    const double sigma2Eps = 0.01;
    const double k = 0.9;
    const double tol = 1e-4;

    int iter = 0;
    for (iter = 1; iter <= maxIter; ++iter) {
        Eigen::VectorXcd fz_local, fbz_local;
        localStep(frames, sigma2Eps, k, fz_global, fbz_global, fz_local, fbz_local);
        Eigen::VectorXcd b_l(2 * numTriangles);
        b_l.head(numTriangles) = fz_local;
        b_l.tail(numTriangles) = fbz_local;

        Eigen::VectorXcd ni = a_l - b_l;
        double norm_ni_sqr = 0.0;
        for (int i = 0; i < W_small.size(); ++i) {
            norm_ni_sqr += W_small(i) * std::norm(ni(i));
        }

        if (std::sqrt(norm_ni_sqr) < tol) {
            break;
        }

        Eigen::VectorXcd niAug(ni.size() + 1);
        niAug.head(ni.size()) = ni;
        niAug(ni.size()) = Complex(0.0, 0.0);

        globalStepATP(MInvJtrans, JWithOnes, c, niAug, norm_ni_sqr, W, c, a_l);
        fz_global = a_l.head(numTriangles);
        fbz_global = a_l.tail(numTriangles);
    }

    result.iterations = iter;
    int numFolds = 0;
    for (int i = 0; i < numTriangles; ++i) {
        if (std::abs(fz_global(i)) < std::abs(fbz_global(i))) {
            ++numFolds;
        }
    }

    bool hasNan = false;
    for (int i = 0; i < c.size(); ++i) {
        if (!std::isfinite(c(i).real()) || !std::isfinite(c(i).imag())) {
            hasNan = true;
            break;
        }
    }

    result.success = (numFolds == 0) && !hasNan;
    result.UVonCones.resize(2 * numCones);
    for (int i = 0; i < numCones; ++i) {
        result.UVonCones(i) = c(i).real();
        result.UVonCones(i + numCones) = c(i).imag();
    }
    return result;
}

NewtonResult runNewton(const Eigen::MatrixXcd& J_fz,
                       const Eigen::MatrixXcd& J_fbz,
                       const Eigen::VectorXd& Area,
                       const Eigen::VectorXd& xInitial,
                       const std::vector<int>& fixedIndices,
                       const Eigen::VectorXd& fixedValues,
                       bool atpSuccess)
{
    NewtonResult result;
    const int n = static_cast<int>(J_fz.cols());
    Eigen::VectorXd x = xInitial;
    Eigen::VectorXd fixedVals = fixedValues;

    int maxIter = atpSuccess ? 500 : 0;
    const bool spdHessian = true;
    const double stepTol = 1e-7;
    const double gradTol = 1e-7;
    const double Etol = 1e-8;
    int EitersCounter = 0;
    const int numConvIters = 5;

    auto energyResult = symDirEnergyGradHess(x, J_fz, J_fbz, Area, fixedIndices, fixedVals, spdHessian);
    double E = energyResult.energy;

    const double scale = optimizeSymDirEnergyByGlobalScaling(energyResult.fz, energyResult.fbz, Area);
    x *= scale;
    fixedVals *= scale;

    int iter = 0;
    for (iter = 1; iter <= maxIter; ++iter) {
        const double Eprev = E;
        energyResult = symDirEnergyGradHess(x, J_fz, J_fbz, Area, fixedIndices, fixedVals, spdHessian);
        E = energyResult.energy;
        const Eigen::VectorXd& G = energyResult.gradient;
        Eigen::MatrixXd H = energyResult.hessian;

        Eigen::LLT<Eigen::MatrixXd> llt(H);
        if (llt.info() != Eigen::Success) {
            break;
        }
        const Eigen::VectorXd d = llt.solve(-G);

        Eigen::VectorXd d_full = Eigen::VectorXd::Zero(2 * n);
        std::vector<bool> isFixed(2 * n, false);
        for (int idx : fixedIndices) {
            if (idx >= 0 && idx < 2 * n) isFixed[idx] = true;
        }
        int j = 0;
        for (int i = 0; i < 2 * n; ++i) {
            if (!isFixed[i]) {
                d_full(i) = d(j++);
            }
        }

        Eigen::VectorXcd d_complex(n);
        for (int i = 0; i < n; ++i) {
            d_complex(i) = Complex(d_full(i), d_full(i + n));
        }
        const Eigen::VectorXcd d_fz = J_fz * d_complex;
        const Eigen::VectorXcd d_fbz = J_fbz * d_complex;

        const double tMax = lineSearchLocalInjectivity(energyResult.fz, energyResult.fbz, d_fz, d_fbz);
        const double t = lineSearchDecreasingEnergy(energyResult.fz, energyResult.fbz, d_fz, d_fbz,
                                                    Area, E, G, d, tMax);

        x += t * d;

        if (G.norm() < gradTol) break;
        if (t * d.norm() < stepTol) break;

        if (std::abs(E - Eprev) < Etol * (E + 1.0)) {
            if (EitersCounter >= numConvIters) break;
            ++EitersCounter;
        } else {
            EitersCounter = 0;
        }
    }

    result.iterations = iter;
    result.finalEnergy = E;
    result.success = (iter > 0) && (iter < maxIter);

    Eigen::VectorXd full(2 * n);
    full.setZero();
    std::vector<bool> isFixed(2 * n, false);
    for (int idx : fixedIndices) {
        if (idx >= 0 && idx < 2 * n) isFixed[idx] = true;
    }
    int j = 0;
    for (int i = 0; i < 2 * n; ++i) {
        if (!isFixed[i]) {
            full(i) = x(j++);
        }
    }
    for (size_t k = 0; k < fixedIndices.size(); ++k) {
        full(fixedIndices[k]) = fixedVals(k);
    }

    result.UVonCones.resize(n, 2);
    for (int i = 0; i < n; ++i) {
        result.UVonCones(i, 0) = full(i);
        result.UVonCones(i, 1) = full(i + n);
    }
    return result;
}

} // namespace FastHGPNumerics
