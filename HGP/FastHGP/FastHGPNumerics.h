#pragma once

#include <complex>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace FastHGPNumerics {

using Complex = std::complex<double>;

struct PerpsResult {
    Eigen::MatrixXcd tc;   // numF x 3 complex perps
    Eigen::VectorXd dblAc; // doubled triangle areas
};

struct SymDirEnergyResult {
    Eigen::VectorXcd fz;
    Eigen::VectorXcd fbz;
    double energy = 0.0;
    Eigen::VectorXd gradient;
    Eigen::MatrixXd hessian;
};

void computeLocalBasis(const Eigen::MatrixXd& v1,
                       const Eigen::MatrixXd& v2,
                       Eigen::MatrixXd& e1,
                       Eigen::MatrixXd& e2);

PerpsResult computePerps(const Eigen::MatrixXd& V,
                         const Eigen::MatrixXi& F,
                         int xlb = 3);

void createJmatrix(const Eigen::SparseMatrix<Complex>& harmonicBasis,
                   const Eigen::MatrixXd& V,
                   const Eigen::MatrixXi& F,
                   Eigen::MatrixXcd& J_fz,
                   Eigen::MatrixXcd& J_fbz,
                   Eigen::VectorXd& Area);

void localStep(const Eigen::VectorXcd& frames,
               double sigma2Eps,
               double k,
               const Eigen::VectorXcd& fz,
               const Eigen::VectorXcd& fbz,
               Eigen::VectorXcd& fz_local,
               Eigen::VectorXcd& fbz_local);

void globalStepMAP(const Eigen::VectorXcd& a_l_in,
                   const Eigen::MatrixXcd& MInvJtrans,
                   const Eigen::MatrixXcd& JWithOnes,
                   Eigen::VectorXcd& c,
                   Eigen::VectorXcd& a_l_out);

void globalStepATP(const Eigen::MatrixXcd& MInvJtrans,
                   const Eigen::MatrixXcd& JWithOnes,
                   const Eigen::VectorXcd& c_i,
                   const Eigen::VectorXcd& ni,
                   double norm_ni_sqr,
                   const Eigen::VectorXd& W,
                   Eigen::VectorXcd& c_i_p_1,
                   Eigen::VectorXcd& a_l);

double symDirEnergyByfzfbz(const Eigen::VectorXcd& fz,
                           const Eigen::VectorXcd& fbz,
                           const Eigen::VectorXd& Area);

SymDirEnergyResult symDirEnergyGradHess(const Eigen::VectorXd& xSmall,
                                        const Eigen::MatrixXcd& J_fz,
                                        const Eigen::MatrixXcd& J_fbz,
                                        const Eigen::VectorXd& Area,
                                        const std::vector<int>& fixedIndices,
                                        const Eigen::VectorXd& fixedValues,
                                        bool spdHessian);

double lineSearchLocalInjectivity(const Eigen::VectorXcd& fz,
                                  const Eigen::VectorXcd& fbz,
                                  const Eigen::VectorXcd& d_fz,
                                  const Eigen::VectorXcd& d_fbz);

double lineSearchDecreasingEnergy(const Eigen::VectorXcd& fz,
                                  const Eigen::VectorXcd& fbz,
                                  const Eigen::VectorXcd& d_fz,
                                  const Eigen::VectorXcd& d_fbz,
                                  const Eigen::VectorXd& Area,
                                  double E,
                                  const Eigen::VectorXd& G,
                                  const Eigen::VectorXd& d,
                                  double tMax);

double optimizeSymDirEnergyByGlobalScaling(const Eigen::VectorXcd& fz,
                                           const Eigen::VectorXcd& fbz,
                                           const Eigen::VectorXd& Area);

void fixFirstCone(const Eigen::VectorXd& initialValue,
                  Eigen::VectorXd& x,
                  std::vector<int>& fixedIndices,
                  Eigen::VectorXd& fixedValues);

Eigen::Vector2d putVertexInKernel(const Eigen::Vector2d& UV,
                                  const Eigen::Matrix2Xd& oneRing);

struct ATPResult {
    bool success = false;
    Eigen::VectorXd UVonCones;
    int iterations = 0;
};

ATPResult ATPForInitialValue(const Eigen::MatrixXcd& J_fz,
                             const Eigen::MatrixXcd& J_fbz,
                             const Eigen::VectorXd& Area,
                             const Eigen::VectorXcd& frames);

struct NewtonResult {
    bool success = false;
    Eigen::MatrixXd UVonCones;
    int iterations = 0;
    double finalEnergy = 0.0;
};

NewtonResult runNewton(const Eigen::MatrixXcd& J_fz,
                       const Eigen::MatrixXcd& J_fbz,
                       const Eigen::VectorXd& Area,
                       const Eigen::VectorXd& xInitial,
                       const std::vector<int>& fixedIndices,
                       const Eigen::VectorXd& fixedValues,
                       bool atpSuccess);

Eigen::VectorXcd expandReducedX(const Eigen::VectorXd& xSmall,
                                int n,
                                const std::vector<int>& fixedIndices,
                                const Eigen::VectorXd& fixedValues);

Eigen::VectorXd packReducedX(const Eigen::VectorXd& fullX,
                             const std::vector<int>& fixedIndices);

} // namespace FastHGPNumerics
