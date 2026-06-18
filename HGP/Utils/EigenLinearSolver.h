#pragma once

#include <memory>
#include <unordered_set>
#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

// Eigen-only sparse solver with selective inverse (for future KKT / harmonic basis).
class EigenLinearSolver {
public:
    explicit EigenLinearSolver(int matrixType = 11, int solver = 0);
    ~EigenLinearSolver();

    void setMatrixType(int matrixType);
    bool init(bool isHighlyIndefinite = false);
    bool preprocess();
    bool reorder();
    bool numericalFactorization();
    bool selectiveInverse();

    bool setMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix);
    bool createPardisoFormatMatrix(const Eigen::SparseMatrix<double>& matrix);
    bool solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& x) const;

    void setSelectiveRows(const std::unordered_set<int>& rows);
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& selectiveInverseMatrix() const;

private:
    int mtype;
    int solver;
    int n;
    bool mIsSymmetric;
    bool mInitialized;

    Eigen::SparseMatrix<double> mA_eigen;
    std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<double>>> mpLU;
    std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> mpLDLT;
    std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>> mpLLT;

    std::unordered_set<int> mSelectiveRows;
    Eigen::SparseMatrix<double, Eigen::RowMajor> mSelectiveInvResult;

    bool doFactorization();
};
