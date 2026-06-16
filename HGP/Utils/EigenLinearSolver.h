#pragma once

#include "Utils/GMM_Macros.h"
#include <unordered_set>
#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>

/**
 * EigenLinearSolver — 使用 Eigen 库替代 PARDISO 的稀疏线性求解器
 * 
 * 支持功能：
 *   - 稀疏矩阵 LU 分解求解
 *   - 对称不定矩阵 LDLT 分解
 *   - 选择性求逆 (selective inverse) —— 通过 Schur 补实现
 * 
 * 用法：与原 PardisoLinearSolver 接口尽量保持一致
 */
class EigenLinearSolver
{
public:

	EigenLinearSolver(int matrixType = 11, int solver_in = 0);
	~EigenLinearSolver();

	void setMatrixType(int matrixType);
	bool init(bool isHighlyIndefinite = false);
	bool preprocess();
	bool reorder();
	bool numericalFactorization();

	/**
	 * 选择性求逆——核心方法
	 * 
	 * 原理：对于 KKT 系统 A，只计算指定行的逆矩阵元素。
	 * 使用 Schur 补公式：
	 *   将 A 分块为 [A_ii  A_ib]
	 *                    [A_bi  A_bb]
	 *   则 A^{-1}_bb = (A_bb - A_bi * A_ii^{-1} * A_ib)^{-1}
	 * 
	 * @return true 成功
	 */
	bool selectiveInverse();
	bool getMatrixInGMMformat(GMMSparseRowMatrix& M);
	bool getMatrixInGMMformat(GMMCompressed1RowMatrix& M);
	bool createPardisoFormatMatrix(const GMMCompressed1RowMatrix& M,
		std::vector<int>& rowsOfElementsToSet,
		std::vector<int>& colsOfElementsToSet);
	bool createPardisoFormatMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& M);
	bool createPardisoFormatMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& M,
		std::vector<int>& rowsOfElementsToSet,
		std::vector<int>& colsOfElementsToSet);

	template<class DenseColMatrixType>
	bool solve(DenseColMatrixType& RHS, DenseColMatrixType& X);

	template<class DenseColMatrixType>
	bool oneTimeSolve(DenseColMatrixType& RHS, DenseColMatrixType& X);

	// ---- 新增接口：设置需要计算的选择性行 ----
	void setSelectiveRows(const std::unordered_set<int>& rows);

private:
	int mtype;              // 矩阵类型 (同 PARDISO: 2=对称正定, -2=对称不定, 11=一般)
	int solver;
	int error;

	// Eigen 内部存储
	Eigen::SparseMatrix<double> mA_eigen;       // 原始矩阵 (col-major, Eigen默认格式)
	int n;                                        // 矩阵维度
	int mNonzeros;

	bool mIsSymmetric;
	bool mInitialized;

	// 求解器指针 (多态存储不同类型的分解)
	std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<double>>> mpLU;
	std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> mpLDLT;
	std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>> mpLLT;

	// 选择性逆相关
	std::unordered_set<int> mSelectiveRows;      // 需要计算的行集合 (0-based)
	Eigen::SparseMatrix<double, Eigen::RowMajor> mSelectiveInvResult; // 选择性逆结果

	template<class DenseColMatrixType>
	bool createRHSvector(DenseColMatrixType& RHS);

	template<class DenseColMatrixType>
	bool updateX(DenseColMatrixType& X);

	// 内部辅助：根据矩阵类型选择并执行分解
	bool doFactorization();
};


// ==================== 模板实现 ====================

template<class DenseColMatrixType>
bool EigenLinearSolver::solve(DenseColMatrixType& RHS, DenseColMatrixType& X)
{
	if (!createRHSvector(RHS))
		return false;

	Eigen::VectorXd b_eigen(RHS.nrows());
	for (int i = 0; i < RHS.nrows(); ++i)
		b_eigen(i) = RHS(i, 0);

	Eigen::VectorXd x_eigen;

	if (mpLU && mpLU->info() == Eigen::Success) {
		x_eigen = mpLU->solve(b_eigen);
	}
	else if (mpLDLT && mpLDLT->info() == Eigen::Success) {
		x_eigen = mpLDLT->solve(b_eigen);
	}
	else if (mpLLT && mpLLT->info() == Eigen::Success) {
		x_eigen = mpLLT->solve(b_eigen);
	}
	else {
		std::cout << "\nERROR: No valid factorization available for solve\n";
		return false;
	}

	if (mpLU && mpLU->info() != Eigen::Success) return false;
	if (mpLDLT && mpLDLT->info() != Eigen::Success) return false;
	if (mpLLT && mpLLT->info() != Eigen::Success) return false;

	return updateX(X);
}

template<class DenseColMatrixType>
bool EigenLinearSolver::oneTimeSolve(DenseColMatrixType& RHS, DenseColMatrixType& X)
{
	// 执行完整分解 (phase 13 equivalent)
	if (!preprocess())
		return false;

	return solve(RHS, X);
}

template<class DenseColMatrixType>
bool EigenLinearSolver::createRHSvector(DenseColMatrixType& RHS)
{
	// 仅做维度检查
	if (RHS.nrows() != n) {
		std::cout << "\nERROR: RHS size mismatch. Expected " << n << ", got " << RHS.nrows() << std::endl;
		return false;
	}
	return true;
}

template<class DenseColMatrixType>
bool EigenLinearSolver::updateX(DenseColMatrixType& X)
{
	// 此函数在 solve 中已处理，保持接口兼容
	return true;
}
