#include "EigenLinearSolver.h"
#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>

// ============================================================
// 构造 / 析构
// ============================================================

EigenLinearSolver::EigenLinearSolver(int matrixType, int solver_in)
	: mtype(matrixType), solver(solver_in)
{
	mIsSymmetric = (mtype == 2) || (mtype == -2);
	mInitialized = false;
	n = 0;
	mNonzeros = 0;
	error = 0;
}

EigenLinearSolver::~EigenLinearSolver()
{
	// 智能指针自动释放
	mpLU.reset();
	mpLDLT.reset();
	mpLLT.reset();
}

// ============================================================
// 基础设置
// ============================================================

void EigenLinearSolver::setMatrixType(int matrixType)
{
	mtype = matrixType;
	mIsSymmetric = (mtype == 2) || (mtype == -2);
}

bool EigenLinearSolver::init(bool isHighlyIndefinite)
{
	// Eigen 不需要显式初始化，这里仅做参数记录
	(void)isHighlyIndefinite;
	mInitialized = true;

	// 预分配求解器（延迟到 preprocess 时真正构建）
	return true;
}

void EigenLinearSolver::setSelectiveRows(const std::unordered_set<int>& rows)
{
	mSelectiveRows = rows;
}

// ============================================================
// 矩阵加载 —— 从 Eigen 稀疏矩阵
// ============================================================

bool EigenLinearSolver::createPardisoFormatMatrix(
	const Eigen::SparseMatrix<double, Eigen::RowMajor>& M)
{
	// 将 RowMajor 转为 ColMajor (Eigen 内部格式)
	mA_eigen = M;  // Eigen 会自动处理转换
	if (!mA_eigen.isCompressed())
		mA_eigen.makeCompressed();

	n = mA_eigen.rows();
	mNonzeros = mA_eigen.nonZeros();

	if (mA_eigen.rows() != mA_eigen.cols()) {
		std::cout << "ERROR: Matrix must be square\n";
		return false;
	}

	return true;
}

bool EigenLinearSolver::createPardisoFormatMatrix(
	const Eigen::SparseMatrix<double, Eigen::RowMajor>& M,
	std::vector<int>& rowsOfElementsToSet,
	std::vector<int>& colsOfElementsToSet)
{
	// 先执行基本加载
	if (!createPardisoFormatMatrix(M))
		return false;

	// 如果需要修改特定元素，在这里可以预处理
	// （当前 PARDISO 接口中此功能用于设置对角线元素等）
	(void)rowsOfElementsToSet;
	(void)colsOfElementsToSet;

	return true;
}

bool EigenLinearSolver::createPardisoFormatMatrix(
	const GMMCompressed1RowMatrix& M,
	std::vector<int>& rowsOfElementsToSet,
	std::vector<int>& colsOfElementsToSet)
{
	// 从 GMM 格式转换为 Eigen 稀疏矩阵
	// GMMCompressed1RowMatrix 使用 1-based CSR 格式
	int nnz = static_cast<int>(M.pr().size());
	int sz = M.nrows();

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(nnz);

	for (int i = 0; i < sz; ++i) {
		for (int j = M.jc()[i] - 1; j < M.jc()[i + 1] - 1; ++j) {
			int row = i;
			int col = static_cast<int>(M.ir_vec()[j]) - 1;  // 1-based -> 0-based
			double val = M.pr()[j];
			triplets.emplace_back(row, col, val);
		}
	}

	mA_eigen.resize(sz, sz);
	mA_eigen.setFromTriplets(triplets.begin(), triplets.end());
	mA_eigen.makeCompressed();

	n = sz;
	mNonzeros = nnz;

	(void)rowsOfElementsToSet;
	(void)colsOfElementsToSet;
	return true;
}

// ============================================================
// 分解阶段
// ============================================================

bool EigenLinearSolver::doFactorization()
{
	// 清理旧的求解器
	mpLU.reset();
	mpLDLT.reset();
	mpLLT.reset();

	try {
		switch (mtype) {
		case 2:  // 对称正定
			mpLLT = std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>>();
			mpLLT->compute(mA_eigen);
			if (mpLLT->info() != Eigen::Success) {
				std::cout << "ERROR: SimplicialLLT factorization failed\n";
				return false;
			}
			break;

		case -2:  // 对称不定 (KKT 系统)
			// 使用 LDLT 处理不定矩阵
			mpLDLT = std::make_unique<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>(mA_eigen);
			if (mpLDLT->info() != Eigen::Success) {
				std::cout << "WARNING: SimplicialLDLT failed, falling back to SparseLU\n";
				// fallback to LU
				mpLDLT.reset();
				mpLU = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<double>>>();
				mpLU->compute(mA_eigen);
				if (mpLU->info() != Eigen::Success) {
					std::cout << "ERROR: SparseLU fallback also failed\n";
					return false;
				}
			}
			break;

		case 11:  // 一般非对称矩阵
		default:
			mpLU = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<double>>>();
			mpLU->compute(mA_eigen);
			if (mpLU->info() != Eigen::Success) {
				std::cout << "ERROR: SparseLU factorization failed\n";
				return false;
			}
			break;
		}
	}
	catch (const std::exception& e) {
		std::cout << "EXCEPTION during factorization: " << e.what() << std::endl;
		return false;
	}

	return true;
}

bool EigenLinearSolver::reorder()
{
	// Eigen 的符号分析在 compute() 中自动完成
	// 此处保持接口兼容：只做初步分析
	return doFactorization();
}

bool EigenLinearSolver::preprocess()
{
	// phase 12 equivalent: 符号分解 + 数值分解
	return doFactorization();
}

bool EigenLinearSolver::numericalFactorization()
{
	// phase 22 equivalent: 仅数值分解
	// Eigen 中无法分离符号/数值阶段，直接调用即可
	return doFactorization();
}

// ============================================================
// 选择性求逆 —— 核心方法
// ============================================================
//
// 对于 KKT 系统的稀疏矩阵 A (n x n)，我们只需要计算 mSelectiveRows
// 中指定的行对应的逆矩阵元素。
//
// 数学原理 (Schur Complement):
//
//   将 A 的行列重排为：
//          [ A_ii   A_ib ]
//   A   =  [ A_bi   A_bb ]
//
//   其中 b = mSelectiveRows 指定的行/列, i = 其余部分
//
//   则 A 的逆在 bb 块为:
//
//   S = A_bb - A_bi * A_ii^{-1} * A_ib    (Schur 补)
//   A^{-1}[b,b] = S^{-1}
//
//   A^{-1}[i,b] = -A_ii^{-1} * A_ib * S^{-1}
//

bool EigenLinearSolver::selectiveInverse()
{
	if (mSelectiveRows.empty()) {
		std::cout << "ERROR: No selective rows specified. Call setSelectiveRows() first.\n";
		return false;
	}

	const int n_total = n;
	int nb = static_cast<int>(mSelectiveRows.size());
	int ni = n_total - nb;

	if (nb >= n_total) {
		std::cout << "ERROR: Selective rows cover entire matrix, use full inverse instead.\n";
		return false;
	}

	// ---- 步骤1: 建立 index 映射 ----
	std::vector<int> idx_b(mSelectiveRows.begin(), mSelectiveRows.end());  // b-indices (sorted)
	std::sort(idx_b.begin(), idx_b.end());

	std::vector<int> global_to_local(n_total, -1);     // 全局索引 -> 局部 b-index
	std::vector<int> global_to_i(n_total, -1);          // 全局索引 -> 局部 i-index

	for (int k = 0; k < nb; ++k)
		global_to_local[idx_b[k]] = k;

	int i_count = 0;
	for (int g = 0; g < n_total; ++g) {
		if (global_to_local[g] == -1) {
			global_to_i[g] = i_count++;
		}
	}

	// ---- 步骤2: 提取子块 A_bb, A_ib, A_bi ----
	// A_bb (nb x nb), 密集存储
	Eigen::MatrixXd A_bb(nb, nb);
	A_bb.setZero();
	{
		const auto& M = mA_eigen;
		for (int k = 0; k < M.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
				int row = it.row();
				int col = it.col();
				int lr = global_to_local[row];
				int lc = global_to_local[col];
				if (lr >= 0 && lc >= 0) {
					A_bb(lr, lc) += it.value();
				}
			}
		}
	}

	// A_ib (ni x nb) — 从 i 行, b 列提取
	Eigen::MatrixXd A_ib(ni, nb);
	A_ib.setZero();
	{
		const auto& M = mA_eigen;
		for (int k = 0; k < M.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
				int row = it.row();
				int col = it.col();
				int ir = global_to_i[row];
				int lc = global_to_local[col];
				if (ir >= 0 && lc >= 0) {
					A_ib(ir, lc) += it.value();
				}
			}
		}
	}

	// A_bi (nb x ni) = A_ib^T
	Eigen::MatrixXd A_bi = A_ib.transpose();

	// ---- 步骤3: 提取 A_ii 并求逆或求解 ----
	// A_ii 是 (ni x ni)，通常很大但稀疏
	// 我们不需要完整的 A_ii^{-1}，只需计算 A_ii^{-1} * A_ib

	Eigen::MatrixXd A_ii_inv_Aib(ni, nb);
	A_ii_inv_Aib.setZero();

	// 方法：用稀疏求解器对 A_ii 的每一列分别求解
	// 构建 A_ii 子矩阵
	std::vector<Eigen::Triplet<double>> ii_triplets;
	{
		const auto& M = mA_eigen;
		for (int k = 0; k < M.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
				int row = it.row();
				int col = it.col();
				int ir = global_to_i[row];
				int ic = global_to_i[col];
				if (ir >= 0 && ic >= 0) {
					ii_triplets.emplace_back(ir, ic, it.value());
				}
			}
		}
	}

	Eigen::SparseMatrix<double> A_ii(ni, ni);
	A_ii.setFromTriplets(ii_triplets.begin(), ii_triplets.end());
	A_ii.makeCompressed();

	// 对 A_ib 的每一列求解 A_ii * x = A_ib(:, j)
	{
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_ii;
		solver_ii.compute(A_ii);
		if (solver_ii.info() != Eigen::Success) {
			// fallback to iterative solver for indefinite systems
			std::cout << "WARNING: SparseLU failed for A_ii, trying BiCGSTAB...\n";

			// 用 Jacobi 预条件的 BiCGSTAB
			Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>,
				Eigen::IncompleteLUT<double>> iter_solver;
			iter_solver.compute(A_ii);
			
			for (int col = 0; col < nb; ++col) {
				Eigen::VectorXd rhs_col = A_ib.col(col);
				Eigen::VectorXd x_col = iter_solver.solve(rhs_col);
				A_ii_inv_Aib.col(col) = x_col;
			}
		}
		else {
			for (int col = 0; col < nb; ++col) {
				Eigen::VectorXd rhs_col = A_ib.col(col);
				Eigen::VectorXd x_col = solver_ii.solve(rhs_col);
				A_ii_inv_Aib.col(col) = x_col;
			}
		}
	}

	// ---- 步骤4: 计算 Schur 补 S = A_bb - A_bi * (A_ii^{-1} * A_ib) ----
	Eigen::MatrixXd S = A_bb - A_bi * A_ii_inv_Aib;

	// ---- 步骤5: 求 S^{-1} ----
	Eigen::MatrixXd S_inv = S.inverse();  // S 很小 (nb x nb), 直接求逆

	// ---- 步骤6: 计算完整的选择性逆结果 ----
	// 结果布局: 与原 PARDISO 相同——返回 CSR 格式的上三角部分
	// 但我们只需要 mSelectiveRows 对应的列 (即 KKT 矩阵右下角的锥点约束列)

	// mSelectiveInvResult 存储格式: 只包含需要的行和列
	// 这里我们按原始 PARDISO 输出格式构造: 上三角部分的 CSR

	// 实际上 FastHGP 只需要读取 selectiveInvSol(row, col) 其中 row 在 mConesAndNearConesMapOfRowsInKKT 中
	// col 在 [conesConstraintsStartRow, conesConstraintsStartRow + mNumDOF) 范围内

	// 我们将结果存储在一个密集矩阵中供 getMatrixInGMMformat 读取
	int result_n = n;
	mSelectiveInvResult.resize(result_n, result_n);
	mSelectiveInvResult.setZero();

	// 填充 bb 块: A^{-1}[b, b] = S^{-1}
	for (int kr = 0; kr < nb; ++kr) {
		for (int kc = 0; kc < nb; ++kc) {
			mSelectiveInvResult.coeffRef(idx_b[kr], idx_b[kc]) = S_inv(kr, kc);
		}
	}

	// 填充 ib 块: A^{-1}[i, b] = -A_ii^{-1} * A_ib * S^{-1}
	Eigen::MatrixXd Aiib_Sinv = A_ii_inv_Aib * S_inv;  // (ni x nb)
	for (int gi = 0; gi < n_total; ++gi) {
		int ir = global_to_i[gi];
		if (ir >= 0) {
			for (int kc = 0; kc < nb; ++kc) {
				mSelectiveInvResult.coeffRef(gi, idx_b[kc]) = -Aiib_Sinv(ir, kc);
			}
		}
	}

	// 注意: bi 块是 ib 块的转置 (对称情况)
	if (mIsSymmetric) {
		for (int gi = 0; gi < n_total; ++gi) {
			int ir = global_to_i[gi];
			if (ir >= 0) {
				for (int kc = 0; kc < nb; ++kc) {
					mSelectiveInvResult.coeffRef(idx_b[kc], gi) = -Aiib_Sinv(ir, kc);
				}
			}
		}
	}

	return true;
}


// ============================================================
// 结果获取
// ============================================================

bool EigenLinearSolver::getMatrixInGMMformat(GMMSparseRowMatrix& M)
{
	M.resize(n, n);
	M.setZero();

	// 从选择性逆结果中填充
	for (int k = 0; k < mSelectiveInvResult.outerSize(); ++k) {
		for (decltype(mSelectiveInvResult)::InnerIterator it(mSelectiveInvResult, k); it; ++it) {
			int r = it.row();
			int c = it.col();
			M(r, c) = it.value();
			if (r != c)
				M(c, r) = it.value();
		}
	}
	return true;
}

bool EigenLinearSolver::getMatrixInGMMformat(GMMCompressed1RowMatrix& M)
{
	if (M.nrows() != static_cast<int>(mSelectiveInvResult.rows()) ||
		M.ncols() != static_cast<int>(mSelectiveInvResult.cols()))
	{
		M.resize(n, n);
	}

	// 转换为 GMM 1-based CSR 格式
	std::vector<std::tuple<int, int, double>> entries;

	for (int k = 0; k < mSelectiveInvResult.outerSize(); ++k) {
		for (decltype(mSelectiveInvResult)::InnerIterator it(mSelectiveInvResult, k); it; ++it) {
			int r2 = it.row(), c2 = it.col();
			entries.emplace_back(r2, c2, it.value());
			if (r2 != c2) {
				entries.emplace_back(c2, r2, it.value());  // 对称填充
			}
		}
	}

	// 按行列排序 (CSR 要求行优先)
	std::sort(entries.begin(), entries.end(),
		[](const auto& a, const auto& b) {
			if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
			return std::get<1>(a) < std::get<1>(b);
		});

	// 构建 CSR 结构
	M.pr().clear();
	M.ir_vec().clear();
	M.jc().clear();
	M.pr().reserve(entries.size());
	M.ir_vec().reserve(entries.size());
	M.jc().reserve(n + 1);

	int current_row = -1;
	for (const auto& entry : entries) {
		int row = std::get<0>(entry);
		int col = std::get<1>(entry);
		double val = std::get<2>(entry);
		while (current_row < row) {
			M.jc().push_back(static_cast<int>(M.pr().size()) + 1);  // 1-based
			current_row++;
		}
		M.ir_vec().push_back(col + 1);  // 1-based column index
		M.pr().push_back(val);
	}
	while (current_row < n - 1) {
		M.jc().push_back(static_cast<int>(M.pr().size()) + 1);
		current_row++;
	}
	M.jc().push_back(static_cast<int>(M.pr().size()) + 1);  // 结束哨兵

	return true;
}
