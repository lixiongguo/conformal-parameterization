#include "EigenLinearSolver.h"
#include <algorithm>
#include <iostream>

EigenLinearSolver::EigenLinearSolver(int matrixType, int solverIn)
    : mtype(matrixType), solver(solverIn), n(0), mIsSymmetric(false), mInitialized(false)
{
    mIsSymmetric = (mtype == 2) || (mtype == -2);
}

EigenLinearSolver::~EigenLinearSolver() = default;

void EigenLinearSolver::setMatrixType(int matrixType)
{
    mtype = matrixType;
    mIsSymmetric = (mtype == 2) || (mtype == -2);
}

bool EigenLinearSolver::init(bool isHighlyIndefinite)
{
    (void)isHighlyIndefinite;
    mInitialized = true;
    return true;
}

void EigenLinearSolver::setSelectiveRows(const std::unordered_set<int>& rows)
{
    mSelectiveRows = rows;
}

const Eigen::SparseMatrix<double, Eigen::RowMajor>& EigenLinearSolver::selectiveInverseMatrix() const
{
    return mSelectiveInvResult;
}

bool EigenLinearSolver::createPardisoFormatMatrix(const Eigen::SparseMatrix<double>& matrix)
{
    Eigen::SparseMatrix<double, Eigen::RowMajor> rowMajor = matrix;
    return setMatrix(rowMajor);
}

bool EigenLinearSolver::setMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix)
{
    mA_eigen = matrix;
    if (!mA_eigen.isCompressed()) {
        mA_eigen.makeCompressed();
    }

    n = static_cast<int>(mA_eigen.rows());
    if (mA_eigen.rows() != mA_eigen.cols()) {
        std::cout << "ERROR: Matrix must be square\n";
        return false;
    }
    return true;
}

bool EigenLinearSolver::doFactorization()
{
    mpLU.reset();
    mpLDLT.reset();
    mpLLT.reset();

    try {
        switch (mtype) {
        case 2:
            mpLLT = std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>>();
            mpLLT->compute(mA_eigen);
            if (mpLLT->info() != Eigen::Success) {
                std::cout << "ERROR: SimplicialLLT factorization failed\n";
                return false;
            }
            break;

        case -2:
            mpLDLT = std::make_unique<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>(mA_eigen);
            if (mpLDLT->info() != Eigen::Success) {
                std::cout << "WARNING: SimplicialLDLT failed, falling back to SparseLU\n";
                mpLDLT.reset();
                mpLU = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<double>>>();
                mpLU->compute(mA_eigen);
                if (mpLU->info() != Eigen::Success) {
                    std::cout << "ERROR: SparseLU fallback also failed\n";
                    return false;
                }
            }
            break;

        case 11:
        default:
            mpLU = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<double>>>();
            mpLU->compute(mA_eigen);
            if (mpLU->info() != Eigen::Success) {
                std::cout << "ERROR: SparseLU factorization failed\n";
                return false;
            }
            break;
        }
    } catch (const std::exception& e) {
        std::cout << "EXCEPTION during factorization: " << e.what() << std::endl;
        return false;
    }

    return true;
}

bool EigenLinearSolver::reorder()
{
    return doFactorization();
}

bool EigenLinearSolver::preprocess()
{
    return doFactorization();
}

bool EigenLinearSolver::numericalFactorization()
{
    return doFactorization();
}

bool EigenLinearSolver::solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& x) const
{
    if (rhs.size() != n) {
        std::cout << "ERROR: RHS size mismatch\n";
        return false;
    }

    if (mpLU && mpLU->info() == Eigen::Success) {
        x = mpLU->solve(rhs);
        return mpLU->info() == Eigen::Success;
    }
    if (mpLDLT && mpLDLT->info() == Eigen::Success) {
        x = mpLDLT->solve(rhs);
        return mpLDLT->info() == Eigen::Success;
    }
    if (mpLLT && mpLLT->info() == Eigen::Success) {
        x = mpLLT->solve(rhs);
        return mpLLT->info() == Eigen::Success;
    }

    std::cout << "ERROR: No valid factorization available for solve\n";
    return false;
}

bool EigenLinearSolver::selectiveInverse()
{
    if (mSelectiveRows.empty()) {
        std::cout << "ERROR: No selective rows specified. Call setSelectiveRows() first.\n";
        return false;
    }

    const int nTotal = n;
    const int nb = static_cast<int>(mSelectiveRows.size());
    const int ni = nTotal - nb;

    if (nb >= nTotal) {
        std::cout << "ERROR: Selective rows cover entire matrix.\n";
        return false;
    }

    std::vector<int> idxB(mSelectiveRows.begin(), mSelectiveRows.end());
    std::sort(idxB.begin(), idxB.end());

    std::vector<int> globalToLocal(nTotal, -1);
    std::vector<int> globalToI(nTotal, -1);

    for (int k = 0; k < nb; ++k) {
        globalToLocal[idxB[k]] = k;
    }

    int iCount = 0;
    for (int g = 0; g < nTotal; ++g) {
        if (globalToLocal[g] == -1) {
            globalToI[g] = iCount++;
        }
    }

    Eigen::MatrixXd A_bb(nb, nb);
    A_bb.setZero();
    for (int k = 0; k < mA_eigen.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mA_eigen, k); it; ++it) {
            const int lr = globalToLocal[it.row()];
            const int lc = globalToLocal[it.col()];
            if (lr >= 0 && lc >= 0) {
                A_bb(lr, lc) += it.value();
            }
        }
    }

    Eigen::MatrixXd A_ib(ni, nb);
    A_ib.setZero();
    for (int k = 0; k < mA_eigen.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mA_eigen, k); it; ++it) {
            const int ir = globalToI[it.row()];
            const int lc = globalToLocal[it.col()];
            if (ir >= 0 && lc >= 0) {
                A_ib(ir, lc) += it.value();
            }
        }
    }

    const Eigen::MatrixXd A_bi = A_ib.transpose();
    Eigen::MatrixXd A_iiInvAib(ni, nb);
    A_iiInvAib.setZero();

    std::vector<Eigen::Triplet<double>> iiTriplets;
    for (int k = 0; k < mA_eigen.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mA_eigen, k); it; ++it) {
            const int ir = globalToI[it.row()];
            const int ic = globalToI[it.col()];
            if (ir >= 0 && ic >= 0) {
                iiTriplets.emplace_back(ir, ic, it.value());
            }
        }
    }

    Eigen::SparseMatrix<double> A_ii(ni, ni);
    A_ii.setFromTriplets(iiTriplets.begin(), iiTriplets.end());
    A_ii.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solverIi;
    solverIi.compute(A_ii);
    if (solverIi.info() != Eigen::Success) {
        std::cout << "WARNING: SparseLU failed for A_ii, trying BiCGSTAB...\n";
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>,
                        Eigen::IncompleteLUT<double>> iterSolver;
        iterSolver.compute(A_ii);
        for (int col = 0; col < nb; ++col) {
            A_iiInvAib.col(col) = iterSolver.solve(A_ib.col(col));
        }
    } else {
        for (int col = 0; col < nb; ++col) {
            A_iiInvAib.col(col) = solverIi.solve(A_ib.col(col));
        }
    }

    const Eigen::MatrixXd S = A_bb - A_bi * A_iiInvAib;
    const Eigen::MatrixXd SInv = S.inverse();

    mSelectiveInvResult.resize(nTotal, nTotal);
    mSelectiveInvResult.setZero();

    for (int kr = 0; kr < nb; ++kr) {
        for (int kc = 0; kc < nb; ++kc) {
            mSelectiveInvResult.coeffRef(idxB[kr], idxB[kc]) = SInv(kr, kc);
        }
    }

    const Eigen::MatrixXd AiibSinv = A_iiInvAib * SInv;
    for (int gi = 0; gi < nTotal; ++gi) {
        const int ir = globalToI[gi];
        if (ir >= 0) {
            for (int kc = 0; kc < nb; ++kc) {
                mSelectiveInvResult.coeffRef(gi, idxB[kc]) = -AiibSinv(ir, kc);
            }
        }
    }

    if (mIsSymmetric) {
        for (int gi = 0; gi < nTotal; ++gi) {
            const int ir = globalToI[gi];
            if (ir >= 0) {
                for (int kc = 0; kc < nb; ++kc) {
                    mSelectiveInvResult.coeffRef(idxB[kc], gi) = -AiibSinv(ir, kc);
                }
            }
        }
    }

    return true;
}
