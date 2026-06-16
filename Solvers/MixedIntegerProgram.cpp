#include "MixedIntegerProgram.h"
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Eigen;

MixedIntegerProgram::MixedIntegerProgram(
    int numVariables,
    std::vector<Constraint> constraints,
    double integerScale)
    : nVars_(numVariables)
    , integerScale_(integerScale)
    , constraints_(std::move(constraints))
{
    nConstraints_ = static_cast<int>(constraints_.size());
    for (int i = 0; i < nConstraints_; ++i) {
        constraints_[static_cast<size_t>(i)].idx = i;
    }
}

void MixedIntegerProgram::setIntegerBounds(int lo, int hi)
{
    integerLo_ = lo;
    integerHi_ = hi;
}

int MixedIntegerProgram::roundToInt(double x)
{
    return static_cast<int>(std::round(x));
}

void MixedIntegerProgram::buildNormalMatrix(SparseMatrix<double>& A) const
{
    A.resize(nVars_, nVars_);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nVars_);

    for (const auto& c : constraints_) {
        if (c.i < 0 || c.j < 0) continue;
        trips.emplace_back(c.i, c.j, -1.0);
        trips.emplace_back(c.j, c.i, -1.0);
        diag(c.i) += 1.0;
        diag(c.j) += 1.0;
    }
    // Pin the first variable to eliminate the Laplacian nullspace (constant shift).
    // Without this, the system L·x = b is singular and the returned solution is
    // defined only up to a global additive constant.
    diag(0) += 1e6;
    for (int i = 0; i < nVars_; ++i) {
        trips.emplace_back(i, i, diag(i) + 1e-8);
    }
    A.setFromTriplets(trips.begin(), trips.end());
}

double MixedIntegerProgram::energy(const VectorXd& variables, const VectorXi& integers) const
{
    double E = 0.0;
    for (const auto& c : constraints_) {
        if (c.i < 0 || c.j < 0) continue;
        const double r = variables(c.i) - variables(c.j) + c.kappa +
                         integerScale_ * static_cast<double>(integers(c.idx));
        E += r * r;
    }
    return E;
}

bool MixedIntegerProgram::solveVariablesGivenIntegers(
    const VectorXi& integers,
    VectorXd& outVariables) const
{
    if (!variableSolverReady_) {
        variableSolver_.compute(normalMatrix_);
        if (variableSolver_.info() != Success) return false;
        variableSolverReady_ = true;
    }

    VectorXd b = VectorXd::Zero(nVars_);
    for (const auto& c : constraints_) {
        if (c.i < 0 || c.j < 0) continue;
        const double z = static_cast<double>(integers(c.idx));
        const double offset = integerScale_ * z + c.kappa;
        b(c.i) -= offset;   //  from ∂/∂x_i  of  (x_i − x_j + κ + s·z)^2
        b(c.j) += offset;
    }
    // Pin variable 0 (its diagonal was boosted in buildNormalMatrix).
    b(0) = 0.0;

    outVariables = variableSolver_.solve(b);
    return variableSolver_.info() == Success && outVariables.size() == nVars_;
}

bool MixedIntegerProgram::optimizeAlternating(VectorXd& variables, VectorXi& integers)
{
    buildNormalMatrix(normalMatrix_);
    variableSolverReady_ = false;
    variableSolver_.compute(normalMatrix_);
    if (variableSolver_.info() != Success) return false;
    variableSolverReady_ = true;

    VectorXd newVariables(nVars_);

    for (int iter = 0; iter < alternatingIters_; ++iter) {
        if (!solveVariablesGivenIntegers(integers, newVariables)) return false;

        for (const auto& c : constraints_) {
            if (c.i < 0 || c.j < 0) continue;
            // Optimal integer given current continuous variables:
            //   z* = round( (x_j − x_i − κ_ij) / scale )
            const double diff = newVariables(c.j) - newVariables(c.i);
            int z = roundToInt((diff - c.kappa) / integerScale_);
            z = std::max(integerLo_, std::min(integerHi_, z));
            integers(c.idx) = z;
        }

        const double change = (newVariables - variables).norm();
        variables = newVariables;
        if (change < 1e-6) break;
    }

    return true;
}

void MixedIntegerProgram::refineIntegersCoordinateDescent(VectorXd& variables, VectorXi& integers)
{
    VectorXd trialVariables(nVars_);

    for (int pass = 0; pass < refinePasses_; ++pass) {
        bool improved = false;

        for (const auto& c : constraints_) {
            if (c.i < 0 || c.j < 0) continue;

            const int saved = integers(c.idx);
            int bestZ = saved;
            double bestE = std::numeric_limits<double>::infinity();

            for (int cand = integerLo_; cand <= integerHi_; ++cand) {
                integers(c.idx) = cand;
                if (!solveVariablesGivenIntegers(integers, trialVariables)) continue;
                const double E = energy(trialVariables, integers);
                if (E < bestE) {
                    bestE = E;
                    bestZ = cand;
                }
            }

            integers(c.idx) = bestZ;
            if (bestZ != saved) improved = true;
        }

        solveVariablesGivenIntegers(integers, variables);
        if (!improved) break;
    }
}

bool MixedIntegerProgram::solve(VectorXd& variables, VectorXi& integers)
{
    if (nVars_ < 1 || variables.size() != nVars_) return false;
    if (integerScale_ == 0.0) return false;
    if (integers.size() != nConstraints_) {
        integers.resize(nConstraints_);
        integers.setZero();
    }

    if (!optimizeAlternating(variables, integers)) return false;
    refineIntegersCoordinateDescent(variables, integers);
    return true;
}
