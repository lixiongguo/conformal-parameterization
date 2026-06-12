#include "CrossFieldIntegerProgram.h"
#include <cmath>
#include <limits>
#include <algorithm>

using namespace Eigen;

CrossFieldIntegerProgram::CrossFieldIntegerProgram(
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

void CrossFieldIntegerProgram::setIntegerBounds(int lo, int hi)
{
    integerLo_ = lo;
    integerHi_ = hi;
}

int CrossFieldIntegerProgram::roundToInt(double x)
{
    return static_cast<int>(std::floor(x + 0.5));
}

void CrossFieldIntegerProgram::buildDifferenceLaplacian(SparseMatrix<double>& L) const
{
    L.resize(nVars_, nVars_);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nVars_);

    for (const auto& c : constraints_) {
        if (c.i < 0 || c.j < 0) continue;
        trips.emplace_back(c.i, c.j, -1.0);
        trips.emplace_back(c.j, c.i, -1.0);
        diag(c.i) += 1.0;
        diag(c.j) += 1.0;
    }
    for (int i = 0; i < nVars_; ++i) {
        trips.emplace_back(i, i, diag(i) + 1e-8);
    }
    L.setFromTriplets(trips.begin(), trips.end());
}

double CrossFieldIntegerProgram::energy(const VectorXd& variables, const VectorXi& integers) const
{
    double E = 0.0;
    for (const auto& c : constraints_) {
        if (c.i < 0 || c.j < 0) continue;
        const double r = variables(c.i) - variables(c.j) +
                         integerScale_ * static_cast<double>(integers(c.idx));
        E += r * r;
    }
    return E;
}

bool CrossFieldIntegerProgram::solveVariablesGivenIntegers(
    const VectorXi& integers,
    VectorXd& outVariables) const
{
    if (!variableSolverReady_) {
        variableSolver_.compute(variableLaplacian_);
        if (variableSolver_.info() != Success) return false;
        variableSolverReady_ = true;
    }

    VectorXd b = VectorXd::Zero(nVars_);
    for (const auto& c : constraints_) {
        if (c.i < 0 || c.j < 0) continue;
        const double z = static_cast<double>(integers(c.idx));
        b(c.i) -= integerScale_ * z;
        b(c.j) += integerScale_ * z;
    }

    outVariables = variableSolver_.solve(b);
    return variableSolver_.info() == Success && outVariables.size() == nVars_;
}

bool CrossFieldIntegerProgram::optimizeAlternating(VectorXd& variables, VectorXi& integers)
{
    buildDifferenceLaplacian(variableLaplacian_);
    variableSolverReady_ = false;
    variableSolver_.compute(variableLaplacian_);
    if (variableSolver_.info() != Success) return false;
    variableSolverReady_ = true;

    VectorXd newVariables(nVars_);

    for (int iter = 0; iter < alternatingIters_; ++iter) {
        if (!solveVariablesGivenIntegers(integers, newVariables)) return false;

        for (const auto& c : constraints_) {
            if (c.i < 0 || c.j < 0) continue;
            const double diff = newVariables(c.j) - newVariables(c.i);
            int z = roundToInt(diff / integerScale_);
            z = std::max(integerLo_, std::min(integerHi_, z));
            integers(c.idx) = z;
        }

        const double change = (newVariables - variables).norm();
        variables = newVariables;
        if (change < 1e-6) break;
    }

    return true;
}

void CrossFieldIntegerProgram::refineIntegersCoordinateDescent(VectorXd& variables, VectorXi& integers)
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

bool CrossFieldIntegerProgram::solve(VectorXd& variables, VectorXi& integers)
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
