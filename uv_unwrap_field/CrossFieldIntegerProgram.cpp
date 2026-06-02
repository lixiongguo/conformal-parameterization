#include "CrossFieldIntegerProgram.h"
#include <cmath>
#include <limits>
#include <algorithm>

using namespace Eigen;

namespace {
constexpr double kQuarterTurn = M_PI / 2.0;
} // namespace

CrossFieldIntegerProgram::CrossFieldIntegerProgram(int numFaces, std::vector<Edge> internalEdges)
    : nF_(numFaces)
    , edges_(std::move(internalEdges))
{
    nE_ = static_cast<int>(edges_.size());
    for (int i = 0; i < nE_; ++i) {
        edges_[static_cast<size_t>(i)].idx = i;
    }
}

void CrossFieldIntegerProgram::setJumpBounds(int lo, int hi)
{
    jumpLo_ = lo;
    jumpHi_ = hi;
}

int CrossFieldIntegerProgram::roundToInt(double x)
{
    return static_cast<int>(std::floor(x + 0.5));
}

void CrossFieldIntegerProgram::buildFaceLaplacian(SparseMatrix<double>& L) const
{
    L.resize(nF_, nF_);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nF_);

    for (const auto& e : edges_) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        trips.emplace_back(e.f1, e.f2, -1.0);
        trips.emplace_back(e.f2, e.f1, -1.0);
        diag(e.f1) += 1.0;
        diag(e.f2) += 1.0;
    }
    for (int i = 0; i < nF_; ++i) {
        trips.emplace_back(i, i, diag(i) + 1e-8);
    }
    L.setFromTriplets(trips.begin(), trips.end());
}

double CrossFieldIntegerProgram::energy(const VectorXd& theta, const VectorXi& jump) const
{
    double E = 0.0;
    for (const auto& e : edges_) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        const double r = theta(e.f1) - theta(e.f2) + kQuarterTurn * static_cast<double>(jump(e.idx));
        E += r * r;
    }
    return E;
}

bool CrossFieldIntegerProgram::solveThetaGivenJumps(
    const VectorXi& jump,
    VectorXd& outTheta) const
{
    if (!faceSolverReady_) {
        faceSolver_.compute(faceLaplacian_);
        if (faceSolver_.info() != Success) return false;
        faceSolverReady_ = true;
    }

    VectorXd b = VectorXd::Zero(nF_);
    for (const auto& e : edges_) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        const double p = static_cast<double>(jump(e.idx));
        b(e.f1) -= kQuarterTurn * p;
        b(e.f2) += kQuarterTurn * p;
    }

    outTheta = faceSolver_.solve(b);
    return faceSolver_.info() == Success && outTheta.size() == nF_;
}

void CrossFieldIntegerProgram::wrapThetaToQuarterTurn(VectorXd& theta) const
{
    for (int i = 0; i < nF_; ++i) {
        theta(i) = std::fmod(theta(i), kQuarterTurn);
        if (theta(i) < 0) theta(i) += kQuarterTurn;
    }
}

bool CrossFieldIntegerProgram::optimizeAlternating(VectorXd& theta, VectorXi& jump)
{
    buildFaceLaplacian(faceLaplacian_);
    faceSolverReady_ = false;
    faceSolver_.compute(faceLaplacian_);
    if (faceSolver_.info() != Success) return false;
    faceSolverReady_ = true;

    VectorXd newTheta(nF_);

    for (int iter = 0; iter < alternatingIters_; ++iter) {
        if (!solveThetaGivenJumps(jump, newTheta)) return false;

        for (const auto& e : edges_) {
            if (e.f1 < 0 || e.f2 < 0) continue;
            const double diff = newTheta(e.f2) - newTheta(e.f1);
            int p = roundToInt(diff / kQuarterTurn);
            p = std::max(jumpLo_, std::min(jumpHi_, p));
            jump(e.idx) = p;
        }

        const double change = (newTheta - theta).norm();
        theta = newTheta;
        if (change < 1e-6) break;
    }

    wrapThetaToQuarterTurn(theta);
    return true;
}

void CrossFieldIntegerProgram::refineJumpsCoordinateDescent(VectorXd& theta, VectorXi& jump)
{
    VectorXd trialTheta(nF_);

    for (int pass = 0; pass < refinePasses_; ++pass) {
        bool improved = false;

        for (const auto& e : edges_) {
            if (e.f1 < 0 || e.f2 < 0) continue;

            const int saved = jump(e.idx);
            int bestP = saved;
            double bestE = std::numeric_limits<double>::infinity();

            for (int cand = jumpLo_; cand <= jumpHi_; ++cand) {
                jump(e.idx) = cand;
                if (!solveThetaGivenJumps(jump, trialTheta)) continue;
                const double E = energy(trialTheta, jump);
                if (E < bestE) {
                    bestE = E;
                    bestP = cand;
                }
            }

            jump(e.idx) = bestP;
            if (bestP != saved) improved = true;
        }

        solveThetaGivenJumps(jump, theta);
        if (!improved) break;
    }

    wrapThetaToQuarterTurn(theta);
}

bool CrossFieldIntegerProgram::solve(VectorXd& theta, VectorXi& jump)
{
    if (nF_ < 1 || theta.size() != nF_) return false;
    if (jump.size() != nE_) {
        jump.resize(nE_);
        jump.setZero();
    }

    if (!optimizeAlternating(theta, jump)) return false;
    refineJumpsCoordinateDescent(theta, jump);
    return true;
}
