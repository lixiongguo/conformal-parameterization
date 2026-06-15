#include <emscripten/bind.h>
#include "Solver.h"
#include "AugmentedLagrangian.h"
#include "MixedIntegerProgram.h"
#include <Eigen/Dense>
#include <vector>

using namespace emscripten;

// ============================================================
//  Solver wrappers: convert Eigen::VectorXd ↔ std::vector<double>
// ============================================================

std::vector<double> solverGetX(Solver& s) {
    return std::vector<double>(s.x.data(), s.x.data() + s.x.size());
}

void solverSetX(Solver& s, const std::vector<double>& xVec) {
    s.x = Eigen::Map<const Eigen::VectorXd>(xVec.data(), xVec.size());
}

void solverSetN(Solver& s, int n) { s.n = n; }
int  solverGetN(const Solver& s)  { return s.n; }

std::vector<double> solverGetObj(const Solver& s) { return s.obj; }

// MeshHandle pointer accessor (raw pointer — use with caution from JS)
MeshHandle* solverGetHandle(Solver& s) { return s.handle; }
void solverSetHandle(Solver& s, MeshHandle* h) { s.handle = h; }

// ============================================================
//  Augmented Lagrangian wrappers
// ============================================================
// solve(x, beta) — need to convert Eigen types
bool alSolve(AugmentedLagrangian& al,
             const std::vector<double>& xInit,
             const std::vector<double>& betaVec,
             std::vector<double>& xOut)
{
    Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(xInit.data(), xInit.size());
    Eigen::VectorXd beta = Eigen::Map<const Eigen::VectorXd>(betaVec.data(), betaVec.size());
    bool ok = al.solve(x, beta);
    xOut.assign(x.data(), x.data() + x.size());
    return ok;
}

// ============================================================
//  MIP wrappers
// ============================================================
bool mipSolve(MixedIntegerProgram& mip,
              const std::vector<double>& vars,
              const std::vector<int>& ints,
              std::vector<double>& outVars,
              std::vector<int>& outInts)
{
    Eigen::VectorXd variables = Eigen::Map<const Eigen::VectorXd>(vars.data(), vars.size());
    Eigen::VectorXi integers(ints.size());
    for (size_t i = 0; i < ints.size(); ++i) integers(i) = ints[i];

    bool ok = mip.solve(variables, integers);

    outVars.assign(variables.data(), variables.data() + variables.size());
    outInts.assign(integers.data(), integers.data() + integers.size());
    return ok;
}

double mipEnergy(MixedIntegerProgram& mip,
                 const std::vector<double>& vars,
                 const std::vector<int>& ints)
{
    Eigen::VectorXd variables = Eigen::Map<const Eigen::VectorXd>(vars.data(), vars.size());
    Eigen::VectorXi integers(ints.size());
    for (size_t i = 0; i < ints.size(); ++i) integers(i) = ints[i];
    return mip.energy(variables, integers);
}

// ============================================================
//  Embind bindings
// ============================================================
EMSCRIPTEN_BINDINGS(SolverModule) {

    class_<MeshHandle>("MeshHandle")
        .constructor<>();

    class_<Solver>("Solver")
        .constructor<>()
        .function("getX",      &solverGetX)
        .function("setX",      &solverSetX)
        .function("getN",      &solverGetN)
        .function("setN",      &solverSetN)
        .function("getObj",    &solverGetObj)
        .function("getHandle", &solverGetHandle, allow_raw_pointers())
        .function("setHandle", &solverSetHandle, allow_raw_pointers())
        .function("gradientDescent", &Solver::gradientDescent)
        .function("newton",    &Solver::newton)
        .function("trustRegion", &Solver::trustRegion)
        .function("lbfgs", emscripten::optional_override(
            [](Solver& s, int m) { s.lbfgs(m); }));

    // --- Augmented Lagrangian ---
    class_<AugmentedLagrangianConfig>("AugmentedLagrangianConfig")
        .constructor<>()
        .property("rhoInit",   &AugmentedLagrangianConfig::rhoInit)
        .property("rhoMax",    &AugmentedLagrangianConfig::rhoMax)
        .property("rhoScale",  &AugmentedLagrangianConfig::rhoScale)
        .property("tolOuter",  &AugmentedLagrangianConfig::tolOuter)
        .property("tolInner",  &AugmentedLagrangianConfig::tolInner)
        .property("maxOuter",  &AugmentedLagrangianConfig::maxOuter)
        .property("maxInner",  &AugmentedLagrangianConfig::maxInner)
        .property("positiveConstraint", &AugmentedLagrangianConfig::positiveConstraint);

    class_<AugmentedLagrangian>("AugmentedLagrangian")
        .function("solve", &alSolve);

    // --- Mixed Integer Program ---
    class_<MixedIntegerProgram::Constraint>("MIPConstraint")
        .constructor<>()
        .property("i",   &MixedIntegerProgram::Constraint::i)
        .property("j",   &MixedIntegerProgram::Constraint::j)
        .property("idx", &MixedIntegerProgram::Constraint::idx);

    register_vector<MixedIntegerProgram::Constraint>("MIPConstraintVector");

    class_<MixedIntegerProgram>("MixedIntegerProgram")
        .constructor<int, std::vector<MixedIntegerProgram::Constraint>, double>()
        .function("solve",  &mipSolve)
        .function("energy", &mipEnergy)
        .function("setIntegerBounds", &MixedIntegerProgram::setIntegerBounds)
        .function("setAlternatingIterations", &MixedIntegerProgram::setAlternatingIterations)
        .function("setRefinePasses", &MixedIntegerProgram::setRefinePasses)
        .function("numVariables", &MixedIntegerProgram::numVariables)
        .function("numConstraints", &MixedIntegerProgram::numConstraints);
}
