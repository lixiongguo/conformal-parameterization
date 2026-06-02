#include "AbelJacobiParameterization.h"
#include <iostream>

AbelJacobiParameterization::AbelJacobiParameterization(Mesh& mesh0)
    : Parameterization(mesh0)
    , aj_(mesh)
{
}

void AbelJacobiParameterization::setSingularities(
    const std::vector<AbelJacobi::SingularPoint>& singularities)
{
    singularities_ = singularities;
}

AbelJacobi::DivisorCheck AbelJacobiParameterization::checkSingularities() const
{
    return aj_.checkDivisor(baseVertex_, singularities_);
}

void AbelJacobiParameterization::parameterize()
{
    if (!aj_.build()) {
        std::cerr << "[AbelJacobiParameterization] build failed, falling back to LSCM.\n";
        Lscm fallback(mesh);
        fallback.parameterize();
        return;
    }

    if (!singularities_.empty()) {
        const auto check = checkSingularities();
        if (!check.poincareHopfOk) {
            std::cerr << "[AbelJacobiParameterization] Warning: Poincare-Hopf constraint not satisfied.\n";
        }
        if (!check.abelJacobiOk) {
            std::cerr << "[AbelJacobiParameterization] Warning: Abel-Jacobi divisor not principal "
                      << "(residual=" << check.latticeResidual << ").\n";
        }
    }

    if (!aj_.quantizeAndParameterize(baseVertex_)) {
        std::cerr << "[AbelJacobiParameterization] quantization failed, falling back to LSCM.\n";
        Lscm fallback(mesh);
        fallback.parameterize();
        return;
    }

    normalize();
}
