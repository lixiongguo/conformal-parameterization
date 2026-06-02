#ifndef ABEL_JACOBI_PARAMETERIZATION_H
#define ABEL_JACOBI_PARAMETERIZATION_H

#include "../Parameterization.h"
#include "AbelJacobi.h"
#include "../uv_unwrap_simple/Lscm.h"
#include <vector>

/**
 * Global conformal parameterization via Abel–Jacobi / period lattice (closed, g >= 1).
 * Falls back to LSCM on disk topology or genus 0.
 */
class AbelJacobiParameterization : public Parameterization {
public:
    AbelJacobiParameterization(Mesh& mesh0);

    void parameterize() override;

    void setBaseVertex(int v) { baseVertex_ = v; }
    void setSingularities(const std::vector<AbelJacobi::SingularPoint>& singularities);

    /** Run Poincaré–Hopf + Abel divisor check before parameterization. */
    AbelJacobi::DivisorCheck checkSingularities() const;

    const AbelJacobi& engine() const { return aj_; }
    AbelJacobi& engine() { return aj_; }

private:
    AbelJacobi aj_;
    int baseVertex_ = 0;
    std::vector<AbelJacobi::SingularPoint> singularities_;
};

#endif
