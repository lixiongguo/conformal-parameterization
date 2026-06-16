#include "BoundedDistortionHarmonicMap.hpp"

#include "../Solvers/Solver.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace bdhm {
namespace {

constexpr double kEps = 1e-12;

struct State {
    VecC phi;
    VecC psi;
};

struct EnergyGradient {
    double energy = 0.0;
    VecC grad_phi;
    VecC grad_psi;
};

void validateInput(const VecC& boundary_polygon, const VecC& energy_samples,
                   const std::vector<PointConstraint>& constraints) {
    if (boundary_polygon.size() < 3) {
        throw BdhmError("boundary_polygon must have at least three vertices");
    }
    if (energy_samples.size() == 0) {
        throw BdhmError("energy_samples must not be empty");
    }
    if (constraints.empty()) {
        throw BdhmError("at least one point constraint is required");
    }
}

Complex linearRowValue(const MatC& M, int row, const VecC& x) {
    Complex out(0.0, 0.0);
    for (int j = 0; j < M.cols(); ++j) {
        out += M(row, j) * x[j];
    }
    return out;
}

void accumulateLinearGradient(const MatC& M, int row, const Complex& grad_y, VecC& grad_x) {
    for (int j = 0; j < M.cols(); ++j) {
        const Complex a = M(row, j);
        const double gr = grad_y.real();
        const double gi = grad_y.imag();
        grad_x[j] += Complex(a.real() * gr + a.imag() * gi,
                             -a.imag() * gr + a.real() * gi);
    }
}

EnergyGradient energyAndGradient(
    const MatC& C_energy,
    const MatC& D_energy,
    const MatC& C_constraints,
    const std::vector<PointConstraint>& constraints,
    const State& state,
    const VecC& reference_phi,
    const Options& options) {
    const double kappa = distortionToKappa(options.distortion_bound);
    EnergyGradient out;
    out.grad_phi = VecC::Zero(state.phi.size());
    out.grad_psi = VecC::Zero(state.psi.size());

    const VecC fz = D_energy * state.phi;
    const VecC gz = D_energy * state.psi;

    for (int i = 0; i < fz.size(); ++i) {
        const double a_abs = std::abs(fz[i]);
        const double b_abs = std::abs(gz[i]);

        const double cone = b_abs - kappa * a_abs;
        if (cone > 0.0) {
            out.energy += options.distortion_penalty * cone * cone;
            Complex grad_a(0.0, 0.0);
            Complex grad_b(0.0, 0.0);
            if (a_abs > kEps) {
                grad_a = options.distortion_penalty * (-2.0 * kappa * cone / a_abs) * fz[i];
            }
            if (b_abs > kEps) {
                grad_b = options.distortion_penalty * (2.0 * cone / b_abs) * gz[i];
            }
            accumulateLinearGradient(D_energy, i, grad_a, out.grad_phi);
            accumulateLinearGradient(D_energy, i, grad_b, out.grad_psi);
        }

        const double jacobian = a_abs * a_abs - b_abs * b_abs;
        const double jacobian_violation = options.min_jacobian - jacobian;
        if (jacobian_violation > 0.0) {
            out.energy += options.jacobian_penalty * jacobian_violation * jacobian_violation;
            const Complex grad_a = options.jacobian_penalty * (-4.0 * jacobian_violation) * fz[i];
            const Complex grad_b = options.jacobian_penalty * (4.0 * jacobian_violation) * gz[i];
            accumulateLinearGradient(D_energy, i, grad_a, out.grad_phi);
            accumulateLinearGradient(D_energy, i, grad_b, out.grad_psi);
        }
    }

    for (int i = 0; i < C_constraints.rows(); ++i) {
        const Complex h = linearRowValue(C_constraints, i, state.phi);
        const Complex g = linearRowValue(C_constraints, i, state.psi);
        const Complex mapped = h + std::conj(g);
        const Complex diff = mapped - constraints[static_cast<std::size_t>(i)].target;
        const double weight = options.position_weight * constraints[static_cast<std::size_t>(i)].weight;

        out.energy += weight * std::norm(diff);
        accumulateLinearGradient(C_constraints, i, 2.0 * weight * diff, out.grad_phi);
        accumulateLinearGradient(C_constraints, i, 2.0 * weight * std::conj(diff), out.grad_psi);
    }

    if (options.anti_holomorphic_weight > 0.0) {
        out.energy += options.anti_holomorphic_weight * state.psi.squaredNorm();
        out.grad_psi += 2.0 * options.anti_holomorphic_weight * state.psi;
    }

    if (options.phi_reference_weight > 0.0) {
        const VecC diff = state.phi - reference_phi;
        out.energy += options.phi_reference_weight * diff.squaredNorm();
        out.grad_phi += 2.0 * options.phi_reference_weight * diff;
    }

    return out;
}

double energyOnly(
    const MatC& C_energy,
    const MatC& D_energy,
    const MatC& C_constraints,
    const std::vector<PointConstraint>& constraints,
    const State& state,
    const VecC& reference_phi,
    const Options& options) {
    return energyAndGradient(C_energy, D_energy, C_constraints, constraints, state, reference_phi, options).energy;
}

VecC constraintsToSources(const std::vector<PointConstraint>& constraints) {
    VecC sources(static_cast<int>(constraints.size()));
    for (int i = 0; i < sources.size(); ++i) {
        sources[i] = constraints[static_cast<std::size_t>(i)].source;
    }
    return sources;
}

Eigen::VectorXd stateToVector(const State& state) {
    const int n = static_cast<int>(state.phi.size());
    Eigen::VectorXd x(4 * n);
    for (int i = 0; i < n; ++i) {
        x[i] = state.phi[i].real();
        x[n + i] = state.phi[i].imag();
        x[2 * n + i] = state.psi[i].real();
        x[3 * n + i] = state.psi[i].imag();
    }
    return x;
}

State vectorToState(const Eigen::VectorXd& x) {
    const int n = static_cast<int>(x.size() / 4);
    State state;
    state.phi.resize(n);
    state.psi.resize(n);
    for (int i = 0; i < n; ++i) {
        state.phi[i] = Complex(x[i], x[n + i]);
        state.psi[i] = Complex(x[2 * n + i], x[3 * n + i]);
    }
    return state;
}

Eigen::VectorXd gradientToVector(const EnergyGradient& eg) {
    const int n = static_cast<int>(eg.grad_phi.size());
    Eigen::VectorXd gradient(4 * n);
    for (int i = 0; i < n; ++i) {
        gradient[i] = eg.grad_phi[i].real();
        gradient[n + i] = eg.grad_phi[i].imag();
        gradient[2 * n + i] = eg.grad_psi[i].real();
        gradient[3 * n + i] = eg.grad_psi[i].imag();
    }
    return gradient;
}

}  // namespace

double distortionToKappa(double distortion_bound) {
    if (distortion_bound < 1.0) {
        throw BdhmError("distortion_bound must be >= 1");
    }
    return (distortion_bound - 1.0) / (distortion_bound + 1.0);
}

VecC evaluateMap(const MatC& C, const VecC& phi, const VecC& psi) {
    return C * phi + (C * psi).conjugate();
}

Stats computeStats(const MatC& D, const VecC& phi, const VecC& psi, double distortion_bound) {
    const double kappa = distortionToKappa(distortion_bound);
    const VecC fz = D * phi;
    const VecC gz = D * psi;

    Stats stats;
    stats.max_distortion = 0.0;
    stats.min_jacobian = std::numeric_limits<double>::infinity();
    stats.max_cone_violation = 0.0;

    for (int i = 0; i < fz.size(); ++i) {
        const double a = std::abs(fz[i]);
        const double b = std::abs(gz[i]);
        const double jacobian = a * a - b * b;
        const double sigma_min = a - b;
        const double distortion = sigma_min > kEps ? (a + b) / sigma_min
                                                   : std::numeric_limits<double>::infinity();

        stats.max_distortion = std::max(stats.max_distortion, distortion);
        stats.min_jacobian = std::min(stats.min_jacobian, jacobian);
        stats.max_cone_violation = std::max(stats.max_cone_violation, std::max(0.0, b - kappa * a));
    }

    return stats;
}

SolveResult solve(
    const VecC& boundary_polygon,
    const VecC& energy_samples,
    const std::vector<PointConstraint>& constraints,
    const Options& options) {
    validateInput(boundary_polygon, energy_samples, constraints);

    const MatC C_energy = cauchy::computeCauchyCoordinates(boundary_polygon, energy_samples);
    const MatC D_energy = cauchy::computeCauchyCoordinateDerivatives(boundary_polygon, energy_samples);
    const VecC constraint_sources = constraintsToSources(constraints);
    const MatC C_constraints = cauchy::computeCauchyCoordinates(boundary_polygon, constraint_sources);

    State state;
    state.phi = boundary_polygon;
    state.psi = VecC::Zero(boundary_polygon.size());
    const VecC reference_phi = state.phi;

    MeshHandle handle;
    handle.computeEnergy = [&](double& energy, const Eigen::VectorXd& x) {
        const State eval_state = vectorToState(x);
        energy = energyOnly(C_energy, D_energy, C_constraints, constraints, eval_state, reference_phi, options);
    };
    handle.computeGradient = [&](Eigen::VectorXd& gradient, const Eigen::VectorXd& x) {
        const State eval_state = vectorToState(x);
        const EnergyGradient eg =
            energyAndGradient(C_energy, D_energy, C_constraints, constraints, eval_state, reference_phi, options);
        gradient = gradientToVector(eg);
    };

    Solver solver;
    solver.handle = &handle;
    solver.x = stateToVector(state);
    solver.n = static_cast<int>(solver.x.size());
    solver.maxIterations = options.max_iterations;
    solver.gradientTolerance = options.gradient_tolerance;
    solver.initialStep = options.initial_step;
    solver.useInitialGuess = true;
    solver.verbose = false;
    if (options.solver_method == SolverMethod::LBFGS) {
        solver.lbfgs(options.lbfgs_history);
    } else {
        solver.gradientDescent();
    }

    state = vectorToState(solver.x);
    const double final_energy = solver.obj.empty() ? std::numeric_limits<double>::infinity() : solver.obj.back();
    const int iterations = solver.obj.empty() ? 0 : static_cast<int>(solver.obj.size()) - 1;

    SolveResult result;
    result.phi = state.phi;
    result.psi = state.psi;
    result.mapped_energy_samples = evaluateMap(C_energy, state.phi, state.psi);
    result.mapped_constraints = evaluateMap(C_constraints, state.phi, state.psi);
    result.stats = computeStats(D_energy, state.phi, state.psi, options.distortion_bound);
    result.final_energy = final_energy;
    result.iterations = iterations;
    return result;
}

NloP2PHarmonicResult nloP2PHarmonic(
    const VecC& boundary_polygon,
    const VecC& energy_samples,
    const std::vector<PointConstraint>& constraints,
    const NloP2PHarmonicInput& options) {
    return solve(boundary_polygon, energy_samples, constraints, options);
}

}  // namespace bdhm
