#include "BoundedDistortionHarmonicMap.hpp"

#include <iostream>
#include <vector>

int main() {
    using namespace bdhm;

    VecC boundary(4);
    boundary << Complex(0.0, 0.0),
                Complex(1.0, 0.0),
                Complex(1.0, 1.0),
                Complex(0.0, 1.0);

    std::vector<Complex> samples;
    for (int iy = 1; iy <= 5; ++iy) {
        for (int ix = 1; ix <= 5; ++ix) {
            samples.emplace_back(ix / 6.0, iy / 6.0);
        }
    }

    VecC energy_samples(static_cast<int>(samples.size()));
    for (int i = 0; i < energy_samples.size(); ++i) {
        energy_samples[i] = samples[static_cast<std::size_t>(i)];
    }

    std::vector<PointConstraint> constraints = {
        {Complex(0.5, 0.5), Complex(0.62, 0.42), 1.0},
        {Complex(0.3, 0.7), Complex(0.28, 0.74), 0.5},
    };

    NloP2PHarmonicInput options;
    options.distortion_bound = 2.0;
    options.min_jacobian = 1e-5;
    options.position_weight = 1e5;
    options.distortion_penalty = 1e3;
    options.jacobian_penalty = 1e3;
    options.anti_holomorphic_weight = 1e-3;
    options.phi_reference_weight = 1e-3;
    options.initial_step = 1e-3;
    options.max_iterations = 1200;

    const NloP2PHarmonicResult result = nloP2PHarmonic(boundary, energy_samples, constraints, options);

    std::cout << "iterations          : " << result.iterations << '\n';
    std::cout << "final energy        : " << result.final_energy << '\n';
    std::cout << "max distortion      : " << result.stats.max_distortion << '\n';
    std::cout << "min jacobian        : " << result.stats.min_jacobian << '\n';
    std::cout << "max cone violation  : " << result.stats.max_cone_violation << '\n';

    double max_constraint_error = 0.0;
    for (int i = 0; i < result.mapped_constraints.size(); ++i) {
        const Complex target = constraints[static_cast<std::size_t>(i)].target;
        const double err = std::abs(result.mapped_constraints[i] - target);
        max_constraint_error = std::max(max_constraint_error, err);
        std::cout << "constraint " << i
                  << " mapped=(" << result.mapped_constraints[i].real() << ", "
                  << result.mapped_constraints[i].imag() << ")"
                  << " target=(" << target.real() << ", " << target.imag() << ")"
                  << " err=" << err << '\n';
    }

    std::cout << "max constraint error: " << max_constraint_error << '\n';

    const bool ok = result.stats.min_jacobian > 0.0 &&
                    result.stats.max_distortion <= options.distortion_bound + 0.25 &&
                    max_constraint_error < 5e-2;
    std::cout << (ok ? "PASS\n" : "FAIL\n");
    return ok ? 0 : 1;
}
