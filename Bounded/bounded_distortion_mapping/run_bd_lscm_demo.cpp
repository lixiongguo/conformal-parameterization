#include "BoundedDistortionMapping.hpp"

#include <Eigen/Core>

#include <iostream>
#include <vector>

int main() {
    using namespace bounded_distortion;

    Eigen::MatrixXd vertices(5, 2);
    vertices << 0.0, 0.0,
                1.0, 0.0,
                1.0, 1.0,
                0.0, 1.0,
                0.5, 0.5;

    Eigen::MatrixXi faces(4, 3);
    faces << 0, 1, 4,
             1, 2, 4,
             2, 3, 4,
             3, 0, 4;

    Eigen::MatrixXd initial_uv = vertices;
    initial_uv(4, 0) = 0.82;
    initial_uv(4, 1) = 0.18;

    std::vector<Anchor> anchors = {
        {0, Eigen::Vector2d(0.0, 0.0)},
        {1, Eigen::Vector2d(1.0, 0.0)},
        {2, Eigen::Vector2d(1.0, 1.0)},
        {3, Eigen::Vector2d(0.0, 1.0)},
    };

    Options options;
    options.distortion_bound = 1.5;
    options.lscm_weight = 1.0;
    options.outer_iterations = 8;
    options.inner_iterations = 400;
    options.initial_step = 1e-2;
    options.distortion_penalty = 5000.0;
    options.positivity_penalty = 5000.0;

    const SolveResult result = solveBoundedDistortionLscm(vertices, faces, initial_uv, anchors, options);

    double lscm_energy = 0.0;
    for (const FaceStats& face : result.faces) {
        lscm_energy += std::norm(face.coeffs.beta);
    }

    std::cout << "BD-LSCM iterations : " << result.iterations << '\n';
    std::cout << "BD-LSCM energy     : " << lscm_energy << '\n';
    std::cout << "max distortion     : " << result.max_distortion << '\n';
    std::cout << "min jacobian       : " << result.min_jacobian << '\n';
    std::cout << "uv:\n" << result.uv << '\n';

    for (std::size_t i = 0; i < result.faces.size(); ++i) {
        std::cout << "face " << i
                  << " |beta|=" << std::abs(result.faces[i].coeffs.beta)
                  << " K=" << result.faces[i].distortion
                  << " J=" << result.faces[i].jacobian
                  << " cone=" << result.faces[i].cone_violation
                  << '\n';
    }

    return result.min_jacobian > 0.0 ? 0 : 1;
}
