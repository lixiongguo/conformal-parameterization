#ifndef GAUSSIAN_CURVATURE_H
#define GAUSSIAN_CURVATURE_H

#include <Eigen/Core>

class Mesh;

namespace geometry {

/**
 * Discrete Gaussian curvature via angle deficit.
 *
 * For an interior vertex v:
 *   K(v) = 2π - Σ_f θ_f(v)
 *
 * For a boundary vertex v:
 *   K(v) = π - Σ_f θ_f(v)
 *
 * This returns the integrated curvature (radians). Use gaussianCurvaturePerArea()
 * if you want a density (radians / area).
 */
Eigen::VectorXd gaussianCurvatureAngleDeficit(const Mesh& mesh);

/**
 * Gaussian curvature density (angle deficit divided by barycentric area).
 *
 * Barycentric area:
 *   A(v) = (1/3) Σ_{f incident to v} area(f)
 *
 * Returns K(v) / max(A(v), eps).
 */
Eigen::VectorXd gaussianCurvaturePerArea(const Mesh& mesh, double eps = 1e-16);

} // namespace geometry

#endif

