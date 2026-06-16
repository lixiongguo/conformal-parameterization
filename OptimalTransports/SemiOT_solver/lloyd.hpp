// MongeAmpere++  (geogram port)
// Original Copyright (C) 2014 Quentin Merigot, CNRS (GPLv3+)
//
// CGAL-free Lloyd / centroidal helpers built on the geogram power diagram and a
// grid source density.

#ifndef MA_LLOYD_HPP
#define MA_LLOYD_HPP

#include <Eigen/Dense>

#include <vector>

#include "../geometry/density_grid.hpp"
#include "../geometry/power_diagram.hpp"
#include "../geometry/quadrature.hpp"

namespace MA {

// Per-cell density mass and (unnormalised) first moment of the coordinates.
template <class Matrix, class Vector>
void first_moment(const otgeo::DensityGrid& source, const Matrix& X,
                  const Vector& weights, Vector& masses, Matrix& centroids) {
    const int N = static_cast<int>(X.rows());
    masses = Vector::Zero(N);
    centroids = Matrix::Zero(N, 2);

    std::vector<double> xy(static_cast<size_t>(N) * 2);
    std::vector<double> w(static_cast<size_t>(N));
    for (int i = 0; i < N; ++i) {
        xy[static_cast<size_t>(i) * 2] = X(i, 0);
        xy[static_cast<size_t>(i) * 2 + 1] = X(i, 1);
        w[static_cast<size_t>(i)] = weights(i);
    }

    const std::vector<otgeo::LaguerreCell> cells =
        otgeo::compute_power_diagram(xy.data(), w.data(), N, source.box());

    for (const otgeo::LaguerreCell& cell : cells) {
        const int idv = cell.site;
        const otgeo::Polygon& p = cell.poly;
        const double area = otgeo::integrate_polygon_centroid(
            p, [&](const otgeo::Vec2& q) { return source(q); });
        const double mx = otgeo::integrate_polygon_centroid(
            p, [&](const otgeo::Vec2& q) { return source(q) * q.x; });
        const double my = otgeo::integrate_polygon_centroid(
            p, [&](const otgeo::Vec2& q) { return source(q) * q.y; });
        masses[idv] += area;
        centroids(idv, 0) += mx;
        centroids(idv, 1) += my;
    }
}

// Centroids of the Laguerre cells (Lloyd relaxation target).
template <class Matrix, class Vector>
void lloyd(const otgeo::DensityGrid& source, const Matrix& X,
           const Vector& weights, Vector& masses, Matrix& centroids) {
    first_moment(source, X, weights, masses, centroids);
    const int N = static_cast<int>(X.rows());
    for (int i = 0; i < N; ++i) {
        if (masses[i] > 0.0) {
            centroids(i, 0) /= masses[i];
            centroids(i, 1) /= masses[i];
        }
    }
}

}  // namespace MA

#endif
