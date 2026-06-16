// MongeAmpere++  (geogram port)
// Original Copyright (C) 2014 Quentin Merigot, CNRS (GPLv3+)
//
// CGAL-free reformulation of the Kantorovich functional for semi-discrete
// optimal transport. The source measure is a density sampled on a regular grid
// (otgeo::DensityGrid); the target is a set of Dirac sites X with prescribed
// masses. The Laguerre (power) cells are computed by the geogram-backed
// geometry layer, and the functional value, gradient and Hessian are obtained
// by integrating the density over each cell and along the shared cell edges.

#ifndef MA_KANTOROVICH_HPP
#define MA_KANTOROVICH_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <vector>

#include "../geometry/density_grid.hpp"
#include "../geometry/power_diagram.hpp"
#include "../geometry/quadrature.hpp"

namespace MA {

// Computes the Kantorovich functional value, its gradient `g` (the Laguerre
// cell masses) and Hessian `h` for the given weights.
//   source  : continuous source density (grid over a box)
//   X       : N x 2 site positions
//   weights : length-N weight vector (the variable of the functional)
// Returns the functional value sum_i (mass_i * w_i) - int rho ||x - x_i||^2.
template <class Matrix, class Vector, class SparseMatrix>
double kantorovich(const otgeo::DensityGrid& source, const Matrix& X,
                   const Vector& weights, Vector& g, SparseMatrix& h) {
    typedef Eigen::Triplet<double> Triplet;

    const int N = static_cast<int>(X.rows());
    g = Vector::Zero(N);

    // Pack sites for the power-diagram backend.
    std::vector<double> xy(static_cast<size_t>(N) * 2);
    std::vector<double> w(static_cast<size_t>(N));
    for (int i = 0; i < N; ++i) {
        xy[static_cast<size_t>(i) * 2] = X(i, 0);
        xy[static_cast<size_t>(i) * 2 + 1] = X(i, 1);
        w[static_cast<size_t>(i)] = weights(i);
    }

    const std::vector<otgeo::LaguerreCell> cells =
        otgeo::compute_power_diagram(xy.data(), w.data(), N, source.box());

    std::vector<Triplet> htri;
    double fval = 0.0;

    for (const otgeo::LaguerreCell& cell : cells) {
        const int idv = cell.site;
        const otgeo::Vec2 xv{X(idv, 0), X(idv, 1)};
        const otgeo::Polygon& p = cell.poly;

        // weighted area (mass) of the cell and transport cost contribution
        const double warea = otgeo::integrate_polygon_centroid(
            p, [&](const otgeo::Vec2& q) { return source(q); });
        const double intg = otgeo::integrate_polygon_ac(
            p, [&](const otgeo::Vec2& q) {
                const otgeo::Vec2 d = q - xv;
                return source(q) * otgeo::dot(d, d);
            });

        g[idv] += warea;
        fval += warea * weights(idv) - intg;

        // Hessian from the shared edges (bisectors with neighbouring sites).
        const std::size_t n = p.size();
        for (std::size_t k = 0; k < n; ++k) {
            const int idw = p.v[k].edge_site;
            if (idw < 0 || idw == idv) continue;
            const otgeo::Vec2& a = p.point(k);
            const otgeo::Vec2& b = p.point((k + 1) % n);
            const double r = otgeo::integrate_segment(
                a, b, [&](const otgeo::Vec2& q) { return source(q); });
            const otgeo::Vec2 xw{X(idw, 0), X(idw, 1)};
            const double d = 2.0 * otgeo::norm(xw - xv);
            if (d <= 0.0) continue;
            htri.push_back(Triplet(idv, idw, -r / d));
            htri.push_back(Triplet(idv, idv, +r / d));
        }
    }

    h = SparseMatrix(N, N);
    h.setFromTriplets(htri.begin(), htri.end());
    h.makeCompressed();
    return fval;
}

}  // namespace MA

#endif
