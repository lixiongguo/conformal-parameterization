#pragma once

// A non-negative density sampled on a regular grid over an axis-aligned box.
// This is the source-measure representation used by the (reformulated)
// semi-discrete OT solver. Replaces the CGAL piecewise-linear-on-triangulation
// density of the old MongeAmpere code.

#include "polygon.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

namespace otgeo {

// Axis-aligned box [xmin,xmax] x [ymin,ymax].
struct Box {
    double xmin = 0.0;
    double ymin = 0.0;
    double xmax = 1.0;
    double ymax = 1.0;

    double width() const { return xmax - xmin; }
    double height() const { return ymax - ymin; }
};

class DensityGrid {
public:
    DensityGrid() = default;
    DensityGrid(Eigen::MatrixXd values, Box box)
        : values_(std::move(values)), box_(box) {}

    const Box& box() const { return box_; }
    int rows() const { return static_cast<int>(values_.rows()); }
    int cols() const { return static_cast<int>(values_.cols()); }

    // Bilinear sample of the (cell-centred) grid at world point (x,y).
    // Points outside the box are clamped to the border value.
    double operator()(const Vec2& q) const {
        const int nr = rows();
        const int nc = cols();
        if (nr == 0 || nc == 0) return 0.0;

        // world -> continuous cell coordinates (cell centres at integer+0.5)
        double fx = (q.x - box_.xmin) / box_.width() * nc - 0.5;
        double fy = (q.y - box_.ymin) / box_.height() * nr - 0.5;
        fx = std::clamp(fx, 0.0, static_cast<double>(nc - 1));
        fy = std::clamp(fy, 0.0, static_cast<double>(nr - 1));

        const int c0 = static_cast<int>(std::floor(fx));
        const int r0 = static_cast<int>(std::floor(fy));
        const int c1 = std::min(c0 + 1, nc - 1);
        const int r1 = std::min(r0 + 1, nr - 1);
        const double tx = fx - c0;
        const double ty = fy - r0;

        const double v00 = values_(r0, c0);
        const double v01 = values_(r0, c1);
        const double v10 = values_(r1, c0);
        const double v11 = values_(r1, c1);
        const double a = v00 * (1.0 - tx) + v01 * tx;
        const double b = v10 * (1.0 - tx) + v11 * tx;
        return a * (1.0 - ty) + b * ty;
    }

    // Total mass (integral over the box).
    double total_mass() const {
        const double cell_area = (box_.width() * box_.height()) /
                                 static_cast<double>(std::max(1, rows() * cols()));
        return values_.sum() * cell_area;
    }

private:
    Eigen::MatrixXd values_;
    Box box_;
};

}  // namespace otgeo
