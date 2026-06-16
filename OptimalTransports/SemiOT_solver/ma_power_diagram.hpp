#ifndef MA_POWER_DIAGRAM_HPP
#define MA_POWER_DIAGRAM_HPP

#include "ma_geometry.hpp"
#include <string>
#include <vector>

namespace MA {

/**
 * Weighted Delaunay / Power diagram cells via Geogram (2D lift).
 * \param xy interleaved site coordinates (n x 2), typically in [0,1]^2
 * \param weights per-site weights (length n)
 * \param n number of sites
 * \param clip_unit_square if true, clip each cell to [0,1]^2
 */
std::vector<Polygon2> compute_power_diagram(
    const double* xy,
    const double* weights,
    int n,
    bool clip_unit_square = true
);

/** Serialize cells to JSON: [[x0,y0,...], ...] */
std::string power_diagram_to_json(const std::vector<Polygon2>& cells);

} // namespace MA

#endif
