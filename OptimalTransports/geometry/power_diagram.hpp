#pragma once

// Weighted Delaunay / power diagram via geogram (2D lift). This replaces the
// CGAL Regular_triangulation_2 + dual construction of the old MongeAmpere code.
//
// For each site we return its Laguerre cell as a convex polygon clipped to the
// domain box. Every polygon edge carries `edge_site`: the id of the neighbour
// site sharing that edge (i.e. the bisector / radical axis between the two
// sites), or otgeo::kBoundaryEdge if the edge lies on the domain boundary. This
// adjacency is what lets the OT solver assemble the Hessian.

#include "density_grid.hpp"  // Box
#include "polygon.hpp"

#include <vector>

namespace otgeo {

struct LaguerreCell {
    int site = -1;     // index of the site owning this cell
    Polygon poly;      // clipped cell; edges tagged with neighbour site ids
};

// Compute the power diagram of `n` weighted sites.
//   xy      : interleaved site coordinates (x0,y0,x1,y1,...), length 2*n
//   weights : per-site weights, length n
//   box     : clipping domain
// Cells with fewer than 3 vertices (empty after clipping) are dropped.
std::vector<LaguerreCell> compute_power_diagram(const double* xy,
                                                const double* weights, int n,
                                                const Box& box);

}  // namespace otgeo
