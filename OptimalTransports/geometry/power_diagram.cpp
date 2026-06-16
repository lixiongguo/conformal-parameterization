#include "power_diagram.hpp"

#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/numeric.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/delaunay_2d.h>
#include <geogram/numerics/predicates.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace otgeo {
namespace {

using namespace GEO;

// Dual of the weighted Delaunay triangulation: builds each Laguerre cell by
// walking the triangles around a site, tagging every cell edge with the
// neighbouring site index. Mirrors the traversal of the previous CGAL code but
// keeps the adjacency information needed for the OT Hessian.
class PowerDiagramBuilder {
public:
    explicit PowerDiagramBuilder(Delaunay& delaunay) : delaunay_(delaunay) {}

    Vec2 dual_vertex(index_t t) const {
        const index_t v1 = delaunay_.cell_to_v()[3 * t];
        const index_t v2 = delaunay_.cell_to_v()[3 * t + 1];
        const index_t v3 = delaunay_.cell_to_v()[3 * t + 2];
        const double* p1 = delaunay_.vertex_ptr(v1);
        const double* p2 = delaunay_.vertex_ptr(v2);
        const double* p3 = delaunay_.vertex_ptr(v3);
        const vec3 q1(p1[0], p1[1], p1[2]);
        const vec3 q2(p2[0], p2[1], p2[2]);
        const vec3 q3(p3[0], p3[1], p3[2]);
        const vec3 c = Geom::triangle_circumcenter(q1, q2, q3);
        return Vec2{c.x, c.y};
    }

    Vec2 infinite_vertex(index_t t, index_t e) const {
        const index_t lv1 = (e + 1) % 3;
        const index_t lv2 = (e + 2) % 3;
        const index_t v1 = delaunay_.cell_to_v()[3 * t + lv1];
        const index_t v2 = delaunay_.cell_to_v()[3 * t + lv2];
        const double* ptr1 = delaunay_.vertex_ptr(v1);
        const double* ptr2 = delaunay_.vertex_ptr(v2);
        GEO::vec2 p1(ptr1[0], ptr1[1]);
        GEO::vec2 p2(ptr2[0], ptr2[1]);
        GEO::vec2 nrm = normalize(p2 - p1);
        nrm = GEO::vec2(nrm.y, -nrm.x);
        const GEO::vec2 r = 0.5 * (p1 + p2) + 100000.0 * nrm;
        return Vec2{r.x, r.y};
    }

    index_t find_vertex(index_t t, index_t v) const {
        for (index_t lv = 0; lv < 3; ++lv) {
            if (delaunay_.cell_to_v()[3 * t + lv] == v) return lv;
        }
        geo_assert_not_reached;
        return 0;
    }

    int site_of(index_t t, index_t lv) const {
        return static_cast<int>(delaunay_.cell_to_v()[3 * t + lv]);
    }

    // Build the (unclipped) Laguerre cell of site `v` reached via local vertex
    // lv of triangle t0. Polygon edges are tagged with neighbour site ids.
    void get_cell(index_t t0, index_t lv, Polygon& cell) const {
        cell.clear();
        const index_t v = delaunay_.cell_to_v()[3 * t0 + lv];
        bool on_border = false;
        index_t t = t0;

        do {
            const index_t e = (lv + 1) % 3;
            const index_t neigh_t = delaunay_.cell_to_cell()[3 * t + e];
            if (neigh_t == NO_INDEX) {
                on_border = true;
                break;
            }
            // edge leaving dual_vertex(t) separates v from local vertex (lv+2)%3
            const int neighbour = site_of(t, (lv + 2) % 3);
            const Vec2 c = dual_vertex(t);
            cell.push(c, neighbour);
            t = neigh_t;
            lv = find_vertex(t, v);
        } while (t != t0);

        if (on_border) {
            cell.clear();
            // first unbounded edge: neighbour is local vertex (lv+2)%3 of t
            cell.push(infinite_vertex(t, (lv + 1) % 3), site_of(t, (lv + 2) % 3));
            for (;;) {
                const int neighbour = site_of(t, (lv + 1) % 3);
                const Vec2 c = dual_vertex(t);
                cell.push(c, neighbour);
                const index_t e = (lv + 2) % 3;
                const index_t neigh_t = delaunay_.cell_to_cell()[3 * t + e];
                if (neigh_t == NO_INDEX) {
                    // closing edge runs off to infinity; treat as boundary
                    cell.push(infinite_vertex(t, e), kBoundaryEdge);
                    break;
                }
                t = neigh_t;
                lv = find_vertex(t, v);
            }
        }
    }

    void build_all(std::vector<LaguerreCell>& out, const Box& box) const {
        const index_t nV = delaunay_.nb_vertices();
        std::vector<bool> visited(static_cast<size_t>(nV), false);
        out.clear();
        out.reserve(static_cast<size_t>(nV));

        for (index_t t = 0; t < delaunay_.nb_cells(); ++t) {
            for (index_t lv = 0; lv < 3; ++lv) {
                const index_t v = delaunay_.cell_to_v()[3 * t + lv];
                if (visited[static_cast<size_t>(v)]) continue;
                visited[static_cast<size_t>(v)] = true;

                Polygon cell;
                get_cell(t, lv, cell);

                Polygon clipped;
                clip_rectangle(cell, box.xmin, box.ymin, box.xmax, box.ymax,
                               clipped);
                if (clipped.size() < 3) continue;

                LaguerreCell lc;
                lc.site = static_cast<int>(v);
                lc.poly = std::move(clipped);
                out.push_back(std::move(lc));
            }
        }
    }

private:
    Delaunay& delaunay_;
};

void ensure_geogram_initialized() {
    static bool done = false;
    if (!done) {
        GEO::initialize();
        done = true;
    }
}

}  // namespace

std::vector<LaguerreCell> compute_power_diagram(const double* xy,
                                                const double* weights, int n,
                                                const Box& box) {
    std::vector<LaguerreCell> empty;
    if (!xy || !weights || n < 1) return empty;

    ensure_geogram_initialized();

    double wMax = weights[0];
    for (int i = 1; i < n; ++i) wMax = std::max(wMax, weights[i]);
    if (wMax < 0.0) wMax = 0.0;

    std::vector<double> verts(static_cast<size_t>(n) * 3);
    for (int i = 0; i < n; ++i) {
        const double w = weights[i];
        const double lift = wMax - w;
        verts[static_cast<size_t>(i) * 3] = xy[i * 2];
        verts[static_cast<size_t>(i) * 3 + 1] = xy[i * 2 + 1];
        verts[static_cast<size_t>(i) * 3 + 2] = lift > 0.0 ? std::sqrt(lift) : 0.0;
    }

    Delaunay_var delaunay = Delaunay::create(3, "BPOW2d");
    if (!delaunay) return empty;

    delaunay->set_vertices(static_cast<index_t>(n), verts.data());

    Delaunay2d* del2d = dynamic_cast<Delaunay2d*>(delaunay.get());
    if (!del2d) return empty;
    if (del2d->has_empty_cells()) return empty;

    PowerDiagramBuilder builder(*delaunay);
    std::vector<LaguerreCell> cells;
    builder.build_all(cells, box);
    return cells;
}

}  // namespace otgeo
