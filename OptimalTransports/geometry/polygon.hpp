#pragma once

// Lightweight 2D convex-polygon primitives used by the geogram-backed geometry
// layer. Replaces the CGAL::Polygon_2 usage of the old MongeAmpere code.
//
// A polygon is an ordered list of vertices. Each vertex additionally stores
// `edge_site`: the id of the neighbouring power-diagram site across the edge
// that STARTS at this vertex (vertex i -> vertex i+1). A value < 0 means the
// edge lies on the clipping domain boundary (no neighbour).

#include <cmath>
#include <cstddef>
#include <vector>

namespace otgeo {

constexpr int kBoundaryEdge = -1;  // edge lies on the domain boundary

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

inline Vec2 operator-(const Vec2& a, const Vec2& b) { return {a.x - b.x, a.y - b.y}; }
inline Vec2 operator+(const Vec2& a, const Vec2& b) { return {a.x + b.x, a.y + b.y}; }
inline Vec2 operator*(double s, const Vec2& a) { return {s * a.x, s * a.y}; }
inline double dot(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }
inline double cross(const Vec2& a, const Vec2& b) { return a.x * b.y - a.y * b.x; }
inline double norm(const Vec2& a) { return std::sqrt(dot(a, a)); }

struct PolyVertex {
    Vec2 p;
    int edge_site = kBoundaryEdge;  // neighbour across edge (this -> next)
};

struct Polygon {
    std::vector<PolyVertex> v;

    void clear() { v.clear(); }
    std::size_t size() const { return v.size(); }
    bool empty() const { return v.empty(); }

    void push(double x, double y, int edge_site = kBoundaryEdge) {
        v.push_back(PolyVertex{Vec2{x, y}, edge_site});
    }
    void push(const Vec2& p, int edge_site = kBoundaryEdge) {
        v.push_back(PolyVertex{p, edge_site});
    }

    const Vec2& point(std::size_t i) const { return v[i].p; }

    // Signed area (positive when vertices are counter-clockwise).
    double signed_area() const {
        const std::size_t n = v.size();
        if (n < 3) return 0.0;
        double a = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const Vec2& p0 = v[i].p;
            const Vec2& p1 = v[(i + 1) % n].p;
            a += cross(p0, p1);
        }
        return 0.5 * a;
    }

    double area() const { return std::fabs(signed_area()); }

    Vec2 centroid() const {
        const std::size_t n = v.size();
        if (n == 0) return Vec2{0.0, 0.0};
        if (n < 3) {
            Vec2 c{0.0, 0.0};
            for (const auto& vv : v) c = c + vv.p;
            return (1.0 / static_cast<double>(n)) * c;
        }
        double a = 0.0;
        Vec2 c{0.0, 0.0};
        for (std::size_t i = 0; i < n; ++i) {
            const Vec2& p0 = v[i].p;
            const Vec2& p1 = v[(i + 1) % n].p;
            const double w = cross(p0, p1);
            a += w;
            c.x += (p0.x + p1.x) * w;
            c.y += (p0.y + p1.y) * w;
        }
        if (std::fabs(a) < 1e-300) return v[0].p;
        const double inv = 1.0 / (3.0 * a);
        return Vec2{c.x * inv, c.y * inv};
    }
};

// Clip a convex polygon against the half-plane to the LEFT of the directed line
// q1 -> q2 (i.e. points p with cross(q2 - q1, p - q1) >= 0 are kept).
// Edges newly created on the clip line are tagged with `clip_edge_site`.
inline void clip_half_plane(const Polygon& in, const Vec2& q1, const Vec2& q2,
                            int clip_edge_site, Polygon& out) {
    out.clear();
    const std::size_t n = in.size();
    if (n == 0) return;

    const Vec2 dir = q2 - q1;
    auto side = [&](const Vec2& p) { return cross(dir, p - q1); };
    const double eps = 0.0;

    auto intersect = [&](const Vec2& a, const Vec2& b) -> Vec2 {
        const double sa = side(a);
        const double sb = side(b);
        const double denom = sa - sb;
        if (std::fabs(denom) < 1e-300) return a;
        const double t = sa / denom;
        return Vec2{a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
    };

    for (std::size_t i = 0; i < n; ++i) {
        const PolyVertex& S = in.v[i];
        const PolyVertex& E = in.v[(i + 1) % n];
        const bool in_S = side(S.p) >= eps;
        const bool in_E = side(E.p) >= eps;

        if (in_E) {
            if (!in_S) {
                // entering: the partial edge lies on the original edge S->E,
                // so it keeps S's neighbour tag.
                out.v.push_back(PolyVertex{intersect(S.p, E.p), S.edge_site});
            }
            // keep E; the edge leaving E is the original edge E->next.
            out.v.push_back(E);
        } else if (in_S) {
            // exiting: the edge leaving this intersection runs along the clip
            // line until we re-enter, so it is a boundary edge.
            out.v.push_back(PolyVertex{intersect(S.p, E.p), clip_edge_site});
        }
    }
}

// Clip against an axis-aligned rectangle [xmin,xmax] x [ymin,ymax].
inline void clip_rectangle(const Polygon& in, double xmin, double ymin,
                           double xmax, double ymax, Polygon& out) {
    Polygon a = in;
    Polygon b;
    // CCW rectangle so that "left of directed edge" == "inside".
    clip_half_plane(a, Vec2{xmin, ymin}, Vec2{xmax, ymin}, kBoundaryEdge, b);
    clip_half_plane(b, Vec2{xmax, ymin}, Vec2{xmax, ymax}, kBoundaryEdge, a);
    clip_half_plane(a, Vec2{xmax, ymax}, Vec2{xmin, ymax}, kBoundaryEdge, b);
    clip_half_plane(b, Vec2{xmin, ymax}, Vec2{xmin, ymin}, kBoundaryEdge, a);
    out = a;
}

}  // namespace otgeo
