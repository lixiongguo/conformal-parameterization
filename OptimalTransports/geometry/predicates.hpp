#pragma once

// Small geometric predicates / helpers shared by the geometry layer.
//
// Robust combinatorial predicates (in-circle, orientation for the Delaunay /
// power triangulation) are handled by geogram inside power_diagram.cpp. The
// helpers here operate on plain doubles and are enough for the polygon-level
// work (line intersections, segment clipping, point location).

#include "polygon.hpp"

#include <algorithm>
#include <cmath>

namespace otgeo {

// Orientation of (b-a, c-a): >0 CCW, <0 CW, 0 collinear.
inline double orient2d(const Vec2& a, const Vec2& b, const Vec2& c) {
    return cross(b - a, c - a);
}

// Intersection of the two lines (p1,p2) and (q1,q2). Assumes they are not
// parallel; returns p1 on (near-)degenerate input.
inline Vec2 line_line_intersection(const Vec2& p1, const Vec2& p2,
                                    const Vec2& q1, const Vec2& q2) {
    const Vec2 d1 = p2 - p1;
    const Vec2 d2 = q2 - q1;
    const double denom = cross(d1, d2);
    if (std::fabs(denom) < 1e-300) return p1;
    const double t = cross(q1 - p1, d2) / denom;
    return Vec2{p1.x + t * d1.x, p1.y + t * d1.y};
}

// Clip segment [a,b] to the rectangle using Liang-Barsky. Returns false if the
// segment lies fully outside; otherwise writes the clipped endpoints.
inline bool clip_segment_rectangle(Vec2 a, Vec2 b, double xmin, double ymin,
                                   double xmax, double ymax, Vec2& ca, Vec2& cb) {
    double t0 = 0.0, t1 = 1.0;
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double p[4] = {-dx, dx, -dy, dy};
    const double q[4] = {a.x - xmin, xmax - a.x, a.y - ymin, ymax - a.y};
    for (int i = 0; i < 4; ++i) {
        if (std::fabs(p[i]) < 1e-300) {
            if (q[i] < 0.0) return false;  // parallel and outside
            continue;
        }
        const double r = q[i] / p[i];
        if (p[i] < 0.0) {
            if (r > t1) return false;
            if (r > t0) t0 = r;
        } else {
            if (r < t0) return false;
            if (r < t1) t1 = r;
        }
    }
    ca = Vec2{a.x + t0 * dx, a.y + t0 * dy};
    cb = Vec2{a.x + t1 * dx, a.y + t1 * dy};
    return true;
}

}  // namespace otgeo
