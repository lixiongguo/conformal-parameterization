#pragma once

// Polygon quadrature on plain doubles. Ported from the CGAL-based
// MA/quadrature.hpp but operating on otgeo::Polygon and scalar functions
// f : Vec2 -> double.

#include "polygon.hpp"

#include <cmath>
#include <cstddef>

namespace otgeo {

inline double triangle_area(const Vec2& a, const Vec2& b, const Vec2& c) {
    return 0.5 * std::fabs(cross(b - a, c - a));
}

// Centroid rule (exact for affine integrands).
template <class F>
double integrate_centroid_tri(const Vec2& a, const Vec2& b, const Vec2& c,
                              const F& f) {
    const Vec2 g{(a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0};
    return triangle_area(a, b, c) * f(g);
}

// Order-3 Albrecht-Collatz rule with 6 points.
template <class F>
double integrate_ac_tri(const Vec2& a, const Vec2& b, const Vec2& c,
                        const F& f) {
    const double _1_2 = 0.5, _1_6 = 1.0 / 6.0, _2_3 = 2.0 / 3.0;
    const double _1_30 = 1.0 / 30.0, _9_30 = 9.0 / 30.0;
    const Vec2 u = b - a;
    const Vec2 v = c - a;
    auto at = [&](double s, double t) { return Vec2{a.x + s * u.x + t * v.x, a.y + s * u.y + t * v.y}; };
    const double r =
        _1_30 * f(at(_1_2, _1_2)) + _1_30 * f(at(_1_2, 0.0)) +
        _1_30 * f(at(0.0, _1_2)) + _9_30 * f(at(_1_6, _2_3)) +
        _9_30 * f(at(_2_3, _1_6)) + _9_30 * f(at(_1_6, _1_6));
    return triangle_area(a, b, c) * r;
}

// Integral of f over the polygon using a triangle fan + centroid rule.
template <class F>
double integrate_polygon_centroid(const Polygon& p, const F& f) {
    const std::size_t n = p.size();
    if (n < 3) return 0.0;
    double r = 0.0;
    for (std::size_t i = 1; i + 1 < n; ++i)
        r += integrate_centroid_tri(p.point(0), p.point(i), p.point(i + 1), f);
    return r;
}

// Integral of f over the polygon using a triangle fan + Albrecht-Collatz rule.
template <class F>
double integrate_polygon_ac(const Polygon& p, const F& f) {
    const std::size_t n = p.size();
    if (n < 3) return 0.0;
    double r = 0.0;
    for (std::size_t i = 1; i + 1 < n; ++i)
        r += integrate_ac_tri(p.point(0), p.point(i), p.point(i + 1), f);
    return r;
}

// Line integral of f along a segment (length * f(midpoint)).
template <class F>
double integrate_segment(const Vec2& a, const Vec2& b, const F& f) {
    const Vec2 mid{0.5 * (a.x + b.x), 0.5 * (a.y + b.y)};
    return norm(b - a) * f(mid);
}

}  // namespace otgeo
