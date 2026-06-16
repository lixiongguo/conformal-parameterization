#include "ma_power_diagram.hpp"

#include "../geometry/power_diagram.hpp"

#include <sstream>
#include <string>

// Compatibility shim: the power diagram now lives in the shared geometry layer
// (otgeo). This keeps the historical MA:: API working for existing callers
// (e.g. tests/test_power_diagram.cpp) while delegating to the geogram backend.

namespace MA {

std::vector<Polygon2> compute_power_diagram(const double* xy,
                                            const double* weights, int n,
                                            bool clip_unit_square) {
    otgeo::Box box{0.0, 0.0, 1.0, 1.0};
    if (!clip_unit_square) {
        // Effectively no clipping: use a large box.
        box = otgeo::Box{-1e6, -1e6, 1e6, 1e6};
    }

    const std::vector<otgeo::LaguerreCell> cells =
        otgeo::compute_power_diagram(xy, weights, n, box);

    std::vector<Polygon2> out;
    out.reserve(cells.size());
    for (const otgeo::LaguerreCell& c : cells) {
        Polygon2 poly;
        poly.xy.reserve(c.poly.size() * 2);
        for (const otgeo::PolyVertex& v : c.poly.v) poly.push(v.p.x, v.p.y);
        out.push_back(std::move(poly));
    }
    return out;
}

std::string power_diagram_to_json(const std::vector<Polygon2>& cells) {
    std::ostringstream os;
    os << '[';
    for (size_t ci = 0; ci < cells.size(); ++ci) {
        if (ci) os << ',';
        os << '[';
        const auto& xy = cells[ci].xy;
        for (size_t k = 0; k < xy.size(); ++k) {
            if (k) os << ',';
            os << xy[k];
        }
        os << ']';
    }
    os << ']';
    return os.str();
}

}  // namespace MA
