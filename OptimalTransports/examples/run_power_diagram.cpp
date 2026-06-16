// Smoke test for the geogram-backed power diagram in the geometry/ layer.

#include "geometry/power_diagram.hpp"

#include <cstdio>

int main() {
    const double xy[] = {0.2, 0.2, 0.8, 0.2, 0.5, 0.8};
    const double w[] = {0.0, 0.0, 0.0};
    const otgeo::Box box{0.0, 0.0, 1.0, 1.0};

    const auto cells = otgeo::compute_power_diagram(xy, w, 3, box);
    if (cells.empty()) {
        std::fprintf(stderr, "run_power_diagram: no cells\n");
        return 1;
    }
    std::printf("cells = %zu\n", cells.size());
    for (const auto& c : cells) {
        std::printf("  site %d : %zu vertices, area = %f\n", c.site,
                    c.poly.size(), c.poly.area());
    }
    return 0;
}
