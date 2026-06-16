#include "MA/ma_power_diagram.hpp"
#include <cstdio>

int main()
{
    const double xy[] = {0.2, 0.2, 0.8, 0.2, 0.5, 0.8};
    const double w[] = {0.0, 0.0, 0.0};
    const auto cells = MA::compute_power_diagram(xy, w, 3, true);
    if (cells.empty()) {
        std::fprintf(stderr, "test_power_diagram: no cells\n");
        return 1;
    }
    const std::string json = MA::power_diagram_to_json(cells);
    std::printf("%s\n", json.c_str());
    return 0;
}
