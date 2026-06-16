#ifndef MA_GEOMETRY_HPP
#define MA_GEOMETRY_HPP

#include <vector>
#include <cmath>

namespace MA {

/** Simple 2D polygon as interleaved coordinates (x0,y0,x1,y1,...). */
struct Polygon2 {
    std::vector<double> xy;

    void clear() { xy.clear(); }
    size_t numVertices() const { return xy.size() / 2; }
    void push(double x, double y) {
        xy.push_back(x);
        xy.push_back(y);
    }
};

} // namespace MA

#endif
