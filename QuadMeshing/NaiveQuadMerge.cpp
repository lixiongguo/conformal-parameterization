#include "NaiveQuadMerge.h"

#include <cstdint>
#include <unordered_map>
#include <utility>

namespace NaiveQuadMerge {
namespace {

struct EdgeUse {
    int face;
    int slot;
};

uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

double edgeLen2(const double* positions, int a, int b) {
    const double* pa = positions + a * 3;
    const double* pb = positions + b * 3;
    const double x = pa[0] - pb[0];
    const double y = pa[1] - pb[1];
    const double z = pa[2] - pb[2];
    return x * x + y * y + z * z;
}

int oppositeVertex(const int* tri, int a, int b) {
    for (int k = 0; k < 3; ++k) {
        const int v = tri[k];
        if (v != a && v != b) return v;
    }
    return -1;
}

} // namespace

Result merge(const double* positions, int nVerts, const int* faces, int nFaces) {
    Result out;
    if (!positions || !faces || nVerts <= 0 || nFaces <= 0) return out;

    StageClock clock(out.stages);
    std::unordered_map<uint64_t, std::vector<EdgeUse>> adj;
    adj.reserve(static_cast<size_t>(nFaces) * 2);

    auto valid = [&](int v) { return v >= 0 && v < nVerts; };

    clock.lap("建立邻接", [&] {
        for (int fi = 0; fi < nFaces; ++fi) {
            const int* tri = faces + fi * 3;
            if (!valid(tri[0]) || !valid(tri[1]) || !valid(tri[2])) continue;
            if (tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0]) continue;
            for (int s = 0; s < 3; ++s) {
                adj[edgeKey(tri[s], tri[(s + 1) % 3])].push_back(EdgeUse{fi, s});
            }
        }
    });

    std::vector<char> used(static_cast<size_t>(nFaces), 0);
    out.quads.reserve(static_cast<size_t>(nFaces) * 2);
    out.leftovers.reserve(static_cast<size_t>(nFaces));

    clock.lap("贪心配对", [&] {
    for (int fi = 0; fi < nFaces; ++fi) {
        if (used[static_cast<size_t>(fi)]) continue;
        const int* tri = faces + fi * 3;
        if (!valid(tri[0]) || !valid(tri[1]) || !valid(tri[2]) ||
            tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0]) {
            used[static_cast<size_t>(fi)] = 1;
            continue;
        }

        int bestSlot = -1;
        int bestNeighbor = -1;
        double bestLen = -1.0;

        for (int s = 0; s < 3; ++s) {
            const int a = tri[s];
            const int b = tri[(s + 1) % 3];
            const double len = edgeLen2(positions, a, b);
            if (!(len > bestLen)) continue;
            const auto it = adj.find(edgeKey(a, b));
            if (it == adj.end()) continue;
            for (const EdgeUse& use : it->second) {
                if (use.face == fi || used[static_cast<size_t>(use.face)]) continue;
                const int* other = faces + use.face * 3;
                const int d = oppositeVertex(other, a, b);
                const int c = tri[(s + 2) % 3];
                if (d < 0 || d == c || !valid(d)) continue;
                bestLen = len;
                bestSlot = s;
                bestNeighbor = use.face;
                break;
            }
        }

        if (bestNeighbor < 0) {
            used[static_cast<size_t>(fi)] = 1;
            out.leftovers.push_back(tri[0]);
            out.leftovers.push_back(tri[1]);
            out.leftovers.push_back(tri[2]);
            continue;
        }

        const int a = tri[bestSlot];
        const int b = tri[(bestSlot + 1) % 3];
        const int c = tri[(bestSlot + 2) % 3];
        const int d = oppositeVertex(faces + bestNeighbor * 3, a, b);
        used[static_cast<size_t>(fi)] = 1;
        used[static_cast<size_t>(bestNeighbor)] = 1;
        out.quads.push_back(a);
        out.quads.push_back(d);
        out.quads.push_back(b);
        out.quads.push_back(c);
    }
    });

    return out;
}

} // namespace NaiveQuadMerge
