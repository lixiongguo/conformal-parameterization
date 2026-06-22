// ============================================================
//  Ruppert Delaunay Refinement — CDT (artem-ogre/CDT) + WASM
// ------------------------------------------------------------
//  Constrained Delaunay triangulation via CDT library;
//  Ruppert refinement inserts Steiner points with CDT::insertVertices.
//
//  Based on CDT (MPL-2.0): https://github.com/artem-ogre/CDT
// ============================================================
#include "CDT.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <emscripten.h>
#include <vector>

namespace {

using namespace CDT;

Triangulation<double>         g_cdt;
std::vector<Edge>             g_constraints;
std::vector<double>           g_outPts;
std::vector<int>              g_outTris;
double                        g_time = 0.0;

inline const V2d<double>& vtx(VertInd i) { return g_cdt.vertices[i]; }

inline double orient(const V2d<double>& a, const V2d<double>& b, const V2d<double>& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

inline bool insideDiametral(const V2d<double>& a, const V2d<double>& b, const V2d<double>& w) {
    const double mx = (a.x + b.x) * 0.5, my = (a.y + b.y) * 0.5;
    const double dx = w.x - mx, dy = w.y - my;
    const double r2 = ((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y)) * 0.25;
    return (dx * dx + dy * dy) < r2 * (1.0 - 1e-9);
}

double minAngleDeg(const V2d<double>& a, const V2d<double>& b, const V2d<double>& c) {
    const double ab = std::hypot(b.x - a.x, b.y - a.y);
    const double bc = std::hypot(c.x - b.x, c.y - b.y);
    const double ca = std::hypot(a.x - c.x, a.y - c.y);
    auto clamp = [](double v) { return v < -1 ? -1 : (v > 1 ? 1 : v); };
    const double A = std::acos(clamp((ab * ab + ca * ca - bc * bc) / (2 * ab * ca + 1e-30)));
    const double B = std::acos(clamp((ab * ab + bc * bc - ca * ca) / (2 * ab * bc + 1e-30)));
    const double C = std::acos(clamp((bc * bc + ca * ca - ab * ab) / (2 * bc * ca + 1e-30)));
    double mn = A < B ? A : B;
    if (C < mn) mn = C;
    return mn * 180.0 / M_PI;
}

double triArea(const V2d<double>& a, const V2d<double>& b, const V2d<double>& c) {
    return std::fabs(orient(a, b, c)) * 0.5;
}

bool circumcenter(const V2d<double>& a, const V2d<double>& b, const V2d<double>& c, V2d<double>& out) {
    const double D = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if (std::fabs(D) < 1e-18) return false;
    const double a2 = a.x * a.x + a.y * a.y;
    const double b2 = b.x * b.x + b.y * b.y;
    const double c2 = c.x * c.x + c.y * c.y;
    out.x = (a2 * (b.y - c.y) + b2 * (c.y - a.y) + c2 * (a.y - b.y)) / D;
    out.y = (a2 * (c.x - b.x) + b2 * (a.x - c.x) + c2 * (b.x - a.x)) / D;
    return true;
}

VertInd insertPoint(const V2d<double>& p) {
    const VertInd before = static_cast<VertInd>(g_cdt.vertices.size());
    g_cdt.insertVertices(std::vector<V2d<double>>{p});
    return before;
}

VertInd splitConstraint(int s) {
    const Edge e = g_constraints[s];
    const VertInd u = e.v1(), v = e.v2();
    const V2d<double> mid{
        (vtx(u).x + vtx(v).x) * 0.5,
        (vtx(u).y + vtx(v).y) * 0.5
    };
    const VertInd midIdx = insertPoint(mid);
    g_constraints[s] = Edge(u, midIdx);
    g_constraints.push_back(Edge(midIdx, v));
    g_cdt.insertEdges({Edge(u, midIdx), Edge(midIdx, v)});
    return midIdx;
}

void runEncroach(std::vector<VertInd>& work) {
    while (!work.empty()) {
        const VertInd w = work.back();
        work.pop_back();
        bool again = true;
        while (again) {
            again = false;
            for (int s = 0; s < static_cast<int>(g_constraints.size()); ++s) {
                const VertInd u = g_constraints[s].v1();
                const VertInd v = g_constraints[s].v2();
                if (w == u || w == v) continue;
                if (insideDiametral(vtx(u), vtx(v), vtx(w))) {
                    const VertInd mid = splitConstraint(s);
                    work.push_back(mid);
                    again = true;
                }
            }
        }
    }
}

bool triIsUsable(const Triangle& t) {
    if (touchesSuperTriangle(t)) return false;
    const VertInd a = t.vertices[0], b = t.vertices[1], c = t.vertices[2];
    return orient(vtx(a), vtx(b), vtx(c)) > 0.0;
}

void collectOutput() {
    g_outPts.clear();
    g_outTris.clear();
    std::vector<int> used(g_cdt.vertices.size(), -1);
    for (const Triangle& t : g_cdt.triangles) {
        if (!triIsUsable(t)) continue;
        for (int k = 0; k < 3; ++k) {
            const VertInd vi = t.vertices[k];
            if (used[vi] < 0) {
                used[vi] = static_cast<int>(g_outPts.size() / 2);
                g_outPts.push_back(vtx(vi).x);
                g_outPts.push_back(vtx(vi).y);
            }
        }
        g_outTris.push_back(used[t.vertices[0]]);
        g_outTris.push_back(used[t.vertices[1]]);
        g_outTris.push_back(used[t.vertices[2]]);
    }
}

void resetAll() {
    g_cdt = Triangulation<double>();
    g_constraints.clear();
    g_outPts.clear();
    g_outTris.clear();
    g_time = 0.0;
}

} // namespace

extern "C" {

EMSCRIPTEN_KEEPALIVE
void dispose() { resetAll(); }

EMSCRIPTEN_KEEPALIVE
int ruppert_run(double* pts, int nCoord, int* segs, int nSegIdx,
                double minAngle, double maxArea, int maxPoints) {
    const auto t0 = std::chrono::steady_clock::now();
    resetAll();

    const int nPts = nCoord / 2;
    const int nSeg = nSegIdx / 2;
    if (nPts < 3 || nSeg < 3) return 1;

    std::vector<V2d<double>> vertices;
    vertices.reserve(nPts);
    for (int i = 0; i < nPts; ++i)
        vertices.emplace_back(V2d<double>{pts[2 * i], pts[2 * i + 1]});

    std::vector<Edge> edges;
    edges.reserve(nSeg);
    for (int s = 0; s < nSeg; ++s) {
        const VertInd a = static_cast<VertInd>(segs[2 * s]);
        const VertInd b = static_cast<VertInd>(segs[2 * s + 1]);
        if (a == b || a >= static_cast<VertInd>(nPts) || b >= static_cast<VertInd>(nPts)) continue;
        edges.emplace_back(a, b);
    }
    if (edges.size() < 3) return 2;

    RemoveDuplicatesAndRemapEdges(vertices, edges);

    try {
        g_cdt.insertVertices(vertices);
        g_cdt.insertEdges(edges);
        g_cdt.eraseOuterTrianglesAndHoles();
    } catch (...) {
        return 3;
    }

    g_constraints = edges;
    const int nInput = static_cast<int>(g_cdt.vertices.size());

    std::vector<VertInd> work;
    for (int i = 0; i < nInput; ++i) work.push_back(static_cast<VertInd>(i));
    runEncroach(work);

    if (maxPoints <= 0) maxPoints = nInput * 8 + 2000;
    int guard = 0;
    const int guardMax = maxPoints * 50 + 100000;

    while (static_cast<int>(g_cdt.vertices.size()) < maxPoints && guard++ < guardMax) {
        int worst = -1;
        double worstKey = minAngle;
        for (TriInd ti = 0; ti < static_cast<TriInd>(g_cdt.triangles.size()); ++ti) {
            const Triangle& t = g_cdt.triangles[ti];
            if (!triIsUsable(t)) continue;
            const V2d<double>& a = vtx(t.vertices[0]);
            const V2d<double>& b = vtx(t.vertices[1]);
            const V2d<double>& c = vtx(t.vertices[2]);
            const double ang = minAngleDeg(a, b, c);
            const bool bigArea = (maxArea > 0.0 && triArea(a, b, c) > maxArea);
            if (ang < minAngle - 1e-9 || bigArea) {
                const double key = bigArea ? -1.0 : ang;
                if (worst < 0 || key < worstKey) {
                    worst = static_cast<int>(ti);
                    worstKey = key;
                }
            }
        }
        if (worst < 0) break;

        const Triangle& wt = g_cdt.triangles[worst];
        V2d<double> cc;
        if (!circumcenter(vtx(wt.vertices[0]), vtx(wt.vertices[1]), vtx(wt.vertices[2]), cc))
            continue;

        bool encroached = false;
        for (int s = 0; s < static_cast<int>(g_constraints.size()); ++s) {
            const VertInd u = g_constraints[s].v1();
            const VertInd v = g_constraints[s].v2();
            if (insideDiametral(vtx(u), vtx(v), cc)) {
                const VertInd mid = splitConstraint(s);
                work.push_back(mid);
                encroached = true;
            }
        }
        if (encroached) {
            runEncroach(work);
            continue;
        }

        const VertInd idx = insertPoint(cc);
        work.push_back(idx);
        runEncroach(work);
    }

    collectOutput();

    const auto t1 = std::chrono::steady_clock::now();
    g_time = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return g_outTris.empty() ? 4 : 0;
}

EMSCRIPTEN_KEEPALIVE int     get_points_size()    { return static_cast<int>(g_outPts.size()); }
EMSCRIPTEN_KEEPALIVE double* get_points()         { return g_outPts.data(); }
EMSCRIPTEN_KEEPALIVE int     get_triangles_size() { return static_cast<int>(g_outTris.size()); }
EMSCRIPTEN_KEEPALIVE int*    get_triangles()      { return g_outTris.data(); }
EMSCRIPTEN_KEEPALIVE double  get_last_time_ms()   { return g_time; }

} // extern "C"
