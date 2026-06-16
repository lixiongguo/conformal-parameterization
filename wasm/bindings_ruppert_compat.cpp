// ============================================================
//  Ruppert Delaunay Refinement — 自包含 C++ 实现 (ccall API)
// ------------------------------------------------------------
//  输入: 平面直线图 (PSLG) —— 边界点 + 边界线段 (闭合环, 可含洞)
//  算法: Bowyer-Watson 增量 Delaunay + Ruppert 加密
//        - 线段被侵占 (端点直径圆内有顶点) 时在中点分裂 → 边界自动恢复为 Delaunay 边
//        - 域内劣质三角形 (最小角 < θ 或 面积 > 上界) 插入外心
//        - 外心若侵占线段则改为分裂线段, 不插入外心 (保证终止 & 不越界)
//  输出: 加密后的 (顶点, 三角形)
//
//  导出 (ccall):
//    int    ruppert_run(double* pts,int nCoord, int* segs,int nSegIdx,
//                       double minAngleDeg, double maxArea, int maxPoints)
//    int    get_points_size();   double* get_points();    // [x0,y0,x1,y1,...]
//    int    get_triangles_size();int*    get_triangles();  // [a0,b0,c0,...]
//    double get_last_time_ms();
//    void   dispose();
// ============================================================
#include <vector>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <emscripten.h>

namespace {

struct Pt { double x, y; };

std::vector<Pt>               g_P;        // 所有顶点 (输入 + 超三角 + 加密点)
std::vector<std::array<int,3>> g_tri;     // 三角形 (顶点索引)
std::vector<char>             g_alive;    // 三角形是否存活
std::vector<signed char>      g_inside;   // 三角形质心是否在多边形内 (-1 未知)
std::vector<std::array<int,2>> g_poly;    // 原始多边形线段 (用于内外判定, 固定不变)
std::vector<std::array<int,2>> g_sub;     // 动态子线段 (会被分裂)
int                           g_superStart = 0;

std::vector<double>           g_outPts;
std::vector<int>              g_outTris;
double                        g_time = 0.0;

// ---------------- 几何谓词 ----------------
inline double orient(const Pt& a, const Pt& b, const Pt& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// >0: d 严格在 a,b,c 外接圆内 (要求 a,b,c 逆时针)
inline double inCircle(const Pt& a, const Pt& b, const Pt& c, const Pt& d) {
    double ax = a.x - d.x, ay = a.y - d.y;
    double bx = b.x - d.x, by = b.y - d.y;
    double cx = c.x - d.x, cy = c.y - d.y;
    double a2 = ax * ax + ay * ay;
    double b2 = bx * bx + by * by;
    double c2 = cx * cx + cy * cy;
    return ax * (by * c2 - b2 * cy)
         - ay * (bx * c2 - b2 * cx)
         + a2 * (bx * cy - by * cx);
}

inline bool isSuper(int i) { return i >= g_superStart && i < g_superStart + 3; }

void addTri(int a, int b, int c) {
    if (orient(g_P[a], g_P[b], g_P[c]) < 0) std::swap(b, c);
    g_tri.push_back({a, b, c});
    g_alive.push_back(1);
    g_inside.push_back(-1);
}

// ---------------- Bowyer-Watson 点插入 ----------------
void insertPoint(int pi) {
    const Pt& p = g_P[pi];
    static std::vector<int> bad;
    bad.clear();
    for (int i = 0; i < (int)g_tri.size(); ++i) {
        if (!g_alive[i]) continue;
        const auto& t = g_tri[i];
        if (inCircle(g_P[t[0]], g_P[t[1]], g_P[t[2]], p) > 0.0) bad.push_back(i);
    }
    if (bad.empty()) return; // 数值退化: 点落在所有外接圆外, 跳过

    // 收集空腔边界边 (仅出现在一个坏三角形中的边)
    static std::vector<std::array<int,2>> edges; // (u,v) 有向, 来自坏三角形
    edges.clear();
    for (int bi : bad) {
        const auto& t = g_tri[bi];
        edges.push_back({t[0], t[1]});
        edges.push_back({t[1], t[2]});
        edges.push_back({t[2], t[0]});
    }
    for (int bi : bad) g_alive[bi] = 0;

    int m = (int)edges.size();
    for (int i = 0; i < m; ++i) {
        int u = edges[i][0], v = edges[i][1];
        if (u < 0) continue; // 已抵消
        bool shared = false;
        for (int j = i + 1; j < m; ++j) {
            if (edges[j][0] < 0) continue;
            // 反向边 (v,u) → 内部边, 抵消
            if (edges[j][0] == v && edges[j][1] == u) {
                edges[j][0] = -1; shared = true; break;
            }
        }
        if (!shared) addTri(pi, u, v);
    }
}

int addPoint(const Pt& p) {
    int idx = (int)g_P.size();
    g_P.push_back(p);
    return idx;
}

// ---------------- 内外判定 (奇偶射线法, 用原始多边形) ----------------
bool pointInPoly(const Pt& p) {
    bool in = false;
    for (const auto& s : g_poly) {
        const Pt& a = g_P[s[0]];
        const Pt& b = g_P[s[1]];
        if (((a.y > p.y) != (b.y > p.y)) &&
            (p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x))
            in = !in;
    }
    return in;
}

bool triInside(int i) {
    if (g_inside[i] >= 0) return g_inside[i] == 1;
    const auto& t = g_tri[i];
    Pt c{ (g_P[t[0]].x + g_P[t[1]].x + g_P[t[2]].x) / 3.0,
          (g_P[t[0]].y + g_P[t[1]].y + g_P[t[2]].y) / 3.0 };
    bool ins = pointInPoly(c);
    g_inside[i] = ins ? 1 : 0;
    return ins;
}

// ---------------- 侵占判定 ----------------
// 点 w 是否在以 (a,b) 为直径的圆内 (严格, 含相对容差防数值死循环)
inline bool insideDiametral(const Pt& a, const Pt& b, const Pt& w) {
    double mx = (a.x + b.x) * 0.5, my = (a.y + b.y) * 0.5;
    double dx = w.x - mx, dy = w.y - my;
    double r2 = ((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y)) * 0.25;
    return (dx * dx + dy * dy) < r2 * (1.0 - 1e-9);
}

// 在中点分裂子线段 s, 返回新中点索引并插入三角剖分
int splitSub(int s) {
    int u = g_sub[s][0], v = g_sub[s][1];
    Pt mid{ (g_P[u].x + g_P[v].x) * 0.5, (g_P[u].y + g_P[v].y) * 0.5 };
    int idx = addPoint(mid);
    insertPoint(idx);
    g_sub[s] = {u, idx};
    g_sub.push_back({idx, v});
    return idx;
}

// 维护不变式: 检查新顶点 (worklist) 是否侵占任何子线段, 若侵占则分裂 (递归)
static std::vector<int> g_work;
void runEncroach() {
    while (!g_work.empty()) {
        int w = g_work.back();
        g_work.pop_back();
        bool again = true;
        while (again) {
            again = false;
            for (int s = 0; s < (int)g_sub.size(); ++s) {
                int u = g_sub[s][0], v = g_sub[s][1];
                if (w == u || w == v) continue;
                if (insideDiametral(g_P[u], g_P[v], g_P[w])) {
                    int mid = splitSub(s);
                    g_work.push_back(mid);
                    again = true; // w 可能仍侵占新的半线段
                }
            }
        }
    }
}

// ---------------- 三角形质量 ----------------
double minAngleDeg(const Pt& a, const Pt& b, const Pt& c) {
    double ab = std::hypot(b.x - a.x, b.y - a.y);
    double bc = std::hypot(c.x - b.x, c.y - b.y);
    double ca = std::hypot(a.x - c.x, a.y - c.y);
    auto clamp = [](double v){ return v < -1 ? -1 : (v > 1 ? 1 : v); };
    double A = std::acos(clamp((ab*ab + ca*ca - bc*bc) / (2*ab*ca + 1e-30)));
    double B = std::acos(clamp((ab*ab + bc*bc - ca*ca) / (2*ab*bc + 1e-30)));
    double C = std::acos(clamp((bc*bc + ca*ca - ab*ab) / (2*bc*ca + 1e-30)));
    double mn = A < B ? A : B; if (C < mn) mn = C;
    return mn * 180.0 / M_PI;
}

double triArea(const Pt& a, const Pt& b, const Pt& c) {
    return std::fabs(orient(a, b, c)) * 0.5;
}

bool circumcenter(const Pt& a, const Pt& b, const Pt& c, Pt& out) {
    double D = 2.0 * (a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y));
    if (std::fabs(D) < 1e-18) return false;
    double a2 = a.x*a.x + a.y*a.y, b2 = b.x*b.x + b.y*b.y, c2 = c.x*c.x + c.y*c.y;
    out.x = (a2*(b.y - c.y) + b2*(c.y - a.y) + c2*(a.y - b.y)) / D;
    out.y = (a2*(c.x - b.x) + b2*(a.x - c.x) + c2*(b.x - a.x)) / D;
    return true;
}

void resetAll() {
    g_P.clear(); g_tri.clear(); g_alive.clear(); g_inside.clear();
    g_poly.clear(); g_sub.clear(); g_work.clear();
    g_outPts.clear(); g_outTris.clear();
    g_superStart = 0;
}

} // namespace

// ============================================================
//                       导出 API
// ============================================================
extern "C" {

EMSCRIPTEN_KEEPALIVE
void dispose() { resetAll(); }

EMSCRIPTEN_KEEPALIVE
int ruppert_run(double* pts, int nCoord, int* segs, int nSegIdx,
                double minAngle, double maxArea, int maxPoints) {
    auto t0 = std::chrono::steady_clock::now();
    resetAll();

    int nPts = nCoord / 2;
    int nSeg = nSegIdx / 2;
    if (nPts < 3 || nSeg < 3) return 1;

    // 1. 读入输入顶点 (带去重: 合并极近点)
    std::vector<int> remap(nPts);
    double bbMinX = 1e300, bbMinY = 1e300, bbMaxX = -1e300, bbMaxY = -1e300;
    for (int i = 0; i < nPts; ++i) {
        double x = pts[2*i], y = pts[2*i+1];
        bbMinX = x < bbMinX ? x : bbMinX; bbMaxX = x > bbMaxX ? x : bbMaxX;
        bbMinY = y < bbMinY ? y : bbMinY; bbMaxY = y > bbMaxY ? y : bbMaxY;
    }
    double diag = std::hypot(bbMaxX - bbMinX, bbMaxY - bbMinY);
    double eps = diag * 1e-7 + 1e-12;
    double eps2 = eps * eps;
    for (int i = 0; i < nPts; ++i) {
        Pt p{ pts[2*i], pts[2*i+1] };
        int found = -1;
        for (int j = 0; j < (int)g_P.size(); ++j) {
            double dx = g_P[j].x - p.x, dy = g_P[j].y - p.y;
            if (dx*dx + dy*dy < eps2) { found = j; break; }
        }
        if (found >= 0) remap[i] = found;
        else { remap[i] = (int)g_P.size(); g_P.push_back(p); }
    }

    // 2. 线段 (去重映射 + 去除退化/重复)
    for (int s = 0; s < nSeg; ++s) {
        int a = remap[segs[2*s]];
        int b = remap[segs[2*s+1]];
        if (a == b) continue;
        g_poly.push_back({a, b});
    }
    g_sub = g_poly;
    if (g_poly.size() < 3) return 2;

    // 3. 超三角形 (包含所有点)
    double cx = (bbMinX + bbMaxX) * 0.5, cy = (bbMinY + bbMaxY) * 0.5;
    double D = (diag > 1e-12 ? diag : 1.0) * 20.0;
    g_superStart = (int)g_P.size();
    g_P.push_back({ cx - D, cy - D });
    g_P.push_back({ cx + D, cy - D });
    g_P.push_back({ cx,      cy + D });
    addTri(g_superStart, g_superStart + 1, g_superStart + 2);

    // 4. 插入所有输入顶点
    int nInput = g_superStart; // 输入点数 (去重后)
    for (int i = 0; i < nInput; ++i) insertPoint(i);

    // 5. 初始边界恢复: 所有输入顶点入队检查侵占
    g_work.clear();
    for (int i = 0; i < nInput; ++i) g_work.push_back(i);
    runEncroach();

    // 6. Ruppert 加密
    if (maxPoints <= 0) maxPoints = nInput * 8 + 2000;
    int guard = 0, guardMax = maxPoints * 50 + 100000;
    while ((int)g_P.size() < maxPoints && guard++ < guardMax) {
        // 找最劣 (角度最小) 的域内劣质三角形
        int worst = -1; double worstAng = minAngle;
        for (int i = 0; i < (int)g_tri.size(); ++i) {
            if (!g_alive[i]) continue;
            const auto& t = g_tri[i];
            if (isSuper(t[0]) || isSuper(t[1]) || isSuper(t[2])) continue;
            if (!triInside(i)) continue;
            const Pt& a = g_P[t[0]]; const Pt& b = g_P[t[1]]; const Pt& c = g_P[t[2]];
            double ang = minAngleDeg(a, b, c);
            bool bigArea = (maxArea > 0.0 && triArea(a, b, c) > maxArea);
            if (ang < minAngle - 1e-9 || bigArea) {
                // 选最小角的三角形优先处理
                double key = bigArea ? -1.0 : ang;
                if (worst < 0 || key < worstAng) { worst = i; worstAng = key; }
            }
        }
        if (worst < 0) break;

        const auto& t = g_tri[worst];
        Pt cc;
        if (!circumcenter(g_P[t[0]], g_P[t[1]], g_P[t[2]], cc)) { g_alive[worst] = 0; continue; }

        // 外心是否侵占子线段?
        bool encroached = false;
        for (int s = 0; s < (int)g_sub.size(); ++s) {
            if (insideDiametral(g_P[g_sub[s][0]], g_P[g_sub[s][1]], cc)) {
                int mid = splitSub(s);
                g_work.push_back(mid);
                encroached = true;
            }
        }
        if (encroached) { runEncroach(); continue; }

        // 否则插入外心
        int idx = addPoint(cc);
        insertPoint(idx);
        g_work.push_back(idx);
        runEncroach();
    }

    // 7. 收集输出 (域内、不含超三角顶点的存活三角形)
    std::vector<int> used(g_P.size(), -1);
    for (int i = 0; i < (int)g_tri.size(); ++i) {
        if (!g_alive[i]) continue;
        const auto& t = g_tri[i];
        if (isSuper(t[0]) || isSuper(t[1]) || isSuper(t[2])) continue;
        if (!triInside(i)) continue;
        for (int k = 0; k < 3; ++k) {
            if (used[t[k]] < 0) {
                used[t[k]] = (int)(g_outPts.size() / 2);
                g_outPts.push_back(g_P[t[k]].x);
                g_outPts.push_back(g_P[t[k]].y);
            }
        }
        g_outTris.push_back(used[t[0]]);
        g_outTris.push_back(used[t[1]]);
        g_outTris.push_back(used[t[2]]);
    }

    auto t1 = std::chrono::steady_clock::now();
    g_time = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return 0;
}

EMSCRIPTEN_KEEPALIVE int     get_points_size()    { return (int)g_outPts.size(); }
EMSCRIPTEN_KEEPALIVE double* get_points()         { return g_outPts.data(); }
EMSCRIPTEN_KEEPALIVE int     get_triangles_size() { return (int)g_outTris.size(); }
EMSCRIPTEN_KEEPALIVE int*    get_triangles()      { return g_outTris.data(); }
EMSCRIPTEN_KEEPALIVE double  get_last_time_ms()   { return g_time; }

} // extern "C"
