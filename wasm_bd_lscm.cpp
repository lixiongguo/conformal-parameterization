/**
 * wasm_bd_lscm.cpp — BD_LSCM WebAssembly 入口 (C 接口)
 *
 * 将 bounded_distortion.h + 内置 LSCM 初始映射编译为 WASM 模块，
 * 通过 export_name="BDLSCMSolver" 导出给 JS 层调用。
 *
 * 使用 C 风格 flat API（与 lscm_solver.wasm 一致）:
 *   solve_bd_lscm(posPtr, posLen, facePtr, faceLen, C, maxIter) -> int
 *   get_bd_uv_result()        -> double*
 *   get_bd_uv_result_size()   -> int
 *   get_bd_last_time_ms()     -> double
 *   bd_dispose()              -> void
 *
 * 编译（独立于 libigl，仅依赖 Eigen）:
 *   emcc wasm_bd_lscm.cpp
 *       -I../conformal-parameterization/deps/Eigen
 *       -I.
 *       -s MODULARIZE=1 -s EXPORT_NAME="BDLSCMSolver"
 *       -s ALLOW_MEMORY_GROWTH=1 -s WASM=1
 *       -std=c++17 -O2
 *       -o ../../assets/wasm/bd_lscm_solver.js
 */

#include <emscripten.h>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <set>

#include "bounded_distortion.h"

// ==================== 全局状态 ====================
static double g_bd_lastTimeMs = 0.0;
static std::vector<double> g_bd_uv_result;

// ==================== 辅助：边界检测 (简单边计数) ====================

/**
 * 找出所有只被一个三角形引用的边（网格边界）。
 * 输入: F = #F x 3  输出: bnd = 按环绕顺序排列的边界顶点索引
 */
static void detect_boundary_loop(
    const bdm::MatrixXi& F,
    bdm::VectorXi& bnd)
{
    std::map<std::pair<int,int>, int> edge_count;
    const int nF = F.rows();

    for (int j = 0; j < nF; ++j) {
        for (int k = 0; k < 3; ++k) {
            int a = F(j, k);
            int b = F(j, (k+1) % 3);
            if (a > b) std::swap(a, b);
            edge_count[{a, b}]++;
        }
    }

    // 收集边界边
    std::vector<std::pair<int,int>> b_edges;
    for (auto& kv : edge_count) {
        if (kv.second == 1) {
            b_edges.push_back(kv.first);
        }
    }

    if (b_edges.empty()) {
        bnd.resize(0);
        return;
    }

    // 构建邻接表并沿边界走一圈
    std::map<int, std::vector<int>> adj;
    for (auto& e : b_edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }

    std::vector<int> loop;
    int curr = b_edges[0].first;
    int prev = -1;
    do {
        loop.push_back(curr);
        int next = -1;
        for (int nbr : adj[curr]) {
            if (nbr != prev) { next = nbr; break; }
        }
        if (next < 0) break;
        prev = curr;
        curr = next;
    } while (curr != loop[0] && (int)loop.size() < (int)b_edges.size() + 1);

    bnd.resize(loop.size());
    for (size_t i = 0; i < loop.size(); ++i) bnd(i) = loop[i];
}

/** 将边界顶点映射到单位圆上 */
static void map_boundary_to_circle(
    const bdm::MatrixXd& V,
    const bdm::VectorXi& bnd,
    bdm::MatrixXd& bnd_uv)
{
    const int nB = bnd.size();
    bnd_uv.resize(nB, 2);
    for (int k = 0; k < nB; ++k) {
        double angle = 2.0 * M_PI * k / nB;
        bnd_uv(k, 0) = std::cos(angle);
        bnd_uv(k, 1) = std::sin(angle);
    }
}

/** 将顶点数组从 C double[] 复制到 Eigen MatrixXd */
static void copy_to_eigen(const double* src, int nV, bdm::MatrixXd& V)
{
    V.resize(nV, 3);
    for (int i = 0; i < nV; ++i) {
        V(i, 0) = src[i * 3];
        V(i, 1) = src[i * 3 + 1];
        V(i, 2) = src[i * 3 + 2];
    }
}

static void copy_to_eigen_int(const int* src, int nF, bdm::MatrixXi& F)
{
    F.resize(nF, 3);
    for (int i = 0; i < nF; ++i) {
        F(i, 0) = src[i * 3];
        F(i, 1) = src[i * 3 + 1];
        F(i, 2) = src[i * 3 + 2];
    }
}

// ==================== 导出的 C 函数 ====================

extern "C" {

/**
 * solve_bd_lscm — 执行 Bounded Distortion LSCM 展开
 *
 * @param posPtr   顶点坐标 [x0,y0,z0, x1,y1,z1, ...] (double)
 * @param posLen   顶点坐标数组长度 (= 3 * numVertices)
 * @param facePtr  面索引   [i0,i1,i2, i0,i1,i2, ...] (int, 0-based)
 * @param faceLen  面数组长度 (= 3 * numFaces)
 * @param C        扭曲上限 (sigma1/sigma2 <= C), 典型值 2~20
 * @param maxIter  最大迭代次数, 典型值 30~50
 * @return 0 成功, -1 失败
 */
EMSCRIPTEN_KEEPALIVE
int solve_bd_lscm(
    double* posPtr, int posLen,
    int*    facePtr, int faceLen,
    double C, int maxIter)
{
    if (!posPtr || !facePtr || posLen < 9 || faceLen < 3) return -1;

    const int nV = posLen / 3;
    const int nF = faceLen / 3;

    auto t0 = std::chrono::high_resolution_clock::now();

    try {
        // 1. 复制数据到 Eigen
        bdm::MatrixXd V;
        bdm::MatrixXi F;
        copy_to_eigen(posPtr, nV, V);
        copy_to_eigen_int(facePtr, nF, F);

        // 2. 检测边界并映射到圆
        bdm::VectorXi bnd;
        detect_boundary_loop(F, bnd);
        if (bnd.size() < 3) return -1;

        bdm::MatrixXd bnd_uv;
        map_boundary_to_circle(V, bnd, bnd_uv);

        // 3. 计算 LSCM 初始映射
        bdm::MatrixXd U_init;
        bdm::lscm_initial_map(V, F, bnd, bnd_uv, U_init);

        // 4. Bounded Distortion 优化
        bdm::BDOptions opts;
        opts.C = C;
        opts.max_iter = maxIter;
        opts.tol = 1e-5;
        opts.verbose = false;
        opts.use_arap = false;  // LSCM 能量

        bdm::BDResult result;
        bdm::bounded_distortion_map(V, F, bnd, bnd_uv, U_init, opts, result);

        // 5. 归一化 UV 到 [0,1]×[0,1]
        const bdm::MatrixXd& U = result.U;
        double min_u = U.col(0).minCoeff();
        double max_u = U.col(0).maxCoeff();
        double min_v = U.col(1).minCoeff();
        double max_v = U.col(1).maxCoeff();

        double range_u = max_u - min_u;
        double range_v = max_v - min_v;
        if (range_u < 1e-10) range_u = 1.0;
        if (range_v < 1e-10) range_v = 1.0;

        g_bd_uv_result.clear();
        g_bd_uv_result.reserve(nV * 2);
        for (int i = 0; i < nV; ++i) {
            g_bd_uv_result.push_back((U(i,0) - min_u) / range_u);
            g_bd_uv_result.push_back((U(i,1) - min_v) / range_v);
        }

    } catch (...) {
        return -1;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    g_bd_lastTimeMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

    return 0;
}

EMSCRIPTEN_KEEPALIVE
double* get_bd_uv_result() {
    return g_bd_uv_result.empty() ? nullptr : g_bd_uv_result.data();
}

EMSCRIPTEN_KEEPALIVE
int get_bd_uv_result_size() {
    return (int)g_bd_uv_result.size();
}

EMSCRIPTEN_KEEPALIVE
double get_bd_last_time_ms() {
    return g_bd_lastTimeMs;
}

EMSCRIPTEN_KEEPALIVE
void bd_dispose() {
    g_bd_uv_result.clear();
    g_bd_uv_result.shrink_to_fit();
}

} // extern "C"
