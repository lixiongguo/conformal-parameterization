// ============================================================
//  uv_unwrap_simple 兼容层 — 所有参数化方法 (ccall API)
// ============================================================
#include "Mesh.h"
#include "MeshIO.h"
#include "Lscm.h"
#include "Scp.h"
#include "AbfPlusPlus.h"
#include "LinAbf.h"
#include "ARAP.h"
#include "Tutte.h"
#include "Solver.h"
#include "CirclePatterns.h"
#include "Cetm.h"
#include "RicciFlow.h"
#include <chrono>
#include <vector>
#include <cmath>
#include <limits>
#include <emscripten.h>

static Mesh*               g_mesh = nullptr;
static std::vector<double> g_uv;
static double              g_last_time = 0.0;
static std::vector<double> g_qc_errors;
static std::vector<double> g_qc_colors;

// ========== 辅助 ==========
static void buildMeshFrom(double* flatPos, int posLen, int* flatFace, int faceLen) {
    int nV = posLen / 3, nF = faceLen / 3;
    MeshData data;
    data.positions.reserve(nV);
    for (int i = 0; i < nV; ++i)
        data.positions.push_back(Eigen::Vector3d(flatPos[i*3], flatPos[i*3+1], flatPos[i*3+2]));
    data.indices.reserve(nF);
    for (int f = 0; f < nF; ++f) {
        data.indices.push_back({Index(flatFace[f*3], -1, -1),
                                Index(flatFace[f*3+1], -1, -1),
                                Index(flatFace[f*3+2], -1, -1)});
    }
    MeshIO::buildMesh(data, *g_mesh);
}

static void extractUV() {
    g_uv.clear();
    for (auto& v : g_mesh->vertices) {
        g_uv.push_back(v.uv.x());
        g_uv.push_back(v.uv.y());
    }
}

static void recordTime(const std::chrono::steady_clock::time_point& t0) {
    auto t1 = std::chrono::steady_clock::now();
    g_last_time = std::chrono::duration<double, std::milli>(t1 - t0).count();
}

static bool hasValidUvSpread() {
    if (!g_mesh || g_mesh->vertices.empty()) return false;
    double minU = std::numeric_limits<double>::infinity(), maxU = -minU;
    double minV = std::numeric_limits<double>::infinity(), maxV = -minV;
    for (const auto& v : g_mesh->vertices) {
        minU = std::min(minU, v.uv.x()); maxU = std::max(maxU, v.uv.x());
        minV = std::min(minV, v.uv.y()); maxV = std::max(maxV, v.uv.y());
    }
    return (maxU - minU) >= 1e-12 || (maxV - minV) >= 1e-12;
}

// ========== 导出 ==========
extern "C" {

EMSCRIPTEN_KEEPALIVE
void dispose() {
    if (g_mesh) { delete g_mesh; g_mesh = nullptr; }
    g_uv.clear(); g_qc_errors.clear(); g_qc_colors.clear();
    g_last_time = 0.0;
}

EMSCRIPTEN_KEEPALIVE
int solve_lscm(double* pos, int posLen, int* face, int faceLen, int /*a1*/, int /*a2*/) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    Lscm lscm(*g_mesh);
    auto t0 = std::chrono::steady_clock::now();
    lscm.parameterize();
    recordTime(t0);
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_scp(double* pos, int posLen, int* face, int faceLen) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    Scp scp(*g_mesh);
    auto t0 = std::chrono::steady_clock::now();
    scp.parameterize();
    recordTime(t0);
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_tutte_circle(double* pos, int posLen, int* face, int faceLen, int w) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    Tutte t(*g_mesh,
            (w == 1) ? TutteBoundary::CIRCLE : TutteBoundary::CIRCLE,
            (w == 1) ? TutteWeight::UNIFORM : TutteWeight::COTAN);
    auto t0 = std::chrono::steady_clock::now();
    t.parameterize();
    recordTime(t0);
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_tutte_square(double* pos, int posLen, int* face, int faceLen, int w) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    Tutte t(*g_mesh, TutteBoundary::SQUARE,
            (w == 1) ? TutteWeight::UNIFORM : TutteWeight::COTAN);
    auto t0 = std::chrono::steady_clock::now();
    t.parameterize();
    recordTime(t0);
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_linabf(double* pos, int posLen, int* face, int faceLen) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    if (g_mesh->boundaries.empty()) return -3;
    LinAbf abf(*g_mesh);
    auto t0 = std::chrono::steady_clock::now();
    abf.parameterize();
    recordTime(t0);
    if (!hasValidUvSpread()) return -2;
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_abfpp(double* pos, int posLen, int* face, int faceLen) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    if (g_mesh->boundaries.empty()) return -3;
    AbfPlusPlus abf(*g_mesh);
    auto t0 = std::chrono::steady_clock::now();
    abf.parameterize();
    recordTime(t0);
    if (!hasValidUvSpread()) return -2;
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_arap(double* pos, int posLen, int* face, int faceLen, int maxIter) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    ARAP arap(*g_mesh, maxIter > 0 ? maxIter : 30);
    auto t0 = std::chrono::steady_clock::now();
    arap.parameterize();
    recordTime(t0);
    extractUV();
    return 0;
}

// ---- Circle Pattern / CETM / Ricci ----
EMSCRIPTEN_KEEPALIVE
int solve_cp(double* pos, int posLen, int* face, int faceLen, int optScheme,
             int* coneIdx, int coneIdxLen, double* coneAngles, int coneAnglesLen) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    auto t0 = std::chrono::steady_clock::now();
    g_mesh->delaunayize();
    CirclePatterns p(*g_mesh, optScheme);
    if (coneIdx && coneAngles && coneIdxLen > 0) {
        std::vector<int> idx(coneIdx, coneIdx + coneIdxLen);
        std::vector<double> angles(coneAngles, coneAngles + std::min(coneIdxLen, coneAnglesLen));
        p.setConeSingulars(idx, angles);
    }
    p.parameterize();
    recordTime(t0);
    if (!hasValidUvSpread()) return -2;
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_cetm(double* pos, int posLen, int* face, int faceLen, int optScheme,
               int* coneIdx, int coneIdxLen, double* coneAngles, int coneAnglesLen) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    auto t0 = std::chrono::steady_clock::now();
    g_mesh->delaunayize();
    Cetm p(*g_mesh, optScheme);
    if (coneIdx && coneAngles && coneIdxLen > 0) {
        std::vector<int> idx(coneIdx, coneIdx + coneIdxLen);
        std::vector<double> angles(coneAngles, coneAngles + std::min(coneIdxLen, coneAnglesLen));
        p.setConeSingulars(idx, angles);
    }
    p.parameterize();
    recordTime(t0);
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_ricci(double* pos, int posLen, int* face, int faceLen, int optScheme,
                int* coneIdx, int coneIdxLen, double* coneAngles, int coneAnglesLen) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    auto t0 = std::chrono::steady_clock::now();
    g_mesh->delaunayize();
    RicciFlow p(*g_mesh, optScheme);
    if (coneIdx && coneAngles && coneIdxLen > 0) {
        std::vector<int> idx(coneIdx, coneIdx + coneIdxLen);
        std::vector<double> angles(coneAngles, coneAngles + std::min(coneIdxLen, coneAnglesLen));
        p.setConeSingulars(idx, angles);
    }
    p.parameterize();
    recordTime(t0);
    extractUV();
    return 0;
}

// ---- UV 结果 ----
EMSCRIPTEN_KEEPALIVE
int get_uv_result_size() { return (int)g_uv.size(); }
EMSCRIPTEN_KEEPALIVE
double* get_uv_result() { return g_uv.data(); }

// ---- 耗时 ----
EMSCRIPTEN_KEEPALIVE
double get_last_time_ms() { return g_last_time; }

// ---- load_mesh_with_uv (QC 误差用) ----
EMSCRIPTEN_KEEPALIVE
int load_mesh_with_uv(double* pos, int posLen, int* face, int faceLen,
                       double* uvFlat, int /*uvLen*/) {
    dispose(); g_mesh = new Mesh();
    buildMeshFrom(pos, posLen, face, faceLen);
    int nV = (int)g_mesh->vertices.size();
    for (int i = 0; i < nV; ++i) {
        g_mesh->vertices[i].uv.x() = uvFlat[i*2];
        g_mesh->vertices[i].uv.y() = uvFlat[i*2+1];
    }
    return 0;
}

// ---- QC 误差 ----
EMSCRIPTEN_KEEPALIVE
int compute_qc_error() {
    if (!g_mesh) return 0;
    g_qc_errors.clear(); g_qc_colors.clear();
    int nF = (int)g_mesh->faces.size();

    for (auto& f : g_mesh->faces) {
        if (f.isBoundary()) continue;
        HalfEdgeIter he = f.he;
        Eigen::Vector2d u0 = he->vertex->uv;
        Eigen::Vector2d u1 = he->next->vertex->uv;
        Eigen::Vector2d u2 = he->next->next->vertex->uv;
        Eigen::Vector2d e1 = u1 - u0, e2 = u2 - u0;
        double jac = e1.x() * e2.y() - e1.y() * e2.x();
        double error;
        if (jac > 1e-12) {
            double s1 = e1.squaredNorm(), s2 = e2.squaredNorm();
            double sMax = std::max(s1, s2), sMin = std::min(s1, s2);
            error = (sMin > 1e-12) ? std::sqrt(sMax / sMin) : 10.0;
            error = std::min(error, 10.0);
        } else { error = 10.0; }
        g_qc_errors.push_back(error);
        double t = std::min((error - 1.0) / 1.0, 1.0);
        g_qc_colors.push_back(t);
        g_qc_colors.push_back(1.0 - t);
        g_qc_colors.push_back(0.15);
    }
    return nF;
}

EMSCRIPTEN_KEEPALIVE
int get_qc_errors_size() { return (int)g_qc_errors.size(); }
EMSCRIPTEN_KEEPALIVE
double* get_qc_errors() { return g_qc_errors.data(); }
EMSCRIPTEN_KEEPALIVE
int get_qc_colors_size() { return (int)g_qc_colors.size(); }
EMSCRIPTEN_KEEPALIVE
double* get_qc_colors() { return g_qc_colors.data(); }

} // extern "C"
