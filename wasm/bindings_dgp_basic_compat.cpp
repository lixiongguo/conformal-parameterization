// ============================================================
//  dgp_basic 兼容层 — 高斯曲率 (ccall API)
// ============================================================
#include "Mesh.h"
#include "MeshIO.h"
#include "GaussianCurvature.h"
#include <chrono>
#include <vector>
#include <emscripten.h>

static Mesh*               g_mesh = nullptr;
static std::vector<double> g_gc;
static double              g_last_time = 0.0;

extern "C" {

EMSCRIPTEN_KEEPALIVE
void dispose() {
    if (g_mesh) { delete g_mesh; g_mesh = nullptr; }
    g_gc.clear();
    g_last_time = 0.0;
}

EMSCRIPTEN_KEEPALIVE
int load_mesh(double* flatPos, int posLen, int* flatFace, int faceLen) {
    dispose();
    g_mesh = new Mesh();

    int nV = posLen / 3;
    int nF = faceLen / 3;

    MeshData data;
    data.positions.reserve(nV);
    for (int i = 0; i < nV; ++i)
        data.positions.push_back(Eigen::Vector3d(flatPos[i*3], flatPos[i*3+1], flatPos[i*3+2]));

    data.indices.reserve(nF);
    for (int f = 0; f < nF; ++f) {
        std::vector<Index> tri;
        tri.push_back(Index(flatFace[f*3], -1, -1));
        tri.push_back(Index(flatFace[f*3+1], -1, -1));
        tri.push_back(Index(flatFace[f*3+2], -1, -1));
        data.indices.push_back(tri);
    }

    return MeshIO::buildMesh(data, *g_mesh) ? 0 : -1;
}

EMSCRIPTEN_KEEPALIVE
int compute_gauss_curvature_per_area() {
    if (!g_mesh) return -1;
    auto t0 = std::chrono::steady_clock::now();

    Eigen::VectorXd K = geometry::gaussianCurvaturePerArea(*g_mesh);
    g_gc.assign(K.data(), K.data() + K.size());

    auto t1 = std::chrono::steady_clock::now();
    g_last_time = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int get_gc_result_size() { return (int)g_gc.size(); }

EMSCRIPTEN_KEEPALIVE
double* get_gc_result() { return g_gc.data(); }

EMSCRIPTEN_KEEPALIVE
double get_last_time_ms() { return g_last_time; }

} // extern "C"
