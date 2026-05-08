/**
 * wasm_tutte_arap.cpp - WASM entry for Tutte (circle/square) & ARAP
 * Exports: solve_tutte_circle, solve_tutte_square, solve_arap
 */

#include <emscripten.h>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <cmath>

#include "Mesh.h"
#include "MeshIO.h"
#include "Tutte.h"
#include "ARAP.h"

// ==================== Globals ====================
static Mesh* g_mesh = nullptr;
static double g_lastTimeMs = 0.0;
static std::vector<double> g_uv_result;

// ==================== Helpers ====================

static bool loadMesh(const double* positions, int posLen,
                     const int* faces, int faceLen)
{
    if (g_mesh) { delete g_mesh; g_mesh = nullptr; }
    if (posLen < 9 || faceLen < 3) return false;

    g_mesh = new Mesh();
    size_t nV = posLen / 3, nF = faceLen / 3;

    std::stringstream ss;
    for (size_t i = 0; i < nV; i++)
        ss << "v " << positions[i*3] << " " << positions[i*3+1] << " " << positions[i*3+2] << "\n";
    for (size_t i = 0; i < nF; i++)
        ss << "f " << faces[i*3]+1 << " " << faces[i*3+1]+1 << " " << faces[i*3+2]+1 << "\n";

    std::string str = ss.str();
    FILE* fp = fopen("/tmp/ta_input.obj", "wb");
    if (!fp) return false;
    fwrite(str.c_str(), 1, str.size(), fp);
    fclose(fp);

    if (!g_mesh->read("/tmp/ta_input.obj")) { delete g_mesh; g_mesh = nullptr; return false; }
    return true;
}

static void extractUV()
{
    g_uv_result.clear();
    if (!g_mesh) return;
    g_uv_result.reserve(g_mesh->vertices.size() * 2);

    double minU = 1e30, maxU = -1e30, minV = 1e30, maxV = -1e30;
    for (auto& v : g_mesh->vertices) {
        if (v.uv.x() < minU) minU = v.uv.x();
        if (v.uv.x() > maxU) maxU = v.uv.x();
        if (v.uv.y() < minV) minV = v.uv.y();
        if (v.uv.y() > maxV) maxV = v.uv.y();
    }
    double rU = std::max(maxU - minU, 1e-10);
    double rV = std::max(maxV - minV, 1e-10);
    for (auto& v : g_mesh->vertices) {
        g_uv_result.push_back((v.uv.x() - minU) / rU);
        g_uv_result.push_back((v.uv.y() - minV) / rV);
    }
}

// ==================== Exports ====================
extern "C" {

EMSCRIPTEN_KEEPALIVE
int solve_tutte_circle(double* pos, int posLen, int* faces, int faceLen) {
    if (!loadMesh(pos, posLen, faces, faceLen)) return -1;
    auto t0 = std::chrono::high_resolution_clock::now();
    g_mesh->delaunayize();
    Tutte p(*g_mesh, TutteBoundary::CIRCLE);
    p.parameterize();
    auto t1 = std::chrono::high_resolution_clock::now();
    g_lastTimeMs = std::chrono::duration<double,std::milli>(t1-t0).count();
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_tutte_square(double* pos, int posLen, int* faces, int faceLen) {
    if (!loadMesh(pos, posLen, faces, faceLen)) return -1;
    auto t0 = std::chrono::high_resolution_clock::now();
    g_mesh->delaunayize();
    Tutte p(*g_mesh, TutteBoundary::SQUARE);
    p.parameterize();
    auto t1 = std::chrono::high_resolution_clock::now();
    g_lastTimeMs = std::chrono::duration<double,std::milli>(t1-t0).count();
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int solve_arap(double* pos, int posLen, int* faces, int faceLen, int maxIter) {
    if (!loadMesh(pos, posLen, faces, faceLen)) return -1;
    if (maxIter < 1) maxIter = 30;
    auto t0 = std::chrono::high_resolution_clock::now();
    g_mesh->delaunayize();
    ARAP p(*g_mesh, maxIter);
    p.parameterize();
    auto t1 = std::chrono::high_resolution_clock::now();
    g_lastTimeMs = std::chrono::duration<double,std::milli>(t1-t0).count();
    extractUV();
    return 0;
}

EMSCRIPTEN_KEEPALIVE double* get_ta_uv_result() { return g_uv_result.empty() ? nullptr : g_uv_result.data(); }
EMSCRIPTEN_KEEPALIVE int get_ta_uv_result_size() { return (int)g_uv_result.size(); }
EMSCRIPTEN_KEEPALIVE double get_ta_last_time_ms() { return g_lastTimeMs; }

EMSCRIPTEN_KEEPALIVE
void ta_dispose() {
    if (g_mesh) { delete g_mesh; g_mesh = nullptr; }
    g_uv_result.clear();
    g_uv_result.shrink_to_fit();
}

} // extern "C"
