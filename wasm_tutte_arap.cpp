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
#include "QcError.h"

// ==================== Globals ====================
static Mesh* g_mesh = nullptr;
static double g_lastTimeMs = 0.0;
static std::vector<double> g_uv_result;
static std::vector<double> g_qc_errors; // per-face QC error [0..N_faces)
static std::vector<double> g_qc_colors; // per-face RGB color [0..3*N_faces)

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
int compute_qc_error() {
    if (!g_mesh) return -1;
    g_qc_errors.clear();
    g_qc_colors.clear();
    int nFaces = 0;

    for (FaceCIter f = g_mesh->faces.begin(); f != g_mesh->faces.end(); f++) {
        if (f->isBoundary()) continue;
        nFaces++;

        std::vector<Eigen::Vector3d> p, q;
        HalfEdgeCIter he = f->he;
        do {
            p.push_back(he->vertex->position);
            q.push_back(Eigen::Vector3d(he->vertex->uv.x(), he->vertex->uv.y(), 0));
            he = he->next;
        } while (he != f->he);

        double err = QuasiConformalError::compute(p, q);
        g_qc_errors.push_back(err);

        Eigen::Vector3d col = QuasiConformalError::color(err);
        g_qc_colors.push_back(col.x());
        g_qc_colors.push_back(col.y());
        g_qc_colors.push_back(col.z());
    }

    return nFaces;
}

// Load mesh WITH existing UV coords (for QC computation after any param method)
EMSCRIPTEN_KEEPALIVE
int load_mesh_with_uv(double* pos, int posLen, int* faces, int faceLen,
                      double* uvs, int uvLen) {
    if (!loadMesh(pos, posLen, faces, faceLen)) return -1;
    int nV = posLen / 3;
    if (uvLen < nV * 2) return -1;
    for (int i = 0; i < nV; i++) {
        g_mesh->vertices[i].uv = Eigen::Vector2d(uvs[i*2], uvs[i*2+1]);
    }
    return 0;
}

EMSCRIPTEN_KEEPALIVE double* get_qc_errors()  { return g_qc_errors.empty() ? nullptr : g_qc_errors.data(); }
EMSCRIPTEN_KEEPALIVE int   get_qc_errors_size()  { return (int)g_qc_errors.size(); }
EMSCRIPTEN_KEEPALIVE double* get_qc_colors()  { return g_qc_colors.empty() ? nullptr : g_qc_colors.data(); }
EMSCRIPTEN_KEEPALIVE int   get_qc_colors_size()  { return (int)g_qc_colors.size(); }

EMSCRIPTEN_KEEPALIVE
void ta_dispose() {
    if (g_mesh) { delete g_mesh; g_mesh = nullptr; }
    g_uv_result.clear();
    g_uv_result.shrink_to_fit();
    g_qc_errors.clear();
    g_qc_colors.clear();
}

} // extern "C"
