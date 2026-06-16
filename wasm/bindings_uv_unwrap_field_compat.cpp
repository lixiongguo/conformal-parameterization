/**
 * uv_unwrap_field WASM wrapper.
 *
 * Pipelines:
 *   QuadCover: step1 principal field -> step2 matching -> step3 quadcover
 *   MIQ:       step1 (optional warm start) -> solve_miq (MIQP cross field + Poisson UV)
 */

#include <emscripten.h>
#include <emscripten/html5.h>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>

#include "Mesh.h"
#include "MeshIO.h"
#include "QuadCover.h"
#include "MIQQuad.h"
#include "HolomorphicOneForm.h"

namespace {

class QuadCoverWasm : public QuadCover {
public:
    using QuadCover::QuadCover;
    using QuadCover::initMeshData;
    using QuadCover::estimatePrincipalCurvature;
    using QuadCover::computeFaceTheta;
    using QuadCover::syncFaceDirsFromTheta;
    using QuadCover::smoothCrossField;
    using QuadCover::computeMatching;
    using QuadCover::computeLayerShift;
    using QuadCover::getCrossField;
    using QuadCover::matching;
    using QuadCover::edgeList;
};

class MIQQuadWasm : public MIQQuad {
public:
    using MIQQuad::MIQQuad;
    using MIQQuad::faceTheta;
    using MIQQuad::edgeJumps;
    using MIQQuad::numFaces;
    using MIQQuad::numEdges;
};

Mesh* g_mesh = nullptr;
QuadCoverWasm* g_qc = nullptr;
MIQQuadWasm* g_miq = nullptr;
bool g_qc_inited = false;

std::vector<double> g_face_dirs;
std::vector<int> g_matching;
std::vector<double> g_face_theta;
std::vector<int> g_edge_jumps;
std::vector<double> g_uv;
double g_last_time_ms = 0.0;
double g_last_miq_energy = 0.0;

bool loadMeshFromArrays(const double* positions, int posLen, const int* faces, int faceLen)
{
    if (g_qc) {
        delete g_qc;
        g_qc = nullptr;
    }
    if (g_miq) {
        delete g_miq;
        g_miq = nullptr;
    }
    if (g_mesh) {
        delete g_mesh;
        g_mesh = nullptr;
    }
    g_qc_inited = false;
    g_face_dirs.clear();
    g_matching.clear();
    g_face_theta.clear();
    g_edge_jumps.clear();
    g_uv.clear();
    g_last_miq_energy = 0.0;

    if (!positions || !faces) return false;
    if (posLen < 9 || faceLen < 3) return false;
    if (posLen % 3 != 0 || faceLen % 3 != 0) return false;

    const size_t nV = static_cast<size_t>(posLen) / 3;
    const size_t nF = static_cast<size_t>(faceLen) / 3;

    std::stringstream ss;
    for (size_t i = 0; i < nV; ++i) {
        ss << "v " << positions[i * 3] << " " << positions[i * 3 + 1] << " "
           << positions[i * 3 + 2] << "\n";
    }
    for (size_t i = 0; i < nF; ++i) {
        ss << "f " << faces[i * 3] + 1 << " " << faces[i * 3 + 1] + 1 << " "
           << faces[i * 3 + 2] + 1 << "\n";
    }

    const std::string obj = ss.str();
    FILE* fp = fopen("/tmp/uv_unwrap_field.obj", "wb");
    if (!fp) return false;
    fwrite(obj.c_str(), 1, obj.size(), fp);
    fclose(fp);

    g_mesh = new Mesh();
    if (!g_mesh->read("/tmp/uv_unwrap_field.obj")) {
        delete g_mesh;
        g_mesh = nullptr;
        return false;
    }

    g_qc = new QuadCoverWasm(*g_mesh);
    g_miq = new MIQQuadWasm(*g_mesh);
    return true;
}

bool ensureQuadCoverInit()
{
    if (!g_qc) return false;
    if (g_qc_inited) return true;
    if (!g_qc->initMeshData()) return false;
    g_qc_inited = true;
    return true;
}

void snapshotFaceDirsFromQC()
{
    if (!g_qc) return;
    const auto& dirs = g_qc->getCrossField();
    const int nF = static_cast<int>(dirs.rows());
    g_face_dirs.resize(static_cast<size_t>(nF) * 3);
    for (int fi = 0; fi < nF; ++fi) {
        g_face_dirs[static_cast<size_t>(fi) * 3] = dirs(fi, 0);
        g_face_dirs[static_cast<size_t>(fi) * 3 + 1] = dirs(fi, 1);
        g_face_dirs[static_cast<size_t>(fi) * 3 + 2] = dirs(fi, 2);
    }
}

void snapshotMatchingFromQC()
{
    if (!g_qc) return;
    const int nE = static_cast<int>(g_qc->matching.size());
    g_matching.resize(static_cast<size_t>(nE));
    for (int ei = 0; ei < nE; ++ei) {
        g_matching[static_cast<size_t>(ei)] = static_cast<int>(g_qc->matching[ei]);
    }
}

void snapshotMIQFields()
{
    if (!g_miq) return;
    const int nF = g_miq->numFaces();
    const int nE = g_miq->numEdges();
    const Eigen::VectorXd& th = g_miq->faceTheta();
    const Eigen::VectorXi& jmp = g_miq->edgeJumps();

    g_face_theta.resize(static_cast<size_t>(nF));
    for (int fi = 0; fi < nF; ++fi) {
        g_face_theta[static_cast<size_t>(fi)] = th(fi);
    }

    g_edge_jumps.resize(static_cast<size_t>(nE));
    for (int ei = 0; ei < nE; ++ei) {
        g_edge_jumps[static_cast<size_t>(ei)] = jmp(ei);
    }

    g_last_miq_energy = g_miq->crossFieldEnergy();
}

void snapshotUV()
{
    if (!g_mesh) return;
    const int nV = static_cast<int>(g_mesh->vertices.size());
    g_uv.resize(static_cast<size_t>(nV) * 2);
    for (int vi = 0; vi < nV; ++vi) {
        g_uv[static_cast<size_t>(vi) * 2] = g_mesh->vertices[vi].uv[0];
        g_uv[static_cast<size_t>(vi) * 2 + 1] = g_mesh->vertices[vi].uv[1];
    }
}
bool hasValidUvSpread()
{
    if (!g_mesh || g_mesh->vertices.empty()) return false;
    double minU = std::numeric_limits<double>::infinity(), maxU = -minU;
    double minV = std::numeric_limits<double>::infinity(), maxV = -minV;
    for (const auto& v : g_mesh->vertices) {
        minU = std::min(minU, v.uv[0]); maxU = std::max(maxU, v.uv[0]);
        minV = std::min(minV, v.uv[1]); maxV = std::max(maxV, v.uv[1]);
    }
    return (maxU - minU) >= 1e-12 || (maxV - minV) >= 1e-12;
}


} // namespace

extern "C" {

EMSCRIPTEN_KEEPALIVE
int load_mesh(double* pos, int posLen, int* faces, int faceLen)
{
    return loadMeshFromArrays(pos, posLen, faces, faceLen) ? 0 : -1;
}

EMSCRIPTEN_KEEPALIVE
int step1_compute_principal_field()
{
    if (!ensureQuadCoverInit()) return -1;
    g_qc->estimatePrincipalCurvature();
    g_qc->computeFaceTheta();
    g_qc->syncFaceDirsFromTheta();
    snapshotFaceDirsFromQC();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int step2_smooth_and_matching(int smoothIters)
{
    if (!ensureQuadCoverInit()) return -1;
    g_qc->computeFaceTheta();
    g_qc->syncFaceDirsFromTheta();
    if (smoothIters > 0) g_qc->smoothCrossField(smoothIters);
    g_qc->computeMatching();
    snapshotFaceDirsFromQC();
    snapshotMatchingFromQC();
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int step3_solve_quadcover(int full)
{
    if (!ensureQuadCoverInit()) return -1;
    const double t0 = emscripten_get_now();
    bool ok = true;
    if (full) ok = g_qc->parameterizeFull();
    else {
        g_qc->parameterize();
        ok = true;
    }
    g_last_time_ms = emscripten_get_now() - t0;
    if (!ok) return -2;
    snapshotUV();
    return 0;
}

/** MIQ: alternating MIQP + coordinate-descent integer refinement + Poisson UV. */
EMSCRIPTEN_KEEPALIVE
int solve_miq(int crossIters, int jumpRefinePasses)
{
    if (!g_mesh || !g_miq) return -1;
    if (crossIters < 1) crossIters = 8;
    if (jumpRefinePasses < 0) jumpRefinePasses = 2;

    g_miq->setCrossFieldIterations(crossIters);
    g_miq->setJumpRefinePasses(jumpRefinePasses);

    const double t0 = emscripten_get_now();
    g_miq->parameterize();
    g_last_time_ms = emscripten_get_now() - t0;

    snapshotMIQFields();
    snapshotUV();
    return 0;
}

/** HOF (Holomorphic One-Form) 鈥?global conformal via 1-form integration. */
EMSCRIPTEN_KEEPALIVE
int solve_hof()
{
    if (!g_mesh) return -1;
    const double t0 = emscripten_get_now();
    HolomorphicOneForm hof(*g_mesh);
    hof.parameterize();
    g_last_time_ms = emscripten_get_now() - t0;
    snapshotUV();
    if (!hasValidUvSpread()) return -2;
    return 0;
}

EMSCRIPTEN_KEEPALIVE
double get_last_time_ms() { return g_last_time_ms; }

EMSCRIPTEN_KEEPALIVE
double get_miq_energy() { return g_last_miq_energy; }

EMSCRIPTEN_KEEPALIVE
double* get_face_dirs()
{
    return g_face_dirs.empty() ? nullptr : g_face_dirs.data();
}

EMSCRIPTEN_KEEPALIVE
int get_face_dirs_size() { return static_cast<int>(g_face_dirs.size()); }

EMSCRIPTEN_KEEPALIVE
int* get_matching()
{
    return g_matching.empty() ? nullptr : g_matching.data();
}

EMSCRIPTEN_KEEPALIVE
int get_matching_size() { return static_cast<int>(g_matching.size()); }

EMSCRIPTEN_KEEPALIVE
double* get_face_theta()
{
    return g_face_theta.empty() ? nullptr : g_face_theta.data();
}

EMSCRIPTEN_KEEPALIVE
int get_face_theta_size() { return static_cast<int>(g_face_theta.size()); }

EMSCRIPTEN_KEEPALIVE
int* get_edge_jumps()
{
    return g_edge_jumps.empty() ? nullptr : g_edge_jumps.data();
}

EMSCRIPTEN_KEEPALIVE
int get_edge_jumps_size() { return static_cast<int>(g_edge_jumps.size()); }

EMSCRIPTEN_KEEPALIVE
double* get_uv_result() { return g_uv.empty() ? nullptr : g_uv.data(); }

EMSCRIPTEN_KEEPALIVE
int get_uv_result_size() { return static_cast<int>(g_uv.size()); }

EMSCRIPTEN_KEEPALIVE
void dispose()
{
    g_face_dirs.clear();
    g_matching.clear();
    g_face_theta.clear();
    g_edge_jumps.clear();
    g_uv.clear();
    g_last_time_ms = 0.0;
    g_last_miq_energy = 0.0;
    g_qc_inited = false;
    if (g_qc) {
        delete g_qc;
        g_qc = nullptr;
    }
    if (g_miq) {
        delete g_miq;
        g_miq = nullptr;
    }
    if (g_mesh) {
        delete g_mesh;
        g_mesh = nullptr;
    }
}

} // extern "C"
