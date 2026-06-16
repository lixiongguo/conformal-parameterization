/**
 * Bounded distortion mesh deformation WASM (ccall API).
 * Wraps bounded_distortion::solveBoundedDistortionMap for interactive handle dragging.
 */

#include "BoundedDistortionMapping.hpp"
#include <vector>
#include <cmath>
#include <emscripten.h>

namespace {

Eigen::MatrixXd g_rest;
Eigen::MatrixXi g_faces;
Eigen::MatrixXd g_uv;
bounded_distortion::Options g_opt;
bool g_ready = false;

double g_max_distortion = 1.0;
double g_min_jacobian = 1.0;
int g_flip_count = 0;
int g_last_iterations = 0;
double g_last_time_ms = 0.0;
int g_bound_on = 1;

void clearState() {
    g_rest.resize(0, 0);
    g_faces.resize(0, 0);
    g_uv.resize(0, 0);
    g_ready = false;
    g_max_distortion = 1.0;
    g_min_jacobian = 1.0;
    g_flip_count = 0;
    g_last_iterations = 0;
    g_last_time_ms = 0.0;
}

void applyBoundOptions() {
    if (!g_bound_on) {
        g_opt.distortion_penalty = 0.0;
        g_opt.positivity_penalty = 0.0;
    }
}

void updateStats(const bounded_distortion::SolveResult& result) {
    g_max_distortion = result.max_distortion;
    g_min_jacobian = result.min_jacobian;
    g_flip_count = 0;
    for (const auto& f : result.faces) {
        if (f.jacobian <= 0.0) ++g_flip_count;
    }
    g_last_iterations = result.iterations;
}

} // namespace

extern "C" {

EMSCRIPTEN_KEEPALIVE
void dispose() { clearState(); }

EMSCRIPTEN_KEEPALIVE
int setup_mesh(double* pos, int posLen, int* faces, int faceLen) {
    clearState();
    if (!pos || !faces) return -1;
    if (posLen < 6 || faceLen < 3) return -1;
    if (posLen % 2 != 0 || faceLen % 3 != 0) return -2;

    const int nV = posLen / 2;
    const int nF = faceLen / 3;
    g_rest.resize(nV, 2);
    for (int i = 0; i < nV; ++i) {
        g_rest(i, 0) = pos[2 * i];
        g_rest(i, 1) = pos[2 * i + 1];
    }
    g_faces.resize(nF, 3);
    for (int f = 0; f < nF; ++f) {
        g_faces(f, 0) = faces[3 * f];
        g_faces(f, 1) = faces[3 * f + 1];
        g_faces(f, 2) = faces[3 * f + 2];
    }
    g_uv = g_rest;
    g_ready = true;
    return 0;
}

EMSCRIPTEN_KEEPALIVE
void set_options(
    double distortion_bound,
    double min_alpha_real,
    double lscm_weight,
    double distortion_penalty,
    double positivity_penalty,
    double reference_weight,
    double smoothness_weight,
    double initial_step,
    int bound_on)
{
    g_bound_on = bound_on ? 1 : 0;
    g_opt.distortion_bound = distortion_bound;
    g_opt.min_alpha_real = min_alpha_real;
    g_opt.lscm_weight = lscm_weight;
    g_opt.distortion_penalty = distortion_penalty;
    g_opt.positivity_penalty = positivity_penalty;
    g_opt.reference_weight = reference_weight;
    g_opt.smoothness_weight = smoothness_weight;
    g_opt.initial_step = initial_step;
    g_opt.outer_iterations = 1;
    g_opt.inner_iterations = 1;
    applyBoundOptions();
}

EMSCRIPTEN_KEEPALIVE
void reset_uv() {
    if (!g_ready) return;
    g_uv = g_rest;
}

EMSCRIPTEN_KEEPALIVE
int optimize_step(int max_steps, int* anchorIdx, int anchorLen, double* anchorTargets) {
    if (!g_ready) return -1;
    if (max_steps < 1) max_steps = 1;

    std::vector<bounded_distortion::Anchor> anchors;
    if (anchorIdx && anchorTargets && anchorLen > 0) {
        anchors.reserve(static_cast<size_t>(anchorLen));
        for (int i = 0; i < anchorLen; ++i) {
            bounded_distortion::Anchor a;
            a.vertex = anchorIdx[i];
            a.target.x() = anchorTargets[2 * i];
            a.target.y() = anchorTargets[2 * i + 1];
            anchors.push_back(a);
        }
    }

    bounded_distortion::Options opts = g_opt;
    opts.outer_iterations = 1;
    opts.inner_iterations = max_steps;
    applyBoundOptions();
    if (!g_bound_on) {
        opts.distortion_penalty = 0.0;
        opts.positivity_penalty = 0.0;
    }

    const double t0 = emscripten_get_now();
    try {
        const bounded_distortion::SolveResult result =
            bounded_distortion::solveBoundedDistortionMap(g_rest, g_faces, g_uv, anchors, opts);
        g_uv = result.uv;
        updateStats(result);
        g_last_time_ms = emscripten_get_now() - t0;
        return 0;
    } catch (const std::exception&) {
        return -2;
    }
}

EMSCRIPTEN_KEEPALIVE
int get_vertex_count() { return g_ready ? static_cast<int>(g_rest.rows()) : 0; }

EMSCRIPTEN_KEEPALIVE
double* get_uv_ptr() { return g_ready ? g_uv.data() : nullptr; }

EMSCRIPTEN_KEEPALIVE
double get_max_distortion() { return g_max_distortion; }

EMSCRIPTEN_KEEPALIVE
double get_min_jacobian() { return g_min_jacobian; }

EMSCRIPTEN_KEEPALIVE
int get_flip_count() { return g_flip_count; }

EMSCRIPTEN_KEEPALIVE
int get_last_iterations() { return g_last_iterations; }

EMSCRIPTEN_KEEPALIVE
double get_last_time_ms() { return g_last_time_ms; }

} // extern "C"
