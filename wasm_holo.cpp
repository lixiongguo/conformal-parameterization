#include <emscripten.h>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <cstdio>
#include <limits>
#include "Mesh.h"
#include "MeshIO.h"
#include "HolomorphicOneForm.h"

static Mesh* g_mesh_ptr = nullptr;
static double g_holo_time_ms = 0.0;
static std::vector<double> g_holo_uv;

static bool load_mesh(const double* p, int pl, const int* f, int fl) {
    if (g_mesh_ptr) { delete g_mesh_ptr; g_mesh_ptr = nullptr; }
    if (pl < 9 || fl < 3) return false;
    g_mesh_ptr = new Mesh();
    size_t nV = pl / 3, nF = fl / 3;
    std::stringstream ss;
    ss << "# Holo\n# " << nV << " v " << nF << " f\n\n";
    for (size_t i = 0; i < nV; i++)
        ss << "v " << p[i*3] << " " << p[i*3+1] << " " << p[i*3+2] << "\n";
    ss << "\n";
    for (size_t i = 0; i < nF; i++)
        ss << "f " << f[i*3]+1 << " " << f[i*3+1]+1 << " " << f[i*3+2]+1 << "\n";
    std::string s = ss.str();
    FILE* fp = fopen("/tmp/input_holo.obj", "wb");
    if (!fp) return false;
    fwrite(s.c_str(), 1, s.size(), fp);
    fclose(fp);
    if (!g_mesh_ptr->read("/tmp/input_holo.obj")) {
        delete g_mesh_ptr; g_mesh_ptr = nullptr; return false;
    }
    return true;
}

extern "C" {

EMSCRIPTEN_KEEPALIVE
int solve_holo(double* pp, int pl, int* fp, int fl) {
    if (!pp || !fp || pl < 9 || fl < 3) return -1;
    if (!load_mesh(pp, pl, fp, fl)) return -1;
    auto t0 = std::chrono::high_resolution_clock::now();
    HolomorphicOneForm hf(*g_mesh_ptr);
    hf.parameterize();
    auto t1 = std::chrono::high_resolution_clock::now();
    g_holo_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    g_holo_uv.clear();
    g_holo_uv.reserve(g_mesh_ptr->vertices.size() * 2);
    double mu = std::numeric_limits<double>::infinity(), Mu = -mu, mv = mu, Mv = -mu;
    for (VertexCIter v = g_mesh_ptr->vertices.begin(); v != g_mesh_ptr->vertices.end(); ++v) {
        double u = v->uv.x(), w = v->uv.y();
        if (u < mu) mu = u; if (u > Mu) Mu = u;
        if (w < mv) mv = w; if (w > Mv) Mv = w;
    }
    double ru = Mu - mu, rv = Mv - mv;
    if (ru < 1e-10) ru = 1.0; if (rv < 1e-10) rv = 1.0;
    for (VertexCIter v = g_mesh_ptr->vertices.begin(); v != g_mesh_ptr->vertices.end(); ++v) {
        g_holo_uv.push_back((v->uv.x() - mu) / ru);
        g_holo_uv.push_back((v->uv.y() - mv) / rv);
    }
    return 0;
}

EMSCRIPTEN_KEEPALIVE double* get_holo_uv_result() { return g_holo_uv.empty() ? nullptr : g_holo_uv.data(); }
EMSCRIPTEN_KEEPALIVE int get_holo_uv_result_size() { return (int)g_holo_uv.size(); }
EMSCRIPTEN_KEEPALIVE double get_holo_last_time_ms() { return g_holo_time_ms; }
EMSCRIPTEN_KEEPALIVE void holo_dispose() {
    if (g_mesh_ptr) { delete g_mesh_ptr; g_mesh_ptr = nullptr; }
    g_holo_uv.clear(); g_holo_uv.shrink_to_fit();
}

}
