/**
 * Abel-Jacobi global parameterization WASM (blog section 4.5).
 *
 * Pipeline: load mesh -> build canonical 1-forms + period lattice -> quantize periods -> UV.
 * Closed mesh with genus >= 1 required.
 */

#include <emscripten.h>
#include <emscripten/html5.h>
#include <cstdio>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>

#include "Mesh.h"
#include "MeshIO.h"
#include "AbelJacobiParameterization.h"

namespace {

Mesh* g_mesh = nullptr;
AbelJacobiParameterization* g_param = nullptr;

double g_last_time_ms = 0.0;
std::vector<double> g_uv;
std::vector<double> g_lattice_flat;

int g_genus = 0;
int g_built = 0;
int g_poincare_ok = 1;
int g_abel_ok = 1;
double g_lattice_residual = 0.0;

bool isClosedMesh(const Mesh& m)
{
    if (!m.boundaries.empty()) return false;
    for (EdgeCIter e = m.edges.begin(); e != m.edges.end(); ++e) {
        if (e->isBoundary()) return false;
    }
    return true;
}

int computeGenusClosed(const Mesh& m)
{
    const int nV = static_cast<int>(m.vertices.size());
    const int nE = static_cast<int>(m.edges.size());
    int nF = 0;
    for (FaceCIter f = m.faces.begin(); f != m.faces.end(); ++f) {
        if (!f->isBoundary()) ++nF;
    }
    const int chi = nV - nE + nF;
    const int g = (2 - chi) / 2;
    return g < 0 ? 0 : g;
}

bool hasValidUvSpread(const Mesh& m)
{
    if (m.vertices.empty()) return false;
    double minU = std::numeric_limits<double>::infinity(), maxU = -minU;
    double minV = std::numeric_limits<double>::infinity(), maxV = -minV;
    for (const auto& v : m.vertices) {
        minU = std::min(minU, v.uv.x()); maxU = std::max(maxU, v.uv.x());
        minV = std::min(minV, v.uv.y()); maxV = std::max(maxV, v.uv.y());
    }
    return (maxU - minU) >= 1e-12 || (maxV - minV) >= 1e-12;
}

bool loadMeshFromArrays(const double* positions, int posLen, const int* faces, int faceLen)
{
    if (g_param) {
        delete g_param;
        g_param = nullptr;
    }
    if (g_mesh) {
        delete g_mesh;
        g_mesh = nullptr;
    }

    g_uv.clear();
    g_lattice_flat.clear();
    g_genus = 0;
    g_built = 0;
    g_poincare_ok = 1;
    g_abel_ok = 1;
    g_lattice_residual = 0.0;

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
    FILE* fp = fopen("/tmp/uv_unwrap_abel_jacobi.obj", "wb");
    if (!fp) return false;
    fwrite(obj.c_str(), 1, obj.size(), fp);
    fclose(fp);

    g_mesh = new Mesh();
    if (!g_mesh->read("/tmp/uv_unwrap_abel_jacobi.obj")) {
        delete g_mesh;
        g_mesh = nullptr;
        return false;
    }

    g_param = new AbelJacobiParameterization(*g_mesh);
    g_genus = isClosedMesh(*g_mesh) ? computeGenusClosed(*g_mesh) : 0;
    return true;
}

void snapshotLattice()
{
    g_lattice_flat.clear();
    if (!g_param || !g_param->engine().isReady()) return;

    const Eigen::MatrixXd& L = g_param->engine().latticeGenerators();
    g_lattice_flat.reserve(static_cast<size_t>(L.size()));
    for (int r = 0; r < L.rows(); ++r) {
        for (int c = 0; c < L.cols(); ++c) {
            g_lattice_flat.push_back(L(r, c));
        }
    }
}

void snapshotUV()
{
    g_uv.clear();
    if (!g_mesh) return;

    double minU = std::numeric_limits<double>::infinity();
    double maxU = -std::numeric_limits<double>::infinity();
    double minV = std::numeric_limits<double>::infinity();
    double maxV = -std::numeric_limits<double>::infinity();

    for (const auto& v : g_mesh->vertices) {
        minU = std::min(minU, v.uv.x());
        maxU = std::max(maxU, v.uv.x());
        minV = std::min(minV, v.uv.y());
        maxV = std::max(maxV, v.uv.y());
    }

    const double rU = std::max(maxU - minU, 1e-10);
    const double rV = std::max(maxV - minV, 1e-10);

    g_uv.reserve(g_mesh->vertices.size() * 2);
    for (const auto& v : g_mesh->vertices) {
        g_uv.push_back((v.uv.x() - minU) / rU);
        g_uv.push_back((v.uv.y() - minV) / rV);
    }
}

} // namespace

extern "C" {

EMSCRIPTEN_KEEPALIVE
int load_mesh(double* pos, int posLen, int* faces, int faceLen)
{
    return loadMeshFromArrays(pos, posLen, faces, faceLen) ? 0 : -1;
}

EMSCRIPTEN_KEEPALIVE
int is_closed_mesh()
{
    if (!g_mesh) return 0;
    return isClosedMesh(*g_mesh) ? 1 : 0;
}

EMSCRIPTEN_KEEPALIVE
int get_genus()
{
    return g_genus;
}

EMSCRIPTEN_KEEPALIVE
int build_abel_jacobi()
{
    if (!g_param) return -1;
    const double t0 = emscripten_get_now();
    const bool ok = g_param->engine().build();
    g_last_time_ms = emscripten_get_now() - t0;
    g_built = ok ? 1 : 0;
    if (ok) snapshotLattice();
    return ok ? 0 : -2;
}

EMSCRIPTEN_KEEPALIVE
int solve_abel_jacobi(
    double* pos,
    int posLen,
    int* faces,
    int faceLen,
    int baseVertex,
    int* singVerts,
    int* singValences,
    int* singDegrees,
    int nSing)
{
    if (!loadMeshFromArrays(pos, posLen, faces, faceLen)) return -1;
    if (!isClosedMesh(*g_mesh) || g_genus < 1) return -4;

    g_param->setBaseVertex(baseVertex);

    std::vector<AbelJacobi::SingularPoint> sings;
    if (singVerts && nSing > 0) {
        sings.reserve(static_cast<size_t>(nSing));
        for (int i = 0; i < nSing; ++i) {
            AbelJacobi::SingularPoint s;
            s.vertexIndex = singVerts[i];
            s.valence = (singValences && singValences[i] > 0) ? singValences[i] : 4;
            s.degree = singDegrees ? singDegrees[i] : (4 - s.valence);
            sings.push_back(s);
        }
        g_param->setSingularities(sings);
    }

    const double t0 = emscripten_get_now();

    const bool builtBefore = g_param->engine().isReady();
    if (!builtBefore) {
        g_built = g_param->engine().build() ? 1 : 0;
    } else {
        g_built = 1;
    }

    if (!g_built) return -2;

    if (!sings.empty()) {
        const auto check = g_param->checkSingularities();
        g_poincare_ok = check.poincareHopfOk ? 1 : 0;
        g_abel_ok = check.abelJacobiOk ? 1 : 0;
        g_lattice_residual = check.latticeResidual;
    }

    g_param->parameterize();

    g_last_time_ms = emscripten_get_now() - t0;

    if (g_param->engine().isReady()) snapshotLattice();
    snapshotUV();

    if (g_uv.empty()) return -3;
    if (!hasValidUvSpread(*g_mesh)) return -5;
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int get_built() { return g_built; }

EMSCRIPTEN_KEEPALIVE
int get_poincare_ok() { return g_poincare_ok; }

EMSCRIPTEN_KEEPALIVE
int get_abel_ok() { return g_abel_ok; }

EMSCRIPTEN_KEEPALIVE
double get_lattice_residual() { return g_lattice_residual; }

EMSCRIPTEN_KEEPALIVE
double* get_lattice_generators()
{
    return g_lattice_flat.empty() ? nullptr : g_lattice_flat.data();
}

EMSCRIPTEN_KEEPALIVE
int get_lattice_generators_size() { return static_cast<int>(g_lattice_flat.size()); }

EMSCRIPTEN_KEEPALIVE
int get_lattice_rows()
{
    if (!g_param || !g_param->engine().isReady()) return 0;
    return g_param->engine().latticeGenerators().rows();
}

EMSCRIPTEN_KEEPALIVE
int get_lattice_cols()
{
    if (!g_param || !g_param->engine().isReady()) return 0;
    return g_param->engine().latticeGenerators().cols();
}

EMSCRIPTEN_KEEPALIVE
double get_last_time_ms() { return g_last_time_ms; }

EMSCRIPTEN_KEEPALIVE
double* get_uv_result() { return g_uv.empty() ? nullptr : g_uv.data(); }

EMSCRIPTEN_KEEPALIVE
int get_uv_result_size() { return static_cast<int>(g_uv.size()); }

EMSCRIPTEN_KEEPALIVE
void dispose()
{
    g_uv.clear();
    g_lattice_flat.clear();
    g_last_time_ms = 0.0;
    g_genus = 0;
    g_built = 0;
    g_poincare_ok = 1;
    g_abel_ok = 1;
    g_lattice_residual = 0.0;
    if (g_param) {
        delete g_param;
        g_param = nullptr;
    }
    if (g_mesh) {
        delete g_mesh;
        g_mesh = nullptr;
    }
}

} // extern "C"
