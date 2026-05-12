/**
 * gauss_curvature_wasm.cpp — 调用 libigl gaussian_curvature，编译为 WASM
 *
 * API:
 *   compute_gauss_curvature(V_ptr, V_rows, F_ptr, F_rows) → int (0=ok, -1=err)
 *   get_gc_result()            → double* (V_rows 个 curvature 值)
 *   get_gc_result_size()       → int
 *   gc_dispose()               → void
 */
#include <emscripten.h>
#include <Eigen/Dense>
#include <igl/gaussian_curvature.h>
#include <vector>

static std::vector<double> g_result;

extern "C" {

EMSCRIPTEN_KEEPALIVE
int compute_gauss_curvature(double* V_ptr, int V_rows, int* F_ptr, int F_rows) {
    if (!V_ptr || !F_ptr || V_rows < 3 || F_rows < 1) return -1;

    try {
        Eigen::Map<Eigen::MatrixXd> V(V_ptr, V_rows, 3);
        Eigen::Map<Eigen::MatrixXi> F(F_ptr, F_rows, 3);

        Eigen::VectorXd K;
        igl::gaussian_curvature(V, F, K);

        g_result.resize(V_rows);
        for (int i = 0; i < V_rows; i++)
            g_result[i] = K(i);

        return 0;
    } catch (...) { return -1; }
}

EMSCRIPTEN_KEEPALIVE double* get_gc_result() { return g_result.empty() ? nullptr : g_result.data(); }
EMSCRIPTEN_KEEPALIVE int get_gc_result_size() { return (int)g_result.size(); }
EMSCRIPTEN_KEEPALIVE void gc_dispose() { g_result.clear(); g_result.shrink_to_fit(); }

}
