/**
 * gauss_curvature_wasm.cpp — 调用 libigl gaussian_curvature，编译为 WASM
 *
 * API:
 *   compute_gauss_curvature(V_ptr, V_rows, F_ptr, F_rows) → int (0=ok, -1=err)
 *   get_gc_result()            → double* (V_rows 个按顶点面积归一化后的 curvature 值)
 *   get_gc_result_size()       → int
 *   gc_dispose()               → void
 */
#include <emscripten.h>
#include <Eigen/Dense>
#include <igl/gaussian_curvature.h>
#include <cmath>
#include <vector>

static std::vector<double> g_result;

extern "C" {

EMSCRIPTEN_KEEPALIVE
int compute_gauss_curvature(double* V_ptr, int V_rows, int* F_ptr, int F_rows) {
    if (!V_ptr || !F_ptr || V_rows < 3 || F_rows < 1) return -1;

    try {
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>> V_map(V_ptr, V_rows, 3);
        Eigen::Map<const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>> F_map(F_ptr, F_rows, 3);
        Eigen::MatrixXd V = V_map;
        Eigen::MatrixXi F = F_map;

        Eigen::VectorXd K;
        igl::gaussian_curvature(V, F, K);

        std::vector<double> vertex_area(V_rows, 0.0);
        for (int i = 0; i < F_rows; i++) {
            const int a = F(i, 0), b = F(i, 1), c = F(i, 2);
            if (a < 0 || a >= V_rows || b < 0 || b >= V_rows || c < 0 || c >= V_rows) return -1;
            const Eigen::Vector3d ab = V.row(b) - V.row(a);
            const Eigen::Vector3d ac = V.row(c) - V.row(a);
            const double area = 0.5 * ab.cross(ac).norm();
            if (!std::isfinite(area) || area <= 0.0) continue;
            const double one_third = area / 3.0;
            vertex_area[a] += one_third;
            vertex_area[b] += one_third;
            vertex_area[c] += one_third;
        }

        g_result.resize(V_rows);
        for (int i = 0; i < V_rows; i++) {
            g_result[i] = vertex_area[i] > 1e-16 ? K(i) / vertex_area[i] : 0.0;
        }

        return 0;
    } catch (...) { return -1; }
}

EMSCRIPTEN_KEEPALIVE double* get_gc_result() { return g_result.empty() ? nullptr : g_result.data(); }
EMSCRIPTEN_KEEPALIVE int get_gc_result_size() { return (int)g_result.size(); }
EMSCRIPTEN_KEEPALIVE void gc_dispose() { g_result.clear(); g_result.shrink_to_fit(); }

}
