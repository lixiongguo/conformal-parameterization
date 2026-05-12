/**
 * principal_curvature_wasm.cpp — 调用 libigl principal_curvature，编译为 WASM
 *
 * API:
 *   compute_principal_curvature(V_ptr, V_rows, F_ptr, F_rows, radius) → int
 *   get_pc_pd1()  → double* (V_rows*3, max direction per vertex)
 *   get_pc_pd2()  → double* (V_rows*3, min direction per vertex)
 *   get_pc_pv1()  → double* (V_rows, max curvature value)
 *   get_pc_pv2()  → double* (V_rows, min curvature value)
 *   get_pc_num_verts() → int
 *   pc_dispose()  → void
 */
#include <emscripten.h>
#include <Eigen/Dense>
#include <igl/principal_curvature.h>
#include <igl/average_onto_faces.h>
#include <vector>

static std::vector<double> g_pd1, g_pd2, g_pv1, g_pv2;
static int g_num_verts = 0;

extern "C" {

EMSCRIPTEN_KEEPALIVE
int compute_principal_curvature(double* V_ptr, int V_rows, int* F_ptr, int F_rows, int radius) {
    if (!V_ptr || !F_ptr || V_rows < 3 || F_rows < 1) return -1;
    if (radius < 1) radius = 5;

    try {
        Eigen::Map<Eigen::MatrixXd> V(V_ptr, V_rows, 3);
        Eigen::Map<Eigen::MatrixXi> F(F_ptr, F_rows, 3);

        Eigen::MatrixXd PD1, PD2;
        Eigen::VectorXd PV1, PV2;
        igl::principal_curvature(V, F, PD1, PD2, PV1, PV2, radius, true);
        // Average onto faces for better visualization
        Eigen::MatrixXd PD1f, PD2f;
        Eigen::VectorXd PV1f, PV2f;
        igl::average_onto_faces(F, PD1, PD1f);
        igl::average_onto_faces(F, PD2, PD2f);
        igl::average_onto_faces(F, PV1, PV1f);
        igl::average_onto_faces(F, PV2, PV2f);

        g_num_verts = V_rows;
        g_pd1.resize(V_rows * 3);
        g_pd2.resize(V_rows * 3);
        g_pv1.resize(V_rows);
        g_pv2.resize(V_rows);
        for (int i = 0; i < V_rows; i++) {
            for (int j = 0; j < 3; j++) {
                g_pd1[i*3+j] = PD1(i,j);
                g_pd2[i*3+j] = PD2(i,j);
            }
            g_pv1[i] = PV1(i);
            g_pv2[i] = PV2(i);
        }
        return 0;
    } catch (...) { return -1; }
}

EMSCRIPTEN_KEEPALIVE double* get_pc_pd1() { return g_pd1.empty() ? nullptr : g_pd1.data(); }
EMSCRIPTEN_KEEPALIVE double* get_pc_pd2() { return g_pd2.empty() ? nullptr : g_pd2.data(); }
EMSCRIPTEN_KEEPALIVE double* get_pc_pv1() { return g_pv1.empty() ? nullptr : g_pv1.data(); }
EMSCRIPTEN_KEEPALIVE double* get_pc_pv2() { return g_pv2.empty() ? nullptr : g_pv2.data(); }
EMSCRIPTEN_KEEPALIVE int get_pc_num_verts() { return g_num_verts; }
EMSCRIPTEN_KEEPALIVE void pc_dispose() { g_pd1.clear(); g_pd2.clear(); g_pv1.clear(); g_pv2.clear(); g_num_verts = 0; }

}
