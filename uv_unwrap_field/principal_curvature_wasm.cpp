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
#include <Eigen/SVD>

// The vendored libigl headers use a few Eigen 3.3+/3.4 names while this
// repository currently vendors Eigen 3.2.x under the eigen-3.4.0 directory.
// Keep the shim local to this WASM target so other builds stay untouched.
#if !EIGEN_VERSION_AT_LEAST(3, 3, 0)
namespace Eigen {
    typedef DenseIndex Index;
    static const int all = 0;

    template <typename MatrixType>
    class CompleteOrthogonalDecomposition {
    public:
        explicit CompleteOrthogonalDecomposition(const MatrixType& matrix) : matrix_(matrix) {}

        template <typename Rhs>
        auto solve(const Rhs& rhs) const {
            return matrix_.jacobiSvd(ComputeThinU | ComputeThinV).solve(rhs);
        }

    private:
        MatrixType matrix_;
    };
}
#endif

#include <igl/principal_curvature.h>
#include <vector>

static std::vector<double> g_pd1, g_pd2, g_pv1, g_pv2;
static int g_num_verts = 0;

extern "C" {

EMSCRIPTEN_KEEPALIVE
int compute_principal_curvature(double* V_ptr, int V_rows, int* F_ptr, int F_rows, int radius) {
    if (!V_ptr || !F_ptr || V_rows < 3 || F_rows < 1) return -1;
    if (radius < 1) radius = 5;

    try {
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>> V_map(V_ptr, V_rows, 3);
        Eigen::Map<const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>> F_map(F_ptr, F_rows, 3);
        Eigen::MatrixXd V = V_map;
        Eigen::MatrixXi F = F_map;

        Eigen::MatrixXd PD1, PD2;
        Eigen::VectorXd PV1, PV2;
        igl::principal_curvature(V, F, PD1, PD2, PV1, PV2, radius, true);

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
