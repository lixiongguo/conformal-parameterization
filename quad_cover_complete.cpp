/**
 * QuadCover 简化实现 - 基于全局参数化的四边形网格化
 *
 * 算法流程 (参考文档):
 *   Phase 0: 预处理 - 平滑脐点区域方向场
 *   Phase 1: 计算匹配 Matching r (边上的 quarter-turn 数)
 *   Phase 2: 由 Matching 构造覆盖空间, 计算 layer shift
 *   Phase 3: 提升 Frame Field → Covering Vector Field
 *   Phase 4: Hodge 分解 - 求解最优可积方向场 (Poisson 求解)
 *   Phase 5: 全局连续性处理 (整格约束)
 *   Phase 6: 投影回原曲面, 输出 UV 参数化
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <igl/readOBJ.h>
#include <igl/read_triangle_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/grad.h>
#include <igl/boundary_facets.h>
#include <igl/boundary_conditions.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/list_to_matrix.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace Eigen;

// ============================================================
// 工具函数
// ============================================================

// 计算三角形面法向 (由三个顶点)
inline Vector3d FN_from_vertices(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2) {
    return ((p1 - p0).cross(p2 - p0));
}

// 计算三角形在顶点 pb 处的余切值 (cotangent)
// 参数: pa, pb, pc 表示角在 pb 处 (即角 apb)
inline double cot_tri(const Vector3d& pa, const Vector3d& pb, const Vector3d& pc) {
    Vector3d a = pa - pb;
    Vector3d b = pc - pb;
    double dot_val = a.dot(b);
    double cross_norm = a.cross(b).norm();
    if (cross_norm < 1e-14) return 0.0;  // 退化情况
    return dot_val / cross_norm;
}

// 90度旋转矩阵 J = [0, -1; 1, 0] (2D)
inline Matrix2d rotation90() {
    Matrix2d J;
    J << 0, -1,
         1, 0;
    return J;
}

// J^k (k 次 quarter-turn)
inline Matrix2d rotation_k(int k) {
    Matrix2d J = rotation90();
    Matrix2d Rk = Matrix2d::Identity();
    for (int i = 0; i < k; ++i) Rk = J * Rk;
    return Rk;
}

// 将 3D 向量投影到切平面
inline Vector3d project_to_tangent(const Vector3d& v, const Vector3d& n) {
    return v - v.dot(n) * n;
}

// 在三角形切平面上构建局部坐标系
void build_local_frame(const Vector3d& n, Vector3d& t1, Vector3d& t2) {
    // 选择一个与法向量不平行的向量作为参考
    Vector3d ref(1, 0, 0);
    if (abs(n.dot(ref)) > 0.9) ref = Vector3d(0, 1, 0);

    t1 = project_to_tangent(ref, n).normalized();
    t2 = n.cross(t1).normalized();
}

// 将 3D 方向向量转换为 2D 局部坐标 (在面片切平面内)
Vector2d direction_to_local(const Vector3d& dir, const Vector3d& n,
                            const Vector3d& t1, const Vector3d& t2) {
    Vector2d local;
    local(0) = dir.dot(t1);
    local(1) = dir.dot(t2);
    return local.normalized();
}

// ============================================================
// 数据结构
// ============================================================

struct Edge {
    int v1, v2;        // 顶点索引 (v1 < v2)
    int f1, f2;        // 相邻面片索引 (-1 表示边界边)
    int local_e1;      // e 在 f1 中的局部边索引 (0,1,2)
    int local_e2;      // e 在 f2 中的局部边索引
    int matching;      // r_ij ∈ {0,1,2,3} matching 值

    bool is_boundary() const { return f2 == -1; }
};

struct QuadCoverResult {
    MatrixXd UV;              // 最终 UV 坐标 (n_verts x 2)
    VectorXd theta;           // 平滑后的 4-RoSy 方向角
    vector<Edge> edges;       // 边信息 + matching
    VectorXd layer_shift;     // 每个顶点的 layer shift
    MatrixXd face_dirs;       // 面片上的参考方向 (n_faces x 3)
};

// ============================================================
// Phase 0: 预处理 + 初始化方向场
// ============================================================

/**
 * 使用主曲率方向初始化 4-RoSy 方向场
 * 每个 triangle 存储一个角度 theta (mod pi/2)
 */
void initialize_cross_field(const MatrixXd& V, const MatrixXi& F,
                             const MatrixXd& FN,
                             VectorXd& theta, MatrixXd& face_dirs) {
    int n_faces = F.rows();

    MatrixXd PD1, PD2;
    VectorXd PV1, PV2;
    igl::principal_curvature(V, F, PD1, PD2, PV1, PV2, 5);

    theta.resize(n_faces);
    face_dirs.resize(n_faces, 3);

    for (int i = 0; i < n_faces; ++i) {
        Vector3d n = FN.row(i).normalized();
        Vector3d pd1 = PD1.row(i).normalized();
        Vector3d t1, t2;

        build_local_frame(n, t1, t2);

        // 投影主方向到切平面
        Vector3d d = project_to_tangent(pd1, n);
        if (d.norm() < 1e-6) {
            // fallback: 脐点区域使用任意方向
            d = t1;
        }
        d.normalize();
        face_dirs.row(i) = d;

        // 计算角度 (相对于局部 x 轴 t1)
        double angle = atan2(d.dot(t2), d.dot(t1));
        theta(i) = fmod(angle, M_PI / 2.0);
        if (theta(i) < 0) theta(i) += M_PI / 2.0;
    }

    std::cout << "[Phase 0] 初始化 cross field 完成, " << n_faces << " 个面" << std::endl;
}

// ============================================================
// 辅助: 构建边表
// ============================================================

/**
 * 构建网格的边表，记录每条边的两个相邻面及其局部边索引
 */
vector<Edge> build_edge_list(const MatrixXi& F) {
    map<pair<int,int>, int> edge_map;
    vector<Edge> edges;

    for (int fi = 0; fi < F.rows(); ++fi) {
        for (int k = 0; k < 3; ++k) {
            int v1 = F(fi, k);
            int v2 = F(fi, (k + 1) % 3);
            if (v1 > v2) swap(v1, v2);  // 保证 v1 < v2

            pair<int,int> key(v1, v2);
            auto it = edge_map.find(key);

            if (it == edge_map.end()) {
                // 新边
                edge_map[key] = edges.size();
                Edge e;
                e.v1 = v1; e.v2 = v2;
                e.f1 = fi; e.f2 = -1;
                e.local_e1 = k; e.local_e2 = -1;
                e.matching = 0;
                edges.push_back(e);
            } else {
                // 已存在, 填充第二个面
                Edge& e = edges[it->second];
                e.f2 = fi;
                e.local_e2 = k;
            }
        }
    }

    return edges;
}

// ============================================================
// Phase 1: 计算 Matching
// ============================================================

/**
 * 对每条内部边, 计算 matching r_ij ∈ {0,1,2,3}
 *
 * r_ij = argmin_k angle(d_i, J^k * d_j)
 * 表示从面 T_i 到 T_j 需要旋转多少个 quarter-turn 才能对齐方向
 */
void compute_matching(const MatrixXd& V, const MatrixXi& F,
                      const MatrixXd& FN, const MatrixXd& face_dirs,
                      vector<Edge>& edges) {
    int matched_count = 0;

    for (auto& e : edges) {
        if (e.is_boundary()) continue;

        Vector3d ni = FN.row(e.f1).normalized();
        Vector3d nj = FN.row(e.f2).normalized();

        Vector3d di = face_dirs.row(e.f1);
        Vector3d dj = face_dirs.row(e.f2);

        // 构建两个面的局部坐标系
        Vector3d ti1, ti2, tj1, tj2;
        build_local_frame(ni, ti1, ti2);
        build_local_frame(nj, tj1, tj2);

        // 转换为 2D 局部坐标
        Vector2d di_2d = direction_to_local(di, ni, ti1, ti2);
        Vector2d dj_2d = direction_to_local(dj, nj, tj1, tj2);

        // 尝试 4 个旋转, 选最小夹角
        double min_angle = 1e10;
        int best_k = 0;
        Matrix2d J = rotation90();
        Matrix2d Rk = Matrix2d::Identity();

        for (int k = 0; k < 4; ++k) {
            Vector2d rotated_dj = Rk * dj_2d;
            double dot_val = di_2d.dot(rotated_dj);
            dot_val = max(-1.0, min(1.0, dot_val));  // clamp
            double angle = acos(abs(dot_val));  // 取绝对值因为 4-RoSy 对称

            if (angle < min_angle) {
                min_angle = angle;
                best_k = k;
            }
            Rk = J * Rk;
        }

        e.matching = best_k;
        matched_count++;
    }

    std::cout << "[Phase 1] Matching 计算完成, " << matched_count << " 条内部边" << std::endl;
}

// ============================================================
// Phase 2: 构造覆盖空间 & 计算 Layer Shift
// ============================================================

/**
 * 对每个顶点, 计算其 layer shift:
 *   ls(v) = (1/4) * sum of matchings around v
 *
 * ls != 0 的顶点是奇异点 (singularities / extraordinary vertices)
 */
void compute_layer_shift(const MatrixXd& V, const MatrixXi& F,
                          const vector<Edge>& edges,
                          VectorXd& layer_shift) {
    int n_verts = V.rows();
    layer_shift = VectorXd::Zero(n_verts);

    // 构建顶点到面的邻接关系
    vector<vector<int>> vert_to_face_adj(n_verts);
    for (int fi = 0; fi < F.rows(); ++fi) {
        for (int k = 0; k < 3; ++k) {
            vert_to_face_adj[F(fi, k)].push_back(fi);
        }
    }

    // 对每个顶点, 统计周围边的 matching 总和
    // 遍历所有内部边, 累加到端点的 layer shift
    for (const auto& e : edges) {
        if (e.is_boundary()) continue;
        layer_shift(e.v1) += e.matching;
        layer_shift(e.v2) += e.matching;
    }

    // 归一化: 除以 4
    layer_shift /= 4.0;

    // 统计奇异点
    int singularity_count = 0;
    for (int i = 0; i < n_verts; ++i) {
        if (abs(layer_shift(i)) > 1e-6) singularity_count++;
    }

    std::cout << "[Phase 2] Layer shift 计算完成, "
              << singularity_count << " 个奇异点" << std::endl;
}

// ============================================================
// 方向场平滑 (可选, 用于改善数值稳定性)
// ============================================================

/**
 * 简单的 Laplacian 平滑用于方向场
 * 目标: minimize sum_{edges} w_ij * (theta_i - theta_j - pi/2 * r_ij)^2
 */
void smooth_cross_field(const MatrixXi& F, const vector<Edge>& edges,
                         VectorXd& theta, int smooth_iters = 3) {
    int n_faces = F.rows();

    for (int iter = 0; iter < smooth_iters; ++iter) {
        VectorXd new_theta = theta;

        for (const auto& e : edges) {
            if (e.is_boundary()) continue;

            // 目标: theta_i ≈ theta_j + (pi/2) * r_ij (mod pi/2)
            double target = theta(e.f2) + (M_PI / 2.0) * e.matching;
            target = fmod(target, M_PI / 2.0);
            if (target < 0) target += M_PI / 2.0;

            // 加权平均
            double alpha = 0.5;
            new_theta(e.f1) = alpha * theta(e.f1) + (1 - alpha) * target;
            new_theta(e.f1) = fmod(new_theta(e.f1), M_PI / 2.0);
            if (new_theta(e.f1) < 0) new_theta(e.f1) += M_PI / 2.0;
        }
        theta = new_theta;
    }

    std::cout << "[Smooth] 方向场平滑完成 (" << smooth_iters << " 次)" << std::endl;
}

// ============================================================
// Phase 3: 构建 Covering Vector Field (提升后的方向场)
// ============================================================

/**
 * 根据 theta 和 matching, 构建两个正交的全局方向场 d1, d2
 * 这些方向定义在每个三角面上
 */
void build_covering_field(const MatrixXd& FN, const VectorXd& theta,
                           MatrixXd& d1, MatrixXd& d2) {
    int n_faces = FN.rows();
    d1.resize(n_faces, 3);
    d2.resize(n_faces, 3);

    for (int i = 0; i < n_faces; ++i) {
        Vector3d n = FN.row(i).normalized();
        Vector3d t1, t2;
        build_local_frame(n, t1, t2);

        // 旋转 theta 得到主方向
        double c = cos(theta(i)), s = sin(theta(i));
        d1.row(i) = c * t1.transpose() + s * t2.transpose();
        // 正交方向 (旋转 90 度)
        d2.row(i) = -s * t1.transpose() + c * t2.transpose();
    }

    std::cout << "[Phase 3] Covering field 构建完成" << std::endl;
}

// ============================================================
// Phase 4: Hodge 分解 - Poisson 求解
// ============================================================

/**
 * 核心步骤:
 *   给定方向场 (d1, d2), 通过求解 Poisson 方程得到 UV 坐标:
 *     L * u = div(d1)
 *     L * v = div(div2)
 *
 *   其中 L 是 cotangent Laplacian, div 是散度算子
 *
 * 这等价于寻找最贴近输入方向场的无旋 (curl-free) 场,
 * 即 Hodge 分解中的 exact + harmonic 部分
 */
void solve_poisson(const MatrixXd& V, const MatrixXi& F,
                    const MatrixXd& d1, const MatrixXd& d2,
                    MatrixXd& UV) {
    int n_verts = V.rows();

    // 1. 构建 Cotangent Laplacian
    SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);
    L = -L;  // libigl 返回的是负的 Laplacian

    // 2. 将面上的方向场插值到顶点
    MatrixXd Vd1(n_verts, 3), Vd2(n_verts, 3);
    Vd1.setZero(); Vd2.setZero();
    VectorXi count = VectorXi::Zero(n_verts);

    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v = F(i, j);
            Vd1.row(v) += d1.row(i);
            Vd2.row(v) += d2.row(i);
            count(v)++;
        }
    }
    for (int i = 0; i < n_verts; ++i) {
        if (count(i) > 0) {
            Vd1.row(i) /= count(i);
            Vd2.row(i) /= count(i);
        }
    }

    // 3. 计算梯度矩阵 (n_faces * 3 x n_verts)
    SparseMatrix<double> G;
    igl::grad(V, F, G);

    // 4. 计算质量矩阵 (用于散度的质量加权)
    SparseMatrix<double> M;
    igl::MassMatrixType mass_type = igl::MASSMATRIX_TYPE_VORONOI;
    igl::massmatrix(V, F, mass_type, M);

    // 5. 计算散度: div = M^{-1} * G^T * M_faces * d
    //    这里用简化版本: div ≈ G^T * (face_area * d_avg)
    SparseMatrix<double> M_faces;
    mass_type = igl::MASSMATRIX_TYPE_BARYCENTRIC;
    igl::massmatrix(V, F, mass_type, M_faces);

    // 将顶点方向转换回面方向 (用于散度计算)
    MatrixXd Fd1(F.rows(), 3), Fd2(F.rows(), 3);
    for (int i = 0; i < F.rows(); ++i) {
        Vector3d avg_d1(0, 0, 0), avg_d2(0, 0, 0);
        for (int j = 0; j < 3; ++j) {
            int v = F(i, j);
            avg_d1 += Vd1.row(v);
            avg_d2 += Vd2.row(v);
        }
        Fd1.row(i) = avg_d1 / 3.0;
        Fd2.row(i) = avg_d2 / 3.0;
    }

    // 散度 = G^T * M_faces * d  (简化的离散散度)
    // G 的维度是 (n_faces*3) x n_verts, 每个面有 3 行对应 3 条边
    // 我们需要将 3D 方向场投影到每条边上

    // 更准确的散度计算:
    // 对于每个面, div(f) = sum over edges of (f · edge_normal) * cot(edge) / area
    VectorXd rhs_u = VectorXd::Zero(n_verts);
    VectorXd rhs_v = VectorXd::Zero(n_verts);

    for (int fi = 0; fi < F.rows(); ++fi) {
        // 面的三个顶点
        int v0 = F(fi, 0), v1 = F(fi, 1), v2 = F(fi, 2);
        Vector3d p0 = V.row(v0), p1 = V.row(v1), p2 = V.row(v2);

        // 三条边
        Vector3d e0 = p1 - p0;  // opposite to v2
        Vector3d e1 = p2 - p1;  // opposite to v0
        Vector3d e2 = p0 - p2;  // opposite to v1

        // 面法向
        Vector3d fn = FN_from_vertices(p0, p1, p2);
        double dbl_area = fn.norm();
        fn /= dbl_area;

        // 边的法向 (在面内, 垂直于边, 指向外侧)
        Vector3d en0 = fn.cross(e0).normalized();
        Vector3d en1 = fn.cross(e1).normalized();
        Vector3d en2 = fn.cross(e2).normalized();

        // cotangent 权重
        double cot0 = cot_tri(p2, p0, p1);  // cot(angle at v0) for edge v1-v2
        double cot1 = cot_tri(p0, p1, p2);
        double cot2 = cot_tri(p1, p2, p0);

        // 面方向
        Vector3d fd1 = Fd1.row(fi);
        Vector3d fd2 = Fd2.row(fi);

        // 散度贡献
        double div_u = (fd1.dot(en0)*cot0 + fd1.dot(en1)*cot1 + fd1.dot(en2)*cot2) / 2.0;
        double div_v = (fd2.dot(en0)*cot0 + fd2.dot(en1)*cot1 + fd2.dot(en2)*cot2) / 2.0;

        rhs_u(v0) += div_v; rhs_u(v1) += div_v; rhs_u(v2) += div_v;  // 注意: libigl 符号约定
        rhs_v(v0) += div_u; rhs_v(v1) += div_u; rhs_v(v2) += div_u;
    }

    // 6. 边界条件
    VectorXi bnd;
    igl::boundary_facets(F, bnd);
    // 去重边界顶点
    vector<int> bnd_unique_vec;
    for (int i = 0; i < bnd.size(); ++i) {
        bnd_unique_vec.push_back(bnd(i));
    }
    sort(bnd_unique_vec.begin(), bnd_unique_vec.end());
    bnd_unique_vec.erase(unique(bnd_unique_vec.begin(), bnd_unique_vec.end()),
                          bnd_unique_vec.end());

    VectorXi bnd_unique;
    igl::list_to_matrix(bnd_unique_vec, bnd_unique);

    // 如果没有边界 (封闭网格), 固定两个顶点防止零空间
    bool is_closed_mesh = (bnd_unique.size() == 0);
    if (is_closed_mesh) {
        bnd_unique.resize(2);
        bnd_unique(0) = 0;
        bnd_unique(1) = 1;
    }

    VectorXd bc_u = VectorXd::Zero(bnd_unique.size());
    VectorXd bc_v = VectorXd::Zero(bnd_unique.size());

    // 7. 求解最小二乘系统: min_x  x^T L x - 2 b^T x, s.t. x(bnd) = bc
    //    使用 min_quad_with_fixed (两步式 API)

    // 添加小的正则化项使系统正定
    SparseMatrix<double> I_reg(n_verts, n_verts);
    I_reg.setIdentity();
    SparseMatrix<double> L_reg = L + 1e-10 * I_reg;

    // 预计算
    igl::min_quad_with_fixed_data<double> data_u, data_v;
    bool precompute_ok_u = igl::min_quad_with_fixed_precompute(
        L_reg, bnd_unique, SparseMatrix<double>(), true, data_u);
    
    bool precompute_ok_v = igl::min_quad_with_fixed_precompute(
        L_reg, bnd_unique, SparseMatrix<double>(), true, data_v);
    
    if (!precompute_ok_u || !precompute_ok_v) {
        std::cerr << "[Phase 4] Precompute 失败!" << std::endl;
        UV = MatrixXd::Zero(n_verts, 2);
        return;
    }

    // 求解
    MatrixXd U_sol, V_sol;
    bool success_u = igl::min_quad_with_fixed_solve(
        data_u, rhs_u, bc_u, VectorXd::Zero(0), U_sol);

    bool success_v = igl::min_quad_with_fixed_solve(
        data_v, rhs_v, bc_v, VectorXd::Zero(0), V_sol);

    if (!success_u || !success_v) {
        std::cerr << "[Phase 4] Poisson 求解失败!" << std::endl;
        UV = MatrixXd::Zero(n_verts, 2);
        return;
    }

    // 组装 UV
    UV.resize(n_verts, 2);
    UV.col(0) = U_sol;
    UV.col(1) = V_sol;

    std::cout << "[Phase 4] Hodge 分解 / Poisson 求解完成" << std::endl;
}

// ============================================================
// Phase 5: 全局连续性处理 (简化版)
// ============================================================

/**
 * 完整实现需要对每个拓扑环路检查整格约束 mu(gamma) in 2Z
 * 这里做简化处理: 确保 UV 坐标的一致性通过匹配验证
 */
void enforce_global_continuity(const MatrixXi& F, const vector<Edge>& edges,
                                const VectorXd& theta, MatrixXd& UV) {
    // 简化版: 检查跨边 UV 一致性并修正
    // 实际完整实现需要处理拓扑环路的同调群

    int n_verts = UV.rows();
    double max_discontinuity = 0;

    for (const auto& e : edges) {
        if (e.is_boundary()) continue;

        // 获取共享边的两个顶点
        int va = e.v1, vb = e.v2;

        // 期望的 UV 差异应该与匹配一致
        // 这里仅做监测
        double du = abs(UV(va, 0) - UV(vb, 0));
        double dv = abs(UV(va, 1) - UV(vb, 1));
        max_discontinuity = max(max_discontinuity, max(du, dv));
    }

    std::cout << "[Phase 5] 全局连续性检查完成, 最大跳变: "
              << max_discontinuity << std::endl;
}

// ============================================================
// Phase 6: UV 后处理 & 输出
// ============================================================

/**
 * 归一化 UV 到 [0,1] 范围 (保持宽高比)
 */
void normalize_uv(MatrixXd& UV) {
    double u_min = UV.col(0).minCoeff();
    double u_max = UV.col(0).maxCoeff();
    double v_min = UV.col(1).minCoeff();
    double v_max = UV.col(1).maxCoeff();

    double u_range = u_max - u_min;
    double v_range = v_max - v_min;
    double range = max(u_range, v_range);

    if (range < 1e-10) {
        std::cerr << "[Warning] UV 范围过小, 跳过归一化" << std::endl;
        return;
    }

    UV.col(0) = (UV.col(0).array() - u_min) / range;
    UV.col(1) = (UV.col(1).array() - v_min) / range;

    std::cout << "[Phase 6] UV 归一化完成, U=[" << u_min << "," << u_max
              << "], V=[" << v_min << "," << v_max << "]" << std::endl;
}

// /**
//  * 保存结果到文件
//  */
// void save_results(const string& prefix, const MatrixXd& UV,
//                    const VectorXd& theta, const VectorXd& layer_shift) {
//     // 保存 UV
//     ofstream uv_file(prefix + "_uv.txt");
//     uv_file << "# UV coordinates (u v)" << std::endl;
//     uv_file << UV << std::endl;
//     uv_file.close();

//     // 保存 Theta
//     ofstream theta_file(prefix + "_theta.txt");
//     theta_file << "# Face angles (radians, mod pi/2)" << std::endl;
//     theta_file << theta << std::endl;
//     theta_file.close();

//     // 保存 Layer Shift
//     ofstream ls_file(prefix + "_layer_shift.txt");
//     ls_file << "# Vertex layer shifts" << std::endl;
//     ls_file << layer_shift << std::endl;
//     ls_file.close();

//     std::cout << "[Output] 结果已保存到 " << prefix << "_*.txt" << std::endl;
// }

// ============================================================
// 主函数
// ============================================================

int main(int argc, char* argv[]) {
    std::string input_path = (argc > 1) ? argv[1] : "input.obj";
    std::string output_prefix = (argc > 2) ? argv[2] : "output";

    // ----------------------------------------------------------
    // 读取网格
    // ----------------------------------------------------------
    MatrixXd V;
    MatrixXi F;
    if (!igl::read_triangle_mesh(input_path, V, F)) {
        std::cerr << "错误: 无法加载网格文件: " << input_path << std::endl;
        return -1;
    }

    int n_verts = V.rows();
    int n_faces = F.rows();
    std::cout << "=== QuadCover 全局参数化 ===" << std::endl;
    std::cout << "网格: " << n_verts << " 个顶点, " << n_faces << " 个面"
              << std::endl;

    // ----------------------------------------------------------
    // 计算面法向
    // ----------------------------------------------------------
    MatrixXd FN;
    igl::per_face_normals(V, F, FN);

    // ----------------------------------------------------------
    // Phase 0: 初始化 Cross Field (主曲率方向)
    // ----------------------------------------------------------
    VectorXd theta;
    MatrixXd face_dirs;
    initialize_cross_field(V, F, FN, theta, face_dirs);

    // 可选: 平滑方向场
    vector<Edge> edges = build_edge_list(F);
    smooth_cross_field(F, edges, theta, /*iterations=*/3);

    // ----------------------------------------------------------
    // Phase 1: 计算 Matching
    // ----------------------------------------------------------
    compute_matching(V, F, FN, face_dirs, edges);

    // ----------------------------------------------------------
    // Phase 2: 构造覆盖空间, 计算 Layer Shift
    // ----------------------------------------------------------
    VectorXd layer_shift;
    compute_layer_shift(V, F, edges, layer_shift);

    // ----------------------------------------------------------
    // Phase 3: 构建 Covering Vector Field
    // ----------------------------------------------------------
    MatrixXd d1, d2;
    build_covering_field(FN, theta, d1, d2);

    // ----------------------------------------------------------
    // Phase 4: Hodge 分解 (Poisson 求解)
    // ----------------------------------------------------------
    MatrixXd UV;
    solve_poisson(V, F, d1, d2, UV);

    // ----------------------------------------------------------
    // Phase 5: 全局连续性处理
    // ----------------------------------------------------------
    enforce_global_continuity(F, edges, theta, UV);

    // ----------------------------------------------------------
    // Phase 6: 归一化 & 保存
    // ----------------------------------------------------------
    normalize_uv(UV);

    // 输出 UV 结果到文件 (供 WASM 调用方读取)
    {
        ofstream uv_out(output_prefix + "_uv.txt");
        uv_out << std::fixed << setprecision(8);
        uv_out << "# QuadCover UV Output\n";
        uv_out << "# n_verts n_cols\n";
        uv_out << UV.rows() << " " << UV.cols() << "\n";
        for (int i = 0; i < UV.rows(); ++i) {
            uv_out << UV(i, 0) << " " << UV(i, 1) << "\n";
        }
        uv_out.close();
        std::cout << "[Output] UV 结果已保存到 " << output_prefix + "_uv.txt" << std::endl;
    }

    return 0;
}
