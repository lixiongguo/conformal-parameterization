#include <emscripten/bind.h>
#include "Mesh.h"
#include "Lscm.h"
#include "BoundedDistortionMapping.hpp"
#include <Eigen/Dense>
#include <vector>

using namespace emscripten;

// ============================================================
//  Eigen ↔ std::vector 转换
// ============================================================
static Eigen::MatrixXd _vecToMatrixXd(const std::vector<double>& data, int rows, int cols) {
    Eigen::MatrixXd m(rows, cols);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
            m(r, c) = data[r * cols + c];
    return m;
}

static std::vector<double> _matrixXdToVec(const Eigen::MatrixXd& m) {
    std::vector<double> v;
    v.reserve(m.rows() * m.cols());
    for (int r = 0; r < m.rows(); ++r)
        for (int c = 0; c < m.cols(); ++c)
            v.push_back(m(r, c));
    return v;
}

static Eigen::MatrixXi _vecToMatrixXi(const std::vector<int>& data, int rows, int cols) {
    Eigen::MatrixXi m(rows, cols);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
            m(r, c) = data[r * cols + c];
    return m;
}

// ============================================================
//  Mesh helpers
// ============================================================
static std::vector<double> _vertexUvs(const Mesh& mesh) {
    std::vector<double> r;
    r.reserve(mesh.vertices.size() * 2);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        r.push_back(v->uv.x());
        r.push_back(v->uv.y());
    }
    return r;
}

static std::vector<double> _vertexPositions2D(const Mesh& mesh) {
    std::vector<double> r;
    r.reserve(mesh.vertices.size() * 2);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        r.push_back(v->position.x());
        r.push_back(v->position.y());
    }
    return r;
}

static std::vector<int> _faceIndices(const Mesh& mesh) {
    std::vector<int> r;
    r.reserve(mesh.faces.size() * 3);
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        HalfEdgeCIter he = f->he;
        for (int k = 0; k < 3; k++) {
            r.push_back(he->vertex->index);
            he = he->next;
        }
    }
    return r;
}

static void _setVertexUvs(Mesh& mesh, const std::vector<double>& uv_flat) {
    int idx = 0;
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->uv.x() = uv_flat[idx++];
        v->uv.y() = uv_flat[idx++];
    }
}

// ============================================================
//  BD-LSCM 管线结果
// ============================================================
struct BDLscmResult {
    std::vector<double> uv_flat;
    double max_distortion;
    double min_jacobian;
    double final_energy;
    int iterations;
    std::vector<double> face_distortions;
    std::vector<double> face_jacobians;
    std::vector<double> face_cones;
};

// 前向声明
static BDLscmResult _solveCore(const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                               const Eigen::MatrixXd&, const std::vector<bounded_distortion::Anchor>&,
                               const bounded_distortion::Options&);

// ============================================================
//  主求解函数: LSCM 初始 + 有界失真优化
//  使用 emscripten::val 直接接收 JS 数组
// ============================================================
static BDLscmResult _solveBDLscmFull(
    val vertices2D_js,
    val faces_js,
    val initial_uv_js,
    val anchor_vertices_js,
    val anchor_targets_js,
    double distortion_bound,
    int outer_iterations,
    int inner_iterations,
    double lscm_weight,
    double distortion_penalty,
    double positivity_penalty,
    double initial_step)
{
    // 提取 JS 数组到 std::vector
    auto vecDbl = [](val arr) {
        std::vector<double> v;
        int len = arr["length"].as<int>();
        v.reserve(len);
        for (int i = 0; i < len; ++i) v.push_back(arr[i].as<double>());
        return v;
    };
    auto vecInt = [](val arr) {
        std::vector<int> v;
        int len = arr["length"].as<int>();
        v.reserve(len);
        for (int i = 0; i < len; ++i) v.push_back(arr[i].as<int>());
        return v;
    };

    std::vector<double> v2D = vecDbl(vertices2D_js);
    std::vector<int>    fi  = vecInt(faces_js);
    std::vector<double> uv0 = vecDbl(initial_uv_js);
    std::vector<int>    av  = vecInt(anchor_vertices_js);
    std::vector<double> at  = vecDbl(anchor_targets_js);

    const int nV = static_cast<int>(v2D.size() / 2);
    const int nF = static_cast<int>(fi.size() / 3);

    // 构建 Eigen 矩阵
    Eigen::MatrixXd vertices = _vecToMatrixXd(v2D, nV, 2);
    Eigen::MatrixXi faces    = _vecToMatrixXi(fi, nF, 3);
    Eigen::MatrixXd init_uv  = _vecToMatrixXd(uv0, nV, 2);

    // 构建锚点
    std::vector<bounded_distortion::Anchor> anchors;
    anchors.reserve(av.size());
    for (size_t i = 0; i < av.size(); ++i) {
        bounded_distortion::Anchor a;
        a.vertex = av[i];
        a.target.x() = at[2 * i];
        a.target.y() = at[2 * i + 1];
        anchors.push_back(a);
    }

    // 配置选项
    bounded_distortion::Options opts;
    opts.distortion_bound    = distortion_bound;
    opts.outer_iterations    = outer_iterations;
    opts.inner_iterations    = inner_iterations;
    opts.lscm_weight         = lscm_weight;
    opts.distortion_penalty  = distortion_penalty;
    opts.positivity_penalty  = positivity_penalty;
    opts.initial_step        = initial_step;

    // 求解
    return _solveCore(vertices, faces, init_uv, anchors, opts);
}

// 核心求解（被 _solveBDLscmFull 和 _runOnMesh 共用）
static BDLscmResult _solveCore(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXd& init_uv,
    const std::vector<bounded_distortion::Anchor>& anchors,
    const bounded_distortion::Options& opts)
{
    bounded_distortion::SolveResult result =
        bounded_distortion::solveBoundedDistortionLscm(vertices, faces, init_uv, anchors, opts);

    BDLscmResult out;
    out.uv_flat         = _matrixXdToVec(result.uv);
    out.max_distortion  = result.max_distortion;
    out.min_jacobian    = result.min_jacobian;
    out.final_energy    = result.final_energy;
    out.iterations      = result.iterations;

    out.face_distortions.reserve(result.faces.size());
    out.face_jacobians.reserve(result.faces.size());
    out.face_cones.reserve(result.faces.size());
    for (const auto& s : result.faces) {
        out.face_distortions.push_back(s.distortion);
        out.face_jacobians.push_back(s.jacobian);
        out.face_cones.push_back(s.cone_violation);
    }
    return out;
}

// ============================================================
//  便捷函数: 在已有 LSCM 的 Mesh 上运行 BD-LSCM
// ============================================================
static BDLscmResult _runOnMesh(
    Mesh& mesh,
    const std::vector<int>& anchor_vertices,
    const std::vector<double>& anchor_targets,
    double distortion_bound)
{
    // 提取网格数据
    std::vector<double> uv_flat  = _vertexUvs(mesh);
    std::vector<double> pos2D    = _vertexPositions2D(mesh);
    std::vector<int>    faces    = _faceIndices(mesh);

    bounded_distortion::Options opts;
    opts.distortion_bound   = distortion_bound;
    opts.lscm_weight        = 1.0;
    opts.outer_iterations   = 8;
    opts.inner_iterations   = 400;
    opts.distortion_penalty = 5000.0;
    opts.positivity_penalty = 5000.0;
    opts.initial_step       = 1e-2;

    // 构建锚点
    std::vector<bounded_distortion::Anchor> anchors;
    anchors.reserve(anchor_vertices.size());
    for (size_t i = 0; i < anchor_vertices.size(); ++i) {
        bounded_distortion::Anchor a;
        a.vertex = anchor_vertices[i];
        a.target.x() = anchor_targets[2 * i];
        a.target.y() = anchor_targets[2 * i + 1];
        anchors.push_back(a);
    }

    // 构建 Eigen 矩阵
    int nV = (int)mesh.vertices.size(), nF = (int)mesh.faces.size();
    Eigen::MatrixXd V = _vecToMatrixXd(pos2D, nV, 2);
    Eigen::MatrixXi F = _vecToMatrixXi(faces, nF, 3);
    Eigen::MatrixXd UV0 = _vecToMatrixXd(uv_flat, nV, 2);

    BDLscmResult result = _solveCore(V, F, UV0, anchors, opts);

    // 回写 UV
    _setVertexUvs(mesh, result.uv_flat);

    return result;
}

// ============================================================
//  Embind 注册
// ============================================================
EMSCRIPTEN_BINDINGS(BD_LSCMModule) {

    // --- 注册 vector 类型 ---
    register_vector<double>("VectorDouble");
    register_vector<int>("VectorInt");

    // --- Mesh (self-contained) ---
    class_<Mesh>("Mesh")
        .constructor<>()
        .function("read",  &Mesh::read)
        .function("write", &Mesh::write)
        .function("getVertexUvs",      &_vertexUvs)
        .function("getVertexPositions2D", &_vertexPositions2D)
        .function("getFaceIndices",    &_faceIndices)
        .function("setVertexUvs",      &_setVertexUvs);

    // --- LSCM ---
    class_<Lscm, base<Parameterization>>("Lscm")
        .constructor<Mesh&>()
        .function("parameterize",  &Lscm::parameterize)
        .function("computeQcError", &Lscm::computeQcError);

    // --- BD-LSCM Result ---
    value_object<BDLscmResult>("BDLscmResult")
        .field("uv_flat",          &BDLscmResult::uv_flat)
        .field("max_distortion",   &BDLscmResult::max_distortion)
        .field("min_jacobian",     &BDLscmResult::min_jacobian)
        .field("final_energy",     &BDLscmResult::final_energy)
        .field("iterations",       &BDLscmResult::iterations)
        .field("face_distortions", &BDLscmResult::face_distortions)
        .field("face_jacobians",   &BDLscmResult::face_jacobians)
        .field("face_cones",       &BDLscmResult::face_cones);

    // --- 求解函数 ---
    function("solveBDLscmFull", &_solveBDLscmFull);
    function("runBDLscmOnMesh", &_runOnMesh);
}
