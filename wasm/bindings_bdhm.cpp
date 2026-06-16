// ============================================================
//  BD-HM 绑定 — 有界失真调和映射 (Bounded Distortion Harmonic Map)
// ============================================================
#include "BoundedDistortionHarmonicMap.hpp"
#include "CauchyCoordinates.hpp"
#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>

using namespace emscripten;

// ========== 辅助函数 ==========
static bdhm::VecC _valToVecC(val arr) {
    int len = arr["length"].as<int>();
    bdhm::VecC v(len / 2);
    for (int i = 0; i < len / 2; ++i)
        v[i] = bdhm::Complex(arr[i * 2].as<double>(), arr[i * 2 + 1].as<double>());
    return v;
}

// ========== 结果结构 ==========
struct BDHMResult {
    std::vector<double> mapped_boundary;       // [re,im,...]
    std::vector<double> mapped_samples;        // [re,im,...]
    std::vector<double> mapped_constraints;    // [re,im,...]
    double max_distortion;
    double min_jacobian;
    double max_cone_violation;
    double final_energy;
    int iterations;
};

static BDHMResult _packResult(const bdhm::SolveResult& r) {
    BDHMResult out;
    for (int i = 0; i < r.phi.size(); ++i) {
        out.mapped_boundary.push_back(r.phi[i].real());
        out.mapped_boundary.push_back(r.phi[i].imag());
    }
    for (int i = 0; i < r.mapped_energy_samples.size(); ++i) {
        out.mapped_samples.push_back(r.mapped_energy_samples[i].real());
        out.mapped_samples.push_back(r.mapped_energy_samples[i].imag());
    }
    for (int i = 0; i < r.mapped_constraints.size(); ++i) {
        out.mapped_constraints.push_back(r.mapped_constraints[i].real());
        out.mapped_constraints.push_back(r.mapped_constraints[i].imag());
    }
    out.max_distortion    = r.stats.max_distortion;
    out.min_jacobian      = r.stats.min_jacobian;
    out.max_cone_violation = r.stats.max_cone_violation;
    out.final_energy      = r.final_energy;
    out.iterations         = r.iterations;
    return out;
}

// ========== 主求解函数 ==========
static BDHMResult _solve(
    val boundary_js,             // [re,im,...] — 目标边界多边形
    val energy_samples_js,       // [re,im,...] — 内部采样点
    val constraints_src_js,      // [re,im,...] — 约束源点
    val constraints_dst_js,      // [re,im,...] — 约束目标点
    val constraints_weight_js,   // [w,...] — 约束权重
    double distortion_bound,
    double min_jacobian,
    double distortion_penalty,
    double jacobian_penalty,
    double position_weight,
    double anti_holomorphic_weight,
    double phi_reference_weight,
    double initial_step,
    int max_iterations)
{
    // 提取数据
    bdhm::VecC boundary      = _valToVecC(boundary_js);
    bdhm::VecC samples       = _valToVecC(energy_samples_js);
    bdhm::VecC src_c         = _valToVecC(constraints_src_js);
    bdhm::VecC dst_c         = _valToVecC(constraints_dst_js);

    int nc = src_c.size();
    int nw = constraints_weight_js["length"].as<int>();
    std::vector<bdhm::PointConstraint> constraints(nc);
    for (int i = 0; i < nc; ++i) {
        constraints[i].source = src_c[i];
        constraints[i].target = dst_c[i];
        constraints[i].weight = (i < nw) ? constraints_weight_js[i].as<double>() : 1.0;
    }

    // 配置
    bdhm::Options opts;
    opts.distortion_bound          = distortion_bound;
    opts.min_jacobian              = min_jacobian;
    opts.distortion_penalty        = distortion_penalty;
    opts.jacobian_penalty          = jacobian_penalty;
    opts.position_weight           = position_weight;
    opts.anti_holomorphic_weight   = anti_holomorphic_weight;
    opts.phi_reference_weight      = phi_reference_weight;
    opts.initial_step              = initial_step;
    opts.max_iterations            = max_iterations;
    opts.gradient_tolerance        = 1e-10;
    opts.solver_method             = bdhm::SolverMethod::GradientDescent;

    // 求解
    bdhm::SolveResult result = bdhm::nloP2PHarmonic(boundary, samples, constraints, opts);
    return _packResult(result);
}

// ========== Cauchy 重心坐标（直接调和映射，无有界失真优化）==========
struct CauchyResult {
    std::vector<double> mapped_samples;     // [re,im,...]
};

static CauchyResult _evaluateCauchy(
    val boundary_js,             // 源边界 [re,im,...]（单位正方形）
    val target_boundary_js,      // 目标边界 [re,im,...]
    val query_points_js)         // 查询点 [re,im,...]
{
    cauchy::VecC cage     = _valToVecC(boundary_js);
    cauchy::VecC target   = _valToVecC(target_boundary_js);
    cauchy::VecC queries  = _valToVecC(query_points_js);

    cauchy::VecC mapped = cauchy::evaluateCauchyMap(cage, target, queries);

    CauchyResult out;
    for (int i = 0; i < mapped.size(); ++i) {
        out.mapped_samples.push_back(mapped[i].real());
        out.mapped_samples.push_back(mapped[i].imag());
    }
    return out;
}

// ============================================================
//  Embind 注册
// ============================================================
EMSCRIPTEN_BINDINGS(BDHMModule) {

    value_object<CauchyResult>("CauchyResult")
        .field("mapped_samples", &CauchyResult::mapped_samples);

    value_object<BDHMResult>("BDHMResult")
        .field("mapped_boundary",    &BDHMResult::mapped_boundary)
        .field("mapped_samples",     &BDHMResult::mapped_samples)
        .field("mapped_constraints", &BDHMResult::mapped_constraints)
        .field("max_distortion",     &BDHMResult::max_distortion)
        .field("min_jacobian",       &BDHMResult::min_jacobian)
        .field("max_cone_violation", &BDHMResult::max_cone_violation)
        .field("final_energy",       &BDHMResult::final_energy)
        .field("iterations",         &BDHMResult::iterations);

    function("solve", &_solve);
    function("evaluateCauchy", &_evaluateCauchy);
}
