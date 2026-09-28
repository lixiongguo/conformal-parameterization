#include "ARAP.h"
#include "Tutte.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>

namespace {

// Cotangent of the angle at vertex `at`, opposite edge (b, c).
double cotanAt(const Eigen::Vector3d& at,
               const Eigen::Vector3d& b,
               const Eigen::Vector3d& c)
{
    const Eigen::Vector3d e1 = b - at;
    const Eigen::Vector3d e2 = c - at;
    const double sinA = e1.cross(e2).norm();
    if (sinA < 1e-14) return 0.0;
    return e1.dot(e2) / sinA;
}

// Stable positive weight (same clamp as Tutte embedding).
double edgeWeight(const Eigen::Vector3d& opp,
                  const Eigen::Vector3d& a,
                  const Eigen::Vector3d& b)
{
    return std::max(cotanAt(opp, a, b), 1e-8);
}

Eigen::Matrix2d rotationFromCovariance(const Eigen::Matrix2d& S)
{
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix2d U = svd.matrixU();
    Eigen::Matrix2d V = svd.matrixV();
    Eigen::Matrix2d R = U * V.transpose();
    if (R.determinant() < 0.0) {
        U.col(1) *= -1.0;
        R = U * V.transpose();
    }
    return R;
}

// 与 Parameterization::normalize() 完全相同的变换（面积重心 + 最大半径归一化），
// 区别是该变换会**同时**作用到 Tutte 初始化快照 initUv 上，从而让
// 「Tutte 初始 UV」与「ARAP 优化后 UV」处在同一坐标系、可以直接叠加对比。
void normalizeMeshAndInit(Mesh& mesh, std::vector<double>& initUv)
{
    double totalArea = 0.0;
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        const Eigen::Vector2d& a(f->he->vertex->uv);
        const Eigen::Vector2d& b(f->he->next->vertex->uv);
        const Eigen::Vector2d& c(f->he->next->next->vertex->uv);
        const Eigen::Vector2d u = b - a;
        const Eigen::Vector2d v = c - a;
        const double area = 0.5 * (u.x() * v.y() - v.x() * u.y());
        center += area * ((a + b + c) / 3.0);
        totalArea += area;
    }
    if (std::abs(totalArea) > 0.0) center /= totalArea;

    double r = 0.0;
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->uv -= center;
        r = std::max(r, v->uv.squaredNorm());
    }
    r = std::sqrt(r);
    if (r <= 0.0) return;
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->uv /= r;
    }

    for (size_t i = 0; i + 1 < initUv.size(); i += 2) {
        initUv[i]     = (initUv[i]     - center.x()) / r;
        initUv[i + 1] = (initUv[i + 1] - center.y()) / r;
    }
}

void collectFaceCorners(FaceCIter f,
                        std::vector<Eigen::Vector3d>& pos3D,
                        std::vector<Eigen::Vector2d>& uvs,
                        std::vector<int>& vIdx)
{
    pos3D.resize(3);
    uvs.resize(3);
    vIdx.resize(3);
    HalfEdgeCIter he = f->he;
    for (int i = 0; i < 3; i++) {
        pos3D[i] = he->vertex->position;
        uvs[i] = he->vertex->uv;
        vIdx[i] = he->vertex->index;
        he = he->next;
    }
}

} // namespace

ARAP::ARAP(Mesh& mesh0, int maxIter, TutteBoundary boundaryShape, TutteWeight weight)
: Parameterization(mesh0),
  m_initBoundary(boundaryShape),
  m_initWeight(weight),
  m_maxIter(maxIter) {}

void ARAP::computeLocalFrames()
{
    m_localRefs.clear();
    m_localRefs.resize(mesh.faces.size());
    m_rotations.resize(mesh.faces.size());

    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        const int fi = f->index;

        std::vector<Eigen::Vector3d> pos3D(3);
        HalfEdgeCIter he = f->he;
        for (int i = 0; i < 3; i++) {
            pos3D[i] = he->vertex->position;
            he = he->next;
        }

        const Eigen::Vector3d e1 = pos3D[1] - pos3D[0];
        const Eigen::Vector3d e2 = pos3D[2] - pos3D[0];

        Eigen::Vector3d xAxis = e1;
        const double e1Len = xAxis.norm();
        if (e1Len < 1e-14) continue;
        xAxis /= e1Len;

        Eigen::Vector3d zAxis = xAxis.cross(e2);
        const double zLen = zAxis.norm();
        if (zLen < 1e-14) continue;
        zAxis /= zLen;

        const Eigen::Vector3d yAxis = zAxis.cross(xAxis);

        m_localRefs[fi].resize(3);
        m_localRefs[fi][0] = Eigen::Vector2d(0, 0);
        m_localRefs[fi][1] = Eigen::Vector2d(e1.dot(xAxis), e1.dot(yAxis));
        m_localRefs[fi][2] = Eigen::Vector2d(e2.dot(xAxis), e2.dot(yAxis));

        m_rotations[fi] = Eigen::Matrix2d::Identity();
    }
}

void ARAP::localStep()
{
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        const int fi = f->index;
        if (m_localRefs[fi].size() != 3) continue;

        std::vector<Eigen::Vector3d> pos3D;
        std::vector<Eigen::Vector2d> uvs;
        std::vector<int> vIdx;
        collectFaceCorners(f, pos3D, uvs, vIdx);

        const auto& x = m_localRefs[fi];
        Eigen::Matrix2d S = Eigen::Matrix2d::Zero();

        const int edges[3][2] = {{1, 2}, {2, 0}, {0, 1}};
        const int opposite[3] = {0, 1, 2};

        for (int e = 0; e < 3; e++) {
            const int j = edges[e][0];
            const int k = edges[e][1];
            const int opp = opposite[e];

            const Eigen::Vector2d du = uvs[j] - uvs[k];
            const Eigen::Vector2d dx = x[j] - x[k];
            const double w = edgeWeight(pos3D[opp], pos3D[j], pos3D[k]);
            S += w * du * dx.transpose();
        }

        m_rotations[fi] = rotationFromCovariance(S);
    }
}

void ARAP::globalStep()
{
    const int N = (int)mesh.vertices.size();

    std::vector<bool> isBoundary(N, false);
    for (HalfEdgeIter heStart : mesh.boundaries) {
        HalfEdgeCIter h = heStart;
        do {
            isBoundary[h->vertex->index] = true;
            h = h->next;
        } while (h != heStart);
    }

    std::vector<int> interiorIdx(N, -1);
    int nInterior = 0;
    for (int i = 0; i < N; i++) {
        if (!isBoundary[i]) interiorIdx[i] = nInterior++;
    }
    if (nInterior == 0) return;

    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd bx = Eigen::VectorXd::Zero(nInterior);
    Eigen::VectorXd by = Eigen::VectorXd::Zero(nInterior);

    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        const int fi = f->index;
        if (m_localRefs[fi].size() != 3) continue;

        const Eigen::Matrix2d& R = m_rotations[fi];
        const auto& x = m_localRefs[fi];

        std::vector<Eigen::Vector3d> pos3D;
        std::vector<Eigen::Vector2d> uvs;
        std::vector<int> vIdx;
        collectFaceCorners(f, pos3D, uvs, vIdx);

        const int edges[3][2] = {{1, 2}, {2, 0}, {0, 1}};
        const int opposite[3] = {0, 1, 2};

        for (int e = 0; e < 3; e++) {
            const int a = edges[e][0];
            const int b = edges[e][1];
            const int opp = opposite[e];

            const double w = edgeWeight(pos3D[opp], pos3D[a], pos3D[b]);
            const Eigen::Vector2d rotDx = R * (x[a] - x[b]);

            const int vi = vIdx[a];
            const int vj = vIdx[b];

            if (!isBoundary[vi]) {
                const int ri = interiorIdx[vi];
                triplets.emplace_back(ri, ri, w);
                bx(ri) += w * rotDx.x();
                by(ri) += w * rotDx.y();

                if (isBoundary[vj]) {
                    bx(ri) += w * mesh.vertices[vj].uv.x();
                    by(ri) += w * mesh.vertices[vj].uv.y();
                } else {
                    triplets.emplace_back(ri, interiorIdx[vj], -w);
                }
            }

            if (!isBoundary[vj]) {
                const int rj = interiorIdx[vj];
                triplets.emplace_back(rj, rj, w);
                bx(rj) -= w * rotDx.x();
                by(rj) -= w * rotDx.y();

                if (isBoundary[vi]) {
                    bx(rj) += w * mesh.vertices[vi].uv.x();
                    by(rj) += w * mesh.vertices[vi].uv.y();
                } else {
                    triplets.emplace_back(rj, interiorIdx[vi], -w);
                }
            }
        }
    }

    Eigen::SparseMatrix<double> L(nInterior, nInterior);
    L.setFromTriplets(triplets.begin(), triplets.end());
    L.makeCompressed();

    // Mild regularization helps near-degenerate cot weights.
    for (int i = 0; i < nInterior; i++) {
        L.coeffRef(i, i) += 1e-10;
    }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Eigen::Success) return;

    const Eigen::VectorXd ux = solver.solve(bx);
    const Eigen::VectorXd uy = solver.solve(by);
    if (solver.info() != Eigen::Success) return;

    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        const int vi = v->index;
        if (!isBoundary[vi]) {
            const int ri = interiorIdx[vi];
            v->uv = Eigen::Vector2d(ux(ri), uy(ri));
        }
    }
}

void ARAP::initTutte()
{
    // ARAP 的前置步骤：Tutte（调和映射）初始化，边界形状/权重由构造参数决定。
    // 关闭 Tutte 内部归一化以保留原始调和解——稍后与 ARAP 结果用同一变换归一化。
    Tutte tutte(mesh, m_initBoundary, m_initWeight);
    tutte.setNormalize(false);
    tutte.parameterize();

    // 快照 Tutte 初始化 UV（扁平 2N），供页面展示「初始 / 优化后」对比
    m_initUv.clear();
    m_initUv.reserve(mesh.vertices.size() * 2);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        m_initUv.push_back(v->uv.x());
        m_initUv.push_back(v->uv.y());
    }
}

// 记录一条收敛指标。调用前 m_rotations 必须是当前 uv 的逐面最优旋转
// （即刚跑完 localStep），此时能量是真正的 ARAP 能量 min_R E(u,R)。
void ARAP::recordIteration()
{
    const int edges[3][2] = {{1, 2}, {2, 0}, {0, 1}};
    const int opposite[3] = {0, 1, 2};

    double energy = 0.0;
    double sumAreaLocal = 0.0;
    double sumAreaUv = 0.0;
    int flipped = 0;

    m_areaRatios.clear();

    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        const int fi = f->index;
        if (m_localRefs[fi].size() != 3) continue;

        std::vector<Eigen::Vector3d> pos3D;
        std::vector<Eigen::Vector2d> uvs;
        std::vector<int> vIdx;
        collectFaceCorners(f, pos3D, uvs, vIdx);

        const auto& x = m_localRefs[fi];
        const Eigen::Matrix2d& R = m_rotations[fi];

        // ARAP 能量：与 global step 目标函数完全同一个量
        for (int e = 0; e < 3; e++) {
            const int a = edges[e][0];
            const int b = edges[e][1];
            const int opp = opposite[e];
            const double w = edgeWeight(pos3D[opp], pos3D[a], pos3D[b]);
            const Eigen::Vector2d du = uvs[a] - uvs[b];
            energy += w * (du - R * (x[a] - x[b])).squaredNorm();
        }

        // 面积畸变与翻面：局部标架 (x0,x1,x2) 恒为正向，故 UV 有向面积为负即折叠
        const Eigen::Vector2d l1 = x[1] - x[0];
        const Eigen::Vector2d l2 = x[2] - x[0];
        const Eigen::Vector2d m1 = uvs[1] - uvs[0];
        const Eigen::Vector2d m2 = uvs[2] - uvs[0];
        const double areaLocal = 0.5 * (l1.x() * l2.y() - l2.x() * l1.y());
        const double areaUv    = 0.5 * (m1.x() * m2.y() - m2.x() * m1.y());

        if (areaUv < 0.0) flipped++;
        if (areaLocal > 1e-12 && std::abs(areaUv) > 1e-15) {
            sumAreaLocal += areaLocal;
            sumAreaUv += std::abs(areaUv);
            m_areaRatios.push_back(std::abs(areaUv) / areaLocal);
        }
    }

    // 面积畸变：先扣掉 UV 的整体缩放（对参数化而言，整体缩放不是畸变），
    // 再统计 mean|ln(r_i / r_global)|；完美等距（含等比例缩放）时为 0。
    double areaDist = 0.0;
    if (!m_areaRatios.empty() && sumAreaLocal > 1e-12 && sumAreaUv > 1e-15) {
        const double rGlobal = sumAreaUv / sumAreaLocal;
        double sum = 0.0;
        for (double r : m_areaRatios) sum += std::abs(std::log(r / rGlobal));
        areaDist = sum / (double)m_areaRatios.size();
    }

    m_iterEnergy.push_back(energy);
    m_iterAreaDist.push_back(areaDist);
    m_iterFlipped.push_back(flipped);
}

void ARAP::parameterize()
{
    if (mesh.boundaries.empty()) return;

    computeLocalFrames();
    initTutte();

    m_iterEnergy.clear();
    m_iterAreaDist.clear();
    m_iterFlipped.clear();

    // 原循环是 localStep(); globalStep();。这里把 local 步提到度量之前：
    //   localStep 只写 m_rotations、不改 uv，所以“先 local 后 global”与
    //   “先 global 后 local”对 uv 的结果完全一致；但这样每轮迭代的 R_f 恰好就是
    //   度量所需要的“当前 uv 的最优旋转”，统计不会额外增加 SVD 开销。
    localStep();
    recordIteration();                 // 迭代 0：Tutte 初始化状态（收益基线）

    for (int iter = 0; iter < m_maxIter; iter++) {
        globalStep();
        localStep();                   // 新 uv 的最优旋转，同时也供下一轮 global 使用
        recordIteration();             // 迭代 iter + 1
    }

    // 归一化：网格 UV 与 Tutte 初始快照共用同一变换，保证两者处于同一坐标系
    normalizeMeshAndInit(mesh, m_initUv);
}
