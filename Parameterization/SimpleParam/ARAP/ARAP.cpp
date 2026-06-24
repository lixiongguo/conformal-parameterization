#include "ARAP.h"
#include "Tutte.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SVD>
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

ARAP::ARAP(Mesh& mesh0, int maxIter)
: Parameterization(mesh0), m_maxIter(maxIter) {}

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
    // ARAP always starts from harmonic map with circle boundary (cot-Laplace).
    Tutte tutte(mesh, TutteBoundary::CIRCLE, TutteWeight::COTAN);
    tutte.parameterize();
}

void ARAP::parameterize()
{
    if (mesh.boundaries.empty()) return;

    computeLocalFrames();
    initTutte();

    for (int iter = 0; iter < m_maxIter; iter++) {
        localStep();
        globalStep();
    }

    normalize();
}
