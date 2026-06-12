#include "BoundedLscm.h"
#include "Lscm.h"
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cassert>

namespace {

constexpr double kPi     = 3.14159265358979323846;
constexpr double kTwoPi  = 6.28318530717958647692;

// cos⁻¹(4/5) ≈ 36.87°  →  tan(θ_max) = 3/4
constexpr double kCosMax = 4.0 / 5.0;
constexpr double kTanMax = 3.0 / 4.0;

// ----------------------------------------------------------------
//  Local-frame coordinates for face f (same as Lscm's helper)
// ----------------------------------------------------------------
void localBasis(const std::vector<Eigen::Vector3d>& pos,
                std::vector<Eigen::Vector2d>& coords)
{
    Eigen::Vector3d x = pos[1] - pos[0];
    Eigen::Vector3d y = pos[2] - pos[0];
    Eigen::Vector3d xhat = x; xhat.normalize();
    Eigen::Vector3d zhat = xhat.cross(y); zhat.normalize();
    Eigen::Vector3d yhat = zhat.cross(xhat); yhat.normalize();
    coords = { Eigen::Vector2d(0, 0),
               Eigen::Vector2d(x.norm(), 0),
               Eigen::Vector2d(y.dot(xhat), y.dot(yhat)) };
}

// ----------------------------------------------------------------
//  Compute α_j, β_j for face f in its *current* UV embedding
//  (no local frame – just canonical coordinates of triangle)
// ----------------------------------------------------------------
void faceAlphaBeta(const std::vector<Eigen::Vector2d>& uvLocal,
                   const std::vector<Eigen::Vector2d>& coord,
                   std::complex<double>& alpha,
                   std::complex<double>& beta)
{
    // Affine map  J = [a b; c d]  maps local-coord → uvLocal
    //  [u1-u0, u2-u0] = J · [p1-p0, p2-p0]
    double p1x = coord[1].x(), p1y = coord[1].y();
    double p2x = coord[2].x(), p2y = coord[2].y();
    double u1x = uvLocal[1].x() - uvLocal[0].x();
    double u1y = uvLocal[1].y() - uvLocal[0].y();
    double u2x = uvLocal[2].x() - uvLocal[0].x();
    double u2y = uvLocal[2].y() - uvLocal[0].y();

    double det = p1x*p2y - p2x*p1y;
    if (std::abs(det) < 1e-14) {
        alpha = std::complex<double>(1, 0);
        beta  = std::complex<double>(0, 0);
        return;
    }
    double invDet = 1.0 / det;
    double a =  invDet * ( u1x*p2y - u2x*p1y);
    double b =  invDet * (-u1x*p2x + u2x*p1x);
    double c =  invDet * ( u1y*p2y - u2y*p1y);
    double d =  invDet * (-u1y*p2x + u2y*p1x);

    // α = (a+d)/2 + i(c-b)/2
    // β = (a-d)/2 + i(b+c)/2
    alpha = std::complex<double>(0.5*(a + d), 0.5*(c - b));
    beta  = std::complex<double>(0.5*(a - d), 0.5*(b + c));
}

// ----------------------------------------------------------------
//  Conformal distortion  K = (|α|+|β|) / (|α|-|β|)
// ----------------------------------------------------------------
double conformalDistortion(const std::complex<double>& alpha,
                           const std::complex<double>& beta)
{
    double an = std::abs(alpha);
    double bn = std::abs(beta);
    if (an <= bn + 1e-12) return 1e10;               // degenerate / flipped
    return (an + bn) / (an - bn);
}

} // anonymous namespace


// ================================================================
//  Constructor
// ================================================================

BoundedLscm::BoundedLscm(Mesh& mesh0, double C, int maxIter)
    : Parameterization(mesh0)
    , m_C(std::max(C, 1.01))
    , m_maxIter(std::max(1, maxIter))
    , m_kappa((m_C - 1.0) / (m_C + 1.0))
{
}

// ================================================================
//  parameterize() – main entry point
// ================================================================

void BoundedLscm::parameterize()
{
    // ---- step 0: fallback for closed meshes ----
    if (mesh.boundaries.empty()) {
        Lscm fallback(mesh);
        fallback.parameterize();
        return;
    }

    // ---- step 1: standard LSCM as initial guess ----
    {
        Lscm init(mesh);
        init.parameterize();
    }

    // ---- step 2: count interior faces ----
    nFaces_ = 0;
    faceIndices_.clear();
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            faceIndices_.push_back(f->index);
            nFaces_++;
        }
    }
    if (nFaces_ == 0) return;

    // ---- step 3: iterative Lipman bounded optimisation ----
    int nVerts = static_cast<int>(mesh.vertices.size());
    Eigen::VectorXd uv(nVerts * 2);       // [u0,v0, u1,v1, ...]

    for (int iter = 0; iter < m_maxIter; iter++) {
        // pack current UV
        for (int i = 0; i < nVerts; i++) {
            uv(2*i)   = mesh.vertices[i].uv.x();
            uv(2*i+1) = mesh.vertices[i].uv.y();
        }

        double maxK = maxConformalDistortion(uv);
        if (maxK <= m_C + 1e-4)
            break;                           // all faces already satisfy bound

        // -- build local frames aligned with current α_j --
        std::vector<std::complex<double>> frameDir(nFaces_);
        buildLocalFrames(uv, frameDir);

        // -- build QP --
        Eigen::SparseMatrix<double> P, Aineq;
        Eigen::VectorXd bineq;
        int nVars;
        if (!buildQP(uv, frameDir, P, Aineq, bineq, nVars))
            break;

        // -- solve --
        Eigen::VectorXd sol;
        if (!solveLscmQP(P, Aineq, bineq, sol))
            break;

        // -- update UV --
        // sol contains [u0,v0, u1,v1, ...] for ALL vertices (pinned included)
        for (int i = 0; i < nVerts; i++) {
            mesh.vertices[i].uv = Eigen::Vector2d(sol(2*i), sol(2*i+1));
        }
        normalize();
    }

    normalize();
}


// ================================================================
//  buildLocalFrames – align each frame so α_j points along ℝ⁺
// ================================================================

void BoundedLscm::buildLocalFrames(
    const Eigen::VectorXd& uv,
    std::vector<std::complex<double>>& frameDir) const
{
    // frameDir[j] = unit complex number e^{i·θ} where θ = arg(α_j_prev)
    // Multiplying by conj(frameDir[j]) rotates α_j so its imaginary
    // part ≈ 0 in the new frame.

    std::vector<std::complex<double>> alpha(nFaces_), beta(nFaces_);
    computeAlphaBeta(uv, alpha, beta);

    for (int j = 0; j < nFaces_; j++) {
        double mag = std::abs(alpha[j]);
        if (mag < 1e-12)
            frameDir[j] = std::complex<double>(1, 0);
        else
            frameDir[j] = alpha[j] / mag;        // e^{i·arg(α)}
    }
}


// ================================================================
//  buildQP
//    Objective: min Σ area_j · |β_j|²   (quadratic in UV)
//    Constraints (per face):
//      1. Half‑plane:  |Im(α_j')|  ≤  kTanMax · Re(α_j')
//                       (α_j' = α_j rotated by conj(frameDir))
//      2. Bounded distortion:  |β_j'|²  ≤  m_kappa² · |α_j_prev|²
//         linearised as  Re(conj(β_prev')·β_j') ≤ m_kappa·|α_j_prev|·|β_j_prev|
//         (tangent line to the circle in the complex plane)
// ================================================================

bool BoundedLscm::buildQP(
    const Eigen::VectorXd& uvPrev,
    const std::vector<std::complex<double>>& frameDir,
    Eigen::SparseMatrix<double>& P,
    Eigen::SparseMatrix<double>& Aineq,
    Eigen::VectorXd& bineq,
    int& nVars)
{
    int nVerts = static_cast<int>(mesh.vertices.size());
    nVars = 2 * nVerts;                        // all vertices are variables

    // ---- collect per-face α_j, β_j (in original coords) ----
    std::vector<std::complex<double>> alphaPrev(nFaces_), betaPrev(nFaces_);
    computeAlphaBeta(uvPrev, alphaPrev, betaPrev);

    // ---- build objective  P  (LSCM energy = Σ area·|β|²) ----
    {
        std::vector<Eigen::Triplet<double>> pTriplets;
        pTriplets.reserve(nFaces_ * 18);  // ~ 6 nnz per row, 3 rows per face

        for (int fj = 0; fj < nFaces_; fj++) {
            int fi = faceIndices_[fj];
            FaceCIter f = mesh.faces.begin() + fi;

            // vertices of this face
            HalfEdgeCIter he = f->he;
            int v0 = he->vertex->index;
            he = he->next;
            int v1 = he->vertex->index;
            he = he->next;
            int v2 = he->vertex->index;

            // local coordinates
            std::vector<Eigen::Vector3d> pos = {
                mesh.vertices[v0].position,
                mesh.vertices[v1].position,
                mesh.vertices[v2].position
            };
            std::vector<Eigen::Vector2d> coord;
            localBasis(pos, coord);

            // Edge vectors in local coords
            Eigen::Vector2d e1 = coord[1] - coord[0];  // v1 - v0
            Eigen::Vector2d e2 = coord[2] - coord[0];  // v2 - v0

            double det = e1.x()*e2.y() - e2.x()*e1.y();
            if (std::abs(det) < 1e-14) continue;
            double invDet = 1.0 / det;

            // β_j = ½[(a-d) + i(b+c)]
            // a = invDet*(u1*e2y - u2*e1y),  d = invDet*(v1*e2x - v2*e1x)  with sign
            // b = invDet*(u2*e1x - u1*e2x),  c = invDet*(v1*e2y - v2*e1y)
            //
            // β_j is linear in UV.  Write |β_j|²·area as quadratic form.
            // β_re = ½(a-d) = ½·invDet·[ u1*e2y - u2*e1y + v1*e2x - v2*e1x ]
            //       (sign of d: -d from (a-d), and d = invDet*(-u1*e2x + u2*e1x)
            //        wait, need to re-derive carefully)
            //
            // Actually from faceAlphaBeta:
            // a = invDet*( u1x*p2y - u2x*p1y)   with p1=e1, p2=e2
            // b = invDet*(-u1x*p2x + u2x*p1x)
            // c = invDet*( u1y*p2y - u2y*p1y)
            // d = invDet*(-u1y*p2x + u2y*p1x)
            //
            // β_re = ½(a-d) = ½·invDet·(u1x*e2y - u2x*e1y + u1y*e2x - u2y*e1x)
            // β_im = ½(b+c) = ½·invDet·(-u1x*e2x + u2x*e1x + u1y*e2y - u2y*e1y)
            //
            // where u1x = uv(v1).x - uv(v0).x, etc.
            //
            // So β_re = Σ_k c_k^re · uv_k,  β_im = Σ_k c_k^im · uv_k
            // and |β|² = (β_re)² + (β_im)² = uv^T·C·uv
            // where C = (c^re)(c^re)^T + (c^im)(c^im)^T

            const double h = 0.5 * invDet;
            const double e2y = e2.y(), e1y = e1.y();
            const double e2x = e2.x(), e1x = e1.x();

            // coefficients for (u0,v0) → -sum of others
            // (u1,v1):
            double c1re = h * e2y;
            double c1im = h * (-e2x);
            // (u2,v2):
            double c2re = h * (-e1y);
            double c2im = h * e1x;
            // (v1,v1):
            double d1re = h * e2x;
            double d1im = h * e2y;
            // (v2,v2):
            double d2re = h * (-e1x);
            double d2im = h * (-e1y);

            // For (u0,v0): c0re = -(c1re+c2re), c0im = -(c1im+c2im)
            //             d0re = -(d1re+d2re), d0im = -(d1im+d2im)

            // area weight
            double area = 0.5 * std::abs(det);

            // Build symmetric block for the 6 variables (u0,v0, u1,v1, u2,v2)
            // P = C_re + C_im  where C_* = c_* · c_*^T
            auto addBlock = [&](int iu, int iv, int ju, int jv,
                                 double cire, double ciim,
                                 double cjre, double cjim) {
                double val = area * (cire*cjre + ciim*cjim);
                pTriplets.emplace_back(iu, ju, val);
                pTriplets.emplace_back(iv, jv, val);
            };

            // Self-blocks (diagonal of upper-left 2×2 for each vertex)
            auto addSelf = [&](int u, int v, double cre, double cim) {
                double val = area * (cre*cre + cim*cim);
                pTriplets.emplace_back(u, u, val);
                pTriplets.emplace_back(v, v, val);
            };

            int ui0 = 2*v0, vi0 = 2*v0+1;
            int ui1 = 2*v1, vi1 = 2*v1+1;
            int ui2 = 2*v2, vi2 = 2*v2+1;

            double c0re = -(c1re + c2re);
            double c0im = -(c1im + c2im);
            double d0re = -(d1re + d2re);
            double d0im = -(d1im + d2im);

            addSelf(ui0, vi0, c0re, c0im);
            addSelf(ui0, vi0, d0re, d0im);
            addSelf(ui1, vi1, c1re, c1im);
            addSelf(ui1, vi1, d1re, d1im);
            addSelf(ui2, vi2, c2re, c2im);
            addSelf(ui2, vi2, d2re, d2im);

            // Cross blocks u-part (c) and v-part (d) are separate.
            // For simplicity we just add both real and imag parts.
            addBlock(ui0, vi0, ui1, vi1, c0re, c0im, c1re, c1im);
            addBlock(ui0, vi0, ui1, vi1, d0re, d0im, d1re, d1im);
            addBlock(ui0, vi0, ui2, vi2, c0re, c0im, c2re, c2im);
            addBlock(ui0, vi0, ui2, vi2, d0re, d0im, d2re, d2im);
            addBlock(ui1, vi1, ui2, vi2, c1re, c1im, c2re, c2im);
            addBlock(ui1, vi1, ui2, vi2, d1re, d1im, d2re, d2im);
        }

        P.resize(nVars, nVars);
        P.setFromTriplets(pTriplets.begin(), pTriplets.end());
    }

    // ---- build inequality constraints  A·uv ≤ b ----
    {
        std::vector<Eigen::Triplet<double>> aTriplets;
        std::vector<double> bVec;
        aTriplets.reserve(nFaces_ * 12);  // up to 3 linear constraints per face

        // Also pin two vertices to fix translation+rotation
        // Use same pins as LSCM: longest boundary diameter
        int pinV0 = -1, pinV1 = -1;
        {
            double maxDist = 0;
            for (auto it1 = mesh.boundaries.begin(); it1 != mesh.boundaries.end(); ++it1) {
                HalfEdgeCIter he1 = *it1;
                do {
                    const auto& p1 = he1->vertex->position;
                    int v1 = he1->vertex->index;
                    for (auto it2 = mesh.boundaries.begin(); it2 != mesh.boundaries.end(); ++it2) {
                        HalfEdgeCIter he2 = *it2;
                        do {
                            double d = (he2->vertex->position - p1).squaredNorm();
                            if (d > maxDist) {
                                maxDist = d;
                                pinV0 = v1;
                                pinV1 = he2->vertex->index;
                            }
                            he2 = he2->next;
                        } while (he2 != *it2);
                    }
                    he1 = he1->next;
                } while (he1 != *it1);
            }
        }

        // Pin constraints: (u0=-0.5, v0=0), (u1=0.5, v1=0)
        if (pinV0 >= 0) {
            aTriplets.emplace_back((int)bVec.size(), 2*pinV0,   1.0);
            bVec.push_back(-0.5);
            aTriplets.emplace_back((int)bVec.size(), 2*pinV0,  -1.0);
            bVec.push_back(0.5);
            aTriplets.emplace_back((int)bVec.size(), 2*pinV0+1, 1.0);
            bVec.push_back(0.0);
            aTriplets.emplace_back((int)bVec.size(), 2*pinV0+1,-1.0);
            bVec.push_back(0.0);
        }
        if (pinV1 >= 0) {
            aTriplets.emplace_back((int)bVec.size(), 2*pinV1,   1.0);
            bVec.push_back(0.5);
            aTriplets.emplace_back((int)bVec.size(), 2*pinV1,  -1.0);
            bVec.push_back(-0.5);
            aTriplets.emplace_back((int)bVec.size(), 2*pinV1+1, 1.0);
            bVec.push_back(0.0);
            aTriplets.emplace_back((int)bVec.size(), 2*pinV1+1,-1.0);
            bVec.push_back(0.0);
        }

        for (int fj = 0; fj < nFaces_; fj++) {
            int fi = faceIndices_[fj];
            FaceCIter f = mesh.faces.begin() + fi;

            HalfEdgeCIter he = f->he;
            int v0 = he->vertex->index;
            he = he->next;
            int v1 = he->vertex->index;
            he = he->next;
            int v2 = he->vertex->index;

            std::vector<Eigen::Vector3d> pos = {
                mesh.vertices[v0].position,
                mesh.vertices[v1].position,
                mesh.vertices[v2].position
            };
            std::vector<Eigen::Vector2d> coord;
            localBasis(pos, coord);

            double det = coord[1].x()*coord[2].y() - coord[2].x()*coord[1].y();
            if (std::abs(det) < 1e-14) continue;
            double invDet = 1.0 / det;

            const double e1x = coord[1].x(), e1y = coord[1].y();
            const double e2x = coord[2].x(), e2y = coord[2].y();

            // α_j, β_j linear in UV (same derivation as above)
            // α_re = ½(a+d) = ½·invDet·(u1x*e2y - u2x*e1y - u1y*e2x + u2y*e1x)
            // α_im = ½(c-b) = ½·invDet·(u1y*e2y - u2y*e1y + u1x*e2x - u2x*e1x)
            //
            // Coefficients for α:
            double a1re =  0.5*invDet*e2y,  a2re = -0.5*invDet*e1y;
            double a1im =  0.5*invDet*e2x,  a2im = -0.5*invDet*e1x;
            double b1re = -0.5*invDet*e2x,  b2re =  0.5*invDet*e1x;
            double b1im =  0.5*invDet*e2y,  b2im = -0.5*invDet*e1y;

            auto a0re = -(a1re + a2re), a0im = -(a1im + a2im);
            auto b0re = -(b1re + b2re), b0im = -(b1im + b2im);

            // Rotate α,β into the local frame
            std::complex<double> fd = frameDir[fj];
            double fdr = fd.real(), fdi = fd.imag();

            // α' = conj(fd) * α
            // α'_re = fdr*α_re + fdi*α_im
            // α'_im = fdr*α_im - fdi*α_re

            // For each vertex k, its contribution to α' is a linear form.
            // We store coeffs for α'_re and α'_im per vertex.
            auto rotCoeff = [&](double cre, double cim) -> std::pair<double,double> {
                return {fdr*cre + fdi*cim, fdr*cim - fdi*cre};
            };

            auto [a0r, a0i] = rotCoeff(a0re, a0im);
            auto [b0r, b0i] = rotCoeff(b0re, b0im);
            auto [a1r, a1i] = rotCoeff(a1re, a1im);
            auto [b1r, b1i] = rotCoeff(b1re, b1im);
            auto [a2r, a2i] = rotCoeff(a2re, a2im);
            auto [b2r, b2i] = rotCoeff(b2re, b2im);

            // --- Half-plane constraint: |Im(α')| ≤ kTanMax · Re(α') ---
            // Split into two linear inequalities:
            //   Im(α') ≤  kTanMax * Re(α')
            //   Im(α') ≥ -kTanMax * Re(α')   →  -Im(α') ≤ kTanMax * Re(α')
            {
                // Im(α') - kTanMax*Re(α') ≤ 0
                int row = static_cast<int>(bVec.size());
                aTriplets.emplace_back(row, 2*v0,   a0i - kTanMax*a0r);
                aTriplets.emplace_back(row, 2*v0+1, b0i - kTanMax*b0r);
                aTriplets.emplace_back(row, 2*v1,   a1i - kTanMax*a1r);
                aTriplets.emplace_back(row, 2*v1+1, b1i - kTanMax*b1r);
                aTriplets.emplace_back(row, 2*v2,   a2i - kTanMax*a2r);
                aTriplets.emplace_back(row, 2*v2+1, b2i - kTanMax*b2r);
                bVec.push_back(0.0);

                // -Im(α') - kTanMax*Re(α') ≤ 0
                row = static_cast<int>(bVec.size());
                aTriplets.emplace_back(row, 2*v0,  -a0i - kTanMax*a0r);
                aTriplets.emplace_back(row, 2*v0+1,-b0i - kTanMax*b0r);
                aTriplets.emplace_back(row, 2*v1,  -a1i - kTanMax*a1r);
                aTriplets.emplace_back(row, 2*v1+1,-b1i - kTanMax*b1r);
                aTriplets.emplace_back(row, 2*v2,  -a2i - kTanMax*a2r);
                aTriplets.emplace_back(row, 2*v2+1,-b2i - kTanMax*b2r);
                bVec.push_back(0.0);
            }

            // --- Bounded-distortion constraint (linearised) ---
            // |β'_j|² ≤ m_kappa² · |α_prev|²
            // Linearise:  Re(conj(β_prev') · β') ≤ m_kappa · |α_prev| · |β_prev|
            {
                double aPrev = std::abs(alphaPrev[fj]);
                double bPrev = std::abs(betaPrev[fj]);
                if (bPrev < 1e-12) continue;      // already satisfied

                std::complex<double> bp = betaPrev[fj];  // original coords
                // Rotate: β' = conj(fd) * β
                std::complex<double> bpRot(fdr*bp.real() + fdi*bp.imag(),
                                           fdr*bp.imag() - fdi*bp.real());
                double bpr = bpRot.real(), bpi = bpRot.imag();
                double bn = std::abs(bpRot);
                if (bn < 1e-12) continue;
                double dirR = bpr / bn, dirI = bpi / bn;  // unit direction

                double rhs = m_kappa * aPrev * bn;

                // Coefficients for β' in rotated frame: β' = conj(fd)*β
                // β = ... from earlier. β coefficients (original):
                // β_re: same as c1re, c2re, d1re, d2re from objective section
                // Actually let me re-derive quickly.
                // β_re = ½(a-d) = ½·invDet·(u1x*e2y - u2x*e1y + u1y*e2x - u2y*e1x)
                // β_im = ½(b+c) = ½·invDet·(-u1x*e2x + u2x*e1x + u1y*e2y - u2y*e1y)
                double hf = 0.5 * invDet;
                double c1reB =  hf * e2y,  c2reB = -hf * e1y;
                double d1reB =  hf * e2x,  d2reB = -hf * e1x;
                double c1imB = -hf * e2x,  c2imB =  hf * e1x;
                double d1imB =  hf * e2y,  d2imB = -hf * e1y;

                auto c0reB = -(c1reB + c2reB), c0imB = -(c1imB + c2imB);
                auto d0reB = -(d1reB + d2reB), d0imB = -(d1imB + d2imB);

                // Rotate: β'_re = fdr*β_re + fdi*β_im
                //        β'_im = fdr*β_im - fdi*β_re
                auto rotB = [&](double cre, double cim) -> std::pair<double,double> {
                    return {fdr*cre + fdi*cim, fdr*cim - fdi*cre};
                };

                auto [br0, bi0] = rotB(c0reB, c0imB);
                auto [dr0, di0] = rotB(d0reB, d0imB);
                auto [br1, bi1] = rotB(c1reB, c1imB);
                auto [dr1, di1] = rotB(d1reB, d1imB);
                auto [br2, bi2] = rotB(c2reB, c2imB);
                auto [dr2, di2] = rotB(d2reB, d2imB);

                // Constraint: dirR·Re(β') + dirI·Im(β') ≤ rhs
                int row = static_cast<int>(bVec.size());
                aTriplets.emplace_back(row, 2*v0,   dirR*br0 + dirI*bi0);
                aTriplets.emplace_back(row, 2*v0+1, dirR*dr0 + dirI*di0);
                aTriplets.emplace_back(row, 2*v1,   dirR*br1 + dirI*bi1);
                aTriplets.emplace_back(row, 2*v1+1, dirR*dr1 + dirI*di1);
                aTriplets.emplace_back(row, 2*v2,   dirR*br2 + dirI*bi2);
                aTriplets.emplace_back(row, 2*v2+1, dirR*dr2 + dirI*di2);
                bVec.push_back(rhs);
            }
        }

        int nIneq = static_cast<int>(bVec.size());
        Aineq.resize(nIneq, nVars);
        Aineq.setFromTriplets(aTriplets.begin(), aTriplets.end());
        bineq = Eigen::VectorXd::Map(bVec.data(), nIneq);
    }

    return true;
}


// ================================================================
//  solveLscmQP — KKT-system approach for  min xᵀPx  s.t. Ax ≤ b
//
//  Uses a simple active-set heuristic: solve once ignoring inequality
//  constraints, then check which are violated and add them as
//  equality constraints.  Repeat at most twice.
// ================================================================

bool BoundedLscm::solveLscmQP(
    const Eigen::SparseMatrix<double>& P,
    const Eigen::SparseMatrix<double>& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& uv)
{
    const int N = static_cast<int>(P.rows());

    // Start with unconstrained minimisation (P is SPD in practice)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(P);
    if (solver.info() != Eigen::Success) return false;

    uv = solver.solve(Eigen::VectorXd::Zero(N));
    if (solver.info() != Eigen::Success) return false;

    // Check inequality violations; collect active set
    Eigen::VectorXd slack = A * uv - b;
    std::vector<int> activeRows;
    for (int i = 0; i < slack.size(); i++) {
        if (slack(i) > 1e-8)
            activeRows.push_back(i);
    }

    // Add active constraints as equalities and re-solve via KKT
    for (int pass = 0; pass < 2 && !activeRows.empty(); pass++) {
        int nActive = static_cast<int>(activeRows.size());
        if (nActive == 0) break;

        // Build sparse A_active (only active rows)
        std::vector<Eigen::Triplet<double>> aActiveTriplets;
        for (int k = 0; k < nActive; k++) {
            int row = activeRows[k];
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, row); it; ++it)
                aActiveTriplets.emplace_back(k, static_cast<int>(it.col()), it.value());
        }
        Eigen::SparseMatrix<double> Aact(nActive, N);
        Aact.setFromTriplets(aActiveTriplets.begin(), aActiveTriplets.end());

        // KKT system:  [ P    Aactᵀ ] [ x ]   [ 0 ]
        //              [ Aact   0   ] [ λ ] = [ b_a ]
        int Ksize = N + nActive;
        std::vector<Eigen::Triplet<double>> kktTriplets;
        kktTriplets.reserve(P.nonZeros() + 2 * Aact.nonZeros());

        // P block
        for (int k = 0; k < P.outerSize(); k++)
            for (Eigen::SparseMatrix<double>::InnerIterator it(P, k); it; ++it)
                kktTriplets.emplace_back(it.row(), it.col(), it.value());

        // Aact block (upper-right) and Aactᵀ (lower-left)
        for (int k = 0; k < Aact.outerSize(); k++)
            for (Eigen::SparseMatrix<double>::InnerIterator it(Aact, k); it; ++it) {
                kktTriplets.emplace_back(it.row(), N + it.col(), it.value());
                kktTriplets.emplace_back(N + it.col(), it.row(), it.value());
            }

        Eigen::SparseMatrix<double> KKT(Ksize, Ksize);
        KKT.setFromTriplets(kktTriplets.begin(), kktTriplets.end());

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> kktSolver;
        kktSolver.compute(KKT);
        if (kktSolver.info() != Eigen::Success) break;

        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(Ksize);
        Eigen::VectorXd bAct(nActive);
        for (int k = 0; k < nActive; k++)
            bAct(k) = b(activeRows[k]);
        rhs.tail(nActive) = bAct;

        uv = kktSolver.solve(rhs).head(N);
        if (kktSolver.info() != Eigen::Success) break;

        // Re-check violations
        slack = A * uv - b;
        activeRows.clear();
        for (int i = 0; i < slack.size(); i++)
            if (slack(i) > 1e-8)
                activeRows.push_back(i);
    }

    return true;
}


// ================================================================
//  computeAlphaBeta
// ================================================================

void BoundedLscm::computeAlphaBeta(
    const Eigen::VectorXd& uv,
    std::vector<std::complex<double>>& alpha,
    std::vector<std::complex<double>>& beta) const
{
    for (int fj = 0; fj < nFaces_; fj++) {
        int fi = faceIndices_[fj];
        FaceCIter f = mesh.faces.begin() + fi;

        HalfEdgeCIter he = f->he;
        int v0 = he->vertex->index;
        he = he->next;
        int v1 = he->vertex->index;
        he = he->next;
        int v2 = he->vertex->index;

        std::vector<Eigen::Vector3d> pos = {
            mesh.vertices[v0].position,
            mesh.vertices[v1].position,
            mesh.vertices[v2].position
        };
        std::vector<Eigen::Vector2d> coord;
        localBasis(pos, coord);

        std::vector<Eigen::Vector2d> uvLoc = {
            Eigen::Vector2d(uv(2*v0), uv(2*v0+1)),
            Eigen::Vector2d(uv(2*v1), uv(2*v1+1)),
            Eigen::Vector2d(uv(2*v2), uv(2*v2+1))
        };

        faceAlphaBeta(uvLoc, coord, alpha[fj], beta[fj]);
    }
}


// ================================================================
//  maxConformalDistortion
// ================================================================

double BoundedLscm::maxConformalDistortion(const Eigen::VectorXd& uv) const
{
    std::vector<std::complex<double>> alpha(nFaces_), beta(nFaces_);
    computeAlphaBeta(uv, alpha, beta);

    double maxK = 1.0;
    for (int j = 0; j < nFaces_; j++) {
        double K = conformalDistortion(alpha[j], beta[j]);
        if (K > maxK) maxK = K;
    }
    return maxK;
}


// ================================================================
//  extractUV (unused — parameterize applies UV directly to mesh)
// ================================================================

void BoundedLscm::extractUV(const Eigen::VectorXd& uvAll)
{
    for (size_t i = 0; i < mesh.vertices.size(); i++) {
        mesh.vertices[i].uv = Eigen::Vector2d(uvAll(2*i), uvAll(2*i+1));
    }
}
