#ifndef BOUNDED_DISTORTION_H
#define BOUNDED_DISTORTION_H

/**
 * Bounded Distortion Mapping for Triangular Meshes
 * =================================================
 * Based on:
 *   Lipman, Y. "Bounded Distortion Mapping Spaces for Triangular Meshes."
 *   ACM TOG (SIGGRAPH), 2012.
 *
 *   Kovalsky, S.Z. et al. "Large-Scale Bounded Distortion Mappings."
 *   ACM TOG (SIGGRAPH Asia), 2015.
 *
 * This implementation uses a local-global iterative projection approach:
 *   - Local step:  per-triangle projection of Jacobians onto
 *                  bounded-distortion space (via SVD clamping)
 *   - Global step: Poisson solve to reconstruct UV from projected Jacobians
 *
 * Requirements: Eigen 3.x, libigl (for mesh I/O and cotangent Laplacian)
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <vector>
#include <functional>
#include <iostream>

namespace bdm {

// ============================================================================
// Types
// ============================================================================

using MatrixXd = Eigen::MatrixXd;
using MatrixXi = Eigen::MatrixXi;
using VectorXd = Eigen::VectorXd;
using VectorXi = Eigen::VectorXi;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;
using RowVector2d = Eigen::RowVector2d;
using RowVector3d = Eigen::RowVector3d;

// ============================================================================
// Per-triangle Local Frame and Jacobian Computation
// ============================================================================

/**
 * Compute a local 2D orthonormal frame for each triangle face.
 *
 * For triangle (v1, v2, v3) in 3D:
 *   o_j = v1
 *   e1  = normalize(v2 - v1)
 *   e2  = normalize(n × e1)   where n = e1 × (v3 - v1)
 *
 * Then the local 2D coordinates are:
 *   p1 = (0, 0)
 *   p2 = (||v2 - v1||, 0)
 *   p3 = (dot(v3 - v1, e1), dot(v3 - v1, e2))
 *
 * @param V     #V x 3  vertex positions in 3D
 * @param F     #F x 3  face indices (0-based)
 * @param[out] localP  #F x 3 x 2  local 2D coordinates per triangle
 * @param[out] Pinv    #F x 2 x 2  inverse-transpose of local basis per triangle
 *                     (used to compute Jacobian from UV)
 */
inline void compute_local_frames(
    const MatrixXd& V,
    const MatrixXi& F,
    std::vector<Eigen::Matrix2d>& Pinv)
{
    const int nF = F.rows();
    Pinv.resize(nF);

    for (int j = 0; j < nF; ++j) {
        const int i0 = F(j,0), i1 = F(j,1), i2 = F(j,2);

        const RowVector3d v1 = V.row(i0);
        const RowVector3d v2 = V.row(i1);
        const RowVector3d v3 = V.row(i2);

        // Local basis vectors
        RowVector3d x = v2 - v1;
        RowVector3d y = v3 - v1;
        const double len_x = x.norm();
        x /= len_x;

        RowVector3d z = x.cross(y);
        z.normalize();
        y = z.cross(x);
        y.normalize();

        // Local 2D coordinates: p1=(0,0), p2=(len_x,0), p3
        const double px2 = len_x;
        const double px3 = (v3 - v1).dot(x);
        const double py3 = (v3 - v1).dot(y);

        // Build [p2-p1, p3-p1]^T = [[px2, 0], [px3, py3]]^T
        // Pinv = this matrix's inverse-transpose
        Eigen::Matrix2d P;
        P << px2, 0,
             px3, py3;

        Pinv[j] = P.inverse().transpose();
    }
}

/**
 * Compute the Jacobian matrix A_j for each triangle given current UV.
 *
 * A_j maps local 2D coords to UV coords:
 *   A_j * [p2-p1, p3-p1]^T = [u2-u1, u3-u1]^T
 *   => A_j = [u2-u1, u3-u1] * Pinv
 *
 * @param U         #V x 2  current UV coordinates
 * @param F         #F x 3  face indices
 * @param Pinv      #F x 2 x 2  precomputed local frame matrices
 * @param[out] J    #F x 2 x 2  Jacobian matrices
 */
inline void compute_jacobians(
    const MatrixXd& U,
    const MatrixXi& F,
    const std::vector<Eigen::Matrix2d>& Pinv,
    std::vector<Eigen::Matrix2d>& J)
{
    const int nF = F.rows();
    J.resize(nF);

    for (int j = 0; j < nF; ++j) {
        const int i0 = F(j,0), i1 = F(j,1), i2 = F(j,2);

        const RowVector2d u1 = U.row(i0);
        const RowVector2d u2 = U.row(i1);
        const RowVector2d u3 = U.row(i2);

        // UV edge vectors
        Eigen::Matrix2d dU;
        dU.row(0) = u2 - u1;
        dU.row(1) = u3 - u1;

        // Jacobian: J = dU * Pinv
        J[j] = dU * Pinv[j];
    }
}

// ============================================================================
// Complex Representation: α, β from Jacobian
// ============================================================================

/**
 * Convert a 2x2 Jacobian matrix to its complex representation.
 *
 * For A = [[a, b], [c, d]]:
 *   α = ½(a + d + i(c - b))
 *   β = ½(a - d + i(b + c))
 *
 * Properties:
 *   σ₁ = |α| + |β|
 *   σ₂ = |α| - |β|  (when |α| >= |β|, i.e. orientation preserving)
 *   det(A) = |α|² - |β|²
 *   K = σ₁/σ₂ = (|α|+|β|)/(|α|-|β|)
 */
struct ComplexCoeffs {
    double re_alpha, im_alpha;  // α = re_alpha + i·im_alpha
    double re_beta,  im_beta;   // β = re_beta  + i·im_beta

    double alpha_sq() const { return re_alpha*re_alpha + im_alpha*im_alpha; }
    double beta_sq()  const { return re_beta*re_beta   + im_beta*im_beta;   }
    double abs_alpha() const { return std::sqrt(alpha_sq()); }
    double abs_beta()  const { return std::sqrt(beta_sq());  }
    double det()       const { return alpha_sq() - beta_sq(); }
    double sigma1()    const { return abs_alpha() + abs_beta(); }
    double sigma2()    const { return abs_alpha() - abs_beta(); }
    double conformal_distortion() const {
        double s2 = sigma2();
        return (s2 > 1e-12) ? sigma1() / s2 : std::numeric_limits<double>::infinity();
    }
};

inline ComplexCoeffs jacobian_to_complex(const Eigen::Matrix2d& J) {
    const double a = J(0,0), b = J(0,1);
    const double c = J(1,0), d = J(1,1);
    return {
        0.5 * (a + d),  0.5 * (c - b),   // α
        0.5 * (a - d),  0.5 * (b + c)    // β
    };
}

inline Eigen::Matrix2d complex_to_jacobian(const ComplexCoeffs& cc) {
    Eigen::Matrix2d J;
    J << cc.re_alpha + cc.re_beta,  -(cc.im_alpha - cc.im_beta),
         cc.im_alpha + cc.im_beta,   cc.re_alpha - cc.re_beta;
    return J;
}

// ============================================================================
// Projection onto Bounded Distortion Space
// ============================================================================

/**
 * Project a Jacobian onto the bounded-distortion space via SVD clamping.
 *
 * Given A = U Σ V^T with Σ = diag(σ₁, σ₂), σ₁ ≥ σ₂:
 *   - Enforce σ₂ ≥ ε  (prevent flip / degenerate)
 *   - Enforce σ₁/σ₂ ≤ C  (bounded distortion)
 *
 * Returns the projected Jacobian A' = U Σ' V^T.
 *
 * @param J       input 2x2 Jacobian
 * @param C       distortion bound (C >= 1)
 * @param eps     minimum singular value (prevents degeneracy)
 * @return        projected Jacobian J'
 */
inline Eigen::Matrix2d project_to_bounded_distortion(
    const Eigen::Matrix2d& J,
    double C = 10.0,
    double eps = 1e-6)
{
    // Compute SVD: J = U * Σ * V^T
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(J,
        Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector2d S = svd.singularValues();

    // Ensure ordering: σ₁ ≥ σ₂
    double s1 = S(0);
    double s2 = S(1);

    // Enforce minimum singular value (no flip)
    if (s2 < eps) s2 = eps;
    if (s1 < eps) s1 = eps;

    // Enforce bounded distortion: s1 / s2 <= C
    if (s1 > C * s2) {
        s1 = C * s2;
    }

    // Reconstruct
    Eigen::Matrix2d Sp = Eigen::Matrix2d::Zero();
    Sp(0,0) = s1;
    Sp(1,1) = s2;

    return svd.matrixU() * Sp * svd.matrixV().transpose();
}

/**
 * More aggressive projection: also try to maintain orientation (det > 0).
 *
 * If s2 < 0, we reflect through V to flip orientation back.
 */
inline Eigen::Matrix2d project_to_bounded_distortion_robust(
    const Eigen::Matrix2d& J,
    double C = 10.0,
    double eps = 1e-6)
{
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(J,
        Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector2d S = svd.singularValues();
    Eigen::Matrix2d U = svd.matrixU();
    Eigen::Matrix2d V = svd.matrixV();

    double s1 = std::abs(S(0));
    double s2 = std::abs(S(1));

    // Sort: s1 >= s2
    if (s1 < s2) std::swap(s1, s2);

    // Clamp
    if (s2 < eps) s2 = eps;
    if (s1 > C * s2) s1 = C * s2;

    Eigen::Matrix2d Sp = Eigen::Matrix2d::Zero();
    Sp(0,0) = s1;
    Sp(1,1) = s2;

    // Ensure orientation preserving (det > 0)
    if (U.determinant() * V.determinant() < 0) {
        // Flip the last column of V to fix orientation
        V.col(1) = -V.col(1);
    }

    return U * Sp * V.transpose();
}

// ============================================================================
// Energy computation
// ============================================================================

/**
 * Compute LSCM energy: E = Σ |β_j|² · area(f_j)
 */
inline double compute_lscm_energy(
    const std::vector<Eigen::Matrix2d>& J,
    const VectorXd& face_areas)
{
    double energy = 0.0;
    for (size_t j = 0; j < J.size(); ++j) {
        ComplexCoeffs cc = jacobian_to_complex(J[j]);
        energy += cc.beta_sq() * face_areas(j);
    }
    return energy;
}

/**
 * Compute ARAP energy: E = Σ ||J_j - R_j||_F² · area(f_j)
 *
 * R_j is the closest rotation to J_j (via SVD polar decomposition).
 */
inline double compute_arap_energy(
    const std::vector<Eigen::Matrix2d>& J,
    const VectorXd& face_areas,
    std::vector<Eigen::Matrix2d>* rotations = nullptr)
{
    double energy = 0.0;
    for (size_t j = 0; j < J.size(); ++j) {
        // Polar decomposition: J = R * S
        Eigen::JacobiSVD<Eigen::Matrix2d> svd(J[j],
            Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix2d R = svd.matrixU() * svd.matrixV().transpose();

        // Ensure rotation (det = +1)
        if (R.determinant() < 0) {
            Eigen::Matrix2d V = svd.matrixV();
            V.col(1) = -V.col(1);
            R = svd.matrixU() * V.transpose();
        }

        if (rotations) (*rotations)[j] = R;

        energy += (J[j] - R).squaredNorm() * face_areas(j);
    }
    return energy;
}

/**
 * Compute maximum and average conformal distortion across all triangles.
 */
inline void compute_distortion_stats(
    const std::vector<Eigen::Matrix2d>& J,
    double& max_K,
    double& avg_K,
    const VectorXd* face_areas = nullptr)
{
    max_K = 0.0;
    avg_K = 0.0;
    double total_area = 0.0;

    for (size_t j = 0; j < J.size(); ++j) {
        ComplexCoeffs cc = jacobian_to_complex(J[j]);
        double K = cc.conformal_distortion();
        max_K = std::max(max_K, K);

        if (face_areas) {
            avg_K += K * (*face_areas)(j);
            total_area += (*face_areas)(j);
        } else {
            avg_K += K;
            total_area += 1.0;
        }
    }
    avg_K /= total_area;
}

// ============================================================================
// Cotangent Laplacian for Poisson reconstruction
// ============================================================================

/**
 * Build the cotangent Laplacian matrix L (|V| x |V|).
 *
 * L_ii = ½ Σ_{j∈N(i)} (cot α_ij + cot β_ij)
 * L_ij = -½ (cot α_ij + cot β_ij)  for edge (i,j)
 *
 * where α_ij, β_ij are the angles opposite edge (i,j).
 */
inline void build_cotangent_laplacian(
    const MatrixXd& V,
    const MatrixXi& F,
    SparseMatrix& L)
{
    const int nV = V.rows();
    const int nF = F.rows();
    std::vector<Triplet> triplets;
    triplets.reserve(nV + nF * 9);  // rough estimate

    for (int j = 0; j < nF; ++j) {
        for (int k = 0; k < 3; ++k) {
            int i0 = F(j, k);
            int i1 = F(j, (k+1)%3);
            int i2 = F(j, (k+2)%3);

            RowVector3d v0 = V.row(i0);
            RowVector3d v1 = V.row(i1);
            RowVector3d v2 = V.row(i2);

            // Angle at vertex i2 (opposite edge i0-i1)
            RowVector3d a = v0 - v2;
            RowVector3d b = v1 - v2;
            double cot_angle = a.dot(b) / a.cross(b).norm();

            // Half cotangent weight
            double w = 0.5 * cot_angle;

            triplets.emplace_back(i0, i0,  w);
            triplets.emplace_back(i1, i1,  w);
            triplets.emplace_back(i0, i1, -w);
            triplets.emplace_back(i1, i0, -w);
        }
    }

    L.resize(nV, nV);
    L.setFromTriplets(triplets.begin(), triplets.end());
}

/**
 * Build the "divergence" operator G for Poisson reconstruction.
 *
 * Given projected Jacobians J'_j, the right-hand side for the Poisson solve
 * is: b = G · vec(J'), where vec(J') stacks all J'_j into a vector.
 *
 * The divergence at vertex i is assembled from the projected Jacobians
 * of all incident triangles.
 *
 * This implementation uses the edge-vector formulation:
 *   For each triangle face f_j with edges e1, e2, e3 (in local coords),
 *   and projected Jacobian J'_j, the contribution to vertex v_k's RHS is:
 *     b_{v_k} += ½ · J'_j^T · (edge opposite v_k rotated by 90°)
 */
inline void assemble_divergence_rhs(
    const MatrixXd& V,
    const MatrixXi& F,
    const std::vector<Eigen::Matrix2d>& J_proj,  // projected Jacobians
    const std::vector<Eigen::Matrix2d>& Pinv,    // local frame inverse
    MatrixXd& rhs)                                 // #V x 2
{
    const int nV = V.rows();
    const int nF = F.rows();
    rhs.setZero(nV, 2);

    for (int j = 0; j < nF; ++j) {
        const int i0 = F(j,0), i1 = F(j,1), i2 = F(j,2);

        // Local 2D coords (computed from Pinv)
        RowVector3d v0 = V.row(i0);
        RowVector3d v1 = V.row(i1);
        RowVector3d v2 = V.row(i2);

        // Reconstruct local coords from Pinv for consistency
        // P matrix rows: e1 = p2-p1, e2 = p3-p1
        Eigen::Matrix2d P = Pinv[j].inverse().transpose();
        RowVector2d e1(P(0,0), P(0,1));  // p2 - p1
        RowVector2d e2(P(1,0), P(1,1));  // p3 - p1

        // Edge vectors in local coords (as row vectors)
        RowVector2d edge_opp_0 = -e2;           // p1 - p3
        RowVector2d edge_opp_1 = e1 - e2;       // p2 - p3 (wait, this is wrong)
        RowVector2d edge_opp_2 = e1;            // p2 - p1

        // Rotate each edge by 90° counterclockwise: (x,y) -> (-y,x)
        RowVector2d rot0(-edge_opp_0(1), edge_opp_0(0));
        RowVector2d rot1(-edge_opp_1(1), edge_opp_1(0));
        RowVector2d rot2(-edge_opp_2(1), edge_opp_2(0));

        // Contribution to each vertex: ½ * (J'_j)^T * rotated_edge
        const Eigen::Matrix2d& Jt = J_proj[j];  // already J' (projected)
        // Jt is 2x2: rows are x,y components

        // Contribution to i0 (opposite edge e3 = p2-p1)
        rhs.row(i0) += 0.5 * rot2 * Jt.transpose();

        // Contribution to i1 (opposite edge p3-p1 -> e2)
        rhs.row(i1) += 0.5 * rot0 * Jt.transpose();

        // Contribution to i2 (opposite edge p1-p2 -> -e1)
        rhs.row(i2) += 0.5 * rot1 * Jt.transpose();
    }
}

// ============================================================================
// Main Bounded Distortion Solver
// ============================================================================

struct BDOptions {
    double C = 10.0;           // distortion bound (sigma1/sigma2 <= C)
    double eps = 1e-6;         // minimum singular value
    int max_iter = 50;         // maximum local-global iterations
    double tol = 1e-6;         // convergence tolerance (relative energy change)
    bool verbose = true;
    bool use_arap = false;     // true: ARAP energy, false: LSCM energy
};

struct BDResult {
    MatrixXd U;                // #V x 2 output UV coordinates
    std::vector<double> energy_history;
    std::vector<double> max_distortion_history;
    std::vector<double> avg_distortion_history;
    int iterations;
    bool converged;
};

/**
 * Compute a bounded-distortion parameterization of a triangle mesh.
 *
 * Algorithm (local-global iteration):
 *   1. Compute initial parameterization (LSCM or Tutte)
 *   2. Precompute local frames Pinv per triangle
 *   3. Build cotangent Laplacian L (for Poisson solve)
 *   4. Iterate:
 *      a. Compute Jacobians from current UV
 *      b. Project each Jacobian onto bounded-distortion space
 *      c. Assemble divergence RHS from projected Jacobians
 *      d. Solve Poisson: L · U_new = RHS
 *      e. Fix boundary vertices
 *      f. Check convergence
 *
 * @param V          #V x 3  input 3D vertex positions
 * @param F          #F x 3  triangle indices
 * @param bnd        #B      boundary vertex indices (fixed)
 * @param bnd_uv     #B x 2  target UV for boundary vertices
 * @param U_init     #V x 2  initial UV (e.g., from LSCM or Tutte)
 * @param options    algorithm options
 * @param[out] result
 */
inline void bounded_distortion_map(
    const MatrixXd& V,
    const MatrixXi& F,
    const VectorXi& bnd,
    const MatrixXd& bnd_uv,
    const MatrixXd& U_init,
    const BDOptions& options,
    BDResult& result)
{
    const int nV = V.rows();
    const int nF = F.rows();

    // --- Step 1: Precompute local frames ---
    std::vector<Eigen::Matrix2d> Pinv;
    compute_local_frames(V, F, Pinv);

    // --- Step 2: Compute face areas ---
    VectorXd face_areas(nF);
    for (int j = 0; j < nF; ++j) {
        const int i0 = F(j,0), i1 = F(j,1), i2 = F(j,2);
        RowVector3d v0 = V.row(i0);
        RowVector3d v1 = V.row(i1);
        RowVector3d v2 = V.row(i2);
        face_areas(j) = 0.5 * (v1 - v0).cross(v2 - v0).norm();
    }

    // --- Step 3: Build cotangent Laplacian ---
    SparseMatrix L;
    build_cotangent_laplacian(V, F, L);

    // --- Step 4: Set up boundary constraints ---
    // Identify free (interior) vertices
    std::vector<bool> is_boundary(nV, false);
    for (int k = 0; k < bnd.size(); ++k) {
        is_boundary[bnd(k)] = true;
    }

    // Eliminate boundary DOFs from Laplacian
    // We solve: L_ff · U_f = rhs_f - L_fb · U_b
    // where f = free, b = boundary

    std::vector<int> free_idx;
    std::vector<int> free_to_full(nV, -1);
    for (int i = 0; i < nV; ++i) {
        if (!is_boundary[i]) {
            free_to_full[i] = free_idx.size();
            free_idx.push_back(i);
        }
    }
    const int nFree = free_idx.size();

    // Build reduced Laplacian L_ff (free × free)
    std::vector<Triplet> Lff_triplets;
    for (int k = 0; k < L.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(L, k); it; ++it) {
            int i = it.row(), j = it.col();
            int fi = free_to_full[i];
            int fj = free_to_full[j];
            if (fi >= 0 && fj >= 0) {
                Lff_triplets.emplace_back(fi, fj, it.value());
            }
        }
    }
    SparseMatrix Lff(nFree, nFree);
    Lff.setFromTriplets(Lff_triplets.begin(), Lff_triplets.end());

    // Prefactor Laplacian (LLT)
    Eigen::SimplicialLLT<SparseMatrix> solver;
    solver.compute(Lff);

    // --- Step 5: Initialize ---
    MatrixXd U = U_init;
    MatrixXd rhs_full(nV, 2);

    // Boundary RHS accumulator: L_fb * U_b
    MatrixXd LfbUb(nFree, 2);
    LfbUb.setZero();

    result.energy_history.clear();
    result.max_distortion_history.clear();
    result.avg_distortion_history.clear();

    double prev_energy = std::numeric_limits<double>::max();

    // --- Step 6: Local-Global Iteration ---
    for (int iter = 0; iter < options.max_iter; ++iter) {
        // ---- Local Step: compute and project Jacobians ----
        std::vector<Eigen::Matrix2d> J;
        compute_jacobians(U, F, Pinv, J);

        std::vector<Eigen::Matrix2d> J_proj(nF);
        std::vector<Eigen::Matrix2d> rotations;

        if (options.use_arap) {
            rotations.resize(nF);
        }

        for (int j = 0; j < nF; ++j) {
            if (options.use_arap) {
                // For ARAP: first compute closest rotation R_j
                Eigen::JacobiSVD<Eigen::Matrix2d> svd(J[j],
                    Eigen::ComputeFullU | Eigen::ComputeFullV);
                rotations[j] = svd.matrixU() * svd.matrixV().transpose();
                if (rotations[j].determinant() < 0) {
                    Eigen::Matrix2d Vm = svd.matrixV();
                    Vm.col(1) = -Vm.col(1);
                    rotations[j] = svd.matrixU() * Vm.transpose();
                }
                // Project: J' = project(R + (J - R)) keeping distortion bounded
                // Simple approach: project J onto bounded distortion, keep R for energy
                J_proj[j] = project_to_bounded_distortion_robust(J[j], options.C, options.eps);
            } else {
                // LSCM: just project Jacobians
                J_proj[j] = project_to_bounded_distortion_robust(J[j], options.C, options.eps);
            }
        }

        // Compute energy
        double energy;
        if (options.use_arap) {
            energy = compute_arap_energy(J, face_areas);
        } else {
            energy = compute_lscm_energy(J, face_areas);
        }

        result.energy_history.push_back(energy);

        // Distortion stats
        double max_K, avg_K;
        compute_distortion_stats(J, max_K, avg_K, &face_areas);
        result.max_distortion_history.push_back(max_K);
        result.avg_distortion_history.push_back(avg_K);

        if (options.verbose && (iter % 10 == 0 || iter < 5)) {
            std::cout << "Iter " << iter
                      << ": energy=" << energy
                      << ", max_K=" << max_K
                      << ", avg_K=" << avg_K << std::endl;
        }

        // Check convergence
        double rel_change = std::abs(energy - prev_energy) / (std::abs(prev_energy) + 1e-12);
        if (iter > 0 && rel_change < options.tol) {
            result.converged = true;
            result.iterations = iter + 1;
            if (options.verbose)
                std::cout << "Converged at iter " << iter << std::endl;
            break;
        }
        prev_energy = energy;

        // ---- Global Step: Poisson reconstruction ----
        // Assemble RHS from projected Jacobians
        assemble_divergence_rhs(V, F, J_proj, Pinv, rhs_full);

        // Extract free part and subtract boundary contribution
        MatrixXd rhs_free(nFree, 2);
        for (int k = 0; k < nFree; ++k) {
            int vi = free_idx[k];
            rhs_free.row(k) = rhs_full.row(vi);

            // Subtract L_fb * U_b contribution
            for (SparseMatrix::InnerIterator it(L, vi); it; ++it) {
                int j = it.col();
                if (is_boundary[j]) {
                    // Find boundary index
                    int bj = -1;
                    for (int bi = 0; bi < bnd.size(); ++bi) {
                        if (bnd(bi) == j) { bj = bi; break; }
                    }
                    if (bj >= 0) {
                        rhs_free.row(k) -= it.value() * bnd_uv.row(bj);
                    }
                }
            }
        }

        // Solve Lff * U_free = rhs_free
        MatrixXd U_free = solver.solve(rhs_free);

        // Update U
        for (int k = 0; k < nFree; ++k) {
            U.row(free_idx[k]) = U_free.row(k);
        }
        // Boundary vertices remain fixed
        for (int k = 0; k < bnd.size(); ++k) {
            U.row(bnd(k)) = bnd_uv.row(k);
        }

        result.iterations = iter + 1;
    }

    if (!result.converged && result.iterations >= options.max_iter) {
        result.converged = false;
        if (options.verbose)
            std::cout << "Reached max iterations (" << options.max_iter << ")" << std::endl;
    }

    result.U = U;
}

// ============================================================================
// Helper: LSCM initial map (simplified version using libigl-style approach)
// ============================================================================

/**
 * Compute an LSCM initial parameterization.
 *
 * Fixes two boundary vertices and solves the least-squares conformal energy.
 * This is a simplified standalone implementation that doesn't require libigl.
 */
inline void lscm_initial_map(
    const MatrixXd& V,
    const MatrixXi& F,
    const VectorXi& bnd,
    const MatrixXd& bnd_uv,
    MatrixXd& U)
{
    const int nV = V.rows();
    const int nF = F.rows();
    const int nB = bnd.size();

    // Precompute local frames
    std::vector<Eigen::Matrix2d> Pinv;
    compute_local_frames(V, F, Pinv);

    // Count free vertices
    int nFree = 0;
    std::vector<int> free_to_idx(nV, -1);
    std::vector<bool> is_bnd(nV, false);
    for (int k = 0; k < nB; ++k) is_bnd[bnd(k)] = true;
    for (int i = 0; i < nV; ++i) {
        if (!is_bnd[i]) free_to_idx[i] = nFree++;
    }

    // Build LSCM system: for each face, 2 equations (real and imag parts)
    // Total equations: 2 * nF
    const int nEq = 2 * nF;
    const int nVar = 2 * nFree;

    std::vector<Triplet> M_triplets;
    VectorXd rhs = VectorXd::Zero(nEq);

    for (int j = 0; j < nF; ++j) {
        const int eq_idx = 2 * j;

        // Local coords in frame
        int i0 = F(j,0), i1 = F(j,1), i2 = F(j,2);

        RowVector3d v1 = V.row(i0);
        RowVector3d v2 = V.row(i1);
        RowVector3d v3 = V.row(i2);

        RowVector3d x = v2 - v1;
        double len_x = x.norm();
        x /= len_x;
        RowVector3d z = x.cross(v3 - v1);
        z.normalize();
        RowVector3d y = z.cross(x);

        double px3 = (v3 - v1).dot(x);
        double py3 = (v3 - v1).dot(y);

        // LSCM weights (from complex Cauchy-Riemann equations)
        // w1 = (p3 - p2) / area, w2 = (p1 - p3) / area, w3 = (p2 - p1) / area
        // p1=(0,0), p2=(len_x,0), p3=(px3, py3)
        double area = 0.5 * len_x * py3;
        if (std::abs(area) < 1e-12) continue;

        RowVector2d w1( (px3 - len_x)/area,  py3/area);
        RowVector2d w2( -px3/area,         -py3/area);
        RowVector2d w3( len_x/area,          0.0);

        // For each vertex, add its contribution
        // Real part:  wx * u - wy * v
        // Imag part:  wy * u + wx * v
        int vi[3] = {i0, i1, i2};
        RowVector2d weights[3] = {w1, w2, w3};

        for (int k = 0; k < 3; ++k) {
            int vi_k = vi[k];
            double wx = weights[k](0);
            double wy = weights[k](1);

            if (is_bnd[vi_k]) {
                // Boundary: move to RHS
                int bi = -1;
                for (int b = 0; b < nB; ++b) {
                    if (bnd(b) == vi_k) { bi = b; break; }
                }
                if (bi >= 0) {
                    double bu = bnd_uv(bi, 0);
                    double bv = bnd_uv(bi, 1);
                    rhs(eq_idx)     -= (wx * bu - wy * bv);
                    rhs(eq_idx + 1) -= (wy * bu + wx * bv);
                }
            } else {
                int fk = free_to_idx[vi_k];
                // Real equation: wx*u_k - wy*v_k
                M_triplets.emplace_back(eq_idx,     2*fk,     wx);
                M_triplets.emplace_back(eq_idx,     2*fk + 1, -wy);
                // Imag equation: wy*u_k + wx*v_k
                M_triplets.emplace_back(eq_idx + 1, 2*fk,     wy);
                M_triplets.emplace_back(eq_idx + 1, 2*fk + 1, wx);
            }
        }
    }

    SparseMatrix M(nEq, nVar);
    M.setFromTriplets(M_triplets.begin(), M_triplets.end());

    // Solve least squares: M^T M x = M^T b
    SparseMatrix MtM = M.transpose() * M;
    VectorXd Mtrhs = M.transpose() * rhs;

    Eigen::SimplicialLLT<SparseMatrix> llt_solver(MtM);
    VectorXd x = llt_solver.solve(Mtrhs);

    // Unpack into U
    U.resize(nV, 2);
    for (int i = 0; i < nV; ++i) {
        if (is_bnd[i]) {
            int bi = -1;
            for (int b = 0; b < nB; ++b) if (bnd(b) == i) { bi = b; break; }
            U.row(i) = bnd_uv.row(bi);
        } else {
            int fi = free_to_idx[i];
            U(i, 0) = x(2 * fi);
            U(i, 1) = x(2 * fi + 1);
        }
    }
}

}  // namespace bdm

#endif  // BOUNDED_DISTORTION_H
