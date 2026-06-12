#ifndef N_ROSY_VECTOR_FIELDS_H
#define N_ROSY_VECTOR_FIELDS_H

#include "Mesh.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <chrono>
#include <vector>

/**
 * Unified N-RoSy vector field computation on triangle meshes.
 *
 * Computes a per-face direction field (nF × 3 matrix, one unit vector per face)
 * using one of several methods:
 *
 *   computeFromCurvature()        — principal curvature direction (Weingarten map)
 *   computeTrivialConnection()    — Levi-Civita connection smoothing (N=1, N-RoSy)
 *
 * The N-RoSy field satisfies:
 *   1. Rotational symmetry: N equivalent directions per face, separated by 2π/N.
 *   2. Smoothness across edges via Levi-Civita parallel transport.
 *   3. User-specified singularities where the total turning is prescribed.
 *
 * Input:  Mesh& (half-edge triangle mesh)
 * Output: per-face unit 3D direction (Eigen::MatrixXd nF×3)
 *         per-face angle in local frame (Eigen::VectorXd nF)
 *         per-edge period jumps (Eigen::VectorXd nE)
 *         per-vertex singularity indices (Eigen::VectorXi nV)
 */
class NRosyVectorFields {
public:
    /**
     * @param mesh  triangle mesh (half-edge, positions in mesh.vertices[i].position)
     * @param N     Rosy number: N=1 = vector field, N=2 = line field,
     *              N=4 = cross field (default)
     */
    NRosyVectorFields(Mesh& mesh, int N = 4);

    // ── Computation ──────────────────────────────────────────

    /** Principal curvature direction per face (Weingarten map). No singularities needed. */
    bool computeFromCurvature();

    /**
     * Trivial connection / N-RoSy field smoothing.
     * Solves: B·φ = 2π·q/N - K, then propagates per-face angles via parallel transport.
     * @param singularities  per-vertex integer indices q[v] (in units of 2π·N round).
     *                       q=0 = regular, q=+1 = source, q=-1 = sink.
     *                       Must satisfy Σ q = N·χ (Gauss-Bonnet).
     */
    bool computeTrivialConnection(const std::vector<int>& singularities);

    /** Auto-generate singularities by placing Q = N·χ at vertices of highest |curvature|. */
    void autoSingularities();

    // ── Output ───────────────────────────────────────────────

    /** Per-face unit direction field (nFaces × 3). */
    const Eigen::MatrixXd& getField() const { return faceField_; }

    /** Per-face angle in local frame (nFaces). */
    const Eigen::VectorXd& getTheta() const { return faceTheta_; }

    /** Per-edge period jumps (nEdges). N-RoSy: jump = round(φ / (2π/N)). */
    const Eigen::VectorXd& getEdgeJumps() const { return edgeJumps_; }

    /** Per-vertex singularity indices (nVertices). */
    const Eigen::VectorXi& getSingularityIndices() const { return singIndices_; }

    /** Wall-clock time of the last computation in ms. */
    double getLastTimeMs() const { return lastTimeMs_; }

    // ── Metadata ─────────────────────────────────────────────

    int numFaces() const { return nFaces_; }
    int numVertices() const { return nVerts_; }
    int numEdges() const { return nEdges_; }
    int N() const { return N_; }

    /** Human-readable name of the last algorithm executed. */
    const std::string& lastAlgorithm() const { return lastAlgorithm_; }

private:
    Mesh& mesh_;
    int N_;
    int nVerts_;
    int nFaces_;
    int nEdges_;

    Eigen::MatrixXd faceField_;    // nF × 3 per-face unit direction
    Eigen::VectorXd faceTheta_;    // nF per-face angle in local frame
    Eigen::VectorXd edgeJumps_;    // nE per-edge period jumps (N-RoSy only)
    Eigen::VectorXi singIndices_;  // nV per-vertex singularity indices
    double lastTimeMs_ = 0.0;
    std::string lastAlgorithm_;

    // ── Internal helpers ─────────────────────────────────────

    int countInteriorFaces() const;
    void buildLocalFrame(FaceCIter f,
                         Eigen::Vector3d& e1,
                         Eigen::Vector3d& e2,
                         Eigen::Vector3d& n) const;
    double vertexCurvature(VertexCIter v) const;
    double transportRotation(HalfEdgeCIter h) const;

    /** Build vertex-edge incidence: B ∈ ℝⁿⱽˣⁿᴱ, B[v][e] = ±1. */
    void buildIncidence(Eigen::SparseMatrix<double>& B) const;

    /** Solve min ||φ||² s.t. B·φ = rhs using KKT system. */
    bool solveKKT(const Eigen::SparseMatrix<double>& B,
                  const Eigen::VectorXd& rhs,
                  Eigen::VectorXd& phi) const;

    /** Propagate per-face angles from connection correction φ. */
    void propagateAngles(const Eigen::VectorXd& phi, double period);
};

#endif
