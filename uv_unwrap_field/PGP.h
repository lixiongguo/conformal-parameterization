#ifndef PGP_H
#define PGP_H

#include "../Parameterization.h"
#include <Eigen/SparseCholesky>
#include <vector>

/**
 * PGP: Periodic Global Parameterization (Ray et al. 2006).
 *
 * Inputs:
 *   1. Triangle mesh — provided via Mesh& in the constructor.
 *   2. 4-RoSy cross field — one unit direction per face (optional;
 *      if omitted, estimated from principal curvature).
 *
 * Output: per-vertex UV (mesh.vertices[i].uv) after parameterize().
 *
 * Algorithm:
 *   1. Extract mesh data (Eigen format), build per-face local frames.
 *   2. Get/estimate cross field directions d_T^u, d_T^v per face.
 *   3. Alternating optimization:
 *      a) Fix integer variables k_T, p_T → solve linear system for (u, v)
 *      b) Fix (u, v) → per-face greedy rounding of k_T, p_T
 *      c) Local ±1 search on p_T to reduce energy
 *   4. Write UV to mesh and normalize.
 */
class PGP : public Parameterization {
public:
    PGP(Mesh& mesh0);

    void parameterize() override;

    /** Set external cross field directions: matrix (nF x 3), one unit vector per face. */
    void setCrossField(const Eigen::MatrixXd& dirs);

    /** Number of alternating optimization iterations (default 3). */
    void setMaxIterations(int iters) { maxIters = iters; }

    /** Enable/disable local search on period jumps (default true). */
    void setLocalSearch(bool enable) { doLocalSearch = enable; }

protected:
    struct InternalEdge {
        int v1, v2;      // vertex indices
        int f1, f2;      // adjacent face indices
        int idx;         // edge index
    };

    // --- Mesh extraction ---
    bool initMeshData();

    // --- Cross field ---
    void computeCrossField();

    // --- Local geometry ---
    void buildLocalFrames();
    void buildLocalFrame(const Eigen::Vector3d& n,
                         Eigen::Vector3d& t1,
                         Eigen::Vector3d& t2);

    // --- Linear system ---
    void assembleLaplacian(Eigen::SparseMatrix<double>& L);
    void assembleRHS(const Eigen::VectorXi& kT,
                     const Eigen::MatrixXi& pT,
                     Eigen::VectorXd& bu,
                     Eigen::VectorXd& bv);
    void solveForUV(const Eigen::VectorXi& kT,
                    const Eigen::MatrixXi& pT,
                    Eigen::VectorXd& u,
                    Eigen::VectorXd& v);
    Eigen::Vector2d computeGradient(const Eigen::VectorXd& u, int f) const;

    // --- Per-face integer optimization ---
    double faceEnergy(int f, const Eigen::Vector2d& gu, const Eigen::Vector2d& gv,
                      int k, const Eigen::Vector2i& p) const;
    void optimizePerFace(const Eigen::VectorXd& u, const Eigen::VectorXd& v,
                         Eigen::VectorXi& kT, Eigen::MatrixXi& pT);
    void localSearchP(const Eigen::VectorXd& u, const Eigen::VectorXd& v,
                      Eigen::VectorXi& kT, Eigen::MatrixXi& pT);

    // --- Output ---
    void writeUV(const Eigen::VectorXd& u, const Eigen::VectorXd& v);

    // --- Mesh data (Eigen format) ---
    int nV, nF;
    Eigen::MatrixXd vertPos;           // vertex positions (nV x 3)
    Eigen::MatrixXi faceIndices;       // face → vertex indices (nF x 3)
    std::vector<InternalEdge> edges;   // internal edges for cotan Laplacian

    // --- Cross field ---
    Eigen::MatrixXd faceDirs;          // per-face direction field (nF x 3)
    std::vector<Vector2d> faceD1;      // per-face d1 in 2D local frame (nF)
    std::vector<Vector2d> faceD2;      // per-face d2 = rot90(d1) in 2D local frame (nF)

    // --- Per-face geometry ---
    Eigen::MatrixXd faceCenters;       // per-face centroid (nF x 3)
    Eigen::MatrixXd faceNormals;       // per-face normal (nF x 3)
    std::vector<Eigen::Matrix2d> faceGrad; // per-face gradient operator (G_T)
    std::vector<double> faceAreas;     // per-face area

    // --- Per-face integer variables ---
    Eigen::VectorXi kVar;              // rotation index k_T ∈ {0,1,2,3}
    Eigen::MatrixXi pVar;              // period jump p_T ∈ Z² (nF x 2)

    // --- Sparse Cholesky solver ---
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    bool solverReady;

    // --- Options ---
    int maxIters;
    bool doLocalSearch;
    bool hasExternalDirs;
};

#endif
