#ifndef ABEL_JACOBI_H
#define ABEL_JACOBI_H

#include "Mesh.h"
#include "TreeCotreeBasis.h"
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <map>
#include <vector>
#include <utility>

/**
 * Discrete Abel–Jacobi machinery for closed triangle meshes (genus g >= 1).
 *
 * Pipeline aligned with blog §4.5:
 *   1. Tree–cotree homology basis {a_i, b_i}
 *   2. Canonical holomorphic 1-forms {φ_j} with ∫_{a_i} φ_j = δ_ij
 *   3. Period lattice Γ spanned by λ_{a_k}, λ_{b_k} ∈ C^g
 *   4. Abel–Jacobi map μ(p) = (∫_{p0}^p φ_1, …, ∫_{p0}^p φ_g) mod Γ
 *   5. Divisor test μ(D) ≡ 0 (mod Γ) and period integerization for seamless UV
 */
class AbelJacobi {
public:
    /** Quad-mesh singularity: degree = 4 - valence (3-valence → +1, 5-valence → -1). */
    struct SingularPoint {
        int vertexIndex = -1;
        int valence = 4;
        int degree = 0;
    };

    struct DivisorCheck {
        bool poincareHopfOk = false;
        bool abelJacobiOk = false;
        /** μ(D) in C^g: row j = (Re, Im) of component j. */
        Eigen::MatrixXd muD;
        /** Nearest integer coefficients on {λ_{a_k}, λ_{b_k}} basis. */
        Eigen::VectorXi latticeCoeffs;
        double latticeResidual = 0.0;
    };

    explicit AbelJacobi(Mesh& mesh);

    /** Extract topology, build canonical basis and period lattice. */
    bool build();

    int genus() const { return genus_; }
    int numEdges() const { return nE_; }
    bool isReady() const { return ready_; }

    /** Poincaré–Hopf: Σ (4-k) n_k = 4χ for quad valences. */
    static bool checkPoincareHopf(int genus, int nBoundary,
                                  const std::vector<SingularPoint>& singularities);

    /** μ(D) = Σ d_i μ(p_i) and membership in Γ. */
    DivisorCheck checkDivisor(int baseVertex,
                              const std::vector<SingularPoint>& singularities) const;

    /** g×2 matrix: row j = (Re μ(p), Im μ(p)) mod fractional part w.r.t. Γ. */
    Eigen::MatrixXd abelJacobiMap(int baseVertex, int targetVertex) const;

    /** Period lattice generators: 2g columns, each a g×2 complex vector flattened as g rows × 2 cols per block. */
    const Eigen::MatrixXd& latticeGenerators() const { return latticeGens_; }

    /** Integrate Σ_j (cRe_j + i cIm_j) φ_j along BFS tree from baseVertex. */
    bool integrateHolomorphic(int baseVertex,
                              const Eigen::VectorXd& coeffRe,
                              const Eigen::VectorXd& coeffIm,
                              std::vector<Eigen::Vector2d>& uv) const;

    /**
     * Heuristic period quantization (small genus): find near-integer periods on
     * homology basis, then integrate to UV.
     */
    bool quantizeAndParameterize(int baseVertex = 0);

private:
    struct EdgeInfo {
        int v0 = 0;
        int v1 = 0;
        int f0 = -1;
        int f1 = -1;
        double length = 0.0;
    };

    struct HoloForm {
        Eigen::VectorXd omega;
        Eigen::VectorXd starOmega;
    };

    void extractMesh();
    int findEdge(int a, int b) const;
    int edgeSign(int from, int to) const;
    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c) const;
    double cotanWeightAtVertex(int vi, int edgeIdx) const;
    void collectHomologyBasis();
    bool solveHarmonicWithPeriods(const std::vector<double>& aPeriodRhs, Eigen::VectorXd& omega) const;
    void computeHodgeConjugate(const Eigen::VectorXd& omega, Eigen::VectorXd& starOmega) const;
    Eigen::Vector2d periodOnCycle(const HoloForm& form,
                                  const topology::TreeCotreeBasis::Cycle& cycle) const;
    void buildLattice();
    Eigen::Vector2d integrateEdge(int from, int to, const HoloForm& form) const;
    Eigen::MatrixXd mapAlongPath(int baseVertex, int targetVertex) const;
    bool nearestLatticePoint(const Eigen::MatrixXd& mu, Eigen::VectorXi& coeffs, double& residual) const;

    Mesh& mesh_;
    Eigen::MatrixXd V_;
    Eigen::VectorXi F_;
    int nV_ = 0;
    int nF_ = 0;
    int nE_ = 0;
    int nB_ = 0;
    int genus_ = 0;
    bool ready_ = false;

    std::vector<EdgeInfo> edges_;
    std::map<std::pair<int, int>, int> edgeMap_;
    std::vector<std::vector<std::pair<int, int>>> vtxInc_;

    std::vector<topology::TreeCotreeBasis::Cycle> aCycles_;
    std::vector<topology::TreeCotreeBasis::Cycle> bCycles_;
    std::vector<HoloForm> canonicalBasis_;
    /** 2g generators, each stored as g×2 (Re, Im per component). */
    Eigen::MatrixXd latticeGens_;
};

#endif
