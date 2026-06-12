#ifndef INCREMENTAL_FLATTENING_H
#define INCREMENTAL_FLATTENING_H

#include "Parameterization.h"
#include <Eigen/SparseCholesky>
#include <map>
#include <vector>

/**
 * Incremental automatic cone selection + global seamless parameterization.
 *
 * Phases:
 *   1. Flattening  – grow flat region Ω, redistribute curvature to cone candidates
 *   2. Rounding    – round cone & homology holonomy to π/2 multiples
 *   3. Cut seams   – Dijkstra paths from cones to boundary (+ tree-cotree for closed)
 *   4. Build M_c   – duplicate seam vertices
 *   5. Seam ARAP   – low-distortion UV with quarter-turn transition constraints
 */
class IncrementalFlattening : public Parameterization {
public:
    explicit IncrementalFlattening(Mesh& mesh0);

    void parameterize() override;

    void setMaxFlatteningIters(int iters) { m_maxFlattenIters = iters; }
    void setMaxArapIters(int iters) { m_maxArapIters = iters; }
    void setEpsilon0(double eps) { m_epsilon0 = eps; }
    void setEpsilonGrowth(double g) { m_epsilonGrowth = g; }

    const std::vector<int>& coneVertices() const { return m_coneVerts; }
    const Eigen::VectorXd& coneCurvatures() const { return m_coneK; }
    int seamEdgeCount() const;
    int cutVertexCount() const;

private:
    struct EdgeInfo {
        int v0 = 0;
        int v1 = 0;
        int f0 = -1;
        int f1 = -1;
        double length = 0.0;
        bool isInterior() const { return f0 >= 0 && f1 >= 0; }
    };

    struct SeamPath {
        std::vector<int> edges;
        std::vector<char> forward;
    };

    struct CutMesh {
        Eigen::MatrixXd cutPos;
        Eigen::VectorXi cutFaces;
        std::vector<int> origVert;
        std::vector<int> seamMate;
        std::vector<int> faceOrig;
        std::vector<int> sheetOfFace;
        std::vector<Eigen::Vector2d> cutUV;
        std::vector<bool> isCutBoundary;
        int nCutVerts = 0;
        int nCutFaces = 0;
    };

    bool extractMeshData();
    void buildCotanLaplacian(Eigen::SparseMatrix<double>& L) const;
    void buildVertexAreas(Eigen::VectorXd& areas) const;
    double cotanWeight(int vi, int vj, int opp) const;

    bool solveConstrainedPhi(const std::vector<int>& constraintVerts,
                             const Eigen::VectorXd& rhs,
                             Eigen::VectorXd& phi) const;

    void flatteningPhase();
    void roundingPhase();
    void globalArapPhase();

    void buildCutSeams();
    bool buildCutMesh();
    void buildTreeCotreeSeamPaths();

    bool identifyIsolatedCones(const std::vector<bool>& inOmega,
                               std::vector<int>& cones) const;

    double roundToHalfPi(double k) const;
    static Eigen::Matrix2d quarterTurnMatrix(int k);

    void collectHomologyCycles(std::vector<std::vector<std::pair<int, int>>>& cycles) const;
    bool solveEdgeOmegas(const Eigen::VectorXd& vertexHolonomy,
                         const std::vector<std::pair<std::vector<std::pair<int, int>>, double>>& extraCycles,
                         Eigen::VectorXd& omega,
                         bool seamOnly) const;

    void computeSeamRotations(const Eigen::VectorXd& omega);
    void solveCrossField(const Eigen::VectorXd& omega, Eigen::VectorXd& theta) const;

    void initCutMeshUV();
    void cutArapLocalStep();
    void enforceSeamRotations();
    void cutArapGlobalStep();
    void projectCutUVToMesh();

    int findEdge(int a, int b) const;
    int edgeSign(int from, int to) const;
    double transportAngle(int fFrom, int fTo, int edgeIdx) const;

    bool dijkstraToBoundary(int source, std::vector<int>& pathEdges, std::vector<char>& pathForward) const;
    bool dijkstraBetween(int source, int target, std::vector<int>& pathEdges, std::vector<char>& pathForward) const;
    bool isMeshBoundaryVertex(int vi) const;

    int nV = 0;
    int nF = 0;
    int nE = 0;
    int nB = 0;
    int genus = 0;

    Eigen::MatrixXd V;
    Eigen::VectorXi F;
    std::vector<EdgeInfo> edges;
    std::map<std::pair<int, int>, int> edgeMap;
    std::vector<std::vector<std::pair<int, int>>> vtxInc;

    Eigen::SparseMatrix<double> m_laplace;
    Eigen::VectorXd m_vertexAreas;
    Eigen::VectorXd m_K0;
    Eigen::VectorXd m_K;
    Eigen::VectorXd m_phi;

    std::vector<int> m_coneVerts;
    Eigen::VectorXd m_coneK;

    std::vector<SeamPath> m_cutPaths;
    std::vector<char> m_isSeamEdge;
    std::vector<Eigen::Matrix2d> m_seamRotation;
    Eigen::VectorXd m_omega;

    CutMesh m_cutMesh;
    std::vector<Eigen::Matrix2d> m_cutRotations;

    int m_maxFlattenIters = 80;
    int m_maxArapIters = 30;
    double m_epsilon0 = -1.0;
    double m_epsilonGrowth = 1.35;
};

#endif
