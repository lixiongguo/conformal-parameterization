#ifndef QUADCOVER_H
#define QUADCOVER_H

#include "Parameterization.h"
#include <Eigen/SparseCholesky>
#include <vector>

/**
 * QuadCover: Global seamless quad parameterization via branch covering.
 *
 * Inputs:
 *   1. Triangle mesh — provided via Mesh& in the constructor (interior faces only).
 *   2. Initial 4-RoSy cross / frame field — one unit tangent direction per face
 *      (nFaces x 3), e.g. principal curvature directions. Pass with setCrossField();
 *      if omitted, principal directions are estimated from the Weingarten map.
 *
 * Output: per-vertex UV (mesh.vertices[i].uv) after parameterize().
 *
 * Algorithm phases (Kälberer et al. 2007):
 *   0. Optional smoothing of the input cross field (umbilic regions)
 *   1. Matching r_ij in {0,1,2,3} per edge
 *   2. Layer shift ls(v) from holonomy; build 4-sheet branch cover
 *   3. Lift frame field to covering vector field
 *   4. Hodge decomposition / Poisson on cover -> phi_tilde
 *   5. Global continuity: round mu_j to 2Z, add harmonic psi
 *   6. Project to base mesh and normalize UV
 */
class QuadCover : public Parameterization {
public:
    QuadCover(Mesh& mesh0);

    void parameterize() override;
    bool parameterizeFull();

    /** Per-face reference direction of the input cross field (rows = faces, cols = 3). */
    void setCrossField(const Eigen::MatrixXd& dirs) { setFaceDirections(dirs); }
    void setFaceDirections(const Eigen::MatrixXd& dirs);
    void setSmoothIterations(int iters);
    void setRequirePureQuads(bool pureQuads);

protected:
    struct Edge { int v1, v2, f1, f2; };

    struct CutPath {
        std::vector<int> edges;
        std::vector<char> forward;
    };

    struct CoverData {
        Eigen::MatrixXd coverPos;
        Eigen::MatrixXi coverFaces;
        Eigen::MatrixXd coverD1;
        Eigen::MatrixXd coverD2;
        Eigen::MatrixXd coverUV;
        std::vector<std::vector<int>> baseToCover;
        std::vector<int> baseSheet0;
        int nCoverVerts = 0;
        int nCoverFaces = 0;
    };

    bool initMeshData();
    bool runFullPipeline();

    void computeVertexNormals();
    void estimatePrincipalCurvature();
    void computeFaceTheta();
    void syncFaceDirsFromTheta();
    void smoothCrossField(int iters);

    void computeMatching();
    void computeLayerShift();

    bool buildBranchCover(CoverData& cover);
    bool integrateOnCover(CoverData& cover);
    bool solveCoverPoisson(const Eigen::MatrixXd& coverPos,
                           const Eigen::MatrixXi& coverFaces,
                           const Eigen::MatrixXd& coverD1,
                           const Eigen::MatrixXd& coverD2,
                           Eigen::MatrixXd& coverUV);

    void buildCutGraphAndPaths();
    void computeMuCoefficients(const CoverData& cover,
                               Eigen::VectorXd& muU,
                               Eigen::VectorXd& muV) const;
    bool buildHarmonicBasis(const CoverData& cover,
                            const Eigen::SparseMatrix<double>& L,
                            std::vector<Eigen::VectorXd>& basisU,
                            std::vector<Eigen::VectorXd>& basisV) const;
    bool enforceGlobalContinuity(CoverData& cover);
    void projectCoverUVToMesh(const CoverData& cover);
    void writeUVToMesh();

    void buildCoverLaplacian(const CoverData& cover,
                             Eigen::SparseMatrix<double>& L,
                             std::vector<std::vector<int>>& adj) const;
    int coverVertexFor(const CoverData& cover, int baseVert) const;
    double computePathPeriod(const CoverData& cover,
                             const CutPath& path,
                             int component) const;
    double roundMu(double value) const;
    int countBoundaryComponents() const;

    void buildLocalFrame(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2);
    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c) const;
    void buildCotLaplacian(Eigen::SparseMatrix<double>& L);
    void buildCoveringField();
    void solvePoisson();

    void computeTransportMatching();
    int signedMatching(int fromFace, int toFace) const;
    int findSharedEdge(int f1, int f2) const;
    void buildDualCyclePath(int cotreeEdgeIdx,
                            const std::vector<int>& parentFace,
                            CutPath& path) const;

    Eigen::MatrixXd vertPos;
    Eigen::VectorXi faces;
    int nVerts, nFaces, nEdges;

    Eigen::MatrixXd vertexNormals;
    Eigen::MatrixXd faceNormals;
    Eigen::MatrixXd faceDirs;
    Eigen::MatrixXd externalFaceDirs;
    Eigen::MatrixXd faceD1, faceD2;
    Eigen::VectorXd faceTheta;
    bool hasExternalFaceDirs;

    Eigen::VectorXi matching;
    Eigen::VectorXd layerShift;
    Eigen::MatrixXd UV;

    std::vector<Edge> edgeList;
    std::vector<CutPath> cutPaths;
    Eigen::VectorXd muU;
    Eigen::VectorXd muV;

    int smoothIters;
    bool requirePureQuads;
};

#endif
