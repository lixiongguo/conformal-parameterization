#ifndef HOLOMORPHIC_ONE_FORM_H
#define HOLOMORPHIC_ONE_FORM_H

#include "CutSeamParameterization.h"
#include "TreeCotreeBasis.h"
#include <Eigen/SparseCholesky>
#include <map>
#include <utility>
#include <vector>

/**
 * Holomorphic 1-form parameterization (Gu & Yau global conformal).
 *
 * Discrete pipeline aligned with doc §4.3–4.5:
 *  - homology basis on the closed mesh
 *  - closedness on faces (dω = 0) and cotangent harmonicity
 *  - period normalization ∫_{e_i} ω_j = δ_ij
 *  - Hodge star inside that harmonic basis (discrete wedge)
 *  - cut along the same cycles, then integrate ω + i *ω
 */
class HolomorphicOneForm : public CutSeamParameterization {
public:
    struct TreeCotree {
        std::vector<std::pair<int, int>> primalTree;
        std::vector<std::pair<int, int>> dualTree;
        std::vector<std::pair<int, int>> homologyEdges;
        std::vector<std::vector<std::pair<int, int>>> cycles;
        int visitedVertices = 0;
        int visitedFaces = 0;
        int vertexCount = 0;
        int faceCount = 0;
        int genus = 0;
    };

    HolomorphicOneForm(Mesh& mesh0);
    void parameterize() override;
    const TreeCotree& treeCotree() const { return treeCotree_; }

private:
    struct EdgeInfo {
        int v0 = 0;
        int v1 = 0;
        int f0 = -1;
        int f1 = -1;
        double length = 0.0;
    };

    void extractMesh();
    int findEdge(int a, int b) const;
    int edgeSign(int from, int to) const;
    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c) const;
    double cotanWeightAtVertex(int vi, int edgeIdx) const;
    void collectBoundaryCycles(std::vector<std::vector<std::pair<int, int>>>& cycles) const;
    void collectTreeCycles(std::vector<std::vector<std::pair<int, int>>>& cycles);
    bool solveHarmonic1Form(const std::vector<std::vector<std::pair<int, int>>>& cycles,
                            const std::vector<double>& periodRhs,
                            Eigen::VectorXd& omega);
    bool computeHodgeConjugate(const std::vector<Eigen::VectorXd>& basis,
                               const Eigen::VectorXd& omega,
                               Eigen::VectorXd& starOmega) const;
    void integrateTree(const Eigen::VectorXd& omega, const Eigen::VectorXd& starOmega);
    double faceWedge(const Eigen::Vector3d& left, const Eigen::Vector3d& right) const;
    double faceStarWedge(int face, const Eigen::Vector3d& left, const Eigen::Vector3d& right) const;
    double triangleArea(int fi) const;
    void faceEdgeValues(int fi, const Eigen::VectorXd& omega, Eigen::Vector3d& vals) const;

    Eigen::MatrixXd V;
    Eigen::VectorXi F;
    int nV = 0;
    int nF = 0;
    int nE = 0;
    int nB = 0;
    int genus = 0;
    int normalizeEdge = 0;

    std::vector<EdgeInfo> edges;
    std::map<std::pair<int, int>, int> edgeMap;
    std::vector<std::vector<std::pair<int, int>>> vtxInc;
    TreeCotree treeCotree_;
};

#endif
