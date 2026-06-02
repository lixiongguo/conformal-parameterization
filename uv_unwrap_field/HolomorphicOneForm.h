#ifndef HOLOMORPHIC_ONE_FORM_H
#define HOLOMORPHIC_ONE_FORM_H

#include "Parameterization.h"
#include <Eigen/SparseCholesky>
#include <map>
#include <vector>

/**
 * Holomorphic 1-form parameterization (Gu & Yau global conformal).
 *
 * Discrete pipeline aligned with doc §4.3–4.5:
 *  - closedness on faces (dω = 0)
 *  - cotangent harmonicity at vertices (Δω = 0)
 *  - period / normalization constraints (genus-aware, web: one UV per vertex)
 *  - Hodge star via discrete wedge / per-face rotation
 *  - integrate ω + i *ω along a spanning tree
 */
class HolomorphicOneForm : public Parameterization {
public:
    HolomorphicOneForm(Mesh& mesh0);
    void parameterize() override;

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
    void collectTreeCycles(std::vector<std::vector<std::pair<int, int>>>& cycles) const;
    bool solveHarmonic1Form(Eigen::VectorXd& omega);
    void computeHodgeConjugate(const Eigen::VectorXd& omega, Eigen::VectorXd& starOmega) const;
    void integrateTree(const Eigen::VectorXd& omega, const Eigen::VectorXd& starOmega);
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
};

#endif
