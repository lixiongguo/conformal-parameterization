#ifndef MIQQUAD_H

#define MIQQUAD_H

#include "GlobalFieldsParameterization.h"

#include "MixedIntegerProgram.h"

#include <Eigen/SparseCholesky>
#include <memory>

/**
 * MIQ (Mixed-Integer Quadrangulation) global parameterization (Bommes et al. 2009).
 *
 * Phase 1 — MixedIntegerProgram on face angles theta and edge jumps p.
 * Phase 2 — Seamless parameterization with cut-graph + integer compatibility.
 */
class MIQQuad : public GlobalFieldsParameterization {

public:

    MIQQuad(Mesh& mesh0);

    void parameterize() override;

    void setCrossFieldIterations(int iters) { crossIters_ = iters; }

    void setJumpRefinePasses(int passes) { jumpRefinePasses_ = passes; }

    void setJumpBounds(int lo, int hi) { jumpLo_ = lo; jumpHi_ = hi; }

    int numFaces() const { return nF; }

    int numEdges() const { return nE; }

    const Eigen::VectorXd& faceTheta() const { return theta; }

    const Eigen::VectorXi& edgeJumps() const { return jump; }

    double crossFieldEnergy() const;

protected:

    struct Edge {
        int v1 = 0;
        int v2 = 0;
        int f1 = -1;
        int f2 = -1;
        int idx = -1;
    };

    // ---- Phase 1 ----

    bool initMeshData();
    void initCrossField();
    bool solveCrossFieldIP();
    std::vector<MixedIntegerProgram::Constraint> buildFaceAdjacencyConstraints() const;
    void wrapCrossFieldAngles();
    void buildTargetDirs();
    void buildCotLaplacian(Eigen::SparseMatrix<double>& L);
    void solvePoisson();

    // ---- Phase 2: seamless parameterization ----

    void buildCutGraph();
    bool buildCutMesh();
    void determineCutRotations();
    bool solveSeamlessParam();
    void buildCotLaplacianCut(Eigen::SparseMatrix<double>& L) const;
    int  findSharedEdge(int f1, int f2) const;
    int  countSingularities() const;

    void buildLocalFrame(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2);
    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);

    // ---- data ----

    Eigen::MatrixXd vertPos;
    Eigen::VectorXi faces;
    int nV = 0;
    int nF = 0;
    int nE = 0;
    std::vector<Edge> edgeList;

    Eigen::MatrixXd faceN;
    Eigen::MatrixXd faceT1;
    Eigen::MatrixXd faceT2;
    Eigen::VectorXd theta;
    Eigen::VectorXi jump;
    Eigen::MatrixXd faceD1;
    Eigen::MatrixXd faceD2;
    Eigen::MatrixXd UV;

    int crossIters_ = 12;
    int jumpRefinePasses_ = 2;
    int jumpLo_ = -4;
    int jumpHi_ = 4;

    std::unique_ptr<MixedIntegerProgram> mixedIntegerProgram_;

    // ---- Phase 2 data ----

    std::vector<int>  cutEdgeMap_;   // edgeList idx → 1 if cut edge, 0 otherwise
    std::vector<int>  cutEdges_;     // list of cut-edge indices into edgeList
    std::vector<int>  cutRot_;       // rotation i ∈ {0,1,2,3} per cut edge
    std::vector<int>  singularVerts_; // vertices with non-zero cross-field index

    // Cut-open mesh  (vertex-duplicated along cut edges)
    int                nCV = 0;
    Eigen::MatrixXd    cutVertPos;   // nCV × 3
    Eigen::VectorXi    cutFaces;     // nF × 3   (indices into cut vertices)
    std::vector<int>   cutParent;    // cut vertex → original vertex (for projection)

    // Integer translations per cut edge  (× 2 components)
    Eigen::VectorXi    cutTx, cutTy;
};

#endif
