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

 * Phase 2 — Poisson parameterization aligned to the optimized 4-RoSy field.

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



    bool initMeshData();

    void initCrossField();

    bool solveCrossFieldIP();

    std::vector<MixedIntegerProgram::Constraint> buildFaceAdjacencyConstraints() const;

    void wrapCrossFieldAngles();

    void buildTargetDirs();

    void buildCotLaplacian(Eigen::SparseMatrix<double>& L);

    void solvePoisson();



    void buildLocalFrame(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2);

    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);



    Eigen::MatrixXd vertPos;

    Eigen::VectorXi faces;

    int nV = 0;

    int nF = 0;

    int nE = 0;

    std::vector<Edge> edgeList;



    Eigen::MatrixXd faceN;

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

};



#endif

