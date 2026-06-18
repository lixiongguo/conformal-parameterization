#pragma once

#include "../HGP/HarmonicParametrization.h"
#include "../HGP/Borders.h"
#include "FastHGPNumerics.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_set>
#include <vector>
#include <string>

class FastHGP : public HarmonicParametrization {
public:
    FastHGP() {
        mMethodName = "FastHGP";
        mHasCones = false;
    }
    ~FastHGP() = default;

    bool run(std::string& objpath, std::string& vfPath) override;

    void setSegSize(int segSize);
    void setFixCot(bool fixCot);

protected:
    void visualize() override;
    void calcDistortion() override;
    void coneAngleDetection(int& numWrongAngles, int& numWrongConeAngles) override;

private:
    int mSegSize = 40;
    bool mFixCot = true;
    int mNumOfBorderVertices = 0;
    int mNumDOF = 0;
    int mSizeOfMatrix = 0;
    bool mCalcFramesFromVecField = false;
    bool mFramesFromFile = false;
    std::string mFramesFilePath;
    bool mATPandNewtonPassed = false;
    bool mATPSuccess = true;
    double mTotalTime = 0.0;

    std::vector<Vertex_handle> mConesAndMetaMap;
    std::vector<int> mIndicesOfMetaVerticesInUVbyGeneralIndex;
    std::vector<bool> mIsMeta;
    std::vector<std::vector<Halfedge_handle>> mMetaVerticesInBorderByHalfEdge;
    std::unordered_set<int> mConesAndNearConesMapOfRowsInKKT;

    Eigen::SparseMatrix<FastHGPNumerics::Complex> mHarmonicBasisEigen;
    Eigen::SparseMatrix<double> mKKtEigen;

    Eigen::MatrixXcd mJfz;
    Eigen::MatrixXcd mJfbz;
    Eigen::VectorXcd mFrames;
    Eigen::MatrixXd mVerticesByHalfedges;
    Eigen::MatrixXi mReducedFaces;
    Eigen::VectorXd mArea;
    Eigen::MatrixXd mReducedSolution;
    std::vector<int> mFixedIndices;
    Eigen::VectorXd mFixedValues;
    std::vector<Facet_handle> mReducedFacets;

    bool initialize(Borders& borders);
    bool loadMesh(std::string& objpath, std::string& vfPath);
    void getSettings();
    void getConesMap();
    void getBordersMapAndSetMetaVertices(Borders& borders);
    void getVerticesThatMustBeMeta(Borders& borders);
    bool prepareReducedMeshData();
    bool runAlgorithm(Borders& borders);

    void constructKKTmatrix(Borders& borders, int& conesConstraintsStartRow);
    void FillLaplacianInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues);
    void updateTermInLaplacianByHalfEdge(Halfedge_handle h, const double& term,
                                         std::vector<Eigen::Triplet<double>>& tripletListValues);
    void updateTermInLaplacianByIndices(int i, int j, const double& term,
                                        std::vector<Eigen::Triplet<double>>& tripletListValues);
    void FillRotationConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues, int& rowInKKT);
    void FillMetaVerticesConstraintsInKKT(Borders& borders, std::vector<Eigen::Triplet<double>>& tripletListValues,
                                          int& rowInKKT);
    void FillConesConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow);
    void SetElementsForPARDISO(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow);

    bool constructHarmonicBasis(int conesConstraintsStartRow);
    bool calculateHarmonicBasisInPARDISO(int conesConstraintsStartRow);
    void createJmatrixInCpp();
    void computeFramesFromVectorFieldInCpp();
    bool loadPrecomputedFramesFromFile();

    bool getATPInitialValue();
    bool getTutteInitialValue(Borders& borders);
    void solveSystemToFindInitialValueEntriesOnNonMainBorder(Borders& borders, int& numMetaVerticesMainBorder,
                                                             Eigen::VectorXd& initialValue);
    void fixFirstConeInCpp(const Eigen::VectorXd& initialValue);

    bool runNewton(int conesConstraintsStartRow, GMMDenseColMatrix& RHS);
    bool getAllUVsByRHS(GMMDenseColMatrix& RHS);

    bool testResult();
    void checkLocationOfFoldoversAndPrint(std::vector<Facet_handle>& flippedTriangles, int& numFoldovers,
                                          int& numFoldsNearCones, int& numFoldsNearBorder);
    void fixCotFoldovers(std::vector<Facet_handle>& flippedTriangles);
    void putVertexInKernelUsingCVX(Vertex_handle v);
    void sendValuesToMatlabReport();
};
