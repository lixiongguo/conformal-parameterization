#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "FastHGP.h"
#include "Utils/EigenLinearSolver.h"
#include "Utils/FramesFile.h"
#include "CGAL/CGAL_Macros.h"
#include <CGAL/Timer.h>
#include <Eigen/SparseLU>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <unordered_set>

#define M_2PI 6.28318530717958647693

using Complex = FastHGPNumerics::Complex;

void FastHGP::setSegSize(int segSize)
{
    if (segSize > 0) {
        mSegSize = segSize;
    }
}

void FastHGP::setFixCot(bool fixCot)
{
    mFixCot = fixCot;
}

bool FastHGP::run(std::string& objpath, std::string& vfPath)
{
    bool res = loadMesh(objpath, vfPath);
    if (!res) {
        std::cout << "Failed to load mesh" << std::endl;
        return false;
    }

    Borders borders(mCgalMesh);
    if (borders.genus() != 0) {
        std::cout << "The code only supports genus 0 for now" << std::endl;
        return false;
    }
    if (mHasCones && borders.numBorders() > 1) {
        std::cout << "The code does not support more than 1 border for models with cones for now" << std::endl;
        return false;
    }
    mNumOfBorderVertices = 0;
    for (int i = 0; i < borders.numBorders(); i++)
        mNumOfBorderVertices += borders.numVertices(i);

    getSettings();

    res = initialize(borders);
    if (!res) return false;

    res = runAlgorithm(borders);
    if (!res) return false;

    updatHalfedgeUVs();

    res = testResult();

    sendValuesToMatlabReport();
    calcDistortion();
    visualize();

    return true;
}

bool FastHGP::initialize(Borders& borders)
{
    if (mHasCones)
        setHalfEdgesInCgal();
    updateIndexOfHEinSystem();
    getBordersMapAndSetMetaVertices(borders);
    bool hasBorder = !mConesAndMetaMap.empty();
    if (!hasBorder && !mHasCones) {
        std::cout << "Model must have cones or border" << std::endl;
        return false;
    }
    if (mHasCones)
        getConesMap();

    mNumDOF = hasBorder ? static_cast<int>(mConesAndMetaMap.size())
                        : static_cast<int>(mConesAndMetaMap.size()) - 1;

    return prepareReducedMeshData();
}

bool FastHGP::loadMesh(std::string& objpath, std::string& vfPath)
{
    MeshBuffer m;
    CGAL::Timer loadObjTimer;
    loadObjTimer.start();
    bool res = Parser::loadOBJ(objpath.c_str(), &mMeshBuffer, &m);
    loadObjTimer.stop();
    if (!res) return false;
    std::cout << "Time to load the mesh: " << loadObjTimer.time() << " seconds.\n";

    CGAL::Timer createCGALmeshTimer;
    createCGALmeshTimer.start();
    std::vector<Point_3> pVec;
    pVec.resize(mMeshBuffer.positions.size());
    for (int i = 0; i < (int)pVec.size(); ++i)
        pVec[i] = Point_3(mMeshBuffer.positions[i][0], mMeshBuffer.positions[i][1], mMeshBuffer.positions[i][2]);
    std::vector<unsigned int> faces;
    faces.resize((int)mMeshBuffer.idx_pos.size());
    for (int i = 0; i < (int)faces.size(); ++i)
        faces[i] = mMeshBuffer.idx_pos[i];
    buildTriangleMesh(pVec, faces, mCgalMesh);
    createCGALmeshTimer.stop();
    std::cout << "Time to create CGAL mesh: " << createCGALmeshTimer.time() << " seconds.\n";

    mHasCones = !mMeshBuffer.cones.empty();

    if (mHasCones) {
        const std::string ffieldPostfix = ".ffield";
        const std::string matPostfix = ".mat";
        const std::string fframesPostfix = ".fframes";
        bool hasFrames = false;

        if (vfPath.size() >= ffieldPostfix.size()
            && vfPath.compare(vfPath.size() - ffieldPostfix.size(), ffieldPostfix.size(), ffieldPostfix) == 0) {
            mCalcFramesFromVecField = true;
            mFramesFromFile = false;
            hasFrames = Parser::loadVectorField(vfPath.c_str(), mCgalMesh);
        } else if ((vfPath.size() >= matPostfix.size()
                    && vfPath.compare(vfPath.size() - matPostfix.size(), matPostfix.size(), matPostfix) == 0)
                   || (vfPath.size() >= fframesPostfix.size()
                       && vfPath.compare(vfPath.size() - fframesPostfix.size(), fframesPostfix.size(), fframesPostfix)
                              == 0)) {
            mCalcFramesFromVecField = false;
            mFramesFromFile = true;
            mFramesFilePath = vfPath;
            hasFrames = true;
        }

        if (!hasFrames) {
            std::cout << "ffield or frames are not provided, cannot use cones" << std::endl;
            mMeshBuffer.cones.clear();
            mHasCones = false;
        } else {
            Parser::setMeshAdditionalData(mCgalMesh, mMeshBuffer);
        }
    }

    return true;
}

void FastHGP::getSettings()
{
    mSegSize = 40;
    mFixCot = true;

    if (const char* envSeg = std::getenv("FASTHGP_SEG_SIZE")) {
        const int v = std::atoi(envSeg);
        if (v > 0) {
            mSegSize = v;
        }
    }
    if (const char* envFix = std::getenv("FASTHGP_FIX_COT")) {
        mFixCot = (std::atoi(envFix) != 0);
    }
}

void FastHGP::getConesMap()
{
    for (int i = 0; i < (int)mMeshBuffer.cones.size(); ++i)
        mConesAndMetaMap.push_back(mCgalMesh.vertex(mMeshBuffer.cones[i].posidx));
}

void FastHGP::getBordersMapAndSetMetaVertices(Borders& borders)
{
    mIndicesOfMetaVerticesInUVbyGeneralIndex.resize(mCgalMesh.size_of_vertices());
    mIsMeta.assign(mCgalMesh.size_of_vertices(), false);
    getVerticesThatMustBeMeta(borders);

    mMetaVerticesInBorderByHalfEdge.resize(borders.numBorders());

    for (int i = 0; i < borders.numBorders(); i++) {
        int borderSegmentSize = mSegSize;
        std::vector<Halfedge_handle>& currBorder = mMetaVerticesInBorderByHalfEdge[i];

        std::vector<Halfedge_handle> borderPath;
        borders.getBorderHalfEdges(borders.vertex(i, 0), borderPath);

        if ((int)borderPath.size() < 3 * mSegSize) {
            borderSegmentSize = (int)borderPath.size() / 3;
        }

        int verticesInSeg = borderSegmentSize;
        for (int j = (int)borderPath.size() - 1; j >= 0; j--) {
            Halfedge_handle HD = borderPath[j];

            if (verticesInSeg >= borderSegmentSize || mIsMeta[HD->vertex()->index()]) {
                currBorder.push_back(HD);
                verticesInSeg = 1;
                mIndicesOfMetaVerticesInUVbyGeneralIndex[HD->vertex()->index()] =
                    static_cast<int>(mConesAndMetaMap.size());
                mConesAndMetaMap.push_back(HD->vertex());
                mIsMeta[HD->vertex()->index()] = true;
            } else {
                verticesInSeg++;
            }
        }
        currBorder.push_back(currBorder[0]);
    }
}

void FastHGP::getVerticesThatMustBeMeta(Borders& borders)
{
    int numBorders = borders.numBorders();
    std::vector<Vertex_handle> firstBridgeVertexInBorder(numBorders, Vertex_handle());

    for (int i = 0; i < borders.numBorders(); i++) {
        for (int j = 0; j < borders.numVertices(i); j++) {
            Vertex_handle v = borders.vertex(i, j);

            if (v->onCut()) {
                mIsMeta[v->index()] = true;
                continue;
            }

            Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
            const Mesh::Halfedge_around_vertex_circulator hEnd = h;

            CGAL_For_all(h, hEnd)
            {
                if (!h->is_border_edge()) {
                    Vertex_handle v1 = h->vertex();
                    Vertex_handle v2 = h->opposite()->vertex();
                    if (v1->is_border() && v2->is_border()) {
                        mIsMeta[v1->index()] = true;
                        mIsMeta[v2->index()] = true;
                        if (firstBridgeVertexInBorder[i] == Vertex_handle())
                            firstBridgeVertexInBorder[i] = v;
                    }
                }
            }
        }
    }

    for (int i = 0; i < borders.numBorders(); i++) {
        if (firstBridgeVertexInBorder[i] == Vertex_handle())
            continue;

        Halfedge_around_vertex_circulator hec = firstBridgeVertexInBorder[i]->vertex_begin();
        Halfedge_around_vertex_circulator hec_end = hec;

        CGAL_For_all(hec, hec_end)
        {
            if (hec->is_border()) break;
        }
        Halfedge_handle he = hec;

        do {
            Halfedge_handle he_next = he->next();
            while (!mIsMeta[he_next->vertex()->index()])
                he_next = he_next->next();
            Halfedge_handle he2 = he_next;
            while (he->vertex() != he2->opposite()->vertex() && he->next()->vertex() != he2->opposite()->vertex()) {
                he = he->next();
                he2 = he2->prev();
            }
            mIsMeta[he2->opposite()->vertex()->index()] = true;
            he = he_next;
        } while (he->vertex() != firstBridgeVertexInBorder[i]);
    }
}

bool FastHGP::prepareReducedMeshData()
{
    mVerticesByHalfedges.resize(mSizeOfSystemVar, 3);
    auto heIt = mCgalMesh.halfedges_begin();
    while (heIt != mCgalMesh.halfedges_end()) {
        int index = mIndexOfHEinSystem[heIt->index()] - 1;
        auto p = heIt->vertex()->point();
        mVerticesByHalfedges(index, 0) = p[0];
        mVerticesByHalfedges(index, 1) = p[1];
        mVerticesByHalfedges(index, 2) = p[2];
        heIt++;
    }

    std::unordered_set<Facet_handle> facesNearCones;

    for (int i = 0; i < (int)mConesAndMetaMap.size(); i++) {
        Vertex_handle v = mConesAndMetaMap[i];
        Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
        const Mesh::Halfedge_around_vertex_circulator hEnd = h;
        CGAL_For_all(h, hEnd)
            if (!h->is_border()) {
                facesNearCones.insert(h->face());
            }
    }

    const int numNearConesFaces = static_cast<int>(facesNearCones.size());
    mReducedFaces.resize(numNearConesFaces, 3);
    mReducedFacets.clear();
    mReducedFacets.reserve(numNearConesFaces);

    int i = 0;
    for (Facet_handle f : facesNearCones) {
        mReducedFacets.push_back(f);
        Halfedge_handle he[3];
        f->getHalfedges(he);
        mReducedFaces(i, 0) = mIndexOfHEinSystem[he[0]->index()];
        mReducedFaces(i, 1) = mIndexOfHEinSystem[he[1]->index()];
        mReducedFaces(i, 2) = mIndexOfHEinSystem[he[2]->index()];

        mConesAndNearConesMapOfRowsInKKT.insert(mIndexOfHEinSystem[he[0]->index()] - 1);
        mConesAndNearConesMapOfRowsInKKT.insert(mIndexOfHEinSystem[he[1]->index()] - 1);
        mConesAndNearConesMapOfRowsInKKT.insert(mIndexOfHEinSystem[he[2]->index()] - 1);
        ++i;
    }

    return true;
}

bool FastHGP::runAlgorithm(Borders& borders)
{
    GMMDenseColMatrix RHS;
    int conesConstraintsStartRow;

    CGAL::Timer totalTime;
    totalTime.start();

    constructKKTmatrix(borders, conesConstraintsStartRow);

    bool res = constructHarmonicBasis(conesConstraintsStartRow);
    if (!res) {
        std::cout << "Failed to construct Harmonic Basis" << std::endl;
        return false;
    }

    bool ATPSucsess = mHasCones ? getATPInitialValue() : getTutteInitialValue(borders);
    bool NewtonSucsess = runNewton(conesConstraintsStartRow, RHS);
    mATPandNewtonPassed = ATPSucsess && NewtonSucsess;

    totalTime.stop();
    mTotalTime = totalTime.time();
    std::cout << "Total time: " << mTotalTime << "\n";

    res = getAllUVsByRHS(RHS);
    if (!res) {
        std::cout << "Failed to Calculate UVs" << std::endl;
        return false;
    }

    return true;
}

void FastHGP::constructKKTmatrix(Borders& borders, int& conesConstraintsStartRow)
{
    std::vector<Eigen::Triplet<double>> tripletListValues;
    tripletListValues.reserve(6 * mCgalMesh.size_of_vertices() + mNumDOF);

    FillLaplacianInKKT(tripletListValues);

    conesConstraintsStartRow = mSizeOfSystemVar;
    if (mHasCones) {
        conesConstraintsStartRow = 2 * mSizeOfSystemVar;
        FillRotationConstraintsInKKT(tripletListValues, conesConstraintsStartRow);
    }

    FillMetaVerticesConstraintsInKKT(borders, tripletListValues, conesConstraintsStartRow);
    mSizeOfMatrix = mHasCones ? conesConstraintsStartRow + 2 * mNumDOF : conesConstraintsStartRow + mNumDOF;

    FillConesConstraintsInKKT(tripletListValues, conesConstraintsStartRow);
    SetElementsForPARDISO(tripletListValues, conesConstraintsStartRow);

    mKKtEigen.resize(mSizeOfMatrix, mSizeOfMatrix);
    mKKtEigen.setFromTriplets(tripletListValues.begin(), tripletListValues.end());
    mKKtEigen.makeCompressed();
}

void FastHGP::FillLaplacianInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues)
{
    for_each_facet(f, mCgalMesh)
    {
        Halfedge_handle h1 = f->halfedge();
        Halfedge_handle h2 = h1->next();
        Halfedge_handle h3 = h2->next();

        double term1 = -0.5 * h1->cot(false);
        double term2 = -0.5 * h2->cot(false);
        double term3 = -0.5 * h3->cot(false);

        updateTermInLaplacianByHalfEdge(h3, term3, tripletListValues);
        updateTermInLaplacianByHalfEdge(h1, term1, tripletListValues);
        updateTermInLaplacianByHalfEdge(h2, term2, tripletListValues);
    }
}

void FastHGP::updateTermInLaplacianByHalfEdge(Halfedge_handle h, const double& term,
                                            std::vector<Eigen::Triplet<double>>& tripletListValues)
{
    const int i = mIndexOfHEinSystem[h->prev()->index()] - 1;
    const int j = mIndexOfHEinSystem[h->index()] - 1;

    if (term != 0) {
        updateTermInLaplacianByIndices(i, j, term, tripletListValues);
        if (mHasCones) {
            updateTermInLaplacianByIndices(mSizeOfSystemVar + i, mSizeOfSystemVar + j, term, tripletListValues);
        }
    }
}

void FastHGP::updateTermInLaplacianByIndices(int i, int j, const double& term,
                                             std::vector<Eigen::Triplet<double>>& tripletListValues)
{
    tripletListValues.push_back(Eigen::Triplet<double>(i, i, -term));
    tripletListValues.push_back(Eigen::Triplet<double>(j, j, -term));
    if (j >= i) {
        tripletListValues.push_back(Eigen::Triplet<double>(i, j, term));
    } else {
        tripletListValues.push_back(Eigen::Triplet<double>(j, i, term));
    }
}

void FastHGP::FillRotationConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues, int& rowInKKT)
{
    int seamSize = (int)mRotationConstraints.nrows();
    double cosAngle[4], sinAngle[4];
    cosAngle[0] = 1; sinAngle[0] = 0;
    cosAngle[1] = 0; sinAngle[1] = -1;
    cosAngle[2] = -1; sinAngle[2] = 0;
    cosAngle[3] = 0; sinAngle[3] = 1;

    for (int i = 0; i < seamSize / 2; ++i) {
        int v0 = mIndexOfHEinSystem[(int)mRotationConstraints(i, 0)] - 1;
        int v1 = mIndexOfHEinSystem[(int)mRotationConstraints(i, 1)] - 1;
        int v3 = mIndexOfHEinSystem[(int)mRotationConstraints(i, 3)] - 1;
        int v4 = mIndexOfHEinSystem[(int)mRotationConstraints(i, 4)] - 1;
        int rot = (int)mRotationConstraints(i, 2);

        int curRow = rowInKKT + 2 * i;

        tripletListValues.push_back(Eigen::Triplet<double>(v0, curRow, 1));
        tripletListValues.push_back(Eigen::Triplet<double>(v1, curRow, -1));
        tripletListValues.push_back(Eigen::Triplet<double>(v4, curRow, -1 * cosAngle[rot]));
        tripletListValues.push_back(Eigen::Triplet<double>(v3, curRow, cosAngle[rot]));
        tripletListValues.push_back(Eigen::Triplet<double>(v4 + mSizeOfSystemVar, curRow, sinAngle[rot]));
        tripletListValues.push_back(Eigen::Triplet<double>(v3 + mSizeOfSystemVar, curRow, -1 * sinAngle[rot]));

        curRow = curRow + 1;

        tripletListValues.push_back(Eigen::Triplet<double>(v0 + mSizeOfSystemVar, curRow, 1));
        tripletListValues.push_back(Eigen::Triplet<double>(v1 + mSizeOfSystemVar, curRow, -1));
        tripletListValues.push_back(Eigen::Triplet<double>(v4, curRow, -1 * sinAngle[rot]));
        tripletListValues.push_back(Eigen::Triplet<double>(v3, curRow, sinAngle[rot]));
        tripletListValues.push_back(Eigen::Triplet<double>(v4 + mSizeOfSystemVar, curRow, -1 * cosAngle[rot]));
        tripletListValues.push_back(Eigen::Triplet<double>(v3 + mSizeOfSystemVar, curRow, cosAngle[rot]));
    }

    rowInKKT = rowInKKT + seamSize;
}

void FastHGP::FillMetaVerticesConstraintsInKKT(Borders& borders,
                                               std::vector<Eigen::Triplet<double>>& tripletListValues, int& rowInKKT)
{
    for (int k = 0; k < borders.numBorders(); k++) {
        std::vector<Halfedge_handle>& currBorder = mMetaVerticesInBorderByHalfEdge[k];
        for (int i = 0; i < (int)currBorder.size() - 1; i++) {
            Halfedge_handle metaVertexHD = currBorder[i];
            Halfedge_handle nextMetaVertexHD = currBorder[i + 1];

            Halfedge_handle h = metaVertexHD;
            double metaEdgeLength = 0.0;
            while (h != nextMetaVertexHD) {
                metaEdgeLength += h->length();
                h = h->prev();
            }

            h = metaVertexHD->prev();
            double lengthFromFirstMeta = 0.0;
            while (h != nextMetaVertexHD) {
                const int firstMetaVertex = mIndexOfHEinSystem[metaVertexHD->index()] - 1;
                const int secondMetaVertex =
                    mIndexOfHEinSystem[nextMetaVertexHD->next()->opposite()->index()] - 1;

                lengthFromFirstMeta += h->next()->length();
                const double t = (1 - lengthFromFirstMeta / metaEdgeLength);

                tripletListValues.push_back(Eigen::Triplet<double>(firstMetaVertex, rowInKKT, t));
                tripletListValues.push_back(Eigen::Triplet<double>(secondMetaVertex, rowInKKT, 1 - t));
                tripletListValues.push_back(Eigen::Triplet<double>(mIndexOfHEinSystem[h->index()] - 1, rowInKKT, -1));

                if (mHasCones) {
                    ++rowInKKT;
                    tripletListValues.push_back(Eigen::Triplet<double>(firstMetaVertex + mSizeOfSystemVar, rowInKKT, t));
                    tripletListValues.push_back(
                        Eigen::Triplet<double>(secondMetaVertex + mSizeOfSystemVar, rowInKKT, 1 - t));
                    tripletListValues.push_back(
                        Eigen::Triplet<double>(mIndexOfHEinSystem[h->index()] - 1 + mSizeOfSystemVar, rowInKKT, -1));
                }

                ++rowInKKT;
                h = h->prev();
            }
        }
    }
}

void FastHGP::FillConesConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues,
                                        int conesConstraintsStartRow)
{
    for (int i = 0; i < mNumDOF; i++) {
        int rowInKKT = conesConstraintsStartRow + i;
        int indexOfConeInSystem =
            mIndexOfHEinSystem[mConesAndMetaMap[i]->halfedge()->index()] - 1;

        tripletListValues.push_back(Eigen::Triplet<double>(indexOfConeInSystem, rowInKKT, 1));
        if (mHasCones)
            tripletListValues.push_back(
                Eigen::Triplet<double>(mSizeOfSystemVar + indexOfConeInSystem, mNumDOF + rowInKKT, 1));
    }
}

void FastHGP::SetElementsForPARDISO(std::vector<Eigen::Triplet<double>>& tripletListValues,
                                    int conesConstraintsStartRow)
{
    for (int i = 0; i < mSizeOfMatrix; i++)
        tripletListValues.push_back(Eigen::Triplet<double>(i, i, 0));

    for (int row : mConesAndNearConesMapOfRowsInKKT) {
        int row2 = row + mSizeOfSystemVar;
        for (int j = 0; j < mNumDOF; j++) {
            int col = j + conesConstraintsStartRow;
            tripletListValues.push_back(Eigen::Triplet<double>(row, col, 0));
            if (mHasCones)
                tripletListValues.push_back(Eigen::Triplet<double>(row2, col, 0));
        }
    }
}

bool FastHGP::constructHarmonicBasis(int conesConstraintsStartRow)
{
    bool res = calculateHarmonicBasisInPARDISO(conesConstraintsStartRow);
    if (!res) return false;

    if (mHasCones && mCalcFramesFromVecField) {
        computeFramesFromVectorFieldInCpp();
    } else if (mHasCones && mFramesFromFile) {
        if (!loadPrecomputedFramesFromFile()) {
            return false;
        }
    }

    createJmatrixInCpp();
    return true;
}

bool FastHGP::calculateHarmonicBasisInPARDISO(int conesConstraintsStartRow)
{
    int mtype = -2;
    EigenLinearSolver eigenSolver(mtype);

    bool res = eigenSolver.init(true);
    if (!res) return false;
    res = eigenSolver.createPardisoFormatMatrix(mKKtEigen);
    if (!res) return false;

    eigenSolver.setSelectiveRows(mConesAndNearConesMapOfRowsInKKT);

    res = eigenSolver.preprocess();
    if (!res) return false;

    res = eigenSolver.selectiveInverse();
    if (!res) return false;

    const Eigen::SparseMatrix<double, Eigen::RowMajor>& selectiveInvSol = eigenSolver.selectiveInverseMatrix();

    std::vector<Eigen::Triplet<Complex>> basisTriplets;
    basisTriplets.reserve(mConesAndNearConesMapOfRowsInKKT.size() * mNumDOF);

    for (int row : mConesAndNearConesMapOfRowsInKKT) {
        for (int j = 0; j < mNumDOF; j++) {
            int col = j + conesConstraintsStartRow;
            double realPart = selectiveInvSol.coeff(row, col);
            double imagPart = mHasCones ? selectiveInvSol.coeff(row + mSizeOfSystemVar, col) : 0.0;
            if (realPart != 0.0 || imagPart != 0.0) {
                basisTriplets.emplace_back(row, j, Complex(realPart, imagPart));
            }
        }
    }

    mHarmonicBasisEigen.resize(mSizeOfSystemVar, mNumDOF);
    mHarmonicBasisEigen.setFromTriplets(basisTriplets.begin(), basisTriplets.end());
    mHarmonicBasisEigen.makeCompressed();

    return true;
}

void FastHGP::createJmatrixInCpp()
{
    FastHGPNumerics::createJmatrix(mHarmonicBasisEigen, mVerticesByHalfedges, mReducedFaces, mJfz, mJfbz, mArea);
}

void FastHGP::computeFramesFromVectorFieldInCpp()
{
    const int numF = static_cast<int>(mReducedFacets.size());
    mFrames.resize(numF);

    for (int f = 0; f < numF; ++f) {
        Facet_handle face = mReducedFacets[f];
        const int he0 = mReducedFaces(f, 0) - 1;
        const int he1 = mReducedFaces(f, 1) - 1;
        const int he2 = mReducedFaces(f, 2) - 1;

        Eigen::Vector3d p0 = mVerticesByHalfedges.row(he0);
        Eigen::Vector3d p1 = mVerticesByHalfedges.row(he1);
        Eigen::Vector3d p2 = mVerticesByHalfedges.row(he2);

        Eigen::Vector3d e1 = p1 - p0;
        Eigen::Vector3d e2 = p2 - p0;
        double e1n = std::max(e1.norm(), 1e-14);
        double e2n = std::max(e2.norm(), 1e-14);
        e1 /= e1n;
        e2 /= e2n;

        Eigen::Vector3d n = e1.cross(e2);
        double nn = std::max(n.norm(), 1e-14);
        n /= nn;

        double vx = face->kv1().x(), vy = face->kv1().y(), vz = face->kv1().z();
        double vn = std::max(std::sqrt(vx * vx + vy * vy + vz * vz), 1e-14);
        vx /= vn; vy /= vn; vz /= vn;

        Eigen::Vector3d v(vx, vy, vz);
        Eigen::Vector3d c1 = e1.cross(v);
        double c1n = c1.norm();
        if (c1n > 1e-14) c1 /= c1n;
        double signSin = c1.dot(n);
        double cosAng = std::max(-1.0, std::min(1.0, v.dot(e1)));
        double sinAng = std::sqrt(std::abs(1.0 - cosAng * cosAng));
        mFrames(f) = Complex(cosAng, signSin * sinAng);
    }
}

bool FastHGP::loadPrecomputedFramesFromFile()
{
    const int numF = static_cast<int>(mReducedFacets.size());
    std::vector<std::complex<double>> frames;
    if (!FramesFile::load(mFramesFilePath, frames)) {
        return false;
    }
    if (static_cast<int>(frames.size()) != numF) {
        std::cout << "Precomputed frames count (" << frames.size() << ") != reduced faces (" << numF
                  << ")\n";
        return false;
    }

    mFrames.resize(numF);
    for (int f = 0; f < numF; ++f) {
        mFrames(f) = frames[static_cast<size_t>(f)];
    }
    return true;
}

bool FastHGP::getATPInitialValue()
{
    FastHGPNumerics::ATPResult atp =
        FastHGPNumerics::ATPForInitialValue(mJfz, mJfbz, mArea, mFrames);
    mATPSuccess = atp.success;
    fixFirstConeInCpp(atp.UVonCones);
    return mATPSuccess;
}

void FastHGP::fixFirstConeInCpp(const Eigen::VectorXd& initialValue)
{
    Eigen::VectorXd x;
    FastHGPNumerics::fixFirstCone(initialValue, x, mFixedIndices, mFixedValues);
    mReducedSolution.resize(mNumDOF, 2);
    Eigen::VectorXd full(2 * mNumDOF);
    full.setZero();
    std::vector<bool> isFixed(2 * mNumDOF, false);
    for (int idx : mFixedIndices) isFixed[idx] = true;
    int j = 0;
    for (int i = 0; i < 2 * mNumDOF; ++i) {
        if (!isFixed[i]) full(i) = x(j++);
    }
    for (size_t k = 0; k < mFixedIndices.size(); ++k) {
        full(mFixedIndices[k]) = mFixedValues(k);
    }
    for (int i = 0; i < mNumDOF; ++i) {
        mReducedSolution(i, 0) = full(i);
        mReducedSolution(i, 1) = full(i + mNumDOF);
    }
}

bool FastHGP::getTutteInitialValue(Borders& borders)
{
    Eigen::VectorXd initialValue(2 * mNumDOF);
    initialValue.setZero();

    Vertex_handle mainBorderVertex = borders.vertex(0, 0);

    std::vector<Halfedge_handle> pathMainBorder;
    borders.getBorderHalfEdges(mainBorderVertex, pathMainBorder);
    int numVerticesMainBorder = (int)pathMainBorder.size();

    double totalBorderLength = 0.0;
    for (int i = 0; i < numVerticesMainBorder; i++)
        totalBorderLength += pathMainBorder[i]->length();
    if (totalBorderLength <= 0.0) return false;

    double totalAngle = 0.0;
    int numMetaVerticesMainBorder = 0;
    double angle = 0.0;
    for (int i = numVerticesMainBorder - 1; i >= 0; --i) {
        totalAngle += angle;
        if (mIsMeta[pathMainBorder[i]->vertex()->index()]) {
            double x = std::cos(totalAngle);
            double y = std::sin(totalAngle);
            int indexInUV = mIndicesOfMetaVerticesInUVbyGeneralIndex[pathMainBorder[i]->vertex()->index()];
            initialValue(indexInUV) = x;
            initialValue(indexInUV + mNumDOF) = y;
            numMetaVerticesMainBorder++;
        }
        angle = 2.0 * M_PI * pathMainBorder[i]->length() / totalBorderLength;
    }

    if (borders.numBorders() != 1) {
        solveSystemToFindInitialValueEntriesOnNonMainBorder(borders, numMetaVerticesMainBorder, initialValue);
    }

    fixFirstConeInCpp(initialValue);
    mATPSuccess = true;
    return true;
}

void FastHGP::solveSystemToFindInitialValueEntriesOnNonMainBorder(Borders& borders,
                                                                  int& numMetaVerticesMainBorder,
                                                                  Eigen::VectorXd& initialValue)
{
    int numEquations = mNumDOF - numMetaVerticesMainBorder;
    Eigen::MatrixXd TutteLaplacianMatrix = Eigen::MatrixXd::Zero(numEquations, numEquations);
    Eigen::MatrixXd TutteRHS = Eigen::MatrixXd::Zero(numEquations, 2);

    std::vector<int> indicesOfNonMainBorderInAllBorder(mNumDOF, -1);
    Vertex_handle v;
    int rowInMat = 0;
    for (int i = 1; i < borders.numBorders(); i++) {
        for (int j = 0; j < borders.numVertices(i); j++) {
            v = borders.vertex(i, j);
            if (mIsMeta[v->index()]) {
                int indexInUV = mIndicesOfMetaVerticesInUVbyGeneralIndex[v->index()];
                indicesOfNonMainBorderInAllBorder[indexInUV] = rowInMat;
                rowInMat++;
            }
        }
    }

    rowInMat = 0;
    for (int i = 1; i < borders.numBorders(); i++) {
        for (int j = 0; j < borders.numVertices(i); j++) {
            v = borders.vertex(i, j);
            if (mIsMeta[v->index()]) {
                Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
                const Mesh::Halfedge_around_vertex_circulator hEnd = h;
                int valance = 0;
                CGAL_For_all(h, hEnd) { ++valance; }
                CGAL_For_all(h, hEnd)
                {
                    double w = 1.0 / (double)valance;
                    Vertex_handle neighborVertex = h->opposite()->vertex();
                    if (neighborVertex->is_border() && mIsMeta[neighborVertex->index()]) {
                        int indexInUV = mIndicesOfMetaVerticesInUVbyGeneralIndex[neighborVertex->index()];
                        TutteLaplacianMatrix(rowInMat, indicesOfNonMainBorderInAllBorder[indexInUV]) -= w;
                    } else {
                        for (int k = 0; k < mNumDOF; k++) {
                            int heRow = mIndexOfHEinSystem[neighborVertex->halfedge()->index()] - 1;
                            double hb = mHarmonicBasisEigen.coeff(heRow, k).real();
                            if (indicesOfNonMainBorderInAllBorder[k] != -1) {
                                TutteLaplacianMatrix(rowInMat, indicesOfNonMainBorderInAllBorder[k]) -= w * hb;
                            } else {
                                TutteRHS(rowInMat, 0) += w * hb * initialValue(k);
                                TutteRHS(rowInMat, 1) += w * hb * initialValue(k + mNumDOF);
                            }
                        }
                    }
                }
                TutteLaplacianMatrix(rowInMat, rowInMat) += 1.0;
                rowInMat++;
            }
        }
    }

    Eigen::MatrixXd initialValueNonMain = TutteLaplacianMatrix.fullPivLu().solve(TutteRHS);

    for (int i = 1; i < borders.numBorders(); i++) {
        for (int j = 0; j < borders.numVertices(i); j++) {
            v = borders.vertex(i, j);
            if (mIsMeta[v->index()]) {
                int indexInUV = mIndicesOfMetaVerticesInUVbyGeneralIndex[v->index()];
                initialValue(indexInUV) = initialValueNonMain(indicesOfNonMainBorderInAllBorder[indexInUV], 0);
                initialValue(indexInUV + mNumDOF) = initialValueNonMain(indicesOfNonMainBorderInAllBorder[indexInUV], 1);
            }
        }
    }
}

bool FastHGP::runNewton(int conesConstraintsStartRow, GMMDenseColMatrix& RHS)
{
    Eigen::VectorXd full(2 * mNumDOF);
    for (int i = 0; i < mNumDOF; ++i) {
        full(i) = mReducedSolution(i, 0);
        full(i + mNumDOF) = mReducedSolution(i, 1);
    }

    Eigen::VectorXd xInitial;
    FastHGPNumerics::fixFirstCone(full, xInitial, mFixedIndices, mFixedValues);

    FastHGPNumerics::NewtonResult newton = FastHGPNumerics::runNewton(
        mJfz, mJfbz, mArea, xInitial, mFixedIndices, mFixedValues, mATPSuccess);

    mReducedSolution = newton.UVonCones;

    gmm::clear(RHS);
    if (mHasCones)
        gmm::resize(RHS, mSizeOfMatrix, 1);
    else
        gmm::resize(RHS, mSizeOfMatrix, 2);

    for (int i = 0; i < mNumDOF; i++) {
        int row = conesConstraintsStartRow + i;
        RHS(row, 0) = mReducedSolution(i, 0);
        if (mHasCones)
            RHS(row + mNumDOF, 0) = mReducedSolution(i, 1);
        else
            RHS(row, 1) = mReducedSolution(i, 1);
    }

    return newton.success;
}

bool FastHGP::getAllUVsByRHS(GMMDenseColMatrix& RHS)
{
    Eigen::SparseMatrix<double> kktSym = mKKtEigen;
    for (int k = 0; k < mKKtEigen.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mKKtEigen, k); it; ++it) {
            if (it.row() != it.col()) {
                kktSym.coeffRef(it.col(), it.row()) = it.value();
            }
        }
    }
    kktSym.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(kktSym);
    if (solver.info() != Eigen::Success) {
        std::cout << "Failed to factorize KKT matrix" << std::endl;
        return false;
    }

    if (mHasCones) {
        Eigen::VectorXd rhs(RHS.nrows());
        for (int i = 0; i < RHS.nrows(); ++i) rhs(i) = RHS(i, 0);
        Eigen::VectorXd X = solver.solve(rhs);
        if (solver.info() != Eigen::Success) return false;

        mUVs.resize(mSizeOfSystemVar, 2);
        for (int i = 0; i < mSizeOfSystemVar; ++i) {
            mUVs(i, 0) = X(i);
            mUVs(i, 1) = X(i + mSizeOfSystemVar);
        }
    } else {
        Eigen::MatrixXd rhs(RHS.nrows(), 2);
        for (int i = 0; i < RHS.nrows(); ++i) {
            rhs(i, 0) = RHS(i, 0);
            rhs(i, 1) = RHS(i, 1);
        }
        Eigen::MatrixXd X = solver.solve(rhs);
        if (solver.info() != Eigen::Success) return false;

        mUVs.resize(mSizeOfSystemVar, 2);
        for (int i = 0; i < mSizeOfSystemVar; ++i) {
            mUVs(i, 0) = X(i, 0);
            mUVs(i, 1) = X(i, 1);
        }
    }

    return true;
}

bool FastHGP::testResult()
{
    int numFoldovers, numFoldsNearCones, numFoldsNearBorder, numWrongAngles, numWrongConeAngles;
    std::vector<Facet_handle> flippedTriangles;
    checkForFoldovers(flippedTriangles);

    if (mFixCot && mATPandNewtonPassed) {
        fixCotFoldovers(flippedTriangles);
        checkForFoldovers(flippedTriangles);
    }

    checkLocationOfFoldoversAndPrint(flippedTriangles, numFoldovers, numFoldsNearCones, numFoldsNearBorder);
    coneAngleDetection(numWrongAngles, numWrongConeAngles);

    if (mATPandNewtonPassed && numFoldovers == 0 && numWrongAngles == 0)
        std::cout << "result: success" << std::endl;
    else if (mATPandNewtonPassed && numFoldsNearCones == 0 && numFoldsNearBorder == 0 && numWrongConeAngles == 0)
        std::cout << "result: partial success" << std::endl;
    else {
        std::cout << "result: fail" << std::endl;
        return false;
    }
    return true;
}

void FastHGP::checkLocationOfFoldoversAndPrint(std::vector<Facet_handle>& flippedTriangles, int& numFoldovers,
                                               int& numFoldsNearCones, int& numFoldsNearBorder)
{
    std::vector<int> flippedTrianglesNearCones;
    std::vector<int> flippedTrianglesNearBorder;

    if (!flippedTriangles.empty()) {
        std::stringstream out;
        for (Facet_handle faceIt : flippedTriangles) {
            out << faceIt->index() << " ";
            Halfedge_handle h = faceIt->halfedge();
            if (faceIt->is_border_face())
                flippedTrianglesNearBorder.push_back(faceIt->index());
            if (h->vertex()->onCut() || h->next()->vertex()->onCut() || h->prev()->vertex()->onCut()) {
                if (h->vertex()->isCone() || h->next()->vertex()->isCone() || h->prev()->vertex()->isCone())
                    flippedTrianglesNearCones.push_back(faceIt->index());
            }
        }
        std::cout << out.str() << std::endl;
    }

    numFoldovers = (int)flippedTriangles.size();
    numFoldsNearCones = (int)flippedTrianglesNearCones.size();
    numFoldsNearBorder = (int)flippedTrianglesNearBorder.size();
}

void FastHGP::fixCotFoldovers(std::vector<Facet_handle>& flippedTriangles)
{
    if (flippedTriangles.empty()) return;
    std::cout << "There are " << flippedTriangles.size() << " local foldovers resulting from cot weights, fixing..."
              << std::endl;
    CGAL::Timer fixCotTime;
    fixCotTime.start();

    for (Facet_handle faceIt : flippedTriangles) {
        Vertex_handle v[3];
        faceIt->getVertices(v);
        for (int j = 0; j < 3; ++j) {
            if (!v[j]->isCone() && !v[j]->is_border() && !v[j]->onCut()) {
                double sum = oneRingAngle(v[j]);
                if (std::abs(sum - M_2PI) > 0.01) {
                    putVertexInKernelUsingCVX(v[j]);
                }
            }
        }
    }

    fixCotTime.stop();
    std::cout << "Fix cot foldovers time: " << fixCotTime.time() << "\n";
}

void FastHGP::putVertexInKernelUsingCVX(Vertex_handle v)
{
    Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
    const Mesh::Halfedge_around_vertex_circulator hEnd = h;
    int numNeighbors = 0;
    CGAL_For_all(h, hEnd) { ++numNeighbors; }

    Eigen::Vector2d UV(v->halfedge()->uv().x(), v->halfedge()->uv().y());
    Eigen::Matrix2Xd oneRing(2, numNeighbors);
    int index = numNeighbors - 1;
    CGAL_For_all(h, hEnd)
    {
        oneRing(0, index) = h->opposite()->uv().x();
        oneRing(1, index) = h->opposite()->uv().y();
        --index;
    }

    Eigen::Vector2d newUV = FastHGPNumerics::putVertexInKernel(UV, oneRing);

    Point_3 UVp = Point_3(newUV.x(), newUV.y(), 0);
    Halfedge_handle he0 = v->halfedge();
    CGAL_For_all(h, hEnd) { h->uv() = UVp; }
    mUVs(mIndexOfHEinSystem[he0->index()] - 1, 0) = newUV.x();
    mUVs(mIndexOfHEinSystem[he0->index()] - 1, 1) = newUV.y();
}

void FastHGP::sendValuesToMatlabReport()
{
    (void)0;
}

void FastHGP::visualize()
{
    // Pure C++ build: no MATLAB seam viewer (see HarmonicParametrization::visualize).
}

void FastHGP::calcDistortion()
{
    double areaWeightedK = 0.0;
    double totalArea = 0.0;
    double minK = 1e30;
    double maxK = 0.0;

    for_each_const_facet(f, mCgalMesh)
    {
        const double k = calcK(f);
        const double area = f->area();
        areaWeightedK += area * k;
        totalArea += area;
        minK = std::min(minK, k);
        maxK = std::max(maxK, k);
    }

    if (totalArea > 0.0) {
        areaWeightedK /= totalArea;
    }
    std::cout << "Distortion k: min=" << minK << " max=" << maxK << " area-weighted mean=" << areaWeightedK
              << "\n";
}

void FastHGP::coneAngleDetection(int& numWrongAngles, int& numWrongConeAngles)
{
    numWrongConeAngles = 0;
    numWrongAngles = 0;

    for (Vertex_iterator v = mCgalMesh.vertices_begin(); v != mCgalMesh.vertices_end(); v++) {
        const double sum = oneRingAngle(v);
        if (v->isCone()) {
            if (std::abs(v->getConeAngle() * M_PI - sum) > 0.01) {
                ++numWrongConeAngles;
                ++numWrongAngles;
            }
        } else if (!v->is_border()) {
            if (std::abs(sum - M_2PI) > 0.01) {
                ++numWrongAngles;
            }
        }
    }

    if (numWrongAngles > 0) {
        std::cout << "Angle check: " << numWrongAngles << " problem vertices ("
                  << numWrongConeAngles << " cones)\n";
    }
}
