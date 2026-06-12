#pragma once

#include "Mesh.h"
#include <Eigen/Core>
#include <complex>
#include <unordered_map>
#include <vector>

// ================================================================
//  ComplexPolyField — N-RoSy field via complex polynomial (N-PV)
//
//  Reference:
//    Diamanti et al., "Designing N-PolyVector Fields with Complex
//    Polynomials", SGP 2014
//
//  N-RoSy (§3.1): q_f = u_f^N, smoothness
//    bar{e}_f^N q_f = bar{e}_g^N q_g
//  Energy: sum w_fg |bar{e}_f^N q_f - bar{e}_g^N q_g|^2
// ================================================================

class ComplexPolyField {
public:
    explicit ComplexPolyField(Mesh& mesh, int N = 4);

    void compute();

    double theta(int faceIndex) const;
    Eigen::Vector3d direction(int faceIndex, int k) const;
    std::vector<Eigen::Vector3d> allDirections(int faceIndex) const;

    /// N-RoSy singularity index at vertex (units of 2*pi/N); 0 = regular, -1 = boundary
    int windingIndex(int vertexIndex) const;

    int numFaces() const { return nFaces_; }
    int numVertices() const { return nVerts_; }
    int N() const { return N_; }

    int localFaceIndex(int globalFaceIndex) const;

private:
    Mesh& mesh_;
    int N_;
    int nFaces_;
    int nVerts_;

    struct LocalFrame {
        Eigen::Vector3d X, Y, N;
    };
    std::vector<LocalFrame> frames_;
    std::unordered_map<int, int> globalToLocalFace_;

    std::vector<std::complex<double>> z_;       // unit |q|=1
    std::vector<std::complex<double>> zRaw_;    // LS solution before normalization

    struct InteriorEdge {
        int fi, fj;
        int va_f, vb_f;
        int va_g, vb_g;
    };
    struct BoundaryEdge {
        int fi, va, vb;
    };
    std::vector<InteriorEdge> interiorEdges_;
    std::vector<BoundaryEdge> boundaryEdges_;

    void buildLocalFrames();
    void buildAdjacency();
    std::complex<double> edgeComplexConjPowN(int fi, int va, int vb) const;
    void generateBoundaryConstraints(
        std::vector<Eigen::Triplet<std::complex<double>>>& trips,
        Eigen::VectorXcd& b,
        int& row) const;
    void generateInteriorEquations(
        std::vector<Eigen::Triplet<std::complex<double>>>& trips,
        int& row) const;
    void generateGaugeConstraint(
        std::vector<Eigen::Triplet<std::complex<double>>>& trips,
        Eigen::VectorXcd& b,
        int& row) const;
};
