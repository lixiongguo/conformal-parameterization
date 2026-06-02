#pragma once

#include "Mesh.h"
#include <Eigen/Core>
#include <complex>
#include <vector>

// ================================================================
//  ComplexPolyField — 4‑RoSy vector field via complex polynomial
//
//  Reference:
//    "Designing N‑PolyVector Fields with Complex Polynomials"
//    Diamanti et al., SGP 2014
//
//  Algorithm:
//    1. Build per‑face local bases {X, Y, N} from triangle geometry
//    2. Boundary constraints: edge‑orthogonal direction = conj(edge)
//       raised to the 4th power
//    3. Internal smoothness: for adjacent faces f,g sharing edge e,
//       e_f^4 · z_f = e_g^4 · z_g  (complex transport equation)
//    4. Solve sparse least‑squares  AᵀA z = Aᵀb
//    5. Recover per‑face angle  θ = arg(z) / 4
//    6. Singularity detection by winding around vertices
// ================================================================

class ComplexPolyField {
public:
    explicit ComplexPolyField(Mesh& mesh, int N = 4);

    // ---- main computation ----
    void compute();

    // ---- per‑face results ----
    /// Per‑face angle in tangent plane (radians, in [0, 2π/N))
    double theta(int faceIndex) const;

    /// 3‑D direction vector for branch k (k = 0..N-1)
    Eigen::Vector3d direction(int faceIndex, int k) const;

    /// All N directions in the face's tangent plane
    std::vector<Eigen::Vector3d> allDirections(int faceIndex) const;

    // ---- singularity detection ----
    /// Winding index around this vertex (±1 = singularity, 0 = regular, -1 = boundary)
    int windingIndex(int vertexIndex) const;

    // ---- accessors ----
    int numFaces()    const { return nFaces_; }
    int numVertices() const { return nVerts_; }
    int N()           const { return N_; }

private:
    Mesh& mesh_;
    int    N_;          // N‑RoSy (default 4)
    int    nFaces_;
    int    nVerts_;

    // Per‑face local bases
    struct LocalFrame { Eigen::Vector3d X, Y, N; };
    std::vector<LocalFrame> frames_;

    // Per‑face complex representation  z = e^{i·N·θ}
    std::vector<std::complex<double>> z_;

    // Edge topology (compactly stored for linear system)
    struct InteriorEdge { int fi, fj, va, vb; };
    struct BoundaryEdge { int fi, va, vb; };
    std::vector<InteriorEdge> interiorEdges_;
    std::vector<BoundaryEdge> boundaryEdges_;

    // ---- helpers ----
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
};
