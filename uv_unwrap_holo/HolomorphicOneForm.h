#ifndef HOLOMORPHIC_ONE_FORM_H
#define HOLOMORPHIC_ONE_FORM_H

#include "Parameterization.h"
#include <Eigen/SparseLU>

/**
 * Holomorphic 1-Form parameterization (genus-0 simplification).
 * 
 * Based on "Global Conformal Parameterization" (Gu & Yau).
 * Computes a harmonic closed 1-form ω, its Hodge-star conjugate *ω,
 * and integrates to UV coordinates.
 * 
 * Key insight: ω + i·*ω is a holomorphic 1-form → integrate → conformal map.
 * 
 * For genus-0: choose one edge as "cut" to break nullspace (|E| equation).
 * System: [closedness |F|-1; harmonicity |V|-1; cut constraint 1] · ω = rhs
 */
class HolomorphicOneForm : public Parameterization {
public:
    HolomorphicOneForm(Mesh& mesh0);
    void parameterize() override;

protected:
    void buildClosedness(Eigen::SparseMatrix<double>& A, int& rowOff);
    void buildHarmonicity(Eigen::SparseMatrix<double>& A, int& rowOff);
    void fixCutConstraint(Eigen::SparseMatrix<double>& A);
    void integrateUV(const Eigen::VectorXd& real, const Eigen::VectorXd& imag);

    // helpers
    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);
    void localFrame(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2);

    // data
    Eigen::MatrixXd V;    // vertices
    Eigen::VectorXi F;    // flattened faces
    int nV, nF, nE;
    struct Edge { int v1, v2, f1, f2; double len; };
    std::vector<Edge> edgeList;
    int cutEdge;   // index of cut edge to fix
    Eigen::MatrixXd UV; // output UV
};

#endif
