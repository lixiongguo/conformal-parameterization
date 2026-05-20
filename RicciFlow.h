#ifndef RICCI_FLOW_H
#define RICCI_FLOW_H

#include "Parameterization.h"
#include "Solver.h"
#include <stack>

/**
 * Discrete Ricci Flow parameterization using circle packing metric.
 * 
 * Theory: The Ricci flow deforms the metric proportionally to curvature:
 *   du_i/dt = K̄_i - K_i
 * where u_i = log(r_i) is the log-radius and K_i is discrete Gaussian curvature.
 * 
 * Edge length in Euclidean circle packing (inversive distance I_ij):
 *   l_ij = sqrt(e^{2u_i} + e^{2u_j} + 2·I_ij·e^{u_i+u_j})
 * 
 * Energy: E(u) = ∫ Σ_i (K_i - K̄_i) du_i  (convex)
 * Gradient: ∇E_i = K_i - K̄_i
 * Hessian: based on cotangent Laplacian
 */
class RicciFlow : public Parameterization {
public:
    // constructor
    RicciFlow(Mesh& mesh0, int optScheme0);

    // parameterize
    void parameterize() override;

protected:
    // initialize inversive distances from original edge lengths
    void initInversiveDistances();
    
    // set target curvature (flat metric: K̄=0 for interior vertices)
    void setTargetCurvature();
    
    // compute energy, gradient and hessian for solver
    void computeEnergy(double& energy, const Eigen::VectorXd& u);
    void computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& u);
    void computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& u);
    
    // edge length from circle packing metric
    double edgeLengthFromMetric(double ui, double uj, double Iij) const;
    
    // compute triangle angle at vertex given three edge lengths
    double triangleAngle(double la, double lb, double lc) const;
    
    // compute all edge lengths and angles from current u
    void computeEdgeLengthsAndAngles();
    
    // optimize radii using solver
    bool optimizeRadii();
    
    // BFS layout to determine UV coordinates
    void performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir,
                           std::unordered_map<int, bool>& visited, std::stack<EdgeCIter>& stack);
    void setUVs();

    // member variables
    std::unordered_map<int, int> index;         // vertex -> optimization variable index
    std::vector<double> inversiveDistance;       // inversive distance per edge (indexed by edge->index)
    Eigen::VectorXd Ktarget;                    // target curvature
    Eigen::VectorXd curvature;                  // current curvature (per optimization variable)
    Eigen::VectorXd edgeLengths;                // scaled edge lengths (for UV layout)
    Eigen::VectorXd halfEdgeAngles;             // angles per half-edge
    Solver solver;
    int OptScheme;
};

#endif
