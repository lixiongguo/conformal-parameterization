#ifndef CIRCLE_PATTERNS_WASM_H
#define CIRCLE_PATTERNS_WASM_H

#include "Parameterization.h"
#include "Utils.h"
#include <stack>

/**
 * Standalone WASM-compatible Circle Patterns parameterization.
 * 
 * This is a complete reimplementation that does NOT depend on Mosek.
 * The angle optimization step (which requires Mosek QP in the original)
 * is replaced with a direct clamped-angle computation.
 * 
 * CONE SINGULARS: Vertices can be marked as cone singulars with a target
 * angle sum (default 2π for flat, any value for cones). The curvature 
 * K_i = 2π - target_sum is concentrated at these vertices.
 * 
 * Algorithm phases:
 *   1. computeAnglesDirect() - edge weights τ_e from original + clamped angles
 *   2. computeRadii()         - convex optimization for face radii ρ = log r
 *   3. setUVs()               - BFS layout from radii + angles
 */
class CirclePatternsWasm : public Parameterization {
public:
    CirclePatternsWasm(Mesh& mesh0, int optScheme0);
    void parameterize() override;
    
    // Cone singular support: set target angle sum for a vertex
    // coneSingularIdx[i] = vertex index, coneSingularAngle[i] = target angle sum
    void setConeSingulars(const std::vector<int>& coneSingularIdx,
                          const std::vector<double>& coneSingularAngle);

protected:
    // Phase 1: compute edge weights τ_e without Mosek
    void computeAnglesDirect();
    
    // Phase 2: convex optimization for face radii
    // Energy E(ρ) = Σ_edges ImLi2Sum(...) + Σ_faces 2π·ρ_f
    void computeEnergy(double& energy, const Eigen::VectorXd& rho);
    void computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& rho);
    void computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& rho);
    void setRadii();
    bool computeRadii();
    
    // Phase 3: UV layout
    void computeAnglesAndEdgeLengths(Eigen::VectorXd& lengths);
    void performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir, 
                           Eigen::VectorXd& lengths,
                           std::unordered_map<int, bool>& visited, 
                           std::stack<EdgeCIter>& stack);
    void setUVs();

    // Member variables
    Eigen::VectorXd angles;          // per half-edge angle
    Eigen::VectorXd thetas;          // per edge weight τ_e
    Eigen::VectorXd radii;           // per face radius
    Solver solver;
    int OptScheme;
    
    // Cone singulars: vertex index → target angle sum (default 2π)
    std::unordered_map<int, double> coneSingulars;
};

#endif
