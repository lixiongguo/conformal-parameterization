#ifndef CIRCLE_PATTERNS_H
#define CIRCLE_PATTERNS_H

#include "Parameterization.h"
#include "Utils.h"
#ifdef __EMSCRIPTEN__
#include "MosekStub.h"
#else
#include "MosekSolver.h"
#endif
#include <stack>
#include <unordered_map>
#include <vector>

class CirclePatterns: public Parameterization {
public:
    // constructor
    CirclePatterns(Mesh& mesh0, int optScheme0);
    
    // parameterize
    void parameterize() override;
    
    // set cone singularities: vertex_index → target_angle_sum (e.g. 4π, 6π)
    void setConeSingulars(const std::vector<int>& coneIdx,
                          const std::vector<double>& coneAngles);
    
protected:
    // sets angle based constraints, bounds and objective function
    void setupAngleOptProblem();
        
    // sets thetas
    void setThetas();
    
    // compute angles
    virtual bool computeAngles();
    
    // computes energy, gradient and hessian
    void computeEnergy(double& energy, const Eigen::VectorXd& rho);
    void computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& rho);
    void computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& rho);

    // sets radii
    void setRadii();
    
    // compute radii
    bool computeRadii();
    
    // computes angles and edge lengths
    void computeAnglesAndEdgeLengths(Eigen::VectorXd& lengths);
    
    // determines position of unfixed face vertex
    void performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir, Eigen::VectorXd& lengths,
                           std::unordered_map<int, bool>& visited, std::stack<EdgeCIter>& stack);
    
    // sets uvs
    void setUVs();
    
    // member variables
    Eigen::VectorXd angles;
    Eigen::VectorXd thetas;
    Eigen::VectorXd radii;
    Eigen::VectorXi eIntIndices;
    int imaginaryHe;
    std::unordered_map<int, double> coneSingulars;  // vertex_index → target_angle_sum
    MosekSolver::Solver mosekSolver;
    Solver solver;
    int OptScheme;
};

#endif 
