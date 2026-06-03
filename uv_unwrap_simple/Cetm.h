#ifndef CETM_H
#define CETM_H

#include "Parameterization.h"
#include "Utils.h"
#include <stack>
#include <unordered_map>
#include <vector>

class Cetm: public Parameterization {
public:
    // constructor
    Cetm(Mesh& mesh0, int optScheme0);

    // parameterize
    void parameterize() override;
    
    // set cone singularities: vertex_index → target_angle_sum (e.g. 4π, 6π)
    void setConeSingulars(const std::vector<int>& coneIdx,
                          const std::vector<double>& coneAngles);
    
protected:
    // sets target theta values
    void setTargetThetas();
    
    // computes energy, gradient and hessian
    void computeEnergy(double& energy, const Eigen::VectorXd& u);
    void computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& u);
    void computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& u);
    
    // sets edge lengths
    void setEdgeLengthsAndAngles();
    
    // computes scale factors
    bool computeScaleFactors();
    
    // determines position of unfixed face vertex
    void performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir, 
                           std::unordered_map<int, bool>& visited, std::stack<EdgeCIter>& stack);
    
    // sets uvs
    void setUVs();
    
    // member variables
    std::unordered_map<int, int> index;
    std::unordered_map<int, double> coneSingulars;  // vertex_index → target_angle_sum
    Eigen::VectorXd thetas;
    Eigen::VectorXd lengths;
    Eigen::VectorXd angles;
    Solver solver;
    int OptScheme;
};

#endif
