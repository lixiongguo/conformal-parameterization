#ifndef CETM_H
#define CETM_H

#include "ConeParameterization.h"
#include "Solver.h"
#include <unordered_map>
#include <vector>

class Cetm: public ConeParameterization {
public:
    // constructor
    Cetm(Mesh& mesh0, int optScheme0);

    // parameterize
    void parameterize() override;
    
    // setConeSingulars() is inherited from ConeParameterization
    
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
    
    // sets uvs
    void setUVs();
    
    // member variables
    std::unordered_map<int, int> index;
    // coneSingulars is inherited from ConeParameterization
    Eigen::VectorXd thetas;
    Eigen::VectorXd lengths;
    Eigen::VectorXd angles;
    Solver solver;
    int OptScheme;
};

#endif
