#ifndef LSCM_H
#define LSCM_H

#include "Parameterization.h"

class Lscm : public Parameterization {
public:
    // constructor
    Lscm(Mesh& mesh0);
    
    // parameterize
    void parameterize() override;

    // 上一次 parameterize() 是否成功完成。
    // 失败时 UV 结果是无效的,调用方应据此返回错误码,而不是继续 extractUV()。
    bool succeeded() const { return lastRunOk; }

private:
    // pin vertices
    void pinVertices();
    
    // checks if vertex is pinned
    bool isPinnedVertex(const int& vIndex, int& shift, Eigen::Vector2d& pinnedPosition) const;
    
    // builds mass matrix; returns false if any row/column index would be out of range
    bool buildMassMatrix(Eigen::SparseMatrix<double>& M, Eigen::VectorXd& b) const;
    
    // set uvs
    void setUvs(const Eigen::VectorXd& x);
    
    // member variables
    std::vector<Eigen::Vector2d> pinnedPositions;
    std::vector<int> pinnedVertices;
    bool lastRunOk = false;
};

#endif
