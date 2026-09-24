#ifndef LSCM_H
#define LSCM_H

#include "Parameterization.h"

class Lscm : public Parameterization {
public:
    // constructor
    Lscm(Mesh& mesh0);
    
    // parameterize
    void parameterize() override;

    // 由调用方指定两个 pin 顶点(顶点序号,不是 2 倍下标)。
    // 仅在 v0/v1 互异且非负时生效;其余情况保持"未设置",parameterize() 会回落到自动选点。
    void setPins(int v0, int v1);

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
