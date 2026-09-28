#ifndef TUTTE_H
#define TUTTE_H

#include "Parameterization.h"

enum class TutteBoundary { CIRCLE, SQUARE, FREE };
enum class TutteWeight { COTAN, UNIFORM };

class Tutte : public Parameterization {
public:
    Tutte(Mesh& mesh0,
          TutteBoundary boundaryShape = TutteBoundary::CIRCLE,
          TutteWeight weight = TutteWeight::COTAN);
    
    void parameterize() override;

    // 关闭内部的 normalize()，保留原始调和解。
    // 供 ARAP 这类上层算法使用：它需要把「Tutte 初始 UV」和「优化后 UV」用
    // 同一个归一化变换放到同一坐标系里做对比（见 ARAP::parameterize）。
    void setNormalize(bool enable) { m_normalize = enable; }

private:
    // fix boundary vertices to circle or square
    void pinBoundary();
    
    // pin two vertices for free-boundary harmonic map (removes rigid DOFs)
    void pinTwoVertices(std::vector<int>& pinVerts);
    
    // find all boundary vertices in order
    void findBoundaryLoop(std::vector<int>& boundaryVerts) const;
    
    // solve Laplacian for vertices not in fixedVerts (cotan or uniform weights)
    void solveLaplacian(const std::vector<int>& fixedVerts);
    
    TutteBoundary m_shape;
    TutteWeight m_weight;
    bool m_normalize = true;
};

#endif
