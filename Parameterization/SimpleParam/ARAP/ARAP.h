#ifndef ARAP_H
#define ARAP_H

#include "Parameterization.h"

class ARAP : public Parameterization {
public:
    ARAP(Mesh& mesh0, int maxIter = 30);
    
    void parameterize() override;

private:
    // compute local 2D reference frame for each triangle
    void computeLocalFrames();
    
    // local step: find best rotation for each triangle
    void localStep();
    
    // global step: solve for vertex positions given rotations
    void globalStep();
    
    // initialise UV with Tutte (circle)
    void initTutte();
    
    // UVs in local triangle reference frame [3 faces x 2 coords]
    std::vector<std::vector<Eigen::Vector2d>> m_localRefs;
    
    // best rotations for each triangle [cos, sin] encoded as Eigen::Matrix2d
    std::vector<Eigen::Matrix2d> m_rotations;
    
    int m_maxIter;
};

#endif
