#ifndef ARAP_H
#define ARAP_H

#include "Parameterization.h"

class ARAP : public Parameterization {
public:
    ARAP(Mesh& mesh0, int maxIter = 30);
    
    void parameterize() override;

    // ---- 逐次迭代的收敛指标（下标 i = 第 i 次迭代后的状态，i=0 为 Tutte 初始化）----
    // energy   : ARAP 能量 E(u) = Σ_f Σ_e w‖(u_i−u_j) − R_f(x_i−x_j)‖²（local/global 两
    //            步都在下降同一个量，故该序列单调不增，可直接作为收敛/收益曲线）
    // areaDist : 面积畸变 mean|ln((A_uv/A_3d) / r_global)|，已扣掉 UV 整体缩放
    //            （几何质量，与能量权重无关；完美等距时为 0）
    // flipped  : UV 中相对面绕序反向（折叠）的面数
    const std::vector<double>& iterationEnergy() const { return m_iterEnergy; }
    const std::vector<double>& iterationAreaDistortion() const { return m_iterAreaDist; }
    const std::vector<int>&    iterationFlipped() const { return m_iterFlipped; }

private:
    // compute local 2D reference frame for each triangle
    void computeLocalFrames();
    
    // local step: find best rotation for each triangle
    void localStep();
    
    // global step: solve for vertex positions given rotations
    void globalStep();
    
    // initialise UV with Tutte (circle)
    void initTutte();

    // 追加一条当前 uv（配合 m_rotations）的收敛指标：能量 / 面积畸变 / 翻面数
    void recordIteration();
    
    // UVs in local triangle reference frame [3 faces x 2 coords]
    std::vector<std::vector<Eigen::Vector2d>> m_localRefs;
    
    // best rotations for each triangle [cos, sin] encoded as Eigen::Matrix2d
    std::vector<Eigen::Matrix2d> m_rotations;

    // 收敛指标历史
    std::vector<double> m_iterEnergy;
    std::vector<double> m_iterAreaDist;
    std::vector<int>    m_iterFlipped;

    // recordIteration 的临时缓冲（复用容量，避免每轮重新分配）
    std::vector<double> m_areaRatios;
    
    int m_maxIter;
};

#endif
