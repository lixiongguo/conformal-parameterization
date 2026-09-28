#ifndef ARAP_H
#define ARAP_H

#include "Parameterization.h"
#include "Tutte.h"

class ARAP : public Parameterization {
public:
    // Tutte 现在是 ARAP 的前置（初始化）步骤，不再是独立算法：
    //   boundaryShape / weight 用来配置这一步的边界形状与调和权重。
    ARAP(Mesh& mesh0, int maxIter = 30,
         TutteBoundary boundaryShape = TutteBoundary::CIRCLE,
         TutteWeight   weight        = TutteWeight::COTAN);

    void parameterize() override;

    // ---- Tutte 初始化 UV ----
    // ARAP 迭代前的 UV（即 Tutte 结果），扁平数组 [u0,v0,u1,v1,...]，长度 = 2 * 顶点数。
    // 与 parameterize() 结束后的 mesh UV 处于**同一坐标系**（同一归一化变换），
    // 因此可以在页面上叠在一起做「初始 / 优化后」对比。
    const std::vector<double>& initialUv() const { return m_initUv; }

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
    
    // initialise UV with Tutte (boundary shape / weight from the constructor)
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

    // Tutte 初始化的配置
    TutteBoundary m_initBoundary;
    TutteWeight   m_initWeight;

    // Tutte 初始化 UV 快照（扁平 2N，与最终 mesh UV 同坐标系）
    std::vector<double> m_initUv;

    int m_maxIter;
};

#endif
