# HGP：MATLAB 脚本与 C++ 对照表

本文档描述 **Harmonic Global Parameterization (BCW17)** 桌面实现 `HGP/` 的迁移状态。  
HGP 与 FastHGP 是**不同算法**：HGP 用 SOCP（CVX/MOSEK）最小化 $\|Lz\|^2$ + Lipman 约束；FastHGP 用子空间 + ATP + Newton（见 `../FastHGP/`）。

相关文件：

| 类型 | 路径 |
|:---|:---|
| 主类 | `HGP/HGP.h`、`HGP/HGP.cpp` |
| 锥点 frame 修复 | `HGP/FrameFixingClass.h`、`.cpp` |
| 公共父类 | `HarmonicParametrization.h`、`.cpp` |
| 网格 / 边界 | `HGP/Parser.*`、`HGP/Borders.*` |
| MATLAB 交换 | `HGP/MatlabInterface.*`、`HGP/MatlabGMMDataExchange.*` |
| **唯一保留的 HGP 求解脚本** | `HGP/HGP_iteration.m`（运行副本）；归档见 `reference_matlab/HGP_iteration.m` |
| 博客解读 | `_posts/.../全局调和参数化.md` |

---

## 1. 架构说明（重要）

HGP **不是**纯 C++ 移植，而是 **C++ 组装约束与数据 ↔ MATLAB 求解** 的混合管线：

```text
HGP::run (C++)
  ├─ 读网格、cut、seam、锥点
  ├─ 组装 W、rotMatrix、Grad、frames、BVFaces … → MATLAB struct HGP.*
  └─ while (迭代)
        ├─ setHarmonic* / computeGradients / frames  (C++)
        ├─ HGP_iteration                          (MATLAB + CVX + MOSEK)  ← 核心
        ├─ 读回 HGP.UV                              (C++)
        └─ FrameFixingClass                         (C++，部分结果仍写 MATLAB)
```

删掉的多份 `.m` 是**辅助脚本**；其逻辑已迁入 `HGP.cpp` / `FrameFixingClass.cpp`，**不**表示算法被砍掉——只要 `HGP_iteration.m` 仍在 MATLAB path 上且 CVX/MOSEK 可用，混合管线可闭合运行。

---

## 2. 主流程对照

| 阶段 | 原 MATLAB / 混合流程 | 现 C++ | 状态 |
|:---|:---|:---|:---|
| 读 OBJ、向量场 / `.mat` frames | `loadMesh` + `load(matLocation)` | `HGP::loadMesh` | ✅ C++ 读 OBJ；`.mat` 仍调 MATLAB `load` |
| 导出 `HGP.V,F,cones,…` | `transferMeshToMatlab` | `HGP::transferMeshToMatlab` | ✅ C++ 写 MATLAB workspace |
| cut 后半边系统索引 | `setHalfEdgesMap` | `HGP::setHalfEdgesMap` | ✅ |
| Seam 旋转约束 | `setRotationsConstraints` | `HGP::setRotationsConstraints` → `HGP.rotMatrix` | ✅ |
| 边界/锥邻域面 `BVFaces` | `setBoundaryFaces` | `HGP::setBoundaryFaces` | ✅ |
| 迭代中更新边长权重 | `updateHalfEdgesMetric` | `HGP::updateHalfEdgesMetric` | ✅ |
| 内部顶点调和约束 | `setHarmonicInternalVerticesConstraints` | 同名 → `HGP.W` | ✅ |
| Seam 顶点调和约束 | `setHarmonicSeamVerticesConstraints` | 同名 → `HGP.WSeam*` | ✅ |
| 梯度算子 $D$ / `Grad` | `computeGrads.m` | `HGP::computeGradientsInCpp` | ✅ 已删 `.m` |
| 向量场 → frames | `getFramesFromVectorField` 等 | `computeFramesFromVectorFieldInCpp` | ✅ |
| 迭代中更新 frames | `setFrames1/2` 等 | `updateFramesFromCurrentFzInCpp` | ✅ |
| **核心 SOCP 优化** | **`HGP_iteration.m`** | `EvalToCout("HGP_iteration")` | ❌ **仍在 MATLAB** |
| 锥点 frame 修复 ABF | `frameFixABF.m` | `FrameFixingClass::extractOneRingAngles` | ✅ 已删 `.m` |
| 锥点局部嵌入 | （原混合脚本） | `FrameFixingClass::setLocalEmbedding` | ✅ C++ 计算，写 `HGP.FrameFix.*` |
| UV 写回半边 | — | `HarmonicParametrization::updatHalfedgeUVs` | ✅ |
| 折叠 / 锥角 / 畸变 | — | `checkForFoldovers`、`coneAngleDetection`、`calcDistortion` | ✅（结果写 MATLAB） |
| 可视化 | `HGP_report` / seam 显示 | `HarmonicParametrization::visualize` | ⚠️ 依赖 MATLAB + `visMatlab` |
| 参数 GUI | `HGP_settings.m` | `HGP::run` 写死 `maxIt=5` | ⚠️ GUI 已删 |

入口：`HGP::run(objPath, vfPath)`。

---

## 3. 已删除 MATLAB 脚本 ↔ C++ 映射

以下脚本已从 `MatlabScripts/HGP/` 移出运行路径，**归档于 `reference_matlab/`**（逻辑已迁入 C++）：

| 原 MATLAB 脚本（典型名称） | C++ 实现 | 说明 |
|:---|:---|:---|
| `computeGrads.m` | `HGP::computeGradientsInCpp` | 输出 `HGP.Grad.m1,m2,S,ti,tj,tk` |
| `getFramesFromVectorField.m` | `HGP::computeFramesFromVectorFieldInCpp` | 首次迭代，`.ffield` 路径 |
| `setFrames1.m` / `setFrames2.m` | `HGP::updateFramesFromCurrentFzInCpp` | 后续迭代，$f_z/|f_z|$ |
| `frameFixABF.m` | `FrameFixingClass::extractOneRingAngles` | CVX 小 QP 改为 C++ 单纯形投影 |
| `HGP_settings.m` | `HGP::run` 内常量 | `maxIt=5`，`useFrameFixing=true`，`visMatlab=0` |
| `HGP_report.m` | `HGP.Result.*` 写 workspace | 无交互 GUI |

---

## 4. 仍保留 / 仍依赖 MATLAB 的部分

### 4.1 `HGP_iteration.m`（不可替代的核心）

| | |
|:---|:---|
| **文件** | `HGP/HGP_iteration.m` |
| **调用** | `HGP::run` → `MatlabInterface::EvalToCout("HGP_iteration")` |
| **求解** | CVX：`minimize norm(HarmonicMat2 * UV)` + 旋转约束 + Lipman 线性化约束 |
| **依赖** | MATLAB Engine、CVX、MOSEK（`MSK_IPAR_PRESOLVE_USE`） |
| **迁移动机** | 纯 C++、无许可证、CI 构建 |
| **难点** | SOCP/QP 规模随系统变量数增长；需 MOSEK C API 或等价求解器 |

### 4.2 其他 MATLAB 钩子

| 调用位置 | MATLAB 用途 |
|:---|:---|
| `loadMesh`（`.mat` frames） | `load(matLocation); HGP.frames=frames` |
| `HGP::run` | `HGP.FrameFix.use = 0/1` |
| `transferMeshToMatlab` | `HGP.F=HGP.F+1`（1-based 索引） |
| `setBoundaryFaces` | `HGP.BVFaces+1; unique` |
| `HarmonicParametrization::visualize` | 读 `visMatlab`，推送 seam 数据 |
| `calcDistortion` / `coneAngleDetection` | 写 `HGP.Result.k`、问题顶点数组 |
| `FrameFixingClass` | 写 `HGP.FrameFix.badOneRing` 等（供调试或下游脚本） |

---

## 5. `HGP_iteration.m` 与 C++ 变量对照

`HGP_iteration.m` 从 workspace 读取的 `HGP.*` 字段，主要由下列 C++ 函数写入：

| `HGP.*` 字段 | 写入函数 |
|:---|:---|
| `V`, `F`, `conesIndices`, `vectorFieldKv1` | `transferMeshToMatlab` |
| `halfEdges` | `setHalfEdgesMap` |
| `rotMatrix` | `setRotationsConstraints` |
| `BVFaces` | `setBoundaryFaces` |
| `W`, `rowsToTake` | `setHarmonicInternalVerticesConstraints` |
| `WSeamX/Y`, `WSeamU1/U2`, `WSeamV1/V2`, `coneIndices`, `rNum` | `setHarmonicSeamVerticesConstraints` |
| `Grad.m1,m2,S,ti,tj,tk` | `computeGradientsInCpp` |
| `frames` | `computeFramesFromVectorFieldInCpp` 或 MATLAB `load` |
| `hN` | `setFramesInMatlab`（首次迭代） |
| `FrameFix.status` | `setFramesInMatlab`（后续迭代） |
| `UV`（输出） | `HGP_iteration.m` 写入 → C++ `GetEngineDenseMatrix("HGP.UV")` |

---

## 6. 父类 `HarmonicParametrization`（HGP 与 FastHGP 共用）

| 方法 | HGP 构建 | FastHGP 构建（`FASTHGP_STANDALONE`） |
|:---|:---|:---|
| `updatHalfedgeUVs` | ✅ C++ | ✅ C++ |
| `checkForFoldovers` | ✅ C++ | ✅ C++ |
| `visualize` | MATLAB | FastHGP override 为空 |
| `calcDistortion` | MATLAB | FastHGP override 打印统计 |
| `coneAngleDetection` | MATLAB | FastHGP override 打印统计 |

HGP 桌面构建**未**定义 `FASTHGP_STANDALONE`，故仍走父类 MATLAB 实现。

---

## 7. 与 FastHGP 的分工

| | **HGP** | **FastHGP** |
|:---|:---|:---|
| 论文 | BCW17 全局调和 + Lipman SOCP | Hefetz et al. 子空间 + 对称 Dirichlet |
| 核心求解 | `HGP_iteration.m` | `FastHGPNumerics`（ATP、Newton） |
| MATLAB | **必需** | 可选（`FASTHGP_STANDALONE`） |
| 典型用途 | 原始高精度 HGP 复现 | 快速桌面 / WASM 简化版 |

二者共享：`Parser`、`Borders`、`HarmonicParametrization` 网格与验证工具。

---

## 8. 完整性结论

| 维度 | 是否完整 |
|:---|:---|
| 原 **C++ + MATLAB 混合** HGP 管线 | ✅ 基本完整（需 MATLAB + CVX + MOSEK + `HGP_iteration.m`） |
| **纯 C++**、无 MATLAB 可运行 | ❌ 不完整（缺 SOCP 核心） |
| 删掉辅助 `.m` 后功能是否缺失 | ✅ 否（逻辑已迁入 C++） |
| 无 GUI / 无报告面板 | ⚠️ `HGP_settings` / `HGP_report` 已移除或弱化 |
| 独立 CMake / 一键构建 | ❌ 未提供（需手动链 CGAL、GMM、MATLAB） |

---

## 9. 后续开发速查

见 [README.md](README.md) 中的待办清单。优先级最高项：**将 `HGP_iteration.m` 迁为 C++ SOCP/QP（MOSEK 或开源替代）**。
