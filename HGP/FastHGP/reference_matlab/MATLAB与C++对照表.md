# FastHGP：MATLAB 脚本与 C++ 对照表

本文档将 `reference_matlab/` 中的脚本，以及历史 MATLAB 混合主程序中的逻辑，映射到当前 C++ 实现。  
数值核心在 `FastHGPNumerics`；网格、KKT、I/O 在 `FastHGP`；选择性矩阵逆在 `EigenLinearSolver`。

相关文件：

| 类型 | 路径 |
|:---|:---|
| MATLAB 参考脚本 | `reference_matlab/*.m` |
| 桌面主流程 | `FastHGP.h` / `FastHGP.cpp` |
| 数值库 | `FastHGPNumerics.h` / `FastHGPNumerics.cpp` |
| 线性求解 | `Utils/EigenLinearSolver.h` / `.cpp` |
| WASM 简化版 | `FastHGPSimple.h` / `.cpp` |
| 共享父类（部分仍调 MATLAB） | `../HGP/HarmonicParametrization.*` |

---

## 1. 主流程对照

历史版本中，主程序在 C++ 与 MATLAB 之间交换 `FastHGP` 结构体；现由 `FastHGP::run` 串联纯 C++ 步骤。

| 阶段 | MATLAB（历史） | C++ | 说明 |
|:---|:---|:---|:---|
| 读网格 + 锥点 | C++ `Parser::loadOBJ` + MATLAB | `FastHGP::loadMesh` | 见 §4 |
| 读向量场 / frames | `HGP.cpp` + `load(matLocation)` | `Parser::loadVectorField`（`.ffield`）；`FramesFile::load`（`.mat` / `.fframes`） |
| 参数 | `FastHGP_settings.m` GUI | `FastHGP::getSettings` | 写死 `segSize=40`, `fixCot=true` |
| meta 顶点 / 边界图 | 历史混合 C++ | `getBordersMapAndSetMetaVertices`、`getVerticesThatMustBeMeta` | 无独立 `.m` |
| 缩减网格 `F_cb` | 历史混合 C++ | `prepareReducedMeshData` | 锥点及邻域三角面 |
| 组装 KKT | 历史混合 C++ | `constructKKTmatrix` 及 `Fill*` 系列 | 拉普拉斯、旋转、meta、锥约束 |
| 调和基 `H` | PARDISO 选择性逆 | `calculateHarmonicBasisInPARDISO` → `EigenLinearSolver::selectiveInverse` | 非 MKL PARDISO，见 §6 |
| `J_fz`, `J_fbz`, `Area` | `createJmatrix.m` | `createJmatrixInCpp` → `FastHGPNumerics::createJmatrix` | |
| frames | `computeFramesFromVectorField`（`HGP.cpp`） | `computeFramesFromVectorFieldInCpp` | 仅 `.ffield` 路径 |
| 初值（有锥） | `ATPForInitialValue.m` | `getATPInitialValue` → `ATPForInitialValue` | |
| 初值（仅边界） | 历史混合 Tutte | `getTutteInitialValue`、`solveSystemToFindInitialValueEntriesOnNonMainBorder` | 无独立 `.m` |
| 固定第一个锥 | `fixFirstCone.m` | `fixFirstConeInCpp` → `fixFirstCone` | |
| Newton 优化 | `Newton.m` | `runNewton` → `FastHGPNumerics::runNewton` | C++ 无 GPU |
| 回代全网格 UV | `general_solve(K,rhs)` 或等价 | `getAllUVsByRHS`（`Eigen::SparseLU`） | |
| cot 折叠修复 | `putVertexInKernel.m` + 循环 | `fixCotFoldovers`、`putVertexInKernelUsingCVX` | 无 CVX |
| 结果检查 | `FastHGP_report.m` | `testResult` | 仅 stdout，无 GUI |
| 可视化 / 畸变 | `HarmonicParametrization` + MATLAB | `visualize`、`calcDistortion` | **仍依赖 MATLAB Engine** |

入口：`FastHGP::run(objPath, vfPath)`。

---

## 2. `reference_matlab/*.m` 逐文件对照

| MATLAB 脚本 | C++ 函数 / 方法 | 所在文件 | 状态 |
|:---|:---|:---|:---|
| `ATPForInitialValue.m` | `FastHGPNumerics::ATPForInitialValue` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `globalStep_MAP.m` | `FastHGPNumerics::globalStepMAP` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `globalStep_ATP.m` | `FastHGPNumerics::globalStepATP` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `localStep.m` | `FastHGPNumerics::localStep` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `Newton.m` | `FastHGPNumerics::runNewton` | `FastHGPNumerics.cpp` | ✅ 已实现（无 `gpuArray`） |
| `symDirEnergyGradHess.m` | `FastHGPNumerics::symDirEnergyGradHess` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `symDirEnergyByfzfbz.m` | `FastHGPNumerics::symDirEnergyByfzfbz` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `lineSearchLocalInjectivity.m` | `FastHGPNumerics::lineSearchLocalInjectivity` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `lineSearchDecreasingEnergy.m` | `FastHGPNumerics::lineSearchDecreasingEnergy` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `optimizeSymDirEnergyByGlobalScaling.m` | `FastHGPNumerics::optimizeSymDirEnergyByGlobalScaling` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `fixFirstCone.m` | `FastHGPNumerics::fixFirstCone` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `createJmatrix.m` | `FastHGPNumerics::createJmatrix` | `FastHGPNumerics.cpp` | ✅ 已实现 |
| `compute_perps.m` | `FastHGPNumerics::computePerps` | `FastHGPNumerics.cpp` | ✅ 已实现（`createJmatrix` 内调用） |
| `compute_local_basis.m` | `FastHGPNumerics::computeLocalBasis` | `FastHGPNumerics.cpp` | ✅ 已实现（`computePerps` 内调用） |
| `putVertexInKernel.m` | `FastHGPNumerics::putVertexInKernel` | `FastHGPNumerics.cpp` | ✅ 已实现（半平面投影，替代 CVX） |
| `general_solve.m` | `EigenLinearSolver::solve` / `getAllUVsByRHS` 中 `SparseLU::solve` | `EigenLinearSolver.cpp`、`FastHGP.cpp` | ✅ 已由 Eigen 替代 `\` |
| `FastHGP_settings.m` | `FastHGP::getSettings` | `FastHGP.cpp` | ⚠️ GUI 移除，参数写死 |
| `FastHGP_report.m` | `testResult` + 空 `sendValuesToMatlabReport` | `FastHGP.cpp` | ⚠️ 仅控制台，无交互报告 |

### 2.1 `Newton.m` 内部调用链（均在 C++ 中）

```text
Newton.m
  ├─ symDirEnergyGradHess.m     → symDirEnergyGradHess()
  ├─ optimizeSymDirEnergyByGlobalScaling.m → optimizeSymDirEnergyByGlobalScaling()
  ├─ lineSearchLocalInjectivity.m          → lineSearchLocalInjectivity()
  └─ lineSearchDecreasingEnergy.m          → lineSearchDecreasingEnergy()
                                      └─ symDirEnergyByfzfbz()
```

桌面封装：`FastHGP::runNewton` 调用 `runNewton`，再组装 RHS 供 `getAllUVsByRHS`。

### 2.2 `ATPForInitialValue.m` 内部调用链

```text
ATPForInitialValue.m
  ├─ globalStep_MAP.m    → globalStepMAP()
  ├─ localStep.m         → localStep()        （每轮迭代）
  └─ globalStep_ATP.m    → globalStepATP()    （每轮迭代）
```

桌面封装：`FastHGP::getATPInitialValue` → `ATPForInitialValue` → `fixFirstConeInCpp`。

### 2.3 `createJmatrix.m` 内部调用链

```text
createJmatrix.m
  └─ compute_perps.m
       └─ compute_local_basis.m
```

桌面封装：`FastHGP::createJmatrixInCpp`（输入 `mHarmonicBasisEigen`、`mVerticesByHalfedges`、`mReducedFaces`）。

---

## 3. 辅助函数（无独立 `.m`，在 C++ 中新增）

| C++ 函数 | 作用 | 对应 MATLAB 习惯用法 |
|:---|:---|:---|
| `FastHGPNumerics::expandReducedX` | 将缩减变量扩回含固定分量的全长向量 | `Newton.m` 末尾 `solution(setdiff(...))` |
| `FastHGPNumerics::packReducedX` | 从全长向量抽出自由分量 | `fixFirstCone` 中 `setdiff` |
| `EigenLinearSolver::selectiveInverse` | 对 KKT 指定行做选择性逆 | 历史 PARDISO 选择性逆求 `HarmonicBasis` |
| `EigenLinearSolver::createPardisoFormatMatrix` | 转为行主序稀疏格式 | 兼容原 PARDISO 数据布局 |

---

## 4. 与 `HGP/` 共享、但 FastHGP 管线用到的逻辑

| 功能 | 原 MATLAB / HGP | FastHGP C++ |
|:---|:---|:---|
| 读 OBJ + 锥点标记 | `Parser::loadOBJ` | 同左 |
| 读 `.ffield` 向量场 | `HGP.cpp` + `Parser::loadVectorField` | `FastHGP::loadMesh` |
| 读 `.mat` 预计算 `frames` | `HGP.cpp`：`load(matLocation); frames` | `FramesFile::load` + `loadPrecomputedFramesFromFile` |
| 边界 / genus | `Borders` | `FastHGP::run`、`initialize` |
| 折叠检测 | `HarmonicParametrization::checkForFoldovers` | 继承使用 |
| 锥角检测 | `HarmonicParametrization::coneAngleDetection` | 继承（**写 MATLAB**） |
| 可视化 seam | `HarmonicParametrization::visualize` | 继承（**读 `visMatlab`**） |
| 畸变 `k` | `HarmonicParametrization::calcDistortion` | 继承（**写 MATLAB**） |

---

## 5. WASM `FastHGPSimple` 与 MATLAB 的关系

`FastHGPSimple` **不调用** `reference_matlab/` 中任何脚本，是独立简化实现：

| 桌面 FastHGP（MATLAB 管线） | WASM `FastHGPSimple` |
|:---|:---|
| KKT + 调和基 | — |
| `ATPForInitialValue.m` | —（改用 LSCM 初值） |
| `Newton.m` + 线搜索 | 梯度下降 + 自有 `localInjectivityStep` / `decreasingEnergyStep` |
| `createJmatrix.m` | 逐三角局部 `invLocal`（`buildTriangleData`） |
| 对称 Dirichlet 能量 | `symmetricDirichletEnergy`（概念相同，实现独立） |

---

## 6. MATLAB 与 C++ 的已知差异

| 项目 | MATLAB | C++ |
|:---|:---|:---|
| Newton 迭代 | 可选 `gpuArray` 加速 | 仅 CPU + Eigen |
| `putVertexInKernel` | 曾用 CVX QP | 半平面投影迭代（`putVertexInKernel`） |
| 调和基求解 | Intel PARDISO 选择性逆 | Eigen `SparseLU` / `SimplicialLDLT`（函数名仍含 `PARDISO`） |
| ATP / Newton 失败标志 | `FastHGP.ATPSucsess` 等写入 struct | `mATPSuccess`、`NewtonResult::success` |
| 报告 / 设置 | GUIDE 图形界面 | 写死参数 + `std::cout` |
| `.mat` frames | `load` + 变量 `frames` | `FramesFile::load`（v5/v7.2；v7.3 请用 `-v7` 或 `.fframes`） |

---

## 7. `FastHGP` 结构体字段 → C++ 成员速查

| MATLAB `FastHGP.*` 字段（常见） | C++ 成员 / 变量 |
|:---|:---|
| `J_fz`, `J_fbz` | `mJfz`, `mJfbz` |
| `Area` | `mArea` |
| `frames` | `mFrames` |
| `HarmonicBasis` | `mHarmonicBasisEigen` |
| KKT 矩阵 | `mKKtEigen` |
| `UVonCones` | `mReducedSolution` |
| `FixedIndices`, `FixedValues` | `mFixedIndices`, `mFixedValues` |
| `F`（缩减面） | `mReducedFaces`, `mReducedFacets` |
| `V`（半边坐标） | `mVerticesByHalfedges` |
| `ATPSucsess` | `mATPSuccess` |
| `numB` / meta 相关 | `mConesAndMetaMap`, `mIsMeta`, `mNumDOF` |

---

## 8. 阅读顺序建议

1. `FastHGP::run` — 总流程  
2. `runAlgorithm` — KKT → 初值 → Newton → 回代  
3. `FastHGPNumerics::ATPForInitialValue` / `runNewton` — 与 `.m` 逐行对照  
4. `reference_matlab/` — 算法细节与论文公式  
5. 博客：`_posts/1.Parameterization/3.几何优化方法/全局调和参数化FastHGP.md`

未实现项与后续待办见 [README.md](README.md)。
