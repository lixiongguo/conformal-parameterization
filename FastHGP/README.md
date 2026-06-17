# FastHGP

快速全局调和参数化（Fast Harmonic Global Parameterization）的桌面版与 WASM 实现。

## 组件

| 类 / 文件 | 平台 | 说明 |
|:---|:---|:---|
| **`FastHGP`**（`FastHGP.h`, `FastHGP.cpp`） | 桌面（CGAL） | 完整管线：KKT 调和基、ATP/Tutte 初值、对称 Dirichlet 能量上的 Newton、cot 折叠修复 |
| **`FastHGPSimple`**（`FastHGPSimple.h`, `.cpp`） | WASM | LSCM 初值 + 对称 Dirichlet 梯度下降 + 单射线搜索 |
| **`FastHGPNumerics`**（`FastHGPNumerics.h`, `.cpp`） | 桌面 | `reference_matlab/*.m` 数值部分的纯 C++ 移植（无 MATLAB） |
| **`Utils/EigenLinearSolver`** | 桌面 | KKT 矩阵的选择性逆，用于构造调和基 |

## 桌面版 `FastHGP` 依赖

- `../HGP/` 中的 CGAL 网格、`Borders`、`Parser`
- Eigen（稠密 + 稀疏）
- GMM：仅在父类 `HarmonicParametrization` 需要处使用（`mUVs`、`mRotationConstraints`）
- **父类中仍有调用 MATLAB Engine 的代码**（见下文 [继承自 HarmonicParametrization 的 MATLAB 钩子](#4-继承自-harmonicparametrization-的-matlab-钩子)）

编译时加入 include 路径 `-I../HGP`（以及 CGAL）。入口：`FastHGP::run(objPath, vfPath)`。

目前 **没有** 独立的 CMake 目标；需手动链接 `FastHGP.cpp`、`FastHGPNumerics.cpp`、`Utils/EigenLinearSolver.cpp` 及所需的 `../HGP/*` 源文件。

## WASM `FastHGPSimple`

通过 `cpp/build_wasm.ps1 -Target uv_unwrap_simple` 或 `wasm/Makefile` 构建（`bindings_uv_unwrap_simple_compat.cpp` 中的 `solve_fasthgp`）。

| 桌面 `FastHGP` | WASM `FastHGPSimple` |
|:---|:---|
| KKT + 调和基 + meta 顶点 | — |
| Seam 旋转约束（锥点） | — |
| ATP / Tutte 初值 + Newton | LSCM 初值 + 梯度下降 |
| Cot 折叠修复（`putVertexInKernel`） | 仅单射线搜索 |
| CGAL `Mesh` | `BaseMesh` |

完整桌面 `FastHGP` **未** 编入 WASM（需要 CGAL 桌面构建）。

## 参数设置（桌面，无 GUI）

默认值在 `FastHGP::getSettings()` 中 **写死**（无配置文件 / 命令行）：

- `segSize = 40`（边界 meta 顶点间距）
- `fixCot = true`（cot 权重引起的局部折叠后，将内点投影修复）

原 MATLAB GUI：`reference_matlab/FastHGP_settings.m`（`segSize`、`fixCot`、`visMatlab`）。

---

## 已实现（桌面核心）

| 功能 | C++ 位置 | MATLAB 参考 |
|:---|:---|:---|
| KKT 组装 + 选择性逆 → 调和基 | `constructKKTmatrix`、`calculateHarmonicBasisInPARDISO`、`EigenLinearSolver` | `general_solve.m`（已由 Eigen 分解替代） |
| `createJmatrix`、缩减网格 `F_cb` | `createJmatrixInCpp`、`prepareReducedMeshData` | `createJmatrix.m`、`compute_perps.m` |
| 从 `.ffield` 向量场计算 frames | `computeFramesFromVectorFieldInCpp` | 与 `HGP.cpp` 中逻辑相同 |
| ATP 初值（有锥点） | `getATPInitialValue` → `FastHGPNumerics::ATPForInitialValue` | `ATPForInitialValue.m` |
| Tutte 初值（仅边界） | `getTutteInitialValue` | 历史混合代码内联 |
| 对称 Dirichlet 上的 Newton | `runNewton` → `FastHGPNumerics::runNewton` | `Newton.m` |
| 能量 / 梯度 / Hessian、线搜索、全局缩放 | `FastHGPNumerics.*` | `symDirEnergy*.m`、`lineSearch*.m`、`optimizeSymDirEnergyByGlobalScaling.m` |
| `fixFirstCone` | `fixFirstConeInCpp` | `fixFirstCone.m` |
| Cot 折叠修复 | `fixCotFoldovers`、`putVertexInKernelUsingCVX` | `putVertexInKernel.m`（CVX 已由 `FastHGPNumerics::putVertexInKernel` 替代） |
| 结果校验 | `testResult`、折叠 / 锥角检测 | `FastHGP_report.m`（仅控制台输出） |

除 GUI / 报告脚本外，`reference_matlab/*.m` 的数值部分已全部迁入 `FastHGPNumerics.cpp`。

---

## 未实现 / 待办

显式占位使用 `NotImplemented.h` 中的 `FASTHGP_NOT_IMPLEMENTED(msg)`。

### 1. 从 `.mat` 加载预计算 frames（阻断不用 `.ffield` 的锥点流程）

| | |
|:---|:---|
| **触发条件** | `vfPath` 以 `.mat` 结尾，且 OBJ 含锥点 |
| **代码位置** | `FastHGP::loadMesh` → 抛出 `FASTHGP_NOT_IMPLEMENTED` |
| **原版实现** | `HGP.cpp`（约 L142–157）：`load(matLocation); HGP.frames = frames;` |
| **临时方案** | 传入 `.ffield` 文件，走 `Parser::loadVectorField` + `computeFramesFromVectorFieldInCpp` |

**实现步骤：**

1. 增加读取器（MAT v7 可用 [matio](https://github.com/tbeu/matio)、HDF5，或将 frames 导出为简单文本/二进制格式）。
2. 读入长度为 `|F_cb|` 的复数向量 `frames`（每个 `mReducedFacets` / `mReducedFaces` 中的三角形一个 frame，顺序与 `computeFramesFromVectorFieldInCpp` 一致）。
3. 设置 `mCalcFramesFromVecField = false`，填入 `mFrames`；在 `constructHarmonicBasis` 中增加分支，当 frames 已从文件加载时跳过 `computeFramesFromVectorFieldInCpp`：

```cpp
// constructHarmonicBasis — 当前仅在 mCalcFramesFromVecField 时计算 frames
if (mHasCones && mCalcFramesFromVecField) {
    computeFramesFromVectorFieldInCpp();
}
// TODO: else if (mHasCones && mFramesLoadedFromFile) { /* mFrames 已就绪 */ }
```

### 2. 设置 GUI / 运行时配置

| | |
|:---|:---|
| **已移除** | `FastHGP_settings.m` 模态 GUI |
| **当前代码** | `getSettings()` 固定 `mSegSize=40`、`mFixCot=true` |
| **待做** | 命令行参数、JSON/INI 配置文件，或在 `run()` 前提供 setter |

### 3. `FastHGP_report` 交互式报告

| | |
|:---|:---|
| **已移除** | `reference_matlab/FastHGP_report.m` |
| **当前代码** | `sendValuesToMatlabReport()` 为空（`FastHGP.cpp` 末尾） |
| **现状** | `testResult()` 向 stdout 打印 `success` / `partial success` / `fail` |
| **待做** | 可选 HTML/JSON 报告，或接入自有查看器 |

### 4. 继承自 `HarmonicParametrization` 的 MATLAB 钩子

`FastHGP::run` 末尾仍调用父类方法，**需要 MATLAB Engine**（除非加守卫或重写）：

| `FastHGP::run` 中的调用 | 父类方法 | MATLAB 用途 |
|:---|:---|:---|
| `visualize()` | `HarmonicParametrization::visualize` | 读取 workspace 变量 `visMatlab`；非零时通过 `MatlabGMMDataExchange` 推送 seam/锥点数据 |
| `calcDistortion()` | `HarmonicParametrization::calcDistortion` | 写入 `FastHGP.Result.k` |
| `coneAngleDetection()` | `HarmonicParametrization::coneAngleDetection` | 角误差时写入问题顶点数组 |

**实现纯 C++、无 MATLAB 的桌面构建：**

- 在 `FastHGP` 中 override 上述三个方法，或
- 在 `HarmonicParametrization.cpp` 中加 `#ifdef FASTHGP_NO_MATLAB` 空实现，或
- 运行前在 MATLAB base workspace 设 `visMatlab = 0`（仅跳过 visualize 的重路径）。

`sendValuesToMatlabReport` 已剥离；上述三处调用尚未处理。

### 5. 桌面构建集成

- `FastHGP/` 下无 `CMakeLists.txt` / VS 工程。
- WASM 构建脚本仅编译 `FastHGPSimple.cpp`。

### 6. 网格 / 算法限制（主动检查，非占位）

| 限制 | 位置 | 提示信息 |
|:---|:---|:---|
| 仅支持 genus 0 | `FastHGP::run` | `"The code only supports genus 0 for now"` |
| 锥点 + 多条边界 | `FastHGP::run` | `"does not support more than 1 border for models with cones"` |
| 必须有边界或锥点 | `initialize` | `"Model must have cones or border"` |
| 有锥点但缺少 `.ffield` / `.mat` | `loadMesh` | 清除锥点，退化为仅边界（`mHasCones = false`） |

### 7. 命名遗留（已实现，但易误解）

| 名称 | 实际实现 |
|:---|:---|
| `calculateHarmonicBasisInPARDISO`、`SetElementsForPARDISO` | **Eigen** `SparseLU` / `SimplicialLDLT` 选择性逆（`EigenLinearSolver.cpp`），**不是** Intel MKL PARDISO |
| `putVertexInKernelUsingCVX` | `FastHGPNumerics::putVertexInKernel`（半平面投影），**无 CVX** |

---

## 后续开发速查

1. **锥点流程不依赖 MATLAB** → 实现 `.mat` frames 读取（§1），或始终使用 `.ffield`。
2. **无头 / CI 构建** → 移除或 override MATLAB 钩子（§4）。
3. **可调 `segSize` / `fixCot`** → 改写 `getSettings()`（§2）。
4. **浏览器与论文管线对齐** → 将桌面管线移植到 WASM，或暴露服务端 `FastHGP::run`（见 WASM 对比表）。
5. **与 MATLAB 对照** → `reference_matlab/` + `HGP.cpp` 中的历史混合代码；数值应与 `FastHGPNumerics` 一致。

## 参考

MATLAB 原版在 `reference_matlab/`。桌面逻辑来自历史 MATLAB 混合实现（`MatlabInterface` / `MatlabGMMDataExchange` 调用已替换为 `FastHGP.cpp` 与 `FastHGPNumerics.cpp` 中的 C++）。

相关博客解读：`_posts/1.Parameterization/3.几何优化方法/全局调和参数化FastHGP.md`。
