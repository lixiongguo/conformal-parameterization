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
- **FastHGP 本体不依赖 MATLAB**（`FASTHGP_STANDALONE` 构建；`visualize` 等已在 `FastHGP` 中重写）

编译：

```bash
cd cpp/conformal-parameterization/FastHGP
cmake -B build -DEIGEN3_INCLUDE_DIR=... -DCGAL_DIR=... -DGMM_INCLUDE_DIR=...
cmake --build build
```

生成静态库 `fasthgp`，定义 `FASTHGP_STANDALONE`（**不链接 MATLAB Engine**）。入口：`FastHGP::run(objPath, vfPath)`。

### 运行时参数

| 方式 | 说明 |
|:---|:---|
| `FastHGP::setSegSize(n)` / `setFixCot(bool)` | 在 `run()` 前调用 |
| 环境变量 `FASTHGP_SEG_SIZE` | 覆盖 meta 顶点间距（默认 40） |
| 环境变量 `FASTHGP_FIX_COT` | `0` 关闭 cot 折叠修复，非 `0` 开启（默认开启） |

### 锥点 frames 输入

| 文件 | 说明 |
|:---|:---|
| `.ffield` | 向量场（与原先相同） |
| `.mat` | MATLAB v5/v7.2，`frames` 变量（`save('-v7',...)`）；不支持 v7.3 HDF5 |
| `.fframes` | 纯文本，每行 `real imag`（推荐无 MATLAB 时使用） |

MATLAB 导出 `.fframes` 示例：

```matlab
writematrix([real(frames), imag(frames)], 'model.fframes', 'Delimiter', ' ');
```

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

显式占位宏 `FASTHGP_NOT_IMPLEMENTED` 仍保留于 `NotImplemented.h`（当前无调用点）。

### 1. ~~从 `.mat` 加载预计算 frames~~ → **已实现**

- 支持 `.mat`（v5/v7.2）与 `.fframes` 文本
- 实现：`Utils/FramesFile.cpp`，`FastHGP::loadPrecomputedFramesFromFile`
- v7.3 HDF5 格式请用 `save('-v7', ...)` 或改用 `.fframes`

### 2. ~~设置 GUI / 运行时配置~~ → **部分实现**

| | |
|:---|:---|
| **已移除** | `FastHGP_settings.m` 模态 GUI |
| **已实现** | `setSegSize` / `setFixCot`；环境变量 `FASTHGP_SEG_SIZE`、`FASTHGP_FIX_COT` |
| **待做** | 命令行解析、JSON/INI 配置文件 |

### 3. `FastHGP_report` 交互式报告

| | |
|:---|:---|
| **已移除** | `reference_matlab/FastHGP_report.m` |
| **当前代码** | `sendValuesToMatlabReport()` 为空（`FastHGP.cpp` 末尾） |
| **现状** | `testResult()` 向 stdout 打印 `success` / `partial success` / `fail` |
| **待做** | 可选 HTML/JSON 报告，或接入自有查看器 |

### 4. ~~继承自 `HarmonicParametrization` 的 MATLAB 钩子~~ → **FastHGP 已脱离**

`FastHGP` 重写 `visualize` / `calcDistortion` / `coneAngleDetection`（纯 C++）。  
CMake 目标 `fasthgp` 定义 `FASTHGP_STANDALONE`，`HarmonicParametrization.cpp` 中不再编译 MATLAB 调用。

仍使用 MATLAB 的是 **`HGP` 类**（未改），与 `FastHGP` 独立。

### 5. ~~桌面构建集成~~ → **已提供 CMake**

- `FastHGP/CMakeLists.txt` → 静态库 `fasthgp`
- WASM 构建脚本仍仅编译 `FastHGPSimple.cpp`

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

1. **MATLAB v7.3 `.mat`** → 用 `save('-v7',...)` 或 `.fframes`。
2. **CLI / 配置文件** → 扩展 `getSettings()` 或增加 `tools/fasthgp_main.cpp`。
3. **交互报告** → 替代 `FastHGP_report.m`（§3）。
4. **浏览器与论文管线对齐** → 移植桌面管线到 WASM 或暴露服务端 `FastHGP::run`。
5. **与 MATLAB 对照** → [MATLAB与C++对照表.md](MATLAB与C++对照表.md)。

## 参考

| 文档 | 内容 |
|:---|:---|
| [MATLAB与C++对照表.md](MATLAB与C++对照表.md) | `reference_matlab/*.m` 与 C++ 函数逐文件对照、主流程与字段映射 |
| `reference_matlab/` | MATLAB 原版脚本 |
| `_posts/1.Parameterization/3.几何优化方法/全局调和参数化FastHGP.md` | 博客算法解读 |

桌面逻辑来自历史 MATLAB 混合实现（`MatlabInterface` / `MatlabGMMDataExchange` 调用已替换为 `FastHGP.cpp` 与 `FastHGPNumerics.cpp` 中的 C++）。
