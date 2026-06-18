# HGP（Harmonic Global Parameterization）

BCW17 全局调和参数化的**桌面混合实现**：C++ 负责网格、约束组装与迭代控制；**核心 SOCP 仍在 MATLAB**（`HGP_iteration.m` + CVX + MOSEK）。

博客解读：`_posts/1.Parameterization/3.几何优化方法/全局调和参数化.md`  
与 FastHGP 的关系：见 [MATLAB与C++对照表.md](MATLAB与C++对照表.md) §7。

---

## 组件

| 文件 | 说明 |
|:---|:---|
| `HGP/HGP.h`, `HGP.cpp` | 主流程：组装 `HGP.*` → 调 `HGP_iteration` |
| `HGP/HGP_iteration.m` | **运行用**求解脚本（CVX SOCP）；与 `reference_matlab/HGP_iteration.m` 同源归档 |
| `reference_matlab/*.m` | 已迁 C++ 的辅助脚本归档（`computeGrads`、`frameFixABF` 等） |
| `HGP/FrameFixingClass.*` | 锥点 frame 修复（`frameFixABF` 已 C++ 化） |
| `HarmonicParametrization.*` | 网格、UV 写回、折叠/畸变检测、可视化 |
| `HGP/Parser.*`, `Borders.*` | OBJ、seam、锥点、边界 |
| `HGP/MatlabInterface.*` | MATLAB Engine 封装 |

---

## 运行依赖

- CGAL、`GMM`（`GMM_Macros.h`）
- **MATLAB**（Engine API）
- **CVX** + **MOSEK**（`HGP_iteration.m` 内 `cvx_solver_settings`）
- MATLAB path 包含 `HGP/HGP_iteration.m` 所在目录

入口：`HGP::run(objPath, vfPath)`。  
`vfPath`：`.ffield` 向量场，或 `.mat` 预计算 `frames`（经 MATLAB `load`）。

默认参数（原 `HGP_settings.m` GUI 已删）：`maxIt = 5`，`useFrameFixing = true`，`visMatlab = 0`（见 `HGP::run`）。

---

## 已实现（C++ 侧）

- 网格加载、cut、系统索引、旋转约束、边界/锥邻域面 `BVFaces`
- 内部 / seam 调和约束矩阵 `W`、`WSeam*`
- 梯度算子 `computeGradientsInCpp`（原 `computeGrads.m`）
- 向量场 frames、迭代中 `updateFramesFromCurrentFzInCpp`
- `FrameFixingClass`（原 `frameFixABF.m` 等）
- 迭代中基于当前 UV 更新边长权重 `updateHalfEdgesMetric`

详见 [MATLAB与C++对照表.md](MATLAB与C++对照表.md)。

---

## 未实现 / 待办

### 1. `HGP_iteration.m` → C++ SOCP（阻塞「纯 C++」）

| | |
|:---|:---|
| **现状** | 每轮 `EvalToCout("HGP_iteration")` |
| **目标** | MOSEK / Clarabel / 自写 QP 接口，输入与 `HGP.*` 相同 |
| **难点** | Lipman 约束 + 调和能量 + 旋转约束的联合 SOCP |

### 2. 脱离 MATLAB 的辅助功能

| 项 | 现状 | 可参考 |
|:---|:---|:---|
| `.mat` frames | MATLAB `load` | `../FastHGP/Utils/FramesFile.cpp` |
| `visualize` / `calcDistortion` / `coneAngleDetection` | 写 MATLAB workspace | `FastHGP` 中的 override |
| `HGP_report` GUI | 已删 | 可选 JSON/HTML 报告 |

### 3. 构建与配置

- 无 `CMakeLists.txt`；需手动链接 CGAL、GMM、MATLAB
- 迭代次数、`useFrameFixing` 写死在 `run()`；可加 setter 或环境变量（比照 FastHGP 的 `FASTHGP_SEG_SIZE`）

### 4. 算法范围（与原版一致）

- 依赖 seam / 锥点、向量场或预计算 frames
- 旋转约束当前为 $\pi/2$ 整数倍（`setRotationsConstraints` 中 `cosAngle[4]`）

---

## 完整性一句话

- **混合架构**（C++ + `HGP_iteration.m` + CVX/MOSEK）：✅ 可视为完整。  
- **纯 C++**：❌ 未完成；最大缺口是 §1。  
- 快速无 MATLAB 路线：使用 **`../FastHGP/`**（不同论文方法，非 HGP 的 SOCP 替代）。

---

## 参考

| 文档 | 内容 |
|:---|:---|
| [reference_matlab/README.md](reference_matlab/README.md) | 已从 git 恢复的 `.m` 归档列表 |
| [MATLAB与C++对照表.md](MATLAB与C++对照表.md) | 逐函数对照、`HGP.*` 字段映射 |
| [代码注释与论文对应.md](代码注释与论文对应.md) | 与博客章节、论文符号对照 |
| `../FastHGP/README.md` | FastHGP 子空间方法（另一套实现） |
