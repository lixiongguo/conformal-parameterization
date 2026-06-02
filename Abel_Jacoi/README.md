# Abel–Jacobi 四边形网格化（C++）

对应博客 §4.5：在闭合三角网格上实现离散 Abel–Jacobi 映射、周期格 $\Gamma$ 与除子可行性检验。

## 文件

| 文件 | 说明 |
|------|------|
| `AbelJacobi.h/cpp` | 核心：同调基、正规全纯 1-形式、格点、$\mu(p)$、$\mu(D)$ |
| `AbelJacobiParameterization.h/cpp` | `Parameterization` 子类，周期取整后输出 UV |
| `abel_jacobi_demo.cpp` | 命令行 demo |

## 算法流程

1. **Tree–Cotree**（`topology/TreeCotreeBasis.h`）构造 $2g$ 条同调环路，前 $g$ 条作 $a$-cycle，后 $g$ 条作 $b$-cycle。
2. **正规基** $\{\varphi_j\}$：对每个 $j$ 求解调和 1-形式，满足 $\int_{a_i}\varphi_j = \delta_{ij}$（离散 cotan Laplacian + 周期约束）。
3. **Hodge 共轭** $\star\varphi_j$，得到全纯 1-形式 $\varphi_j + i\star\varphi_j$。
4. **格点** $\Gamma$：列向量为 $\lambda_{a_k}, \lambda_{b_k} \in \mathbb{C}^g$（沿环路积分）。
5. **Abel–Jacobi 映射**：沿 BFS 树从 $p_0$ 积分到 $p$，得 $\mu(p) \in \mathbb{C}^g$。
6. **除子检验**：$\mu(D)=\sum d_i\mu(p_i)$，判断是否在 $\Gamma$ 中（最近格点 + 残差）。
7. **参数化**：启发式 period quantization 后沿树积分，写入 `mesh.vertices[i].uv`。

## 限制

- 当前仅支持**无边界闭合网格**且 **genus ≥ 1**。
- genus 0 或 build 失败时自动 **fallback 到 LSCM**。
- 周期整数化为小亏格启发式，非完整 MIQ/MIP。

## WASM 构建（前端 `uv_unwrap_global_cross_fields.html`）

```powershell
cd cpp/build
.\abel-jacobi\build_wasm_abel_jacobi.ps1
# 或
.\build_wasm_all.ps1 -Targets abel_jacobi
```

输出：`assets/wasm/uv_unwrap_abel_jacobi.js` + `.wasm`

导出模块名：`UvUnwrapAbelJacobiSolver`  
主入口：`solve_abel_jacobi(pos, faces, baseVertex, …)`

## 编译 demo

在 `cpp/` 目录下（需 Eigen 3.4）：

```powershell
g++ -std=c++17 -O2 `
  -I conformal-parameterization `
  -I deps/eigen-3.4.0 `
  conformal-parameterization/Abel_Jacoi/abel_jacobi_demo.cpp `
  conformal-parameterization/Abel_Jacoi/AbelJacobi.cpp `
  conformal-parameterization/Abel_Jacoi/AbelJacobiParameterization.cpp `
  conformal-parameterization/Mesh.cpp `
  conformal-parameterization/MeshIO.cpp `
  conformal-parameterization/Parameterization.cpp `
  conformal-parameterization/Vertex.cpp `
  conformal-parameterization/Edge.cpp `
  conformal-parameterization/Face.cpp `
  conformal-parameterization/HalfEdge.cpp `
  conformal-parameterization/Solver.cpp `
  conformal-parameterization/QcError.cpp `
  conformal-parameterization/uv_unwrap_simple/Lscm.cpp `
  -o abel_jacobi_demo
```

```bash
./abel_jacobi_demo torus.obj 0
```

## API 示例

```cpp
#include "AbelJacobiParameterization.h"

Mesh mesh;
mesh.read("torus.obj");

AbelJacobiParameterization param(mesh);
param.setBaseVertex(0);

// 可选：指定四边形奇异点（degree = 4 - valence）
std::vector<AbelJacobi::SingularPoint> sings;
sings.push_back({42, 3, +1});  // 3-valence
sings.push_back({17, 5, -1});  // 5-valence
param.setSingularities(sings);

auto check = param.checkSingularities();
// check.poincareHopfOk, check.abelJacobiOk, check.latticeResidual

param.parameterize();
```

## 与文档公式对应

| 公式 | 代码 |
|------|------|
| $\int_{a_i}\varphi_j = \delta_{ij}$ | `solveHarmonicWithPeriods` |
| $\lambda_\gamma = (\int_\gamma\varphi_1,\ldots)$ | `periodOnCycle` |
| $\Gamma = \{\sum s_k\lambda_{a_k}+t_k\lambda_{b_k}\}$ | `buildLattice` |
| $\mu(p) = (\int_{p_0}^p \varphi_1,\ldots) \bmod \Gamma$ | `abelJacobiMap` |
| $\mu(D) = \sum \nu_p(f)\mu(p) = 0$ | `checkDivisor` |
