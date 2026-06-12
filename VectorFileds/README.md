# VectorFields — 向量场计算

对应文章 [`2017-10-01-向量场的计算`](../../_posts/1.Parameterization/4.全局参数化与四边形网格化/2017-10-01-向量场的计算.md) 中的三种上游算法。

| 方法 | 可执行文件 | 说明 |
|:---|:---|:---|
| 联络修正 (N=1) | `trivial_connection` | 约束 LS 求连接修正 φ，再传播向量场 |
| 联络修正 (N-RoSy) | `nrosy_trivial_connection` | 同上，holonomy 模 2π/N |
| 复多项式 (N-RoSy) | `complex_poly_field` | bar{e}^N q 稀疏线性 LS（Diamanti SGP 2014） |
| 度量驱动 | — | 见 `cpp/conformal-parameterization/RicciFlow.h` |

## 构建 (CMake)

```powershell
cd cpp/VectorFileds
cmake -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

## 运行示例

```powershell
.\build\Release\complex_poly_field.exe ..\..\assets\Models2\torus_F7680.obj out_torus --n 4
.\build\Release\nrosy_trivial_connection.exe mesh.obj out_tc --n 4 --auto-singularity
```

## WebAssembly

`vector_field_unified_wasm.cpp` 由 `cpp/build/build_wasm_vector_field.ps1` 编译，导出 `compute_4rosy_field`（复多项式）与 `compute_lc_propagated_field`（纯 LC 传播，无奇异点修正）。
