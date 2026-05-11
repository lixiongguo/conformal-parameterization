# WASM 编译指南

本文档记录将 C++ 参数化算法编译为 WebAssembly 的流程及遇到的问题。

## 环境

- Emscripten: `cpp/emsdk/`（内置 Node / Python）
- C++ 标准: C++17, 优化 `-O2`
- 依赖:
  - Eigen 3（header-only，`cpp/deps/Eigen/`）
  - libigl（`cpp/Libigl-Discrete-Geometry/libigl/include/`，仅 principal_curvature 需要）
  - glm（`cpp/deps/glm/`）
  - mosek（`cpp/deps/mosek/`）

## 一键编译（macOS / Linux）

```bash
# 1. 激活 emsdk
source ../../emsdk/emsdk_env.sh

# 2. 进入目录
cd cpp/conformal-parameterization

# 3. 构建全部 11 个目标
make all

# 或单独构建某个目标
make lscm
make tutte_arap
make cetm
make circle_patterns
make holo
make miq
make quadcover
make ricci
make bff
make bd_lscm
make principal_curvature

# 清理
make clean

# 环境检查
make check
```

## 编译目标

| 目标 | EXPORT_NAME | 输出文件 | 入口文件 |
|---|---|---|---|
| LSCM | `LSCMSolver` | `lscm_solver.js` | `wasm_main.cpp` |
| Tutte+ARAP | `TutteARAPSolver` | `tutte_arap_solver.js` | `wasm_tutte_arap.cpp` |
| CETM | `CETMSolver` | `cetm_solver.js` | `wasm_cetm.cpp` |
| Circle Patterns | `CPSolver` | `cp_solver.js` | `wasm_circle_patterns.cpp` |
| Holomorphic 1-Form | `HoloSolver` | `holo_solver.js` | `wasm_holo.cpp` |
| MIQ Quad | `MIQSolver` | `miq_solver.js` | `wasm_miq.cpp` |
| QuadCover | `QuadCoverSolver` | `quadcover_solver.js` | `wasm_quadcover.cpp` |
| Ricci Flow | `RicciFlowSolver` | `ricci_solver.js` | `wasm_ricci.cpp` |
| BFF | `BFFSolver` | `bff_solver.js` | `wasm_bff.cpp` |
| BD-LSCM | `BDLSCMSolver` | `bd_lscm_solver.js` | `wasm_bd_lscm.cpp` |
| Principal Curvature | `PrincipalCurvatureSolver` | `principal_curvature.js` | `principal_curvature_wasm.cpp` |

## 输出文件

编译后 `assets/wasm/` 目录结构：

```
assets/wasm/
├── lscm_solver.js
├── lscm_solver.wasm
├── tutte_arap_solver.js
├── tutte_arap_solver.wasm
├── cetm_solver.js
├── cetm_solver.wasm
├── cp_solver.js
├── cp_solver.wasm
├── holo_solver.js
├── holo_solver.wasm
├── miq_solver.js
├── miq_solver.wasm
├── quadcover_solver.js
├── quadcover_solver.wasm
├── ricci_solver.js
├── ricci_solver.wasm
├── bff_solver.js
├── bff_solver.js
├── bd_lscm_solver.js
├── bd_lscm_solver.wasm
├── principal_curvature.js
└── principal_curvature.wasm
```

## Makefile 变量覆盖

如需自定义路径，可在 make 命令行覆盖：

```bash
make all EMSDK_ROOT=../../emsdk EMXX=em++ EIGEN_INC=../../deps/Eigen OUT_DIR=../../../assets/wasm
```

或在 Makefile 顶部直接修改默认值。

## 关键编译参数

| 参数 | 含义 |
|---|---|
| `MODULARIZE=1` | 生成工厂函数，JS 侧 `await Factory()` 获取实例 |
| `EXPORT_NAME='XXX'` | 模块全局变量名 |
| `EXPORTED_RUNTIME_METHODS` | 暴露给 JS 的方法（`ccall`, `cwrap`, `getValue`, `setValue`） |
| `EXPORTED_FUNCTIONS` | 需要导出的 C 函数（含 `_malloc`, `_free`） |
| `ALLOW_MEMORY_GROWTH=1` | WASM 堆动态增长 |
| `WASM=1` | 生成 `.wasm`（非 asm.js） |

## JS ↔ WASM 数据交换流程

```
JS 侧                              WASM 侧
──────────────────────────────     ─────────────────
Float64Array 数据
  │
  ├─ _malloc(n * 8) ───────────→ 分配堆内存
  ├─ setValue(ptr, val) ───────→ 逐元素写入
  │
  ├─ ccall('solve_xxx', ...) ───→ C 函数处理
  │                                  ├─ 加载网格
  │                                  ├─ 执行参数化
  │                                  ├─ 计时
  │                                  └─ UV → 全局缓冲
  │
  ├─ ccall('get_uv_size') ──────→ 返回数据长度
  ├─ ccall('get_uv_result') ────→ 返回数据指针
  ├─ getValue(ptr + i*8) ───────→ 逐元素读取
  │
  ├─ _free(ptr) ────────────────→ 释放内存
  └─ ccall('xxx_dispose') ──────→ 清理全局状态
```

## C 导出函数模板

```cpp
extern "C" {

EMSCRIPTEN_KEEPALIVE
int solve_xxx(double* pos, int posLen, int* faces, int faceLen) {
    // 1. 加载网格到全局变量
    // 2. 执行参数化
    // 3. 提取 UV 到全局浮点缓冲
    return 0;
}

EMSCRIPTEN_KEEPALIVE double* get_xxx_uv_result() { return g_uv.data(); }
EMSCRIPTEN_KEEPALIVE int     get_xxx_uv_result_size() { return g_uv.size(); }
EMSCRIPTEN_KEEPALIVE double  get_xxx_last_time_ms() { return g_time; }
EMSCRIPTEN_KEEPALIVE void    xxx_dispose() { delete g_mesh; g_uv.clear(); }

}
```

**要点**：
- `EMSCRIPTEN_KEEPALIVE` 防止 LTO 优化掉导出函数
- 全局变量避免栈上返回悬空指针
- `_malloc` 分配的内存必须由 JS 侧 `_free` 释放

## 遇到的问题 & 解决方案

### 1. `_malloc` / `_free` 导出失败

```
error: undefined exported symbol: "_malloc" in EXPORTED_RUNTIME_METHODS
error: undefined exported symbol: "_free" in EXPORTED_RUNTIME_METHODS
```

**原因**：新版本 Emscripten（3.x+）中 `_malloc`/`_free` 不属于 `EXPORTED_RUNTIME_METHODS`。

**解决**：移到 `EXPORTED_FUNCTIONS` 中（Makefile 中已正确处理）。

### 2. `getValue` / `setValue` 不可用

JS 侧报 `wasmModule.getValue is not a function`。

**解决**：在 `EXPORTED_RUNTIME_METHODS` 中显式声明：
```bash
-s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','getValue','setValue']"
```

### 3. macOS 上 .ps1 / .bat 脚本不可用

旧脚本为 Windows/PowerShell 编写，macOS 无法直接运行。已统一改为 Makefile。

## 文件组织结构

```
cpp/conformal-parameterization/
├── Makefile              # ← 统一构建脚本（本文件）
├── WASM_BUILD_GUIDE.md   # ← 本文档
│
├── wasm_main.cpp          # LSCM 入口
├── wasm_tutte_arap.cpp   # Tutte+ARAP 入口
├── wasm_cetm.cpp         # CETM 入口
├── wasm_circle_patterns.cpp  # Circle Patterns 入口
├── wasm_holo.cpp         # Holomorphic 1-Form 入口
├── wasm_miq.cpp          # MIQ Quad 入口
├── wasm_quadcover.cpp    # QuadCover 入口
├── wasm_ricci.cpp        # Ricci Flow 入口
├── wasm_bff.cpp          # BFF 入口
├── wasm_bd_lscm.cpp      # BD-LSCM 入口
├── principal_curvature_wasm.cpp  # Principal Curvature 入口
│
├── Lscm.cpp / .h         # 各算法实现
├── Cetm.cpp / .h
├── ...
├── Mesh.cpp / .h          # 公共网格数据结构
├── Parameterization.cpp/.h
├── bounded_distortion.h   # BD-LSCM 算法
└── ...
```
