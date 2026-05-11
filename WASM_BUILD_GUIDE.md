# WASM 编译与测试指南

## 目录结构

```
cpp/
├── conformal-parameterization/     # 参数化算法源码
│   ├── build_wasm_all.ps1          # Windows 统一构建脚本（PowerShell）
│   ├── Makefile                    # macOS/Linux 统一构建脚本（make）
│   ├── wasm_*.cpp                  # 各算法 WASM 入口（C 接口）
│   ├── *.h / *.cpp                 # 算法实现 + Half-Edge 网格结构
│   └── bounded_distortion.h       # BD-LSCM 算法（header-only，依赖 Eigen）
├── deps/                           # 第三方依赖
│   ├── eigen-3.4.0/                # Eigen 3.4.0（线性代数）
│   ├── glm/                        # OpenGL Mathematics（向量/矩阵运算）
│   └── mosek/                      # MOSEK 优化库（部分算法可选）
├── Libigl-Discrete-Geometry/       # libigl（几何处理库）
│   └── libigl/include/igl/         # libigl header-only 部分
└── emsdk/                          # Emscripten SDK（C++ → WASM 编译工具链）
```

## WASM 模块列表（共 11 个）

| # | Target | 输出文件（js/wasm） | 入口源文件 | 算法依赖 |
|---|--------|---------------------|-----------|---------|
| 1 | `lscm` | `lscm_solver.js/wasm` | `wasm_main.cpp` | Lscm + Solver + QcError + Mesh |
| 2 | `tutte_arap` | `tutte_arap_solver.js/wasm` | `wasm_tutte_arap.cpp` | Tutte + ARAP + QcError + Mesh |
| 3 | `cetm` | `cetm_solver.js/wasm` | `wasm_cetm.cpp` | Cetm + Solver + QcError + Mesh |
| 4 | `circle_patterns` | `cp_solver.js/wasm` | `wasm_circle_patterns.cpp` | CirclePatternsWasm + Solver + QcError + Mesh |
| 5 | `holo` | `holo_solver.js/wasm` | `wasm_holo.cpp` | HolomorphicOneForm + Mesh |
| 6 | `miq` | `miq_solver.js/wasm` | `wasm_miq.cpp` | MIQQuad + Mesh |
| 7 | `quadcover` | `quadcover_solver.js/wasm` | `wasm_quadcover.cpp` | QuadCover + Mesh |
| 8 | `ricci` | `ricci_solver.js/wasm` | `wasm_ricci.cpp` | RicciFlow + Solver + QcError + Mesh |
| 9 | `bff` | `bff_solver.js/wasm` | `wasm_bff.cpp` | BFF + Mesh |
| 10 | `bd_lscm` | `bd_lscm_solver.js/wasm` | `wasm_bd_lscm.cpp` | bounded_distortion.h（独立，不依赖 Mesh） |
| 11 | `principal_curvature` | `principal_curvature.js/wasm` | `principal_curvature_wasm.cpp` | **libigl**（`igl::principal_curvature`） |

### 公共依赖

- `MESH_SRCS`：`Mesh.cpp MeshIO.cpp Parameterization.cpp Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp`
- `SOLVER_SRCS`：`Solver.cpp QcError.cpp`
- 所有目标均引用 `<Eigen/Core>`（路径：`../deps/eigen-3.4.0`）

## 编译命令

### macOS / Linux（Makefile）

```bash
# 1. 激活 Emscripten 环境
source ../emsdk/emsdk_env.sh

# 2. 检查环境
make check

# 3. 编译全部 11 个目标
make all

# 4. 编译单个目标
make lscm
make tutte_arap
make cetm
# ... 等
```

### Windows（PowerShell）

```powershell
# 进入源码目录
cd cpp\conformal-parameterization

# 编译全部 11 个目标
.\build_wasm_all.ps1

# 编译指定目标
.\build_wasm_all.ps1 -Targets lscm,tutte_arap

# 跳过 libigl 相关目标
.\build_wasm_all.ps1 -SkipLibigl

# 跳过 BD-LSCM
.\build_wasm_all.ps1 -SkipBd

# 同时跳过两者
.\build_wasm_all.ps1 -SkipLibigl -SkipBd
```

### 编译参数说明

| 参数 | 值 | 说明 |
|------|-----|------|
| 编译器 | `em++` | Emscripten C++ 编译器 |
| C++ 标准 | `-std=c++17` | C++17 |
| 优化 | `-O2 -flto` | O2 优化 + 链接时优化 |
| Eigen | `-I../deps/eigen-3.4.0` | 线性代数库 |
| GLM | `-I../deps/glm` | 向量/矩阵库 |
| MOSEK | `-I../deps/mosek` | 优化库（部分算法） |
| libigl | `-I../Libigl-Discrete-Geometry/libigl/include` | 几何处理（仅 principal_curvature） |
| MODULARIZE | `1` | 模块化输出，通过工厂函数加载 |
| WASM | `1` | 输出 .wasm 二进制 |
| ALLOW_MEMORY_GROWTH | `1` | 允许动态增长内存 |
| ENVIRONMENT | `web` | 目标环境为浏览器 |
| FORCE_FILESYSTEM | `0` | 不强制文件系统支持 |

## Emscripten 导出函数（`_` 前缀规则）

每个 WASM 目标通过 `EXPORTED_FUNCTIONS` 指定对外暴露的 C 函数列表。

**关键规则**：`extern "C"` 声明的函数在 `EXPORTED_FUNCTIONS` 中必须加 `_` 前缀，但 JavaScript 端通过 `ccall()` 调用时**不加** `_` 前缀。

### 示例：LSCM 目标

**C++ 端（wasm_main.cpp）**：
```cpp
extern "C" {
    EMSCRIPTEN_KEEPALIVE int solve_lscm(double* posPtr, int posLen, ...);
    EMSCRIPTEN_KEEPALIVE double* get_uv_result();
    EMSCRIPTEN_KEEPALIVE int get_uv_result_size();
    EMSCRIPTEN_KEEPALIVE double get_last_time_ms();
    EMSCRIPTEN_KEEPALIVE void dispose();
}
```

**EXPORTED_FUNCTIONS（efile / JSON）**：
```json
["_malloc","_free","_solve_lscm","_get_uv_result","_get_uv_result_size","_get_last_time_ms","_dispose"]
```

**JavaScript 调用**：
```javascript
const solver = await LSCMSolver();  // MODULARIZE=1 工厂函数
solver.ccall('solve_lscm', 'number',
  ['number','number','number','number','number','number'],
  [posPtr, posLen, facePtr, faceLen, anchor0, anchor1]);
```

## 本地测试（Python HTTP Server）

### 启动服务器

```powershell
# 从项目根目录启动（确保 assets/ 路径可访问）
python -m http.server 8080 --directory "c:\path\to\lixiongguo.github.io"
```

### 访问测试页面

浏览器打开：
```
http://localhost:8080/uv-unwrap.html
```

### 测试页面功能

`uv-unwrap.html` 是完整的 UV 展开演示页面，提供：

- **12 种算法**切换（LSCM / BD-LSCM / Tutte / ARAP / CETM / Ricci / CP / QuadCover / BFF / MIQ / Holo）
- **3 组模型库**（`assets/Models/`、`assets/Models2/`、`assets/Models3/`）
- **参数控制**：BD 上限 C、ARAP 迭代次数、CETM/Ricci/CP 优化方法
- **纹理显示**：棋盘格贴图 / 共形扭曲热力图
- **性能统计**：计算耗时（ms）、翻转面数
- 左侧 3D 模型预览（Three.js OrbitControls）
- 右侧 UV 网格线框展示

### WASM 引擎预加载

页面启动时自动加载所有 WASM 模块：
```javascript
Promise.all([
    initWASMSolver(),      // LSCM
    initBDWASMSolver(),    // BD-LSCM
    initTAWASMSolver(),    // Tutte + ARAP
    initCETMWASMSolver(),  // CETM
    initRicciWASMSolver(), // Ricci Flow
    initCPWASMSolver(),    // Circle Patterns
    initQCWASMSolver(),    // QuadCover
    initBFFWASMSolver(),   // BFF
    initMIQWASolver(),     // MIQ
    initHoloWASMSolver()   // Holomorphic 1-Form
]);
```

## WASM 调用流程

```mermaid
sequenceDiagram
    participant JS as JavaScript
    participant W as WASM Module
    participant C as C++ (wasm_*.cpp)
    participant A as Algorithm

    JS->>W: LSCMSolver() 工厂函数
    W-->>JS: Module 实例

    JS->>JS: 构建 Float64Array(pos) + Int32Array(faces)
    JS->>W: _malloc(posBytes)
    W-->>JS: posPtr
    JS->>W: HEAPF64.set(flatPositions, posPtr/8)
    JS->>W: _malloc(faceBytes)
    W-->>JS: facePtr
    JS->>W: HEAP32.set(flatFaces, facePtr/4)

    JS->>W: ccall('solve_lscm', [posPtr,posLen,facePtr,faceLen,a0,a1])
    W->>C: solve_lscm()
    C->>A: Mesh::read() → parameterize(LSCM)
    A-->>C: UV 坐标存入 Mesh
    C-->>W: return 0

    JS->>W: ccall('get_uv_result')
    W-->>JS: uvPtr (HEAPF64 指针)
    JS->>W: ccall('get_uv_result_size')
    W-->>JS: uvSize
    JS->>JS: new Float64Array(HEAPF64.buffer, uvPtr, uvSize)

    JS->>W: _free(posPtr); _free(facePtr);
    JS->>W: ccall('dispose')
```

## 常见问题

### 1. `emcc: error: undefined exported symbol`

**原因**：`EXPORTED_FUNCTIONS` 中的函数名缺少 `_` 前缀。

**解决**：所有 `extern "C"` 函数在 `EXPORTED_FUNCTIONS` 中均需加 `_` 前缀，如 `_solve_lscm` 而非 `solve_lscm`。

### 2. `wasm-ld: error: ... undefined symbol`

**原因**：缺少源文件或 include 路径不正确。

**解决**：检查 `deps/` 目录结构是否完整，include 路径是否正确。

### 3. WASM 文件 MIME 类型错误

**现象**：浏览器报 `WebAssembly.instantiateStreaming` failed。

**解决**：Python 3.7+ 已内置 `application/wasm` MIME 类型。如使用其他服务器，确保配置了正确的 Content-Type。

## 版本记录

| 日期 | 变更 |
|------|------|
| 2026-05-11 | 整理依赖到 `cpp/deps/`，统一 Mac Makefile 和 Windows PowerShell 构建脚本 |
| 2026-05-11 | 添加 `principal_curvature` 和 `bd_lscm` 两个新 WASM 目标 |
| 2026-05-11 | 修复 `EXPORTED_FUNCTIONS` 中函数名缺 `_` 前缀的编译错误 |
| 2026-05-11 | 创建本指南文档 |
