# WASM 编译指南

本文档记录将 C++ 参数化算法编译为 WebAssembly 的流程及遇到的问题。

## 环境

- Emscripten: `cpp/emsdk/`（内置 Node 22.16 / Python 3.13）
- C++ 标准: C++17, 优化 `-O2`
- 依赖: Eigen 3（header-only，`deps/Eigen/`）

## 编译目标

| 模块 | EXPORT_NAME | 输出文件 | 入口 |
|---|---|---|---|
| LSCM | `LSCMSolver` | `lscm_solver.js` | `wasm_main.cpp` |
| BD-LSCM | `BDLSCMSolver` | `bd_lscm_solver.js` | `BoundedDistortion/wasm_bd_lscm.cpp` |
| Tutte+ARAP | `TutteARAPSolver` | `tutte_arap_solver.js` | `wasm_tutte_arap.cpp` |

## 一键编译

```powershell
cd cpp/conformal-parameterization
powershell -ExecutionPolicy Bypass -File build_wasm_full.ps1
```

## 关键编译参数

| 参数 | 含义 |
|---|---|
| `MODULARIZE=1` | 生成工厂函数，JS 侧 `await Factory()` 获取实例 |
| `EXPORT_NAME='XXX'` | 模块全局变量名 |
| `EXPORTED_RUNTIME_METHODS` | 暴露给 JS 的方法（`ccall`, `cwrap`, `getValue`, `setValue`） |
| `EXPORTED_FUNCTIONS` | 需要导出的 C 函数（含 `_malloc`, `_free`） |
| `ALLOW_MEMORY_GROWTH=1` | WASM 堆动态增长 |
| `WASM=1` | 生成 `.wasm`（非 asm.js） |

## 遇到的问题 & 解决方案

### 1. `_malloc` / `_free` 导出失败

```
error: undefined exported symbol: "_malloc" in EXPORTED_RUNTIME_METHODS
error: undefined exported symbol: "_free" in EXPORTED_RUNTIME_METHODS
```

**原因**：新版本 Emscripten（3.x+）中 `_malloc`/`_free` 不属于 `EXPORTED_RUNTIME_METHODS`。

**解决**：移到 `EXPORTED_FUNCTIONS` 中：

```diff
- -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','_malloc','_free']"
+ -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','getValue','setValue']"
+ -s "EXPORTED_FUNCTIONS=['_malloc','_free','_solve_xxx',...]"
```

> 若使用 `--bind`（Embind），`_malloc`/`_free` 会自动可用。

### 2. `getValue` / `setValue` 不可用

JS 侧报 `wasmModule.getValue is not a function`。

**解决**：在 `EXPORTED_RUNTIME_METHODS` 中显式声明：
```bash
-s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','getValue','setValue']"
```

### 3. PowerShell 中文编码问题

```
The string is missing the terminator: ".
```

**原因**：`.ps1` 中的中文字符编码不一致（未保存为 UTF-8 BOM）。

**解决**：构建脚本中统一使用英文输出，或确保文件编码为 UTF-8 with BOM。

## JS ↔ WASM 数据交换流程

```
JS 侧                              WASM 侧
──────────────────────────────     ─────────────────
Float64Array 数据
  │
  ├─ _malloc(n * 8) ────────────→ 分配堆内存
  ├─ setValue(ptr, val) ────────→ 逐元素写入
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

## 输出文件

编译后 `assets/wasm/` 目录结构：

```
assets/wasm/
├── lscm_solver.js          # LSCM    JS 胶水代码
├── lscm_solver.wasm        # LSCM    WASM 二进制
├── bd_lscm_solver.js       # BD-LSCM JS 胶水代码
├── bd_lscm_solver.wasm     # BD-LSCM WASM 二进制
├── tutte_arap_solver.js    # Tutte+ARAP JS 胶水代码
└── tutte_arap_solver.wasm  # Tutte+ARAP WASM 二进制
```
