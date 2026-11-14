#!/usr/bin/env bash
# 编译 Ruppert 三角化 WASM（自包含，无 Eigen 依赖）。
#
#   ./build_ruppert.sh
#
# 产物：../../assets/wasm/ruppert.js + ruppert.wasm（供 Triangularizition.html 使用）
set -euo pipefail

cd "$(dirname "$0")"

EMSDK_DIR="${EMSDK_DIR:-/Users/lgxgeogo/Desktop/MyDoc/lixiongguo.github.io/cpp/emsdk}"
if [[ -f "$EMSDK_DIR/emsdk_env.sh" ]]; then
  # shellcheck disable=SC1091
  source "$EMSDK_DIR/emsdk_env.sh" >/dev/null 2>&1
fi

if ! command -v em++ >/dev/null 2>&1; then
  echo "em++ 不可用，请先激活 emscripten（或设置 EMSDK_DIR）" >&2
  exit 1
fi

OUT_DIR="../../assets/wasm"
mkdir -p "$OUT_DIR"

echo "[wasm] em++ -> $OUT_DIR/ruppert.js"
em++ -O3 -std=c++17 \
  bindings_ruppert_compat.cpp \
  -o "$OUT_DIR/ruppert.js" \
  -sMODULARIZE=1 -sALLOW_MEMORY_GROWTH=1 -sWASM=1 -sSTACK_SIZE=1048576 \
  -sEXPORTED_FUNCTIONS="['_malloc','_free']" \
  -sEXPORTED_RUNTIME_METHODS="['ccall','cwrap','getValue','setValue','HEAPF64','HEAP32']" \
  -sEXPORT_NAME="RuppertModule"

echo "[wasm] done"
ls -lh "$OUT_DIR/ruppert.js" "$OUT_DIR/ruppert.wasm"
