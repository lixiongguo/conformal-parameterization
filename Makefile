# LSCM Solver WASM Build Script
# 
# 使用 Emscripten 编译 C++ LSCM 求解器为 WebAssembly
#
# 前置条件:
#   1. 安装 Emscripten SDK (emsdk)
#      - 下载: https://emscripten.org/docs/getting_started/downloads.html
#      - 安装: emsdk install latest
#                   emsdk activate latest
#                   emsdk env
#
#   2. 下载 Eigen 库
#      - 下载地址: https://eigen.tuxfamily.org/index.php?title=Main_Page
#      - 解压到本目录的 eigen/ 子目录
#      - 确认文件存在: eigen/Eigen/Dense
#
# 使用方法:
#   Windows (PowerShell): .\Makefile
#   或手动执行:
    # em++ lscm_solver.cpp -o ../assets/wasm/lscm_solver.js \
    #     -O3 -std=c++17 \
    #     -ID:\third_party\eigen-3.4.0 \
    #     --bind \
    #     -s MODULARIZE=1 \
    #     -s EXPORT_NAME="LCMSolver" \
    #     -s ALLOW_MEMORY_GROWTH=1 \
    #     -s TOTAL_MEMORY=256MB \
    #     -s WASM=1
#
#   Linux/macOS:
#     em++ lscm_solver.cpp -o ../assets/wasm/lscm_solver.js \
#         -O3 -std=c++17 \
#         -I../eigen \
#         --bind \
#         -s MODULARIZE=1 \
#         -s EXPORT_NAME="LCMSolver" \
#         -s ALLOW_MEMORY_GROWTH=1 \
#         -s TOTAL_MEMORY=256MB \
#         -s WASM=1

# ========== 配置区域 ==========

# Eigen 库路径 (根据实际情况修改)
# EIGEN_PATH := E:/Dev/libs/eigen-3.4.0
EIGEN_PATH := ./eigen-3.4.0
# 或者使用相对路径 (如果 Eigen 解压在同一父目录)
# EIGEN_PATH := ../eigen

# 输出目录
OUTPUT_DIR := ../assets/wasm

# 编译目标
TARGET_JS := $(OUTPUT_DIR)/lscm_solver.js
TARGET_WASM := $(OUTPUT_DIR)/lscm_solver.wasm

# 源文件
SRC := conformal-parameterization/wasm_bindings.cpp \
       conformal-parameterization/wasm_mesh.cpp \
       conformal-parameterization/MeshIO.cpp \
       conformal-parameterization/Vertex.cpp \
       conformal-parameterization/Edge.cpp \
       conformal-parameterization/Face.cpp \
       conformal-parameterization/HalfEdge.cpp \
       conformal-parameterization/Lscm.cpp \
       conformal-parameterization/Cetm.cpp \
       conformal-parameterization/Solver.cpp \
       conformal-parameterization/Parameterization.cpp \
       conformal-parameterization/QcError.cpp

# 编译选项
CXX := em++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra
INC := -I$(EIGEN_PATH)
LDFLAGS := --bind \
           -s MODULARIZE=1 \
           -s EXPORT_NAME="LSCMSolver" \
           -s ALLOW_MEMORY_GROWTH=1 \
           -s WASM=1 \
           -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap', 'getValue', 'setValue']" \
           -s EXPORTED_FUNCTIONS="['_malloc', '_free']"

# ========== 构建目标 ==========

.PHONY: all build clean watch help

all: build

build: $(TARGET_JS)
	@echo ""
	@echo "构建成功!"
	@echo "生成文件:"
	@echo "  - $(TARGET_JS)"
	@echo "  - $(TARGET_WASM)"
	@echo ""
	@echo "请确保这些文件已提交到仓库并部署到 GitHub Pages。"

$(TARGET_JS): $(SRC)
	@echo "正在编译 LSCM 求解器为 WebAssembly..."
	@echo "源文件: $(SRC)"
	@echo "Eigen 路径: $(EIGEN_PATH)"
	@echo ""
	$(CXX) $(CXXFLAGS) $(INC) $(SRC) -o $(TARGET_JS) $(LDFLAGS)
	@echo ""
	@echo "编译完成!"

# 清理生成的文件
clean:
	@echo "清理生成的文件..."
	-rm -f $(TARGET_JS) $(TARGET_WASM)
	@echo "清理完成。"

# 首次检查依赖
check:
	@echo "检查编译环境..."
	@which em++ > /dev/null 2>&1 && echo "[OK] Emscripten (em++) 已安装" || echo "[FAIL] Emscripten 未安装，请先安装: https://emscripten.org/docs/getting_started/downloads.html"
	@test -d "$(EIGEN_PATH)" && echo "[OK] Eigen 库已找到: $(EIGEN_PATH)" || echo "[FAIL] Eigen 库未找到: $(EIGEN_PATH)"
	@echo ""
	@echo "如需修改 Eigen 路径，请编辑本文件的 EIGEN_PATH 变量。"

# 帮助信息
help:
	@echo "LSCM Solver WASM 编译脚本"
	@echo ""
	@echo "使用方法:"
	@echo "  make          编译 WASM 模块"
	@echo "  make check    检查编译环境"
	@echo "  make clean    清理生成的文件"
	@echo "  make help     显示本帮助信息"
	@echo ""
	@echo "编译前提:"
	@echo "  1. Emscripten SDK: https://emscripten.org"
	@echo "  2. Eigen 库: https://eigen.tuxfamily.org"
	@echo ""
	@echo "首次编译前请运行 'make check' 检查环境配置。"
