# ============================================================
#  Conformal Parameterization — Multi-module WASM Build System
# ============================================================
#
#  Usage:
#    make build-base-mesh       Build BaseMesh WASM module
#    make build-cutseam-mesh    Build CutSeamMesh WASM module
#    make build-solver          Build Solver WASM module
#    make build-lscm            Build LSCM (Lscm + Scp) WASM module
#    make build-abf             Build ABF (Abf++ + LinAbf) WASM module
#    make build-arap            Build ARAP (ARAP + Tutte) WASM module
#    make build-all             Build all 6 modules
#    make clean                 Remove all generated files
#    make check                 Verify toolchain
#    make help                  Show this help

# ========== Configuration ==========
EIGEN_PATH := ../../deps/eigen-3.4.0
OUTPUT_DIR := ../../../assets/wasm

CXX       := em++
CXXFLAGS  := -O3 -std=c++17 -Wall -Wextra

# Common Emscripten link flags
EM_FLAGS  := --bind \
             -s MODULARIZE=1 \
             -s ALLOW_MEMORY_GROWTH=1 \
             -s WASM=1 \
             -s EXPORTED_RUNTIME_METHODS="['ccall','cwrap','getValue','setValue']"

# Compat flags (no --bind, EMSCRIPTEN_KEEPALIVE for ccall API)
EM_FLAGS_COMPAT := -s MODULARIZE=1 \
             -s ALLOW_MEMORY_GROWTH=1 \
             -s WASM=1 \
             -s EXPORTED_FUNCTIONS="['_malloc','_free']" \
             -s EXPORTED_RUNTIME_METHODS="['ccall','cwrap','getValue','setValue']"

# ========== Source file groups ==========
BASEMESH_SRCS := ../BaseMesh/Mesh.cpp \
                 ../BaseMesh/MeshIO.cpp \
                 ../BaseMesh/Vertex.cpp \
                 ../BaseMesh/Edge.cpp \
                 ../BaseMesh/Face.cpp \
                 ../BaseMesh/HalfEdge.cpp \
                 ../BaseMesh/GaussianCurvature.cpp \
                 ../BaseMesh/QcError.cpp \
                 ../Parameterization/Parameterization.cpp

CUTSEAM_SRCS  := $(BASEMESH_SRCS) ../CutSeamMesh/CutSeamMesh.cpp

SOLVER_SRCS   := ../Solvers/Solver.cpp \
                 ../Solvers/AugmentedLagrangian.cpp \
                 ../Solvers/MixedIntegerProgram.cpp

LSCM_SRCS     := $(BASEMESH_SRCS) \
                 ../Parameterization/SimpleParam/LSCM/Lscm.cpp \
                 ../Parameterization/SimpleParam/LSCM/Scp.cpp

ABF_SRCS      := $(BASEMESH_SRCS) $(SOLVER_SRCS) \
                 ../Parameterization/SimpleParam/ABF/AbfPlusPlus.cpp \
                 ../Parameterization/SimpleParam/ABF/LinAbf.cpp \
                 ../Parameterization/SimpleParam/LSCM/Lscm.cpp

ARAP_SRCS     := $(BASEMESH_SRCS) $(SOLVER_SRCS) \
                 ../Parameterization/SimpleParam/ARAP/ARAP.cpp \
                 ../Parameterization/SimpleParam/ARAP/Tutte.cpp

BD_LSCM_SRCS  := $(BASEMESH_SRCS) \
                 ../Parameterization/SimpleParam/LSCM/Lscm.cpp \
                 ../../bounded_distortion_mapping/BoundedDistortionMapping.cpp

# ========== Compat 模块 (ccall API) ==========
DGP_COMPAT_SRCS     := $(BASEMESH_SRCS)
SIMPLE_COMPAT_SRCS  := $(BASEMESH_SRCS) $(SOLVER_SRCS) \
                       ../Parameterization/SimpleParam/LSCM/Lscm.cpp \
                       ../Parameterization/SimpleParam/LSCM/Scp.cpp \
                       ../Parameterization/SimpleParam/ABF/AbfPlusPlus.cpp \
                       ../Parameterization/SimpleParam/ABF/LinAbf.cpp \
                       ../Parameterization/SimpleParam/ARAP/ARAP.cpp \
                       ../Parameterization/SimpleParam/ARAP/Tutte.cpp

SIMPLE_COMPAT_INC := -I$(EIGEN_PATH) -I../BaseMesh -I../Parameterization \
                     -I../Parameterization/SimpleParam/LSCM \
                     -I../Parameterization/SimpleParam/ABF \
                     -I../Parameterization/SimpleParam/ARAP \
                     -I../Solvers -I../Solvers/Mosek

# ========== BD-HM (Bounded Distortion Harmonic Map) ==========
BDHM_SRCS      := ../../bdhm/BoundedDistortionHarmonicMap.cpp \
                  ../../cauchy_coordinates/CauchyCoordinates.cpp \
                  ../../Solvers/Lbfgs.cpp \
                  ../../Solvers/Solver.cpp
BDHM_INC       := -I$(EIGEN_PATH) -I../../bdhm -I../../cauchy_coordinates -I../../Solvers -I../../Solvers/Mosek -I../BaseMesh

# ========== Include paths per module ==========
BASEMESH_INC  := -I$(EIGEN_PATH) -I../BaseMesh -I../Parameterization
CUTSEAM_INC   := -I$(EIGEN_PATH) -I../BaseMesh -I../CutSeamMesh -I../Parameterization
SOLVER_INC    := -I$(EIGEN_PATH) -I../BaseMesh -I../Solvers -I../Solvers/Mosek
LSCM_INC      := -I$(EIGEN_PATH) -I../BaseMesh -I../Parameterization -I../Parameterization/SimpleParam/LSCM
ABF_INC       := -I$(EIGEN_PATH) -I../BaseMesh -I../Parameterization -I../Parameterization/SimpleParam/ABF -I../Parameterization/SimpleParam/LSCM -I../Solvers -I../Solvers/Mosek
ARAP_INC      := -I$(EIGEN_PATH) -I../BaseMesh -I../Parameterization -I../Parameterization/SimpleParam/ARAP -I../Solvers -I../Solvers/Mosek
BD_LSCM_INC   := -I$(EIGEN_PATH) -I../BaseMesh -I../Parameterization -I../Parameterization/SimpleParam/LSCM -I../../bounded_distortion_mapping

# ========== Build targets ==========
.PHONY: build-base-mesh build-cutseam-mesh build-solver \
        build-lscm build-abf build-arap build-bd-lscm \
        build-dgp-compat build-simple-compat build-bdhm \
        build-all clean check help

build-all: build-base-mesh build-cutseam-mesh build-solver build-lscm build-abf build-arap build-bd-lscm build-dgp-compat build-simple-compat build-bdhm

# --- BaseMesh ---
build-base-mesh: $(OUTPUT_DIR)/base_mesh.js
$(OUTPUT_DIR)/base_mesh.js: bindings_base_mesh.cpp $(BASEMESH_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[BaseMesh] Compiling..."
	$(CXX) $(CXXFLAGS) $(BASEMESH_INC) bindings_base_mesh.cpp $(BASEMESH_SRCS) \
	       -o $@ $(EM_FLAGS) -s EXPORT_NAME="BaseMeshModule"
	@echo "[BaseMesh] Done → $@"

# --- CutSeamMesh ---
build-cutseam-mesh: $(OUTPUT_DIR)/cutseam_mesh.js
$(OUTPUT_DIR)/cutseam_mesh.js: bindings_cutseam_mesh.cpp $(CUTSEAM_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[CutSeamMesh] Compiling..."
	$(CXX) $(CXXFLAGS) $(CUTSEAM_INC) bindings_cutseam_mesh.cpp $(CUTSEAM_SRCS) \
	       -o $@ $(EM_FLAGS) -s EXPORT_NAME="CutSeamMeshModule"
	@echo "[CutSeamMesh] Done → $@"

# --- Solver ---
build-solver: $(OUTPUT_DIR)/solver.js
$(OUTPUT_DIR)/solver.js: bindings_solver.cpp $(SOLVER_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[Solver] Compiling..."
	$(CXX) $(CXXFLAGS) $(SOLVER_INC) bindings_solver.cpp $(SOLVER_SRCS) \
	       -o $@ $(EM_FLAGS) -s EXPORT_NAME="SolverModule"
	@echo "[Solver] Done → $@"

# --- LSCM ---
build-lscm: $(OUTPUT_DIR)/lscm.js
$(OUTPUT_DIR)/lscm.js: bindings_lscm.cpp $(LSCM_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[LSCM] Compiling..."
	$(CXX) $(CXXFLAGS) $(LSCM_INC) bindings_lscm.cpp $(LSCM_SRCS) \
	       -o $@ $(EM_FLAGS) -s EXPORT_NAME="LSCMModule"
	@echo "[LSCM] Done → $@"

# --- ABF ---
build-abf: $(OUTPUT_DIR)/abf.js
$(OUTPUT_DIR)/abf.js: bindings_abf.cpp $(ABF_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[ABF] Compiling..."
	$(CXX) $(CXXFLAGS) $(ABF_INC) bindings_abf.cpp $(ABF_SRCS) \
	       -o $@ $(EM_FLAGS) -s EXPORT_NAME="ABFModule"
	@echo "[ABF] Done → $@"

# --- ARAP ---
build-arap: $(OUTPUT_DIR)/arap.js
$(OUTPUT_DIR)/arap.js: bindings_arap.cpp $(ARAP_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[ARAP] Compiling..."
	$(CXX) $(CXXFLAGS) $(ARAP_INC) bindings_arap.cpp $(ARAP_SRCS) \
	       -o $@ $(EM_FLAGS) -s EXPORT_NAME="ARAPModule"
	@echo "[ARAP] Done → $@"

# --- BD-LSCM (LSCM + Bounded Distortion Mapping) ---
build-bd-lscm: $(OUTPUT_DIR)/bd_lscm.js
$(OUTPUT_DIR)/bd_lscm.js: bindings_bd_lscm.cpp $(BD_LSCM_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[BD-LSCM] Compiling..."
	$(CXX) $(CXXFLAGS) $(BD_LSCM_INC) bindings_bd_lscm.cpp $(BD_LSCM_SRCS) \
	       -o $@ $(EM_FLAGS) -s EXPORT_NAME="BD_LSCMModule" \
	       -s FORCE_FILESYSTEM=1 \
	       -s EXPORTED_RUNTIME_METHODS="['ccall','cwrap','getValue','setValue','FS']"
	@echo "[BD-LSCM] Done → $@"

# --- DGP Basic 兼容层 (高斯曲率, ccall API) ---
build-dgp-compat: $(OUTPUT_DIR)/dgp_basic.js
$(OUTPUT_DIR)/dgp_basic.js: bindings_dgp_basic_compat.cpp $(DGP_COMPAT_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[DGP-Compat] Compiling..."
	$(CXX) $(CXXFLAGS) $(BASEMESH_INC) bindings_dgp_basic_compat.cpp $(DGP_COMPAT_SRCS) \
	       -o $@ $(EM_FLAGS_COMPAT) -s EXPORT_NAME="DgpBasicSolver"
	@echo "[DGP-Compat] Done → $@"

# --- Simple 兼容层 (所有参数化方法, ccall API) ---
build-simple-compat: $(OUTPUT_DIR)/uv_unwrap_simple.js
$(OUTPUT_DIR)/uv_unwrap_simple.js: bindings_uv_unwrap_simple_compat.cpp $(SIMPLE_COMPAT_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[Simple-Compat] Compiling..."
	$(CXX) $(CXXFLAGS) $(SIMPLE_COMPAT_INC) bindings_uv_unwrap_simple_compat.cpp $(SIMPLE_COMPAT_SRCS) \
	       -o $@ $(EM_FLAGS_COMPAT) -s EXPORT_NAME="UvUnwrapSimpleSolver"
	@echo "[Simple-Compat] Done → $@"

# --- BD-HM (Bounded Distortion Harmonic Map) ---
build-bdhm: $(OUTPUT_DIR)/bdhm.js
$(OUTPUT_DIR)/bdhm.js: bindings_bdhm.cpp $(BDHM_SRCS)
	@mkdir -p $(OUTPUT_DIR)
	@echo "[BDHM] Compiling..."
	$(CXX) $(CXXFLAGS) $(BDHM_INC) bindings_bdhm.cpp $(BDHM_SRCS) \
	       -o $@ $(EM_FLAGS) -s EXPORT_NAME="BDHMModule"
	@echo "[BDHM] Done → $@"

# ========== Utility ==========
clean:
	@echo "Cleaning WASM outputs..."
	-rm -f $(OUTPUT_DIR)/base_mesh.js  $(OUTPUT_DIR)/base_mesh.wasm
	-rm -f $(OUTPUT_DIR)/cutseam_mesh.js $(OUTPUT_DIR)/cutseam_mesh.wasm
	-rm -f $(OUTPUT_DIR)/solver.js  $(OUTPUT_DIR)/solver.wasm
	-rm -f $(OUTPUT_DIR)/lscm.js   $(OUTPUT_DIR)/lscm.wasm
	-rm -f $(OUTPUT_DIR)/abf.js    $(OUTPUT_DIR)/abf.wasm
	-rm -f $(OUTPUT_DIR)/arap.js   $(OUTPUT_DIR)/arap.wasm
	-rm -f $(OUTPUT_DIR)/bd_lscm.js $(OUTPUT_DIR)/bd_lscm.wasm
	-rm -f $(OUTPUT_DIR)/dgp_basic.js $(OUTPUT_DIR)/dgp_basic.wasm
	-rm -f $(OUTPUT_DIR)/uv_unwrap_simple.js $(OUTPUT_DIR)/uv_unwrap_simple.wasm
	-rm -f $(OUTPUT_DIR)/bdhm.js $(OUTPUT_DIR)/bdhm.wasm
	@echo "Done."

check:
	@echo "=== Toolchain Check ==="
	@which em++ > /dev/null 2>&1 && echo "[OK] Emscripten (em++)" || echo "[FAIL] Emscripten not found"
	@test -f "$(EIGEN_PATH)/Eigen/Dense" && echo "[OK] Eigen $(EIGEN_PATH)" || echo "[FAIL] Eigen not found at $(EIGEN_PATH)"
	@test -f "../BaseMesh/Mesh.h" && echo "[OK] BaseMesh sources" || echo "[FAIL] BaseMesh sources"
	@test -f "../Solvers/Solver.h" && echo "[OK] Solver sources" || echo "[FAIL] Solver sources"
	@test -f "../Parameterization/Parameterization.h" && echo "[OK] Parameterization sources"
	@test -f "../CutSeamMesh/CutSeamMesh.h" && echo "[OK] CutSeamMesh sources"
	@echo ""

help:
	@echo "Conformal Parameterization — WASM Build System"
	@echo ""
	@echo "Targets:"
	@echo "  make build-base-mesh       Build BaseMesh module"
	@echo "  make build-cutseam-mesh    Build CutSeamMesh module"
	@echo "  make build-solver          Build Solver module"
	@echo "  make build-lscm            Build LSCM module (Lscm + Scp)"
	@echo "  make build-abf             Build ABF module (Abf++ + LinAbf)"
	@echo "  make build-arap            Build ARAP module (ARAP + Tutte)"
	@echo "  make build-bd-lscm         Build BD-LSCM module (LSCM + Bounded Distortion)"
	@echo "  make build-dgp-compat      Build dgp_basic compat (Gauss curvature, ccall)"
	@echo "  make build-simple-compat   Build uv_unwrap_simple compat (all methods, ccall)"
	@echo "  make build-bdhm            Build BD-HM module (Bounded Distortion Harmonic Map)"
	@echo "  make build-all             Build all 10 modules"
	@echo "  make clean                 Remove outputs"
	@echo "  make check                 Verify toolchain"
	@echo ""
