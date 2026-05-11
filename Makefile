############################################################################
# Unified WASM Makefile for Conformal Parameterization Algorithms
# Usage (from this directory):
#   make              — show help
#   make check        — check emsdk environment
#   make lscm         — build LSCM only
#   make all          — build all 11 targets
#   make clean        — remove all build artifacts
#
# Prerequisites:
#   - emsdk environment activated (source ../../emsdk/emsdk_env.sh)
#   - Eigen at $(EIGEN_INC)   (default: ../deps/eigen-3.4.0)
#   - libigl at $(LIBIGL_INC) (default: ../Libigl-Discrete-Geometry/libigl/include)
#   - Output: $(OUT_DIR) (default: ../../assets/wasm)
############################################################################

# ==================== Toolchain ====================
# Assumes emsdk_env.sh has been sourced (emcc/em++ in PATH)
EMCC ?= emcc
EMXX ?= em++

# ==================== Dependencies ====================
# cpp/deps/eigen-3.4.0/Eigen/Core  ←  #include <Eigen/Core> resolves here
EIGEN_INC   ?= ../deps/eigen-3.4.0
LIBIGL_INC  ?= ../Libigl-Discrete-Geometry/libigl/include
GLM_INC     ?= ../deps/glm
MOSEK_INC   ?= ../deps/mosek

# ==================== Output ====================
OUT_DIR  ?= ../../assets/wasm
BUILD_DIR ?= build

# ==================== Common Flags ====================
CXXFLAGS_COMMON := -std=c++17 -O2 -flto
INC_FLAGS      := -I$(EIGEN_INC) -I. -I$(GLM_INC) -I$(MOSEK_INC)

# ==================== Common Source Files ====================
MESH_SRCS   := Mesh.cpp MeshIO.cpp Parameterization.cpp Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp
SOLVER_SRCS := Solver.cpp QcError.cpp

# ==================== Header Files ====================
HDRS := $(wildcard *.h)

############################################################################
# Help
############################################################################
.PHONY: all clean check help \
       lscm tutte_arap cetm circle_patterns holo miq quadcover ricci bff bd_lscm principal_curvature

help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  all               Build all 11 WASM modules"
	@echo "  lscm              LSCM (Least Squares Conformal Maps)"
	@echo "  tutte_arap        Tutte Embedding + ARAP"
	@echo "  cetm              CETM (Conformal Energy Minimization)"
	@echo "  circle_patterns   Circle Patterns"
	@echo "  holo              Holomorphic 1-Form"
	@echo "  miq               MIQ Quad"
	@echo "  quadcover         QuadCover"
	@echo "  ricci             Ricci Flow"
	@echo "  bff               BFF (Boundary First Flattening)"
	@echo "  bd_lscm           Bounded Distortion LSCM"
	@echo "  principal_curvature  Principal Curvature (libigl)"
	@echo ""
	@echo "Utilities:"
	@echo "  check   Check emsdk environment"
	@echo "  clean   Remove all build artifacts"
	@echo ""
	@echo "Output directory: $(OUT_DIR)"

############################################################################
# Environment Check
############################################################################
check:
	@echo "=== WASM Build Environment Check ==="
	@echo "EMSDK in PATH: $$(which em++ 2>/dev/null || echo 'NOT FOUND — run: source ../../emsdk/emsdk_env.sh')"
	@echo "EIGEN_INC:   $(EIGEN_INC)"
	@echo "LIBIGL_INC:  $(LIBIGL_INC)"
	@echo "OUT_DIR:     $(OUT_DIR)"
	@if [ -d "$(EIGEN_INC)" ]; then echo "[OK] Eigen found"; else echo "[FAIL] Eigen not found at $(EIGEN_INC)"; fi
	@if [ -d "$(LIBIGL_INC)" ]; then echo "[OK] libigl found"; else echo "[WARN] libigl not found (needed for principal_curvature)"; fi
	@mkdir -p "$(OUT_DIR)"
	@echo ""
	@echo "Tip: source ../../emsdk/emsdk_env.sh  (activate emsdk first)"

############################################################################
# How efiles work:
#   printf '["_fn1","_fn2"]\n' > build/XXX.efile
#   -s EXPORTED_FUNCTIONS=@build/XXX.efile
# The file must contain a valid JSON array of strings.
############################################################################

# ==================== Helper: write efile ====================
# $(call write_efile,TARGET,func1 func2 ...)
define write_efile
	@mkdir -p $(BUILD_DIR)
	$(file > $(BUILD_DIR)/$(1).efile)
	$(foreach f,$(2),$(file >> $(BUILD_DIR)/$(1).efile,"$(f)",))
	$(file >> $(BUILD_DIR)/$(1).efile,$(endl))
endef

############################################################################
# 1. LSCM
# Exports: solve_lscm, get_uv_result, get_uv_result_size,
#          get_last_time_ms, dispose
############################################################################
LSCM_SRCS   := wasm_main.cpp Lscm.cpp $(SOLVER_SRCS) $(MESH_SRCS)
LSCM_EXPORT := LSCMSolver
LSCM_FUNCS  := _malloc _free solve_lscm get_uv_result get_uv_result_size get_last_time_ms dispose
LSCM_OUT    := $(OUT_DIR)/lscm_solver.js

lscm: $(LSCM_OUT)

$(LSCM_OUT): $(LSCM_SRCS) $(HDRS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_lscm","get_uv_result","get_uv_result_size","get_last_time_ms","dispose"]\n' > $(BUILD_DIR)/lscm.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(LSCM_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/lscm.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(LSCM_SRCS) -o $@
	@echo "[DONE] LSCM → $@"

############################################################################
# 2. Tutte + ARAP
# Exports: solve_tutte_circle, solve_tutte_square, solve_arap,
#          get_ta_uv_result, get_ta_uv_result_size, get_ta_last_time_ms,
#          compute_qc_error, load_mesh_with_uv,
#          get_qc_errors, get_qc_errors_size,
#          get_qc_colors, get_qc_colors_size, ta_dispose
############################################################################
TUTTE_SRCS   := wasm_tutte_arap.cpp Tutte.cpp ARAP.cpp QcError.cpp $(MESH_SRCS)
TUTTE_EXPORT := TutteARAPSolver
TUTTE_FUNCS  := _malloc _free solve_tutte_circle solve_tutte_square solve_arap get_ta_uv_result get_ta_uv_result_size get_ta_last_time_ms compute_qc_error load_mesh_with_uv get_qc_errors get_qc_errors_size get_qc_colors get_qc_colors_size ta_dispose
TUTTE_OUT    := $(OUT_DIR)/tutte_arap_solver.js

tutte_arap: $(TUTTE_OUT)

$(TUTTE_OUT): $(TUTTE_SRCS) $(HDRS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_tutte_circle","solve_tutte_square","solve_arap","get_ta_uv_result","get_ta_uv_result_size","get_ta_last_time_ms","compute_qc_error","load_mesh_with_uv","get_qc_errors","get_qc_errors_size","get_qc_colors","get_qc_colors_size","ta_dispose"]\n' > $(BUILD_DIR)/tutte.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(TUTTE_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/tutte.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(TUTTE_SRCS) -o $@
	@echo "[DONE] Tutte+ARAP → $@"

############################################################################
# 3. CETM
# Exports: solve_cetm, get_cetm_uv_result, get_cetm_uv_result_size,
#          get_cetm_last_time_ms, cetm_dispose
############################################################################
CETM_SRCS   := wasm_cetm.cpp Cetm.cpp $(SOLVER_SRCS) $(MESH_SRCS)
CETM_EXPORT := CETMSolver
CETM_FUNCS  := _malloc _free solve_cetm get_cetm_uv_result get_cetm_uv_result_size get_cetm_last_time_ms cetm_dispose
CETM_OUT    := $(OUT_DIR)/cetm_solver.js

cetm: $(CETM_OUT)

$(CETM_OUT): $(CETM_SRCS) $(HDRS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_cetm","get_cetm_uv_result","get_cetm_uv_result_size","get_cetm_last_time_ms","cetm_dispose"]\n' > $(BUILD_DIR)/cetm.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(CETM_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/cetm.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(CETM_SRCS) -o $@
	@echo "[DONE] CETM → $@"

############################################################################
# 4. Circle Patterns
# Exports: solve_cp, get_cp_uv_result, get_cp_uv_result_size,
#          get_cp_last_time_ms, cp_dispose
############################################################################
CP_SRCS   := wasm_circle_patterns.cpp CirclePatternsWasm.cpp $(SOLVER_SRCS) $(MESH_SRCS)
CP_EXPORT := CPSolver
CP_FUNCS  := _malloc _free solve_cp get_cp_uv_result get_cp_uv_result_size get_cp_last_time_ms cp_dispose
CP_OUT    := $(OUT_DIR)/cp_solver.js

circle_patterns: $(CP_OUT)

$(CP_OUT): $(CP_SRCS) $(HDRS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_cp","get_cp_uv_result","get_cp_uv_result_size","get_cp_last_time_ms","cp_dispose"]\n' > $(BUILD_DIR)/cp.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(CP_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/cp.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(CP_SRCS) -o $@
	@echo "[DONE] Circle Patterns → $@"

############################################################################
# 5. Holomorphic 1-Form
# Exports: solve_holo, get_holo_uv_result, get_holo_uv_result_size,
#          get_holo_last_time_ms, holo_dispose
############################################################################
HOLO_SRCS   := wasm_holo.cpp HolomorphicOneForm.cpp $(MESH_SRCS)
HOLO_EXPORT := HoloSolver
HOLO_FUNCS  := _malloc _free solve_holo get_holo_uv_result get_holo_uv_result_size get_holo_last_time_ms holo_dispose
HOLO_OUT    := $(OUT_DIR)/holo_solver.js

holo: $(HOLO_OUT)

$(HOLO_OUT): $(HOLO_SRCS) $(HDRS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_holo","get_holo_uv_result","get_holo_uv_result_size","get_holo_last_time_ms","holo_dispose"]\n' > $(BUILD_DIR)/holo.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(HOLO_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/holo.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(HOLO_SRCS) -o $@
	@echo "[DONE] Holomorphic 1-Form → $@"

############################################################################
# 6. MIQ Quad
# Exports: solve_miq, get_miq_uv_result, get_miq_uv_result_size,
#          get_miq_last_time_ms, miq_dispose
############################################################################
MIQ_SRCS   := wasm_miq.cpp MIQQuad.cpp $(MESH_SRCS)
MIQ_EXPORT := MIQSolver
MIQ_FUNCS  := _malloc _free solve_miq get_miq_uv_result get_miq_uv_result_size get_miq_last_time_ms miq_dispose
MIQ_OUT    := $(OUT_DIR)/miq_solver.js

miq: $(MIQ_OUT)

$(MIQ_OUT): $(MIQ_SRCS) $(HDRS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_miq","get_miq_uv_result","get_miq_uv_result_size","get_miq_last_time_ms","miq_dispose"]\n' > $(BUILD_DIR)/miq.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(MIQ_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/miq.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(MIQ_SRCS) -o $@
	@echo "[DONE] MIQ Quad → $@"

############################################################################
# 7. QuadCover
# Exports: solve_qc, get_qc_uv_result, get_qc_uv_result_size,
#          get_qc_last_time_ms, qc_dispose
############################################################################
QC_SRCS   := wasm_quadcover.cpp QuadCover.cpp $(MESH_SRCS)
QC_EXPORT := QuadCoverSolver
QC_FUNCS  := _malloc _free solve_qc get_qc_uv_result get_qc_uv_result_size get_qc_last_time_ms qc_dispose
QC_OUT    := $(OUT_DIR)/quadcover_solver.js

quadcover: $(QC_OUT)

$(QC_OUT): $(QC_SRCS) $(HDRS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_qc","get_qc_uv_result","get_qc_uv_result_size","get_qc_last_time_ms","qc_dispose"]\n' > $(BUILD_DIR)/qc.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(QC_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/qc.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(QC_SRCS) -o $@
	@echo "[DONE] QuadCover → $@"

############################################################################
# 8. Ricci Flow
# Exports: solve_ricci, get_ricci_uv_result, get_ricci_uv_result_size,
#          get_ricci_last_time_ms, ricci_dispose
############################################################################
RICCI_SRCS   := wasm_ricci.cpp RicciFlow.cpp $(SOLVER_SRCS) $(MESH_SRCS)
RICCI_EXPORT := RicciFlowSolver
RICCI_FUNCS  := _malloc _free solve_ricci get_ricci_uv_result get_ricci_uv_result_size get_ricci_last_time_ms ricci_dispose
RICCI_OUT    := $(OUT_DIR)/ricci_solver.js

ricci: $(RICCI_OUT)

$(RICCI_OUT): $(RICCI_SRCS) $(HDRS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_ricci","get_ricci_uv_result","get_ricci_uv_result_size","get_ricci_last_time_ms","ricci_dispose"]\n' > $(BUILD_DIR)/ricci.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(RICCI_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/ricci.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(RICCI_SRCS) -o $@
	@echo "[DONE] Ricci Flow → $@"

############################################################################
# 9. BFF (Boundary First Flattening)
# Exports: solve_bff, get_bff_uv_result, get_bff_uv_result_size,
#          get_bff_last_time_ms, bff_dispose
############################################################################
BFF_SRCS   := wasm_bff.cpp BFF.cpp $(MESH_SRCS)
BFF_EXPORT := BFFSolver
BFF_FUNCS  := _malloc _free solve_bff get_bff_uv_result get_bff_uv_result_size get_bff_last_time_ms bff_dispose
BFF_OUT    := $(OUT_DIR)/bff_solver.js

bff: $(BFF_OUT)

$(BFF_OUT): $(BFF_SRCS) $(HDRS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_bff","get_bff_uv_result","get_bff_uv_result_size","get_bff_last_time_ms","bff_dispose"]\n' > $(BUILD_DIR)/bff.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(BFF_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/bff.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(BFF_SRCS) -o $@
	@echo "[DONE] BFF → $@"

############################################################################
# 10. BD-LSCM (Bounded Distortion LSCM)
# Exports: solve_bd_lscm, get_bd_uv_result, get_bd_uv_result_size,
#          get_bd_last_time_ms, bd_dispose
############################################################################
BD_SRCS   := wasm_bd_lscm.cpp bounded_distortion.h
BD_EXPORT := BDLSCMSolver
BD_FUNCS  := _malloc _free solve_bd_lscm get_bd_uv_result get_bd_uv_result_size get_bd_last_time_ms bd_dispose
BD_OUT    := $(OUT_DIR)/bd_lscm_solver.js

bd_lscm: $(BD_OUT)

$(BD_OUT): $(BD_SRCS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","solve_bd_lscm","get_bd_uv_result","get_bd_uv_result_size","get_bd_last_time_ms","bd_dispose"]\n' > $(BUILD_DIR)/bd.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(INC_FLAGS) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(BD_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/bd.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  $(BD_SRCS) -o $@
	@echo "[DONE] BD-LSCM → $@"

############################################################################
# 11. Principal Curvature (depends on libigl)
# Exports: compute_principal_curvature, get_result_buffer_size
############################################################################
PC_SRCS   := principal_curvature_wasm.cpp
PC_EXPORT := PrincipalCurvatureSolver
PC_FUNCS  := _malloc _free compute_principal_curvature get_result_buffer_size
PC_OUT    := $(OUT_DIR)/principal_curvature.js
PC_INC    := -I$(EIGEN_INC) -I$(LIBIGL_INC) -I.

principal_curvature: $(PC_OUT)

$(PC_OUT): $(PC_SRCS)
	@mkdir -p $(OUT_DIR) $(BUILD_DIR)
	@printf '["_malloc","_free","compute_principal_curvature","get_result_buffer_size"]\n' > $(BUILD_DIR)/pc.efile
	$(EMXX) $(CXXFLAGS_COMMON) $(PC_INC) \
	  -s MODULARIZE=1 -s EXPORT_NAME="$(PC_EXPORT)" \
	  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 \
	  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" \
	  -s EXPORTED_FUNCTIONS=@$(BUILD_DIR)/pc.efile \
	  -s FORCE_FILESYSTEM=0 -s ENVIRONMENT='web' \
	  --bind \
	  $(PC_SRCS) -o $@
	@echo "[DONE] Principal Curvature → $@"

############################################################################
# Build All
############################################################################
all: lscm tutte_arap cetm circle_patterns holo miq quadcover ricci bff bd_lscm principal_curvature
	@echo ""
	@echo "=== All 11 WASM modules built successfully ==="
	@ls -lh $(OUT_DIR)/*.js $(OUT_DIR)/*.wasm 2>/dev/null || true

############################################################################
# Clean
############################################################################
clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(BUILD_DIR)
	@rm -f $(OUT_DIR)/lscm_solver.js $(OUT_DIR)/lscm_solver.wasm
	@rm -f $(OUT_DIR)/tutte_arap_solver.js $(OUT_DIR)/tutte_arap_solver.wasm
	@rm -f $(OUT_DIR)/cetm_solver.js $(OUT_DIR)/cetm_solver.wasm
	@rm -f $(OUT_DIR)/cp_solver.js $(OUT_DIR)/cp_solver.wasm
	@rm -f $(OUT_DIR)/holo_solver.js $(OUT_DIR)/holo_solver.wasm
	@rm -f $(OUT_DIR)/miq_solver.js $(OUT_DIR)/miq_solver.wasm
	@rm -f $(OUT_DIR)/quadcover_solver.js $(OUT_DIR)/quadcover_solver.wasm
	@rm -f $(OUT_DIR)/ricci_solver.js $(OUT_DIR)/ricci_solver.wasm
	@rm -f $(OUT_DIR)/bff_solver.js $(OUT_DIR)/bff_solver.wasm
	@rm -f $(OUT_DIR)/bd_lscm_solver.js $(OUT_DIR)/bd_lscm_solver.wasm
	@rm -f $(OUT_DIR)/principal_curvature.js $(OUT_DIR)/principal_curvature.wasm
	@echo "[DONE] Clean complete."

############################################################################
# Install emsdk (helper target)
############################################################################
.PHONY: install_emsdk

install_emsdk:
	@echo "Cloning emsdk..."
	@git clone https://github.com/emscripten-core/emsdk.git ../../emsdk
	cd ../../emsdk && ./emsdk install latest && ./emsdk activate latest
	@echo ""
	@echo "To activate in current shell:  source ../../emsdk/emsdk_env.sh"
