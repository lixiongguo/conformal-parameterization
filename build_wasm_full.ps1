# ============================================================
# Full WASM build: LSCM + Tutte/ARAP
# ============================================================
$emsdkPath = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\emsdk"
. "$emsdkPath\emsdk_env.ps1"

$srcDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\conformal-parameterization"
$outDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\assets\wasm"
if (!(Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir -Force | Out-Null }
Set-Location $srcDir

Write-Host "===================================="

# ---- LSCM ----
Write-Host "[1/2] Building LSCM..."
emcc wasm_main.cpp Mesh.cpp MeshIO.cpp Lscm.cpp Parameterization.cpp QcError.cpp Solver.cpp Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp `
  -I./deps -I./deps/Eigen -I. `
  -s MODULARIZE=1 -s "EXPORT_NAME='LSCMSolver'" `
  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" `
  -s "EXPORTED_FUNCTIONS=['_malloc','_free','_solve_lscm','_get_uv_result','_get_uv_result_size','_get_last_time_ms','_dispose']" `
  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 `
  -std=c++17 -O2 `
  -o "$outDir\lscm_solver.js"
if ($LASTEXITCODE -ne 0) { Write-Host "[FAIL] LSCM"; exit 1 }
Write-Host "  [OK] lscm_solver.js"

# ---- Tutte + ARAP ----
Write-Host "[2/2] Building Tutte+ARAP..."
emcc wasm_tutte_arap.cpp Tutte.cpp ARAP.cpp Mesh.cpp MeshIO.cpp Parameterization.cpp QcError.cpp Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp `
  -I./deps -I./deps/Eigen -I. `
  -s MODULARIZE=1 -s "EXPORT_NAME='TutteARAPSolver'" `
  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','getValue','setValue']" `
  -s "EXPORTED_FUNCTIONS=['_malloc','_free','_solve_tutte_circle','_solve_tutte_square','_solve_arap','_compute_qc_error','_load_mesh_with_uv','_get_ta_uv_result','_get_ta_uv_result_size','_get_ta_last_time_ms','_get_qc_errors','_get_qc_errors_size','_get_qc_colors','_get_qc_colors_size','_ta_dispose']" `
  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 `
  -std=c++17 -O2 `
  -o "$outDir\tutte_arap_solver.js"
if ($LASTEXITCODE -ne 0) { Write-Host "[FAIL] Tutte+ARAP"; exit 1 }
Write-Host "  [OK] tutte_arap_solver.js"

Write-Host ""
Write-Host "Verifying output..."
$files = @("$outDir\lscm_solver.js","$outDir\lscm_solver.wasm",
           "$outDir\tutte_arap_solver.js","$outDir\tutte_arap_solver.wasm")
foreach ($f in $files) {
    if (Test-Path $f) { Write-Host "  [OK] $f" }
    else { Write-Host "  [MISS] $f" }
}
Write-Host "Done!"
