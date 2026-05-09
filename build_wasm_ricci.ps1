# ============================================================
# Ricci Flow WASM build: 离散 Ricci 流参数化
# ============================================================
$emsdkPath = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\emsdk"
. "$emsdkPath\emsdk_env.ps1"

$srcDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\conformal-parameterization"
$outDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\assets\wasm"
if (!(Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir -Force | Out-Null }
Set-Location $srcDir

Write-Host "===================================="
Write-Host "Building Ricci Flow Parameterization"
Write-Host "===================================="

# ---- Ricci Flow ----
# 优化方法: 0=梯度下降, 1=牛顿法(推荐), 3=LBFGS
Write-Host "[1/1] Building RicciFlow Solver..."
emcc wasm_ricci.cpp RicciFlow.cpp Mesh.cpp MeshIO.cpp Solver.cpp Parameterization.cpp QcError.cpp Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp `
  -I./deps -I./deps/Eigen -I. `
  -s MODULARIZE=1 -s "EXPORT_NAME='RicciFlowSolver'" `
  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" `
  -s "EXPORTED_FUNCTIONS=['_malloc','_free','_solve_ricci','_get_ricci_uv_result','_get_ricci_uv_result_size','_get_ricci_last_time_ms','_ricci_dispose']" `
  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 `
  -std=c++17 -O2 `
  -o "$outDir\ricci_solver.js"

if ($LASTEXITCODE -ne 0) { 
    Write-Host "[FAIL] RicciFlow build failed" -ForegroundColor Red
    exit 1 
}

Write-Host "  [OK] ricci_solver.js" -ForegroundColor Green

Write-Host ""
Write-Host "Verifying output..."
$files = @("$outDir\ricci_solver.js", "$outDir\ricci_solver.wasm")
foreach ($f in $files) {
    if (Test-Path $f) { Write-Host "  [OK] $f" -ForegroundColor Green }
    else { Write-Host "  [MISS] $f" -ForegroundColor Red }
}
Write-Host ""
Write-Host "Done! Ricci Flow Solver ready." -ForegroundColor Green
