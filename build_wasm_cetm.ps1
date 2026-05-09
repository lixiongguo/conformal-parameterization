# ============================================================
# CETM WASM build: 共形因子参数化 (Conformal Factor Method)
# ============================================================
$emsdkPath = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\emsdk"
. "$emsdkPath\emsdk_env.ps1"

$srcDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\conformal-parameterization"
$outDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\assets\wasm"
if (!(Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir -Force | Out-Null }
Set-Location $srcDir

Write-Host "===================================="
Write-Host "Building CETM Conformal Parameterization"
Write-Host "===================================="

# ---- CETM (共形因子方法) ----
# 优化方法: 0=梯度下降, 1=牛顿法(推荐), 3=LBFGS
Write-Host "[1/1] Building CETM Solver..."
emcc wasm_cetm.cpp Cetm.cpp Mesh.cpp MeshIO.cpp Solver.cpp Parameterization.cpp QcError.cpp Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp `
  -I./deps -I./deps/Eigen -I. `
  -s MODULARIZE=1 -s "EXPORT_NAME='CETMSolver'" `
  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" `
  -s "EXPORTED_FUNCTIONS=['_malloc','_free','_solve_cetm','_get_cetm_uv_result','_get_cetm_uv_result_size','_get_cetm_last_time_ms','_cetm_dispose']" `
  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 `
  -std=c++17 -O2 `
  -o "$outDir\cetm_solver.js"

if ($LASTEXITCODE -ne 0) { 
    Write-Host "[FAIL] CETM build failed" -ForegroundColor Red
    exit 1 
}

Write-Host "  [OK] cetm_solver.js" -ForegroundColor Green

Write-Host ""
Write-Host "Verifying output..."
$files = @("$outDir\cetm_solver.js", "$outDir\cetm_solver.wasm")
foreach ($f in $files) {
    if (Test-Path $f) { Write-Host "  [OK] $f" -ForegroundColor Green }
    else { Write-Host "  [MISS] $f" -ForegroundColor Red }
}
Write-Host ""
Write-Host "Done! CETM Solver ready." -ForegroundColor Green
