# ============================================================
# QuadCover WASM build
# ============================================================
$emsdkPath = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\emsdk"
. "$emsdkPath\emsdk_env.ps1"

$srcDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\conformal-parameterization"
$outDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\assets\wasm"
if (!(Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir -Force | Out-Null }
Set-Location $srcDir

Write-Host "Building QuadCover..."
emcc wasm_quadcover.cpp QuadCover.cpp Mesh.cpp MeshIO.cpp Parameterization.cpp Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp `
  -I./deps -I./deps/Eigen -I. `
  -s MODULARIZE=1 -s "EXPORT_NAME='QuadCoverSolver'" `
  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" `
  -s "EXPORTED_FUNCTIONS=['_malloc','_free','_solve_qc','_get_qc_uv_result','_get_qc_uv_result_size','_get_qc_last_time_ms','_qc_dispose']" `
  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 `
  -std=c++17 -O2 `
  -o "$outDir\quadcover_solver.js"

if($LASTEXITCODE -ne 0){Write-Host "[FAIL]" -ForegroundColor Red;exit 1}
Write-Host "[OK] quadcover_solver.js" -ForegroundColor Green
Write-Host "Done!" -ForegroundColor Green
