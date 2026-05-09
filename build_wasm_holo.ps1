$emsdkPath="c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\emsdk"
."$emsdkPath\emsdk_env.ps1"
$srcDir="c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\conformal-parameterization"
$outDir="c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\assets\wasm"
if(!(Test-Path $outDir)){New-Item -ItemType Directory -Path $outDir -Force|Out-Null}
Set-Location $srcDir
Write-Host "Building Holomorphic 1-Form..."
emcc wasm_holo.cpp HolomorphicOneForm.cpp Mesh.cpp MeshIO.cpp Parameterization.cpp Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp `
  -I./deps -I./deps/Eigen -I. `
  -s MODULARIZE=1 -s "EXPORT_NAME='HoloSolver'" `
  -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','getValue','setValue']" `
  -s "EXPORTED_FUNCTIONS=['_malloc','_free','_solve_holo','_get_holo_uv_result','_get_holo_uv_result_size','_get_holo_last_time_ms','_holo_dispose']" `
  -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 -std=c++17 -O2 -o "$outDir\holo_solver.js"
if($LASTEXITCODE -ne 0){Write-Host "[FAIL]" -ForegroundColor Red;exit 1}
Write-Host "[OK] holo_solver.js" -ForegroundColor Green
