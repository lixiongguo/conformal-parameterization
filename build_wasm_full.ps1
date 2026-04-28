# build_wasm_full.ps1 - 完整编译脚本 (激活环境 + 编译)
# 使用方法: powershell -ExecutionPolicy Bypass -File build_wasm_full.ps1

Write-Host "==================================="
Write-Host "Conformal Parameterization WASM Build"
Write-Host "==================================="

# 激活 Emscripten 环境
$emsdkPath = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\emsdk"
. "$emsdkPath\emsdk_env.ps1"

# 检查 emcc 是否可用
try {
    $emccVersion = emcc --version 2>&1 | Select-Object -First 1
    Write-Host "[OK] 找到 emcc: $emccVersion"
} catch {
    Write-Host "[错误] 未找到 emcc 编译器，请检查 emsdk 安装"
    exit 1
}

# 设置路径
$srcDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\conformal-parameterization"
$outDir = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\assets\wasm"

# 创建输出目录
if (!(Test-Path $outDir)) {
    New-Item -ItemType Directory -Path $outDir -Force | Out-Null
    Write-Host "[创建] 输出目录: $outDir"
}

# 切换到源码目录
Set-Location $srcDir

Write-Host ""
Write-Host "[1/2] 编译 WASM 模块..."

# 编译命令
$emccCmd = @(
    "wasm_main.cpp",
    "Mesh.cpp",
    "MeshIO.cpp",
    "Lscm.cpp",
    "Parameterization.cpp",
    "QcError.cpp",
    "Solver.cpp",
    "Utils.cpp",
    "Vertex.cpp",
    "Edge.cpp",
    "Face.cpp",
    "HalfEdge.cpp",
    "-I./deps",
    "-I./deps/Eigen",
    "-I.",
    "-s", "MODULARIZE=1",
    "-s", "EXPORT_NAME='LSCMSolver'",
    "-s", "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','_malloc','_free']",
    "-s", "ALLOW_MEMORY_GROWTH=1",
    "-s", "WASM=1",
    "--bind",
    "-std=c++17",
    "-O2",
    "-o", "$outDir\lscm_solver.js"
)

$emccArgs = $emccCmd -join " "
Write-Host "执行: emcc $emccArgs"

emcc @emccCmd

if ($LASTEXITCODE -ne 0) {
    Write-Host "[失败] 编译出错 (退出码: $LASTEXITCODE)"
    exit 1
}

Write-Host ""
Write-Host "[2/2] 检查输出文件..."

if (Test-Path "$outDir\lscm_solver.js") {
    Write-Host "  [OK] 已生成: $outDir\lscm_solver.js"
} else {
    Write-Host "  [警告] 未找到 lscm_solver.js"
}

if (Test-Path "$outDir\lscm_solver.wasm") {
    Write-Host "  [OK] 已生成: $outDir\lscm_solver.wasm"
} else {
    Write-Host "  [警告] 未找到 lscm_solver.wasm"
}

Write-Host ""
Write-Host "==================================="
Write-Host "编译完成！"
Write-Host "WASM 文件位于: $outDir"
Write-Host "==================================="
Write-Host ""
Write-Host "接下来请在浏览器中测试 uv-unwrap.html"
