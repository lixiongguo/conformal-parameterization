@echo off
REM build_wasm_simple.bat - 简单编译脚本

echo ===================================
echo Conformal Parameterization WASM Build
echo ===================================

REM 激活 Emscripten 环境
call "%~dp0..\..\cpp\emsdk\emsdk_env.bat"

REM 检查 emcc
where emcc >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo [错误] 未找到 emcc 编译器
    exit /b 1
)

echo [OK] 找到 emcc 编译器
emcc --version | findstr /i "emcc"

REM 设置路径
set SRC_DIR=%~dp0
set OUT_DIR=%SRC_DIR%..\..\assets\wasm

REM 创建输出目录
if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"

echo.
echo [1/2] 编译 WASM 模块...

cd /d "%SRC_DIR%"

emcc ^
    wasm_main.cpp ^
    Mesh.cpp ^
    MeshIO.cpp ^
    Lscm.cpp ^
    Parameterization.cpp ^
    QcError.cpp ^
    Solver.cpp ^
    Utils.cpp ^
    Vertex.cpp ^
    Edge.cpp ^
    Face.cpp ^
    HalfEdge.cpp ^
    -I./deps ^
    -I./deps/Eigen ^
    -I. ^
    -s MODULARIZE=1 ^
    -s EXPORT_NAME="LSCMSolver" ^
    -s "EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','_malloc','_free']" ^
    -s ALLOW_MEMORY_GROWTH=1 ^
    -s WASM=1 ^
    --bind ^
    -std=c++17 ^
    -O2 ^
    -o "%OUT_DIR%\lscm_solver.js"

if %ERRORLEVEL% NEQ 0 (
    echo [失败] 编译出错
    exit /b 1
)

echo.
echo [2/2] 检查输出文件...

if exist "%OUT_DIR%\lscm_solver.js" (
    echo [OK] 已生成: %OUT_DIR%\lscm_solver.js
) else (
    echo [警告] 未找到 lscm_solver.js
)

if exist "%OUT_DIR%\lscm_solver.wasm" (
    echo [OK] 已生成: %OUT_DIR%\lscm_solver.wasm
) else (
    echo [警告] 未找到 lscm_solver.wasm
)

echo.
echo ===================================
echo 编译完成！
echo WASM 文件位于: %OUT_DIR%
echo ===================================
pause
