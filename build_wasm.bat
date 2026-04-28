@echo off
REM build_wasm.bat - 编译 conformal-parameterization 为 WebAssembly
REM 使用前请先运行 emsdk_env.bat 设置环境

echo ===================================
echo Conformal Parameterization WASM Build
echo ===================================

REM 检查 emscripten 编译器
where emcc >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo [错误] 未找到 emcc 编译器
    echo 请先运行: emsdk_env.bat
    exit /b 1
)

REM 设置路径
set SRC_DIR=%~dp0
set OUT_DIR=%SRC_DIR%..\..\assets\wasm
set EMSCRIPTEN_FLAGS=-s MODULARIZE=1 -s EXPORT_NAME="LSCMSolver" -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 --bind

REM 创建输出目录
if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"

echo.
echo [1/2] 编译 WASM 模块...
emcc ^
    -std=c++17 ^
    -O2 ^
    -I"%SRC_DIR%deps" ^
    -I"%SRC_DIR%deps\Eigen" ^
    -I"%SRC_DIR%." ^
    "%SRC_DIR%wasm_main.cpp" ^
    "%SRC_DIR%Mesh.cpp" ^
    "%SRC_DIR%MeshIO.cpp" ^
    "%SRC_DIR%Lscm.cpp" ^
    "%SRC_DIR%Parameterization.cpp" ^
    "%SRC_DIR%QcError.cpp" ^
    "%SRC_DIR%Solver.cpp" ^
    "%SRC_DIR%Utils.cpp" ^
    "%SRC_DIR%Vertex.cpp" ^
    "%SRC_DIR%Edge.cpp" ^
    "%SRC_DIR%Face.cpp" ^
    "%SRC_DIR%HalfEdge.cpp" ^
    %EMSCRIPTEN_FLAGS% ^
    -s "EXPORTED_RUNTIME_METHODS=['ccall', 'cwrap', 'setValue', 'getValue', 'UTF8ToString', 'stringToUTF8', 'lengthBytesUTF8', '_malloc', '_free']" ^
    -s FILESYSTEM=1 ^
    -o "%OUT_DIR%\lscm_solver.js"

if %ERRORLEVEL% NEQ 0 (
    echo [失败] 编译出错，请检查错误信息
    exit /b 1
)

echo.
echo [2/2] 复制依赖文件...
if exist "%OUT_DIR%\lscm_solver.wasm" (
    echo   已生成: %OUT_DIR%\lscm_solver.js
    echo   已生成: %OUT_DIR%\lscm_solver.wasm
) else (
    echo   [警告] 未找到 .wasm 文件，检查编译输出
)

echo.
echo ===================================
echo 编译完成！
echo WASM 文件位于: %OUT_DIR%
echo ===================================
pause
