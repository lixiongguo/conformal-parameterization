@echo off
REM Principal Curvature WASM 编译脚本
REM 使用 Emscripten + libigl + Eigen 编译 principal_curvature_wasm.cpp

setlocal

REM 设置路径 (使用绝对路径)
set "EMSSDK_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\emsdk"
set "EIGEN_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\eigen-3.4.0"
set "LIBIGL_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\Libigl-Discrete-Geometry\libigl\include"
set "LIBIGL_SRC=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\Libigl-Discrete-Geometry\libigl"
set "SRC_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp"
set "OUTPUT_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\assets\wasm"

REM 激活 EMSDK 环境
call "%EMSSDK_DIR%\emsdk_env.bat" >nul 2>&1

echo ============================================
echo  Principal Curvature WASM 编译
echo ============================================
echo.
echo 源文件: %SRC_DIR%\principal_curvature_wasm.cpp
echo Eigen:   %EIGEN_DIR%
echo libigl:  %LIBIGL_DIR%
echo 输出到:  %OUTPUT_DIR%
echo.

REM 创建输出目录
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

REM 执行编译
em++ "%SRC_DIR%\principal_curvature_wasm.cpp" ^
    -o "%OUTPUT_DIR%\principal_curvature.js" ^
    -O2 -std=c++17 ^
    -I"%EIGEN_DIR%" ^
    -I"%LIBIGL_DIR%" ^
    -s MODULARIZE=1 ^
    -s EXPORT_NAME="PrincipalCurvatureSolver" ^
    -s ALLOW_MEMORY_GROWTH=1 ^
    -s TOTAL_MEMORY=512MB ^
    -s WASM=1 ^
    -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap', 'getValue', 'setValue', 'writeArrayToMemory', 'HEAPF32']" ^
    -s EXPORTED_FUNCTIONS="['_malloc', '_free', '_compute_principal_curvature', '_get_result_buffer_size']" ^
    -s FILESYSTEM=0 ^
    -s EXIT_RUNTIME=1 ^
    -s "ENVIRONMENT=web" ^
    -fexceptions

if %ERRORLEVEL% EQU 0 (
    echo.
    echo ============================================
    echo  编译成功!
    echo ============================================
    echo  输出文件:
    echo    %OUTPUT_DIR%\principal_curvature.js
    echo    %OUTPUT_DIR%\principal_curvature.wasm
) else (
    echo.
    echo ============================================
    echo  编译失败! 请检查错误信息。
    echo ============================================
)

endlocal
