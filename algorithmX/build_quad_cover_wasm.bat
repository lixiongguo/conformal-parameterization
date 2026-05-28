@echo off
REM QuadCover WASM 编译脚本
REM 使用 Emscripten 将 quad_cover_complete.cpp 编译为 WebAssembly

setlocal

REM 设置路径 (使用绝对路径)
set "EMSSDK_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\emsdk"
set "EIGEN_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\eigen-3.4.0"
set "LIBIGL_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\Libigl-Discrete-Geometry\libigl\include"
set "SRC_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\Libigl-Discrete-Geometry"
set "OUTPUT_DIR=c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\assets\wasm"

REM 激活 EMSDK 环境
call "%EMSSDK_DIR%\emsdk_env.bat" >nul 2>&1

echo ============================================
echo  QuadCover WASM 编译
echo ============================================
echo.
echo 源文件: %SRC_DIR%\quad_cover_complete.cpp
echo Eigen:   %EIGEN_DIR%
echo libigl:  %LIBIGL_DIR%
echo 输出到:  %OUTPUT_DIR%
echo.

REM 创建输出目录
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

REM 执行编译
em++ "%SRC_DIR%\quad_cover_complete.cpp" ^
    -o "%OUTPUT_DIR%\quad_cover.js" ^
    -O2 -std=c++17 ^
    -I"%EIGEN_DIR%" ^
    -I"%LIBIGL_DIR%" ^
    --bind ^
    -s MODULARIZE=1 ^
    -s EXPORT_NAME="QuadCoverSolver" ^
    -s ALLOW_MEMORY_GROWTH=1 ^
    -s TOTAL_MEMORY=512MB ^
    -s WASM=1 ^
    -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap', 'getValue', 'setValue', 'FS', 'writeArrayToMemory']" ^
    -s EXPORTED_FUNCTIONS="['_malloc', '_free', '_main']" ^
    -s FILESYSTEM=1 ^
    -s EXIT_RUNTIME=1 ^
    -s "ENVIRONMENT=web"

if %ERRORLEVEL% EQU 0 (
    echo.
    echo ============================================
    echo  编译成功!
    echo ============================================
    echo  输出文件:
    echo    %OUTPUT_DIR%\quad_cover.js
    echo    %OUTPUT_DIR%\quad_cover.wasm
) else (
    echo.
    echo ============================================
    echo  编译失败! 请检查错误信息。
    echo ============================================
)

endlocal
