# ============================================================
#  conformal-parameterization WASM 构建（无需 GNU make）
#  等价于同目录 Makefile；用法：
#    .\build.ps1                         # 默认 = build-simple-compat
#    .\build.ps1 -Target build-simple-compat
#    .\build.ps1 -Target build-page-all
# ============================================================
param(
    [string]$Target = "build-simple-compat"
)

$ErrorActionPreference = "Stop"

$WasmDir  = Split-Path -Parent $MyInvocation.MyCommand.Path   # .../conformal-parameterization/wasm
$PkgRoot  = Split-Path -Parent $WasmDir                       # .../conformal-parameterization
$CppRoot  = Split-Path -Parent $PkgRoot                       # .../cpp
$RepoRoot = Split-Path -Parent $CppRoot                       # repo 根目录
$OutDir   = Join-Path $RepoRoot "assets\wasm"
$Eigen    = Join-Path $CppRoot  "deps\eigen-3.4.0"

$EmsdkEnv = Join-Path $CppRoot "emsdk\emsdk_env.ps1"
if (Test-Path $EmsdkEnv) { $env:EMSDK_QUIET = "1"; . $EmsdkEnv }

# 注意：PowerShell 下直接调用 `em++` 可能命中 em++.ps1（内部 Remove-Item 在部分环境会抛
# `Cannot find 'Env:\_PYTHON_SYSCONFIGDATA_NAME'`）。这里显式指向 em++.bat 规避该问题。
$Emxx = Join-Path $CppRoot "emsdk\upstream\emscripten\em++.bat"
if (-not (Test-Path $Emxx)) {
    if ($env:EMSDK) { $Emxx = Join-Path $env:EMSDK "upstream\emscripten\em++.bat" }
}
if (-not (Test-Path $Emxx)) {
    throw "找不到 em++.bat。请先在 cpp\emsdk 执行 .\emsdk activate latest"
}
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

function P([string]$rel) { Join-Path $PkgRoot $rel }
# 注意：不要命名为 H —— PowerShell 里 H 是 Get-History 的别名，会覆盖同名函数
function Hp([string]$rel) { Join-Path $CppRoot $rel }

# ---------- 源文件组（对应 Makefile） ----------
$BASEMESH = @(
    "BaseMesh/Mesh.cpp","BaseMesh/MeshIO.cpp","BaseMesh/Vertex.cpp","BaseMesh/Edge.cpp",
    "BaseMesh/Face.cpp","BaseMesh/HalfEdge.cpp","BaseMesh/GaussianCurvature.cpp",
    "BaseMesh/QcError.cpp","Parameterization/Parameterization.cpp") | ForEach-Object { P $_ }

$SOLVER = @(
    "Solvers/Solver.cpp","Solvers/AugmentedLagrangian.cpp","Solvers/MixedIntegerProgram.cpp"
) | ForEach-Object { P $_ }

$LSCM = @(
    "Parameterization/SimpleParam/LSCM/Lscm.cpp",
    "Parameterization/SimpleParam/LSCM/Scp.cpp") | ForEach-Object { P $_ }

$CUTSEAM = @( (P "CutSeamMesh/CutSeamMesh.cpp") )

$CONE = @(
    "Parameterization/ConeParam/ConeParameterization.cpp",
    "Parameterization/ConeParam/Conformal/CirclePatterns.cpp",
    "Parameterization/ConeParam/Conformal/Cetm.cpp",
    "Parameterization/ConeParam/RicciFlow/RicciFlow.cpp") | ForEach-Object { P $_ }

# ---------- include 路径（对应 Makefile 的 *_INC） ----------
$IncBase   = @("-I$Eigen", "-I$(P 'BaseMesh')", "-I$(P 'Parameterization')")

$IncSimple = $IncBase + @(
    "-I$(P 'Parameterization/SimpleParam/LSCM')",
    "-I$(P 'Parameterization/SimpleParam/ABF')",
    "-I$(P 'Parameterization/SimpleParam/ARAP')",
    "-I$(Hp 'HGP/FastHGP')",
    "-I$(P 'CutSeamMesh')",
    "-I$(P 'Parameterization/ConeParam')",
    "-I$(P 'Parameterization/ConeParam/Conformal')",
    "-I$(P 'Parameterization/ConeParam/RicciFlow')",
    "-I$(P 'Parameterization/ConeParam/IncrementalFlattenning')",
    "-I$(P 'Solvers')","-I$(P 'Solvers/Mosek')")

$IncField  = $IncBase + @(
    "-I$(P 'CutSeamMesh')",
    "-I$(P 'Parameterization/GlobalFieldsParam')",
    "-I$(P 'Parameterization/CutSeamParam')",
    "-I$(P 'Parameterization/CutSeamParam/HoloOneForm')",
    "-I$(P 'Parameterization/SimpleParam/LSCM')",
    "-I$(P 'VectorFileds')",
    "-I$(P 'VectorFileds/PrincipalCurvatureFields')",
    "-I$(P 'Solvers')","-I$(P 'Solvers/Mosek')")

$IncAbel   = $IncBase + @(
    "-I$(P 'CutSeamMesh')",
    "-I$(P 'Parameterization/CutSeamParam')",
    "-I$(P 'Parameterization/CutSeamParam/Abel_Jacoi')",
    "-I$(P 'Parameterization/SimpleParam/LSCM')")

$IncBdLscm = $IncBase + @(
    "-I$(P 'Parameterization/SimpleParam/LSCM')",
    "-I$(P 'Bounded/bounded_distortion_mapping')")

# ---------- 链接参数（对应 EM_FLAGS / EM_FLAGS_COMPAT） ----------
# -DEIGEN_STACK_ALLOCATION_LIMIT=0 + STACK_SIZE 是修复 "memory access out of bounds" 的关键
$Common = @("-O3","-std=c++17","-Wall","-Wextra","-DEIGEN_STACK_ALLOCATION_LIMIT=0")

$Compat = @(
    "-s","MODULARIZE=1","-s","ALLOW_MEMORY_GROWTH=1","-s","WASM=1","-s","STACK_SIZE=1048576",
    "-sEXPORTED_FUNCTIONS=['_malloc','_free']",
    "-sEXPORTED_RUNTIME_METHODS=['ccall','cwrap','getValue','setValue']")

$Embind = @(
    "--bind","-s","MODULARIZE=1","-s","ALLOW_MEMORY_GROWTH=1","-s","WASM=1","-s","STACK_SIZE=1048576",
    "-s","FORCE_FILESYSTEM=1",
    "-sEXPORTED_RUNTIME_METHODS=['ccall','cwrap','getValue','setValue','FS']")

function Invoke-Build {
    param([string]$Name,[string]$OutJs,[string[]]$Sources,[string[]]$Includes,
          [string[]]$Flags,[string]$ExportName)
    Write-Host "[$Name] em++ -> $OutJs" -ForegroundColor Cyan
    & $Emxx @Common @Includes @Sources -o $OutJs @Flags "-sEXPORT_NAME=$ExportName"
    if ($LASTEXITCODE -ne 0) { throw "[$Name] 编译失败 (exit $LASTEXITCODE)" }
    Write-Host "[$Name] OK" -ForegroundColor Green
}

function Build-Simple {
    Invoke-Build -Name "uv_unwrap_simple" -OutJs (Join-Path $OutDir "uv_unwrap_simple.js") `
        -Includes $IncSimple -Flags $Compat -ExportName "UvUnwrapSimpleSolver" -Sources (@(
            (P "wasm/bindings_uv_unwrap_simple_compat.cpp")) + $BASEMESH + $SOLVER + $LSCM + @(
            (P "Parameterization/SimpleParam/ABF/AbfPlusPlus.cpp"),
            (P "Parameterization/SimpleParam/ABF/LinAbf.cpp"),
            (P "Parameterization/SimpleParam/ARAP/ARAP.cpp"),
            (P "Parameterization/SimpleParam/ARAP/Tutte.cpp"),
            (Hp "HGP/FastHGP/FastHGPSimple.cpp"),
            (P "CutSeamMesh/CutSeamMesh.cpp"),
            (P "Parameterization/ConeParam/IncrementalFlattenning/IncrementalFlattening.cpp")) + $CONE)
}

function Build-Dgp {
    Invoke-Build -Name "dgp_basic" -OutJs (Join-Path $OutDir "dgp_basic.js") `
        -Includes $IncBase -Flags $Compat -ExportName "DgpBasicSolver" -Sources (@(
            (P "wasm/bindings_dgp_basic_compat.cpp")) + $BASEMESH)
}

function Build-Field {
    Invoke-Build -Name "uv_unwrap_field" -OutJs (Join-Path $OutDir "uv_unwrap_field.js") `
        -Includes $IncField -Flags $Compat -ExportName "UvUnwrapFieldSolver" -Sources (@(
            (P "wasm/bindings_uv_unwrap_field_compat.cpp")) + $BASEMESH + $CUTSEAM + $SOLVER + $LSCM + @(
            (P "Parameterization/GlobalFieldsParam/GlobalFieldsParameterization.cpp"),
            (P "Parameterization/GlobalFieldsParam/QuadCover.cpp"),
            (P "Parameterization/GlobalFieldsParam/MIQQuad.cpp"),
            (P "Parameterization/CutSeamParam/CutSeamParameterization.cpp"),
            (P "Parameterization/CutSeamParam/HoloOneForm/HolomorphicOneForm.cpp"),
            (P "VectorFileds/NRosyVectorFields.cpp"),
            (P "VectorFileds/PrincipalCurvatureFields/PrincipalCurvatureField.cpp")))
}

function Build-Abel {
    Invoke-Build -Name "uv_unwrap_abel_jacobi" -OutJs (Join-Path $OutDir "uv_unwrap_abel_jacobi.js") `
        -Includes $IncAbel -Flags $Compat -ExportName "UvUnwrapAbelJacobiSolver" -Sources (@(
            (P "wasm/bindings_uv_unwrap_abel_jacobi_compat.cpp")) + $BASEMESH + $CUTSEAM + $LSCM + @(
            (P "Parameterization/CutSeamParam/CutSeamParameterization.cpp"),
            (P "Parameterization/CutSeamParam/Abel_Jacoi/AbelJacobi.cpp"),
            (P "Parameterization/CutSeamParam/Abel_Jacoi/AbelJacobiParameterization.cpp")))
}

function Build-BdLscm {
    Invoke-Build -Name "bd_lscm" -OutJs (Join-Path $OutDir "bd_lscm.js") `
        -Includes $IncBdLscm -Flags $Embind -ExportName "BD_LSCMModule" -Sources (@(
            (P "wasm/bindings_bd_lscm.cpp")) + $BASEMESH + $LSCM + @(
            (P "Bounded/bounded_distortion_mapping/BoundedDistortionMapping.cpp")))
}

function Build-Ruppert {
    Invoke-Build -Name "ruppert" -OutJs (Join-Path $OutDir "ruppert.js") `
        -Includes @() -Flags $Compat -ExportName "RuppertModule" -Sources @(
            (P "wasm/bindings_ruppert_compat.cpp"))
}

switch ($Target) {
    "build-simple-compat" { Build-Simple }
    "build-dgp-compat"    { Build-Dgp }
    "build-field-compat"  { Build-Field }
    "build-abel-compat"   { Build-Abel }
    "build-bd-lscm"       { Build-BdLscm }
    "build-ruppert"       { Build-Ruppert }
    "build-page-all"      { Build-Simple; Build-Dgp; Build-Field; Build-Abel; Build-BdLscm }
    "page-all"            { Build-Simple; Build-Dgp; Build-Field; Build-Abel; Build-BdLscm }
    "page"                { Build-Simple; Build-Dgp; Build-Field; Build-Abel; Build-BdLscm }
    default               { throw "未知 Target: $Target" }
}
Write-Host "全部完成 → $OutDir" -ForegroundColor Green
