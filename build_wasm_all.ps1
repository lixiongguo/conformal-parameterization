# ============================================================
# Unified WASM Build Script for Conformal Parameterization
# Mirrors Makefile (11 targets). Run from this directory.
#
# Usage:
#   .\build_wasm_all.ps1                    # build all 11 targets
#   .\build_wasm_all.ps1 -Targets lscm      # build LSCM only
#   .\build_wasm_all.ps1 -Targets lscm,tutte_arap  # multiple targets
#   .\build_wasm_all.ps1 -Skip libigl        # skip libigl-dependent targets
# ============================================================
param(
    [string[]]$Targets = @(),           # specific targets to build (empty = all)
    [switch]$SkipLibigl = $false,       # skip principal_curvature (needs libigl)
    [switch]$SkipBd = $false,           # skip bd_lscm
    [switch]$Clean = $false             # clean before building
)

$ErrorActionPreference = "Stop"

# ==================== Paths ====================
$emsdkPath     = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\emsdk"
$srcDir        = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\cpp\conformal-parameterization"
$outDir        = "c:\Users\lixio\OneDrive\Desktop\MyDoc\lixiongguo.github.io\assets\wasm"
$buildDir      = "$srcDir\build"
$eigenInc      = "$srcDir\..\deps\eigen-3.4.0"
$libiglInc     = "$srcDir\..\Libigl-Discrete-Geometry\libigl\include"
$glmInc        = "$srcDir\..\deps\glm"
$mosekInc      = "$srcDir\..\deps\mosek"

# ==================== Activate Emscripten ====================
if (!(Test-Path "$emsdkPath\emsdk_env.ps1")) {
    Write-Host "[ERROR] emsdk not found at: $emsdkPath" -ForegroundColor Red
    exit 1
}
Push-Location $emsdkPath
. ".\emsdk_env.ps1" | Out-Null
Pop-Location

# Verify em++
$emxx = Get-Command "em++.bat" -ErrorAction SilentlyContinue
if (!$emxx) { $emxx = Get-Command "em++" -ErrorAction SilentlyContinue }
if (!$emxx) {
    Write-Host "[ERROR] em++ not found. Please run emsdk_env.ps1 first." -ForegroundColor Red
    exit 1
}
Write-Host "[OK] em++ found: $($emxx.Source)" -ForegroundColor Green

# ==================== Check dependencies ====================
$checks = @{
    "Eigen"      = $eigenInc
    "glm"        = $glmInc
    "mosek"      = $mosekInc
    "libigl"     = $libiglInc
}
foreach ($name in $checks.Keys) {
    $path = $checks[$name]
    if (Test-Path $path) {
        Write-Host "[OK] $name found: $path" -ForegroundColor Green
    } else {
        $color = if ($name -eq "libigl" -or $name -eq "mosek") { "Yellow" } else { "Red" }
        Write-Host "[WARN] $name NOT found: $path" -ForegroundColor $color
        if ($color -eq "Red") { exit 1 }
    }
}

Set-Location $srcDir

# ==================== Shared Settings ====================
$cxxFlags     = "-std=c++17 -O2 -flto"
$incFlags     = "-I$eigenInc -I. -I$glmInc -I$mosekInc"
$meshSrcs     = "Mesh.cpp MeshIO.cpp Parameterization.cpp Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp"
$solverSrcs   = "Solver.cpp QcError.cpp"
$allHdrs      = Get-ChildItem -Path $srcDir -Filter "*.h" | ForEach-Object { $_.Name }

$commonEmFlags = @(
    "-s MODULARIZE=1"
    "-s ALLOW_MEMORY_GROWTH=1"
    "-s WASM=1"
    '-s "EXPORTED_RUNTIME_METHODS=[''ccall'',''cwrap'',''UTF8ToString'',''getValue'',''setValue'']"'
    "-s FORCE_FILESYSTEM=0"
    '-s ENVIRONMENT="web"'
)

# ==================== Ensure output directories ====================
if ($Clean) {
    Write-Host "Cleaning previous builds..." -ForegroundColor Yellow
    Remove-Item -Recurse -Force $buildDir -ErrorAction SilentlyContinue
}
New-Item -ItemType Directory -Path $outDir  -Force | Out-Null
New-Item -ItemType Directory -Path $buildDir -Force | Out-Null

# ==================== Helper: write efile (JSON array for EXPORTED_FUNCTIONS) ====================
function Write-Efile($name, $funcs) {
    $json = "[" + (($funcs | ForEach-Object { """$_""" }) -join ",") + "]"
    $json | Out-File -FilePath "$buildDir\$name.efile" -Encoding ascii -NoNewline
    return "@$buildDir\$name.efile"
}

# ==================== Helper: build one target ====================
function Build-Target($name, $srcs, $exportName, $funcs, $output, $extraFlags, $customInc) {
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host "Building: $name" -ForegroundColor Cyan
    Write-Host "  => $output" -ForegroundColor Cyan
    Write-Host "========================================" -ForegroundColor Cyan

    $efile = Write-Efile $name $funcs
    $useInc = if ($customInc) { $customInc } else { $incFlags }

    $allFlags = @(
        $cxxFlags
        $useInc
        $commonEmFlags
        "-s EXPORT_NAME=""$exportName"""
        "-s EXPORTED_FUNCTIONS=$efile"
    )
    if ($extraFlags) { $allFlags += $extraFlags }

    $cmdArgs = $allFlags + ($srcs -split '\s+') + @("-o", $output)

    $cmdLine = "em++ " + ($cmdArgs -join " ")
    Write-Host "[CMD] $cmdLine" -ForegroundColor DarkGray

    $proc = Start-Process -FilePath "em++.bat" -ArgumentList $cmdArgs -NoNewWindow -Wait -PassThru
    if ($proc.ExitCode -ne 0) {
        Write-Host "[FAIL] $name — exit code $($proc.ExitCode)" -ForegroundColor Red
        return $false
    }
    Write-Host "[DONE] $name -> $output" -ForegroundColor Green
    return $true
}

# ==================== Build Targets ====================
$allOk = $true
$totalTime = 0

function Build-One($name, $condition, $scriptBlock) {
    if (-not $condition) { return $true }
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    $result = & $scriptBlock
    $sw.Stop()
    $script:totalTime += $sw.Elapsed.TotalSeconds
    if (-not $result) { $script:allOk = $false }
    return $result
}

$shouldBuild = ($Targets.Count -eq 0)
function ShouldBuild($name) {
    if ($shouldBuild) { return $true }
    return $Targets -contains $name
}

# ------------------------------------------------------------------
# 1. LSCM
# ------------------------------------------------------------------
Build-One "LSCM" (ShouldBuild "lscm") {
    Build-Target "lscm" `
        "wasm_main.cpp Lscm.cpp $solverSrcs $meshSrcs" `
        "LSCMSolver" `
        @("_malloc","_free","_solve_lscm","_get_uv_result","_get_uv_result_size","_get_last_time_ms","_dispose") `
        "$outDir\lscm_solver.js"
}

# ------------------------------------------------------------------
# 2. Tutte + ARAP
# ------------------------------------------------------------------
Build-One "Tutte+ARAP" (ShouldBuild "tutte_arap") {
    Build-Target "tutte_arap" `
        "wasm_tutte_arap.cpp Tutte.cpp ARAP.cpp QcError.cpp $meshSrcs" `
        "TutteARAPSolver" `
        @("_malloc","_free","_solve_tutte_circle","_solve_tutte_square","_solve_arap","_get_ta_uv_result","_get_ta_uv_result_size","_get_ta_last_time_ms","_compute_qc_error","_load_mesh_with_uv","_get_qc_errors","_get_qc_errors_size","_get_qc_colors","_get_qc_colors_size","_ta_dispose") `
        "$outDir\tutte_arap_solver.js"
}

# ------------------------------------------------------------------
# 3. CETM
# ------------------------------------------------------------------
Build-One "CETM" (ShouldBuild "cetm") {
    Build-Target "cetm" `
        "wasm_cetm.cpp Cetm.cpp $solverSrcs $meshSrcs" `
        "CETMSolver" `
        @("_malloc","_free","_solve_cetm","_get_cetm_uv_result","_get_cetm_uv_result_size","_get_cetm_last_time_ms","_cetm_dispose") `
        "$outDir\cetm_solver.js"
}

# ------------------------------------------------------------------
# 4. Circle Patterns
# ------------------------------------------------------------------
Build-One "Circle Patterns" (ShouldBuild "circle_patterns") {
    Build-Target "circle_patterns" `
        "wasm_circle_patterns.cpp CirclePatternsWasm.cpp $solverSrcs $meshSrcs" `
        "CPSolver" `
        @("_malloc","_free","_solve_cp","_get_cp_uv_result","_get_cp_uv_result_size","_get_cp_last_time_ms","_cp_dispose") `
        "$outDir\cp_solver.js"
}

# ------------------------------------------------------------------
# 5. Holomorphic 1-Form
# ------------------------------------------------------------------
Build-One "Holo" (ShouldBuild "holo") {
    Build-Target "holo" `
        "wasm_holo.cpp HolomorphicOneForm.cpp $meshSrcs" `
        "HoloSolver" `
        @("_malloc","_free","_solve_holo","_get_holo_uv_result","_get_holo_uv_result_size","_get_holo_last_time_ms","_holo_dispose") `
        "$outDir\holo_solver.js"
}

# ------------------------------------------------------------------
# 6. MIQ Quad
# ------------------------------------------------------------------
Build-One "MIQ" (ShouldBuild "miq") {
    Build-Target "miq" `
        "wasm_miq.cpp MIQQuad.cpp $meshSrcs" `
        "MIQSolver" `
        @("_malloc","_free","_solve_miq","_get_miq_uv_result","_get_miq_uv_result_size","_get_miq_last_time_ms","_miq_dispose") `
        "$outDir\miq_solver.js"
}

# ------------------------------------------------------------------
# 7. QuadCover
# ------------------------------------------------------------------
Build-One "QuadCover" (ShouldBuild "quadcover") {
    Build-Target "quadcover" `
        "wasm_quadcover.cpp QuadCover.cpp $meshSrcs" `
        "QuadCoverSolver" `
        @("_malloc","_free","_solve_qc","_get_qc_uv_result","_get_qc_uv_result_size","_get_qc_last_time_ms","_qc_dispose") `
        "$outDir\quadcover_solver.js"
}

# ------------------------------------------------------------------
# 8. Ricci Flow
# ------------------------------------------------------------------
Build-One "Ricci" (ShouldBuild "ricci") {
    Build-Target "ricci" `
        "wasm_ricci.cpp RicciFlow.cpp $solverSrcs $meshSrcs" `
        "RicciFlowSolver" `
        @("_malloc","_free","_solve_ricci","_get_ricci_uv_result","_get_ricci_uv_result_size","_get_ricci_last_time_ms","_ricci_dispose") `
        "$outDir\ricci_solver.js"
}

# ------------------------------------------------------------------
# 9. BFF
# ------------------------------------------------------------------
Build-One "BFF" (ShouldBuild "bff") {
    Build-Target "bff" `
        "wasm_bff.cpp BFF.cpp $meshSrcs" `
        "BFFSolver" `
        @("_malloc","_free","_solve_bff","_get_bff_uv_result","_get_bff_uv_result_size","_get_bff_last_time_ms","_bff_dispose") `
        "$outDir\bff_solver.js"
}

# ------------------------------------------------------------------
# 10. BD-LSCM (Bounded Distortion LSCM)
# ------------------------------------------------------------------
Build-One "BD-LSCM" ((ShouldBuild "bd_lscm") -and (-not $SkipBd)) {
    Build-Target "bd_lscm" `
        "wasm_bd_lscm.cpp" `
        "BDLSCMSolver" `
        @("_malloc","_free","_solve_bd_lscm","_get_bd_uv_result","_get_bd_uv_result_size","_get_bd_last_time_ms","_bd_dispose") `
        "$outDir\bd_lscm_solver.js"
}

# ------------------------------------------------------------------
# 11. Principal Curvature (requires libigl)
# ------------------------------------------------------------------
Build-One "Principal Curvature" ((ShouldBuild "principal_curvature") -and (-not $SkipLibigl) -and (Test-Path $libiglInc)) {
    $pcIncFlags = "-I$eigenInc -I$libiglInc -I."
    Build-Target "principal_curvature" `
        "principal_curvature_wasm.cpp" `
        "PrincipalCurvatureSolver" `
        @("_malloc","_free","_compute_principal_curvature","_get_result_buffer_size") `
        "$outDir\principal_curvature.js" `
        @("--bind") `
        $pcIncFlags
}

# ==================== Summary ====================
Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
if ($allOk) {
    Write-Host "BUILD COMPLETE" -ForegroundColor Green
    Write-Host "Total time: $([math]::Round($totalTime, 1))s" -ForegroundColor Green
} else {
    Write-Host "BUILD FINISHED WITH ERRORS" -ForegroundColor Red
}
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# List output files
Write-Host "Output files in ${outDir}:" -ForegroundColor Yellow
Get-ChildItem -Path $outDir -Filter "*.js"  | ForEach-Object { 
    $size = "{0,8:N1} KB" -f ($_.Length / 1024)
    Write-Host "  $size  $($_.Name)" -ForegroundColor White
}
Get-ChildItem -Path $outDir -Filter "*.wasm" | ForEach-Object { 
    $size = "{0,8:N1} KB" -f ($_.Length / 1024)
    Write-Host "  $size  $($_.Name)" -ForegroundColor DarkGray
}

if (-not $allOk) { exit 1 }
