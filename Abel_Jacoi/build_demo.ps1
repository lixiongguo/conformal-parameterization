# Build Abel-Jacobi demo (native g++)
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$srcRoot = Split-Path -Parent $here
$cppRoot = Split-Path -Parent $srcRoot
$eigenInc = Join-Path $cppRoot "deps\eigen-3.4.0"

$sources = @(
    (Join-Path $here "abel_jacobi_demo.cpp"),
    (Join-Path $here "AbelJacobi.cpp"),
    (Join-Path $here "AbelJacobiParameterization.cpp"),
    (Join-Path $srcRoot "Mesh.cpp"),
    (Join-Path $srcRoot "MeshIO.cpp"),
    (Join-Path $srcRoot "Parameterization.cpp"),
    (Join-Path $srcRoot "Vertex.cpp"),
    (Join-Path $srcRoot "Edge.cpp"),
    (Join-Path $srcRoot "Face.cpp"),
    (Join-Path $srcRoot "HalfEdge.cpp"),
    (Join-Path $srcRoot "Solver.cpp"),
    (Join-Path $srcRoot "QcError.cpp"),
    (Join-Path $srcRoot "uv_unwrap_simple\Lscm.cpp")
)

$out = Join-Path $here "abel_jacobi_demo.exe"
$args = @(
    "-std=c++17", "-O2",
    "-I$srcRoot",
    "-I$here",
    "-I$eigenInc"
) + $sources + @("-o", $out)

Write-Host "Building $out ..."
& g++ @args
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Done: $out"
