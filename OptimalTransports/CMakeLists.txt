cmake_minimum_required(VERSION 3.16)
project(OptimalTransports LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------------------
# Dependencies: Eigen3 (headers) and geogram (replaces the former CGAL backend)
# ---------------------------------------------------------------------------
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  if(DEFINED EIGEN3_INCLUDE_DIR)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR
      "Eigen3 not found. Install it or pass -DEIGEN3_INCLUDE_DIR=/path/to/eigen")
  endif()
endif()

# geogram: prefer its CMake package; otherwise fall back to GEOGRAM_ROOT
# (a directory containing include/ and lib/).
find_package(Geogram CONFIG QUIET)
if(TARGET Geogram::geogram)
  set(GEOGRAM_TARGET Geogram::geogram)
elseif(TARGET geogram)
  set(GEOGRAM_TARGET geogram)
else()
  if(NOT DEFINED GEOGRAM_ROOT)
    message(FATAL_ERROR
      "geogram not found. Pass -DGeogram_DIR=... or -DGEOGRAM_ROOT=/path/to/geogram")
  endif()
  add_library(geogram_ext INTERFACE)
  target_include_directories(geogram_ext INTERFACE
    "${GEOGRAM_ROOT}/include"
    "${GEOGRAM_ROOT}/include/geogram1")
  find_library(GEOGRAM_LIB NAMES geogram geogram_full
    PATHS "${GEOGRAM_ROOT}/lib" REQUIRED)
  target_link_libraries(geogram_ext INTERFACE "${GEOGRAM_LIB}")
  set(GEOGRAM_TARGET geogram_ext)
endif()

set(OT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

# ---------------------------------------------------------------------------
# geometry/ : geogram-backed primitives shared by the solvers (no CGAL)
# ---------------------------------------------------------------------------
add_library(ot_geometry STATIC
  geometry/power_diagram.cpp
  SemiOT_solver/ma_power_diagram.cpp   # MA:: compatibility shim
)
target_include_directories(ot_geometry PUBLIC
  ${OT_ROOT} ${OT_ROOT}/geometry ${OT_ROOT}/SemiOT_solver)
target_link_libraries(ot_geometry PUBLIC Eigen3::Eigen ${GEOGRAM_TARGET})

# ---------------------------------------------------------------------------
# Common input/output interface (header-only): ot_problem.hpp
# ---------------------------------------------------------------------------
add_library(ot_interface INTERFACE)
target_include_directories(ot_interface INTERFACE ${OT_ROOT})
target_link_libraries(ot_interface INTERFACE Eigen3::Eigen)

# Semi-discrete solver (header-only) + its adapter, on top of the geometry lib.
add_library(ot_semi INTERFACE)
target_include_directories(ot_semi INTERFACE ${OT_ROOT} ${OT_ROOT}/SemiOT_solver)
target_link_libraries(ot_semi INTERFACE ot_interface ot_geometry)

# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------
add_executable(run_power_diagram examples/run_power_diagram.cpp)
target_link_libraries(run_power_diagram PRIVATE ot_geometry)

add_executable(run_semi_ot examples/run_semi_ot.cpp)
target_link_libraries(run_semi_ot PRIVATE ot_semi)

# Unified example calling ALL three solvers through the common ot::Solver
# interface. Needs the extra backend dependencies (stb headers for AHT,
# surface_mesh + utils for gridOT), so it is opt-in.
option(OT_BUILD_ALL_SOLVERS "Build the example that drives all three solvers" OFF)
if(OT_BUILD_ALL_SOLVERS)
  add_executable(run_all_solvers
    examples/run_all_solvers.cpp
    utils/stb_impl.cpp                 # AHT image IO impl (stb)
    gridOT_solver/otsolver_2dgrid.cpp  # grid solver
    gridOT_solver/transport_map.cpp    # grid transport map
  )
  target_include_directories(run_all_solvers PRIVATE
    ${OT_ROOT}
    ${OT_ROOT}/BenamouOT_solver
    ${OT_ROOT}/gridOT_solver
    ${OT_ROOT}/SemiOT_solver
    ${OT_ROOT}/utils)              # stb_image.h / stb_image_write.h live here
  target_link_libraries(run_all_solvers PRIVATE ot_semi)
endif()

# ---------------------------------------------------------------------------
# Notes:
#  * BenamouOT_solver (aht_ot) is Eigen-only and header-based; its adapter is
#    BenamouOT_solver/ot_adapter.hpp. Building it additionally needs the stb
#    image headers used by utils/image_ops.hpp.
#  * gridOT_solver (otmap) depends on surface_mesh + utils; wire it here once
#    those sources are present in this tree.
#  * The legacy MongeAmpere/, AHT_OT/, include/ and tests/ trees still contain
#    the old CGAL code and are superseded by this build.
# ---------------------------------------------------------------------------
