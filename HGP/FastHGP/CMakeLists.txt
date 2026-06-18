cmake_minimum_required(VERSION 3.16)
project(FastHGP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(FASTHGP_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(HGP_ROOT ${FASTHGP_ROOT}/../HGP)

# ---------------------------------------------------------------------------
# Eigen3
# ---------------------------------------------------------------------------
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  if(DEFINED EIGEN3_INCLUDE_DIR)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
  elseif(EXISTS "${FASTHGP_ROOT}/../../deps/eigen-3.4.0/Eigen")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${FASTHGP_ROOT}/../../deps/eigen-3.4.0")
  elseif(EXISTS "${FASTHGP_ROOT}/../../deps/eigen-3.3.9/Eigen")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${FASTHGP_ROOT}/../../deps/eigen-3.3.9")
  else()
    message(FATAL_ERROR
      "Eigen3 not found. Install it or pass -DEIGEN3_INCLUDE_DIR=/path/to/eigen")
  endif()
endif()

# ---------------------------------------------------------------------------
# CGAL (headers + libs); override with -DCGAL_DIR=...
# ---------------------------------------------------------------------------
find_package(CGAL QUIET)
if(NOT CGAL_FOUND)
  message(WARNING
    "CGAL not found via find_package. Set CGAL_DIR or CGAL_INCLUDE_DIRS manually.")
endif()

# ---------------------------------------------------------------------------
# GMM (header-only template library)
# ---------------------------------------------------------------------------
if(DEFINED GMM_INCLUDE_DIR)
  set(GMM_INC "${GMM_INCLUDE_DIR}")
elseif(EXISTS "${FASTHGP_ROOT}/../../deps/gmm/include")
  set(GMM_INC "${FASTHGP_ROOT}/../../deps/gmm/include")
elseif(EXISTS "${HGP_ROOT}/../deps/gmm/include")
  set(GMM_INC "${HGP_ROOT}/../deps/gmm/include")
else()
  message(WARNING "GMM include path not found; set -DGMM_INCLUDE_DIR=...")
  set(GMM_INC "")
endif()

# ---------------------------------------------------------------------------
# fasthgp — pure C++ (no MATLAB Engine)
# ---------------------------------------------------------------------------
add_library(fasthgp STATIC
  FastHGP.cpp
  FastHGPNumerics.cpp
  Utils/EigenLinearSolver.cpp
  Utils/FramesFile.cpp
  ${HGP_ROOT}/HarmonicParametrization.cpp
  ${HGP_ROOT}/Borders.cpp
  ${HGP_ROOT}/Parser.cpp
)

target_compile_definitions(fasthgp PUBLIC FASTHGP_STANDALONE)

target_include_directories(fasthgp PUBLIC
  ${FASTHGP_ROOT}
  ${HGP_ROOT}
  ${GMM_INC}
)

target_link_libraries(fasthgp PUBLIC Eigen3::Eigen)

if(CGAL_FOUND)
  target_link_libraries(fasthgp PUBLIC CGAL::CGAL)
elseif(DEFINED CGAL_INCLUDE_DIRS)
  target_include_directories(fasthgp PUBLIC ${CGAL_INCLUDE_DIRS})
endif()

# Optional CLI driver (uncomment when you add tools/fasthgp_main.cpp)
# add_executable(fasthgp_run tools/fasthgp_main.cpp)
# target_link_libraries(fasthgp_run PRIVATE fasthgp)

message(STATUS "FastHGP: built static library 'fasthgp' (FASTHGP_STANDALONE, no MATLAB)")
