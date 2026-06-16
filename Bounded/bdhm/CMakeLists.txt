cmake_minimum_required(VERSION 3.16)
project(bounded_distortion_harmonic_map LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(
    ${CMAKE_CURRENT_SOURCE_DIR}/../cauchy_coordinates
    ${CMAKE_CURRENT_BINARY_DIR}/cauchy_coordinates_build
)

add_library(bdhm STATIC
    BoundedDistortionHarmonicMap.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Solvers/Solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Solvers/Lbfgs.cpp
)

target_include_directories(bdhm
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen-3.4.0
        ${CMAKE_CURRENT_SOURCE_DIR}/../Solvers
        ${CMAKE_CURRENT_SOURCE_DIR}/../Solvers/Mosek
        ${CMAKE_CURRENT_SOURCE_DIR}/../conformal-parameterization/BaseMesh
)

target_link_libraries(bdhm
    PUBLIC cauchy_coordinates
)

add_executable(run_bdhm_demo
    run_bdhm_demo.cpp
)

target_link_libraries(run_bdhm_demo
    PRIVATE bdhm
)
