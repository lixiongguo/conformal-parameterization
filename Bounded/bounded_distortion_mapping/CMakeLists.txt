cmake_minimum_required(VERSION 3.16)
project(bounded_distortion_mapping LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(bounded_distortion_mapping STATIC
    BoundedDistortionMapping.cpp
)

target_include_directories(bounded_distortion_mapping
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen-3.4.0
)

add_executable(run_bounded_distortion_demo
    run_bounded_distortion_demo.cpp
)

target_link_libraries(run_bounded_distortion_demo
    PRIVATE bounded_distortion_mapping
)

add_executable(run_bd_lscm_demo
    run_bd_lscm_demo.cpp
)

target_link_libraries(run_bd_lscm_demo
    PRIVATE bounded_distortion_mapping
)
