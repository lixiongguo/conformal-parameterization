cmake_minimum_required(VERSION 3.16)
project(cauchy_coordinates LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(cauchy_coordinates STATIC
    CauchyCoordinates.cpp
)

target_include_directories(cauchy_coordinates
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen-3.4.0
)

add_executable(run_cauchy_coordinates_demo
    run_cauchy_coordinates_demo.cpp
)

target_link_libraries(run_cauchy_coordinates_demo
    PRIVATE cauchy_coordinates
)
