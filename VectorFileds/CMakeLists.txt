cmake_minimum_required(VERSION 3.16)
project(VectorFields LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MESH_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../conformal-parameterization")
set(EIGEN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen-3.3.9")
if(NOT EXISTS "${EIGEN_ROOT}")
  set(EIGEN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen-3.4.0")
endif()

set(MESH_SOURCES
  "${MESH_ROOT}/Mesh.cpp"
  "${MESH_ROOT}/MeshIO.cpp"
  "${MESH_ROOT}/Parameterization.cpp"
  "${MESH_ROOT}/geometry/QcError.cpp"
  "${MESH_ROOT}/Vertex.cpp"
  "${MESH_ROOT}/Edge.cpp"
  "${MESH_ROOT}/Face.cpp"
  "${MESH_ROOT}/HalfEdge.cpp"
  "${MESH_ROOT}/uv_unwrap_simple/Lscm.cpp"
)

set(MESH_INCLUDES
  "${MESH_ROOT}"
  "${MESH_ROOT}/geometry"
  "${MESH_ROOT}/uv_unwrap_simple"
  "${EIGEN_ROOT}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

add_library(vector_field_core
  ComplexPolyField.cpp
)
target_include_directories(vector_field_core PUBLIC ${MESH_INCLUDES})
target_compile_definitions(vector_field_core PRIVATE _USE_MATH_DEFINES)
if(MSVC)
  target_compile_options(vector_field_core PRIVATE /utf-8 /Zc:__cplusplus)
endif()

function(add_vf_executable target main_source)
  add_executable(${target} ${main_source} ${MESH_SOURCES})
  target_link_libraries(${target} PRIVATE vector_field_core)
  target_include_directories(${target} PRIVATE ${MESH_INCLUDES})
  target_compile_definitions(${target} PRIVATE _USE_MATH_DEFINES)
  if(MSVC)
    target_compile_options(${target} PRIVATE /utf-8 /Zc:__cplusplus)
  endif()
endfunction()

add_vf_executable(trivial_connection trivial_connection.cpp)
add_vf_executable(nrosy_trivial_connection nrosy_trivial_connection.cpp)
add_vf_executable(complex_poly_field complex_poly_field_main.cpp)
