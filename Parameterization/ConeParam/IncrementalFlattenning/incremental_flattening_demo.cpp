/**
 * Demo: incremental flattening global parameterization with automatic cone selection.
 *
 * Build (from repo cpp/):
 *   g++ -std=c++17 -O2 -I conformal-parameterization -I deps/eigen-3.3.9 \
 *       conformal-parameterization/IncrementalFlattenning/incremental_flattening_demo.cpp \
 *       conformal-parameterization/IncrementalFlattenning/IncrementalFlattening.cpp \
 *       conformal-parameterization/Mesh.cpp MeshIO.cpp Parameterization.cpp \
 *       conformal-parameterization/Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp \
 *       conformal-parameterization/geometry/GaussianCurvature.cpp \
 *       conformal-parameterization/geometry/QcError.cpp \
 *       conformal-parameterization/uv_unwrap_simple/Tutte.cpp \
 *       -o incremental_flattening_demo
 */
#include "IncrementalFlattening.h"
#include <iostream>

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj>\n";
        return 1;
    }

    Mesh mesh;
    if (!mesh.read(argv[1])) {
        std::cerr << "Failed to read mesh: " << argv[1] << "\n";
        return 1;
    }

    IncrementalFlattening param(mesh);
    param.parameterize();

    std::cout << "QC error = " << param.computeQcError() << "\n";
    std::cout << "Cone vertices:";
    const auto& cones = param.coneVertices();
    const auto& cK = param.coneCurvatures();
    for (size_t i = 0; i < cones.size(); ++i) {
        std::cout << " " << cones[i];
        if (static_cast<int>(i) < cK.size()) {
            std::cout << "(K=" << cK(static_cast<int>(i)) << ")";
        }
    }
    std::cout << "\n";
    std::cout << "Seam edges: " << param.seamEdgeCount()
              << ", cut vertices: " << param.cutVertexCount() << "\n";

    const std::string out = "incremental_flattening_uv.obj";
    if (mesh.write(out)) {
        std::cout << "Wrote " << out << "\n";
    }
    return 0;
}
