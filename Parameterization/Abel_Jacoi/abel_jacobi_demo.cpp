/**
 * Standalone demo: Abel–Jacobi map and divisor check on a closed mesh.
 *
 * Build (from repo cpp/):
 *   g++ -std=c++17 -O2 -I conformal-parameterization -I deps/eigen-3.4.0 \
 *       conformal-parameterization/Abel_Jacoi/abel_jacobi_demo.cpp \
 *       conformal-parameterization/Abel_Jacoi/AbelJacobi.cpp \
 *       conformal-parameterization/Abel_Jacoi/AbelJacobiParameterization.cpp \
 *       conformal-parameterization/Mesh.cpp MeshIO.cpp Parameterization.cpp \
 *       conformal-parameterization/Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp \
 *       conformal-parameterization/Solver.cpp QcError.cpp \
 *       conformal-parameterization/uv_unwrap_simple/Lscm.cpp \
 *       -o abel_jacobi_demo
 */
#include "AbelJacobiParameterization.h"
#include "MeshIO.h"
#include <iostream>

static bool isClosed(const Mesh& m)
{
    if (!m.boundaries.empty()) return false;
    for (EdgeCIter e = m.edges.begin(); e != m.edges.end(); ++e) {
        if (e->isBoundary()) return false;
    }
    return true;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [base_vertex]\n";
        return 1;
    }

    Mesh mesh;
    if (!mesh.read(argv[1])) {
        std::cerr << "Failed to read mesh: " << argv[1] << "\n";
        return 1;
    }

    if (!isClosed(mesh)) {
        std::cerr << "Mesh must be closed (no boundary) for Abel–Jacobi.\n";
        return 1;
    }

    AbelJacobiParameterization param(mesh);
    if (argc >= 3) {
        param.setBaseVertex(std::atoi(argv[2]));
    }

    AbelJacobi aj(mesh);
    if (!aj.build()) {
        std::cerr << "Abel–Jacobi build failed.\n";
        return 1;
    }

    std::cout << "Genus g = " << aj.genus() << "\n";
    std::cout << "Lattice generators (g x 4g):\n" << aj.latticeGenerators() << "\n";

    const int base = (argc >= 3) ? std::atoi(argv[2]) : 0;
    const int nV = static_cast<int>(mesh.vertices.size());
    for (int vi = 0; vi < std::min(nV, 5); ++vi) {
        const Eigen::MatrixXd mu = aj.abelJacobiMap(base, vi);
        std::cout << "mu(" << vi << ") = " << mu.transpose() << "\n";
    }

    param.parameterize();
    std::cout << "QC error = " << param.computeQcError() << "\n";

    const std::string out = "abel_jacobi_uv.obj";
    if (mesh.write(out)) {
        std::cout << "Wrote " << out << "\n";
    }
    return 0;
}
