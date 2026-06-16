#include <emscripten/bind.h>
#include "Mesh.h"
#include "Lscm.h"
#include "Scp.h"
#include <Eigen/Dense>
#include <vector>

using namespace emscripten;

// ============================================================
//  Mesh helpers (self-contained for this module)
// ============================================================
std::vector<double> _lscm_vertexPositions(const Mesh& mesh) {
    std::vector<double> r;
    r.reserve(mesh.vertices.size() * 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        r.push_back(v->position.x());
        r.push_back(v->position.y());
        r.push_back(v->position.z());
    }
    return r;
}

std::vector<double> _lscm_vertexUvs(const Mesh& mesh) {
    std::vector<double> r;
    r.reserve(mesh.vertices.size() * 2);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        r.push_back(v->uv.x());
        r.push_back(v->uv.y());
    }
    return r;
}

std::vector<int> _lscm_faceIndices(const Mesh& mesh) {
    std::vector<int> r;
    r.reserve(mesh.faces.size() * 3);
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        HalfEdgeCIter he = f->he;
        for (int k = 0; k < 3; k++) {
            r.push_back(he->vertex->index);
            he = he->next;
        }
    }
    return r;
}

EMSCRIPTEN_BINDINGS(LSCMModule) {

    // Base Mesh (self-contained)
    class_<Mesh>("Mesh")
        .constructor<>()
        .function("read", &Mesh::read)
        .function("write", &Mesh::write)
        .function("getVertexPositions", &_lscm_vertexPositions)
        .function("getVertexUvs", &_lscm_vertexUvs)
        .function("getFaceIndices", &_lscm_faceIndices);

    // LSCM (Least Squares Conformal Maps)
    class_<Lscm, base<Parameterization>>("Lscm")
        .constructor<Mesh&>()
        .function("parameterize", &Lscm::parameterize)
        .function("computeQcError", &Lscm::computeQcError);

    // SCP (Spectral Conformal Parameterization)
    class_<Scp, base<Parameterization>>("Scp")
        .constructor<Mesh&>()
        .function("parameterize", &Scp::parameterize)
        .function("computeQcError", &Scp::computeQcError);
}
