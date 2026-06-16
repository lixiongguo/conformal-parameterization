#include <emscripten/bind.h>
#include "Mesh.h"
#include "ARAP.h"
#include "Tutte.h"
#include <Eigen/Dense>
#include <vector>

using namespace emscripten;

// ============================================================
//  Mesh helpers (self-contained for this module)
// ============================================================
std::vector<double> _arap_vertexPositions(const Mesh& mesh) {
    std::vector<double> r;
    r.reserve(mesh.vertices.size() * 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        r.push_back(v->position.x());
        r.push_back(v->position.y());
        r.push_back(v->position.z());
    }
    return r;
}

std::vector<double> _arap_vertexUvs(const Mesh& mesh) {
    std::vector<double> r;
    r.reserve(mesh.vertices.size() * 2);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        r.push_back(v->uv.x());
        r.push_back(v->uv.y());
    }
    return r;
}

std::vector<int> _arap_faceIndices(const Mesh& mesh) {
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

EMSCRIPTEN_BINDINGS(ARAPModule) {

    // Base Mesh (self-contained)
    class_<Mesh>("Mesh")
        .constructor<>()
        .function("read", &Mesh::read)
        .function("write", &Mesh::write)
        .function("getVertexPositions", &_arap_vertexPositions)
        .function("getVertexUvs", &_arap_vertexUvs)
        .function("getFaceIndices", &_arap_faceIndices);

    // ARAP (As-Rigid-As-Possible)
    class_<ARAP, base<Parameterization>>("ARAP")
        .constructor<Mesh&>()
        .constructor<Mesh&, int>()
        .function("parameterize", &ARAP::parameterize)
        .function("computeQcError", &ARAP::computeQcError);

    // Tutte embedding
    class_<Tutte, base<Parameterization>>("Tutte")
        .constructor<Mesh&>()
        .constructor<Mesh&, TutteBoundary, TutteWeight>()
        .function("parameterize", &Tutte::parameterize)
        .function("computeQcError", &Tutte::computeQcError);

    // Enum helpers for Tutte parameters (passed as int from JS)
    // TutteBoundary::CIRCLE = 0, TutteBoundary::SQUARE = 1
    // TutteWeight::COTAN = 0, TutteWeight::UNIFORM = 1
    enum_<TutteBoundary>("TutteBoundary")
        .value("CIRCLE", TutteBoundary::CIRCLE)
        .value("SQUARE", TutteBoundary::SQUARE);

    enum_<TutteWeight>("TutteWeight")
        .value("COTAN", TutteWeight::COTAN)
        .value("UNIFORM", TutteWeight::UNIFORM);
}
