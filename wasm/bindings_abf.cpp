#include <emscripten/bind.h>
#include "Mesh.h"
#include "AbfPlusPlus.h"
#include "LinAbf.h"
#include <Eigen/Dense>
#include <vector>

using namespace emscripten;

// ============================================================
//  Mesh helpers (self-contained for this module)
// ============================================================
std::vector<double> _abf_vertexPositions(const Mesh& mesh) {
    std::vector<double> r;
    r.reserve(mesh.vertices.size() * 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        r.push_back(v->position.x());
        r.push_back(v->position.y());
        r.push_back(v->position.z());
    }
    return r;
}

std::vector<double> _abf_vertexUvs(const Mesh& mesh) {
    std::vector<double> r;
    r.reserve(mesh.vertices.size() * 2);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        r.push_back(v->uv.x());
        r.push_back(v->uv.y());
    }
    return r;
}

std::vector<int> _abf_faceIndices(const Mesh& mesh) {
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

EMSCRIPTEN_BINDINGS(ABFModule) {

    // Base Mesh (self-contained)
    class_<Mesh>("Mesh")
        .constructor<>()
        .function("read", &Mesh::read)
        .function("write", &Mesh::write)
        .function("getVertexPositions", &_abf_vertexPositions)
        .function("getVertexUvs", &_abf_vertexUvs)
        .function("getFaceIndices", &_abf_faceIndices);

    // ABF++ (Angle Based Flattening++)
    class_<AbfPlusPlus, base<Parameterization>>("AbfPlusPlus")
        .constructor<Mesh&>()
        .constructor<Mesh&, int>()
        .function("parameterize", &AbfPlusPlus::parameterize)
        .function("computeQcError", &AbfPlusPlus::computeQcError);

    // Linear ABF
    class_<LinAbf, base<Parameterization>>("LinAbf")
        .constructor<Mesh&>()
        .constructor<Mesh&, int>()
        .function("parameterize", &LinAbf::parameterize)
        .function("computeQcError", &LinAbf::computeQcError);
}
