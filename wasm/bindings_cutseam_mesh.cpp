#include <emscripten/bind.h>
#include "Mesh.h"
#include "CutSeamMesh.h"
#include "TreeCotreeBasis.h"
#include <vector>

using namespace emscripten;

// --- Mesh helpers (self-contained for this module) ---
std::vector<double> getVertexPositions_csm(const Mesh& mesh) {
    std::vector<double> result;
    result.reserve(mesh.vertices.size() * 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        result.push_back(v->position.x());
        result.push_back(v->position.y());
        result.push_back(v->position.z());
    }
    return result;
}

std::vector<double> getVertexUvs_csm(const Mesh& mesh) {
    std::vector<double> result;
    result.reserve(mesh.vertices.size() * 2);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        result.push_back(v->uv.x());
        result.push_back(v->uv.y());
    }
    return result;
}

std::vector<int> getFaceIndices_csm(const Mesh& mesh) {
    std::vector<int> result;
    result.reserve(mesh.faces.size() * 3);
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        HalfEdgeCIter he = f->he;
        for (int k = 0; k < 3; k++) {
            result.push_back(he->vertex->index);
            he = he->next;
        }
    }
    return result;
}

// --- CutSeamMesh binding helpers ---
std::vector<int> originalVertexForVertexVec(const CutSeamMesh& cm) {
    const auto& v = cm.originalVertexForVertex();
    return std::vector<int>(v.begin(), v.end());
}

std::vector<int> originalFaceForFaceVec(const CutSeamMesh& cm) {
    const auto& v = cm.originalFaceForFace();
    return std::vector<int>(v.begin(), v.end());
}

std::vector<int> cutBoundaryEdgesVec(const CutSeamMesh& cm) {
    const auto& v = cm.cutBoundaryEdges();
    std::vector<int> result;
    result.reserve(v.size());
    for (bool b : v) result.push_back(b ? 1 : 0);
    return result;
}

// seam vertex pairs as flat int array [a0,b0, a1,b1, ...]
std::vector<int> seamVertexPairsFlat(const CutSeamMesh& cm) {
    const auto& pairs = cm.seamVertexPairs();
    std::vector<int> result;
    result.reserve(pairs.size() * 2);
    for (const auto& p : pairs) {
        result.push_back(p.first);
        result.push_back(p.second);
    }
    return result;
}

// --- TreeCotreeBasis ---
// buildClosedMeshBasis: return flat array of edgeIndex, sign pairs
std::vector<int> treeCotreeBasisVec(const Mesh& mesh) {
    using namespace topology;
    auto cycles = TreeCotreeBasis::buildClosedMeshBasis(mesh);
    std::vector<int> result;
    for (const auto& cyc : cycles) {
        result.push_back(-1); // cycle separator
        for (const auto& se : cyc) {
            result.push_back(se.first);  // edgeIndex
            result.push_back(se.second); // sign
        }
    }
    return result;
}

bool treeCotreeBuildCutSeam(const Mesh& mesh, CutSeamMesh& cutMesh) {
    return topology::TreeCotreeBasis::buildCutSeamMesh(mesh, cutMesh);
}

EMSCRIPTEN_BINDINGS(CutSeamMeshModule) {
    // Base Mesh (self-contained)
    class_<Mesh>("Mesh")
        .constructor<>()
        .function("read", &Mesh::read)
        .function("write", &Mesh::write)
        .function("delaunayize", &Mesh::delaunayize)
        .function("meanEdgeLength", &Mesh::meanEdgeLength)
        .function("getVertexPositions", &getVertexPositions_csm)
        .function("getVertexUvs", &getVertexUvs_csm)
        .function("getFaceIndices", &getFaceIndices_csm);

    // CutSeamMesh extends Mesh
    class_<CutSeamMesh, base<Mesh>>("CutSeamMesh")
        .constructor<>()
        .constructor<const Mesh&>()
        .function("buildFromSeamEdges", &CutSeamMesh::buildFromSeamEdges)
        .function("buildFromSeamVertexPairs", &CutSeamMesh::buildFromSeamVertexPairs)
        .function("originalVertexIndex", &CutSeamMesh::originalVertexIndex)
        .function("originalFaceIndex", &CutSeamMesh::originalFaceIndex)
        .function("isCutBoundaryEdge",
                  select_overload<bool(int)const>(&CutSeamMesh::isCutBoundaryEdge))
        .function("originalVertexForVertex", &originalVertexForVertexVec)
        .function("originalFaceForFace", &originalFaceForFaceVec)
        .function("cutBoundaryEdges", &cutBoundaryEdgesVec)
        .function("seamVertexPairs", &seamVertexPairsFlat);

    // TreeCotreeBasis static methods
    function("treeCotreeBasis", &treeCotreeBasisVec);
    function("treeCotreeBuildCutSeam", &treeCotreeBuildCutSeam);
}
