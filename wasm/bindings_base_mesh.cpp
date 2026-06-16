#include <emscripten/bind.h>
#include "Mesh.h"
#include "GaussianCurvature.h"
#include "QcError.h"
#include <Eigen/Dense>

using namespace emscripten;

// Mesh vertex/edge data helpers for JS interop
std::vector<double> getVertexPositions(const Mesh& mesh) {
    std::vector<double> result;
    result.reserve(mesh.vertices.size() * 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        result.push_back(v->position.x());
        result.push_back(v->position.y());
        result.push_back(v->position.z());
    }
    return result;
}

std::vector<double> getVertexUvs(const Mesh& mesh) {
    std::vector<double> result;
    result.reserve(mesh.vertices.size() * 2);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        result.push_back(v->uv.x());
        result.push_back(v->uv.y());
    }
    return result;
}

std::vector<int> getFaceIndices(const Mesh& mesh) {
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

int getVertexCount(const Mesh& mesh) { return (int)mesh.vertices.size(); }
int getFaceCount(const Mesh& mesh)   { return (int)mesh.faces.size(); }
int getEdgeCount(const Mesh& mesh)   { return (int)mesh.edges.size(); }

// Gaussian curvature → vector of doubles
std::vector<double> gaussianCurvatureAngleDeficitVec(const Mesh& mesh) {
    Eigen::VectorXd k = geometry::gaussianCurvatureAngleDeficit(mesh);
    return std::vector<double>(k.data(), k.data() + k.size());
}

std::vector<double> gaussianCurvaturePerAreaVec(const Mesh& mesh, double eps) {
    Eigen::VectorXd k = geometry::gaussianCurvaturePerArea(mesh, eps);
    return std::vector<double>(k.data(), k.data() + k.size());
}

EMSCRIPTEN_BINDINGS(BaseMeshModule) {
    class_<Mesh>("Mesh")
        .constructor<>()
        .function("read", &Mesh::read)
        .function("write", &Mesh::write)
        .function("delaunayize", &Mesh::delaunayize)
        .function("meanEdgeLength", &Mesh::meanEdgeLength)
        .function("getVertexPositions", &getVertexPositions)
        .function("getVertexUvs", &getVertexUvs)
        .function("getFaceIndices", &getFaceIndices)
        .function("getVertexCount", &getVertexCount)
        .function("getFaceCount", &getFaceCount)
        .function("getEdgeCount", &getEdgeCount);

    function("gaussianCurvatureAngleDeficit", &gaussianCurvatureAngleDeficitVec);
    function("gaussianCurvaturePerArea",       &gaussianCurvaturePerAreaVec);

    class_<QuasiConformalError>("QuasiConformalError")
        .class_function("compute", emscripten::optional_override(
            [](const std::vector<double>& p1, const std::vector<double>& p2,
               const std::vector<double>& p3,
               const std::vector<double>& q1, const std::vector<double>& q2,
               const std::vector<double>& q3) {
                std::vector<Eigen::Vector3d> p = {
                    Eigen::Vector3d(p1[0], p1[1], p1[2]),
                    Eigen::Vector3d(p2[0], p2[1], p2[2]),
                    Eigen::Vector3d(p3[0], p3[1], p3[2])};
                std::vector<Eigen::Vector3d> q = {
                    Eigen::Vector3d(q1[0], q1[1], q1[2]),
                    Eigen::Vector3d(q2[0], q2[1], q2[2]),
                    Eigen::Vector3d(q3[0], q3[1], q3[2])};
                return QuasiConformalError::compute(p, q);
            }))
        .class_function("color", emscripten::optional_override(
            [](double qc) {
                Eigen::Vector3d c = QuasiConformalError::color(qc);
                return std::vector<double>{c.x(), c.y(), c.z()};
            }));
}
