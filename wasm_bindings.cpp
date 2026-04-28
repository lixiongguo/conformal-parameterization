/**
 * WebAssembly bindings for conformal-parameterization library
 * Exposes LSCM and CETM algorithms to JavaScript
 */

#include <emscripten.h>
#include <emscripten/bind.h>
#include <vector>
#include <string>
#include <memory>
#include <cstring>

#include "Mesh.h"
#include "Lscm.h"
#include "Cetm.h"

using namespace emscripten;

// Global mesh instance managed by WASM module
static std::unique_ptr<Mesh> g_mesh = nullptr;
static std::unique_ptr<Parameterization> g_param = nullptr;

class ConformalParam {
public:
    ConformalParam() {}

    // Initialize mesh from vertex positions and face indices (from JavaScript/Three.js)
    bool initMesh(const std::vector<float>& positions, const std::vector<int>& faces) {
        if (positions.size() == 0 || faces.size() == 0) return false;
        if (positions.size() % 3 != 0 || faces.size() % 3 != 0) return false;

        size_t nVertices = positions.size() / 3;
        size_t nFaces = faces.size() / 3;

        g_mesh.reset(new Mesh());

        // Build OBJ content in memory and parse it
        std::string objContent;
        objContent.reserve(nVertices * 30 + nFaces * 40);

        // Write vertices
        for (size_t i = 0; i < nVertices; i++) {
            objContent += "v ";
            objContent += std::to_string(positions[i * 3]);
            objContent += " ";
            objContent += std::to_string(positions[i * 3 + 1]);
            objContent += " ";
            objContent += std::to_string(positions[i * 3 + 2]);
            objContent += "\n";
        }

        // Write faces (1-indexed as per OBJ format)
        for (size_t i = 0; i < nFaces; i++) {
            objContent += "f ";
            objContent += std::to_string(faces[i * 3] + 1);
            objContent += " ";
            objContent += std::to_string(faces[i * 3 + 1] + 1);
            objContent += " ";
            objContent += std::to_string(faces[i * 3 + 2] + 1);
            objContent += "\n";
        }

        // Write to a virtual file that Emscripten can read
        const char* filename = "/tmp/input_mesh.obj";
        FILE* f = fopen(filename, "w");
        if (!f) return false;
        fputs(objContent.c_str(), f);
        fclose(f);

        bool success = g_mesh->read(std::string(filename));
        return success;
    }

    // Initialize mesh from OBJ string (alternative method)
    bool initMeshFromObj(const std::string& objStr) {
        g_mesh.reset(new Mesh());

        const char* filename = "/tmp/input_mesh.obj";
        FILE* f = fopen(filename, "w");
        if (!f) return false;
        fputs(objStr.c_str(), f);
        fclose(f);

        return g_mesh->read(std::string(filename));
    }

    // Run LSCM parameterization
    bool computeLSCM() {
        if (!g_mesh) return false;

        g_param.reset(new Lscm(*g_mesh));
        g_param->parameterize();

        return true;
    }

    // Run CETM parameterization with specified optimization scheme:
    // 0 = gradient descent, 1 = newton, 2 = lbfgs
    bool computeCETM(int optScheme) {
        if (!g_mesh) return false;

        try {
            g_param.reset(new Cetm(*g_mesh, optScheme));
            g_param->parameterize();
        } catch (...) {
            return false;
        }

        return true;
    }

    // Get UV coordinates as flat array [u0, v0, u1, v1, ...]
    std::vector<float> getUVs() {
        std::vector<float> uvs;
        if (!g_mesh) return uvs;

        uvs.reserve(g_mesh->vertices.size() * 2);
        for (VertexCIter v = g_mesh->vertices.begin(); v != g_mesh->vertices.end(); v++) {
            uvs.push_back(static_cast<float>(v->uv.x()));
            uvs.push_back(static_cast<float>(v->uv.y()));
        }

        return uvs;
    }

    // Get number of vertices
    int getNumVertices() {
        if (!g_mesh) return 0;
        return static_cast<int>(g_mesh->vertices.size());
    }

    // Get number of faces
    int getNumFaces() {
        if (!g_mesh) return 0;
        return static_cast<int>(g_mesh->faces.size());
    }

    // Get quasi-conformal error after parameterization
    double getQCError() {
        if (!g_param) return -1.0;
        return g_param->computeQcError();
    }

    // Free resources
    void cleanup() {
        g_param.reset();
        g_mesh.reset();
    }
};

// Register C++ classes and functions for JavaScript access
EMSCRIPTEN_BINDINGS(conformal_parameterization_module) {
    class_<ConformalParam>("ConformalParam")
        .constructor<>()
        .function("initMesh", &ConformalParam::initMesh)
        .function("initMeshFromObj", &ConformalParam::initMeshFromObj)
        .function("computeLSCM", &ConformalParam::computeLSCM)
        .function("computeCETM", &ConformalParam::computeCETM)
        .function("getUVs", &ConformalParam::getUVs)
        .function("getNumVertices", &ConformalParam::getNumVertices)
        .function("getNumFaces", &ConformalParam::getNumFaces)
        .function("getQCError", &ConformalParam::getQCError)
        .function("cleanup", &ConformalParam::cleanup)
        ;

    // Register vector types for interop
    register_vector<float>("VectorFloat");
    register_vector<int>("VectorInt");
}
