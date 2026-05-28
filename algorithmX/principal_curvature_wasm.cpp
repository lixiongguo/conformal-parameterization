/**
 * Principal Curvature WASM Module
 * 
 * Uses libigl's principal_curvature to compute per-face principal
 * curvature directions and values. Compiled to WebAssembly via Emscripten.
 * 
 * API (exposed via EMSCRIPTEN_BINDINGS / cwrap):
 *   compute_principal_curvature(V_ptr, V_rows, F_ptr, F_rows)
 *   -> { PD1, PD2, PV1, PV2, numFaces }
 * 
 *   - V: vertices as flat float array [x0,y0,z0, x1,y1,z1, ...]
 *   - F: faces as flat int array [v0,v1,v2, ...]
 *   Returns pointer to result struct with all arrays
 */

#include <emscripten.h>
#include <emscripten/bind.h>

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>

#include <igl/principal_curvature.h>
#include <igl/average_onto_faces.h>

// ============================================================
// Result structure stored in WASM memory
// Layout: header + data arrays
// Header (all int32):
//   [0]  = magic number (0xDEADBEEF)
//   [1]  = numVertices
//   [2]  = numFaces
//   [3]  = offset to PD1_face (float array, numFaces*3)
//   [4]  = offset to PD2_face (float array, numFaces*3)
//   [5]  = offset to PV1_face (float array, numFaces*1)
//   [6]  = offset to PV2_face (float array, numFaces*1)
//   [7]  = offset to PD1_vert (float array, numVert*3)
//   [8]  = offset to PD2_vert (float array, numVert*3)
//   [9]  = offset to PV1_vert (float array, numVert*1)
//   [10] = offset to PV2_vert (float array, numVert*1)
//   Total header size: 44 bytes (11 ints)
// ============================================================

struct CurvatureResult {
    int numVertices;
    int numFaces;
    
    // Per-face principal directions and curvatures
    Eigen::MatrixXd PD1_face; // numFaces x 3
    Eigen::MatrixXd PD2_face; // numFaces x 3
    Eigen::VectorXd PV1_face; // numFaces x 1
    Eigen::VectorXd PV2_face; // numFaces x 1
    
    // Per-vertex principal directions and curvatures
    Eigen::MatrixXd PD1_vert; // numV x 3
    Eigen::MatrixXd PD2_vert; // numV x 3
    Eigen::VectorXd PV1_vert; // numV x 1
    Eigen::VectorXd PV2_vert; // numV x 1
    
    // Serialized buffer
    std::vector<float> buffer;
};

extern "C" {

/**
 * Compute principal curvature for a triangle mesh.
 * 
 * @param V_ptr  Pointer to vertex positions (float array, V_rows*3 elements, row-major)
 * @param V_rows  Number of vertices
 * @param F_ptr  Pointer to face indices (int array, F_rows*3 elements, row-major)
 * @param F_rows  Number of faces
 * @return  Pointer to serialized result buffer in WASM memory.
 *          Caller must call free_result(ptr) when done.
 */
EMSCRIPTEN_KEEPALIVE
int* compute_principal_curvature(
    float* V_ptr,
    int V_rows,
    int* F_ptr,
    int F_rows
) {
    static CurvatureResult result;
    
    try {
        // ---- Build libigl matrices from raw pointers ----
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>> 
            V_map(V_ptr, V_rows, 3);
        
        Eigen::MatrixXd V = V_map.cast<double>();
        
        Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>> 
            F_map(F_ptr, F_rows, 3);
        
        Eigen::MatrixXi F = F_map;
        
        result.numVertices = V_rows;
        result.numFaces = F_rows;
        
        // ---- Compute principal curvature at vertices ----
        Eigen::MatrixXd PD1_vert, PD2_vert;
        Eigen::VectorXd PV1_vert, PV2_vert;
        
        igl::principal_curvature(
            V, F,
            PD1_vert, PD2_vert, PV1_vert, PV2_vert,
            5   // radius (average edge length multiplier)
            , true  // use K-ring neighborhood
        );
        
        result.PD1_vert = PD1_vert;
        result.PD2_vert = PD2_vert;
        result.PV1_vert = PV1_vert;
        result.PV2_vert = PV2_vert;
        
        // ---- Average onto faces ----
        Eigen::MatrixXd PD1_face, PD2_face;
        Eigen::VectorXd PV1_face, PV2_face;
        
        igl::average_onto_faces(F, PD1_vert, PD1_face);
        igl::average_onto_faces(F, PD2_vert, PD2_face);
        igl::average_onto_faces(F, PV1_vert, PV1_face);
        igl::average_onto_faces(F, PV2_vert, PV2_face);
        
        result.PD1_face = PD1_face;
        result.PD2_face = PD2_face;
        result.PV1_face = PV1_face;
        result.PV2_face = PV2_face;
        
        // ---- Serialize to flat buffer ----
        // Header: 11 ints
        int hdrSize = 11;
        int faceDataSize = F_rows * 3;       // PD1_face, PD2_face
        int vertDataSize = V_rows * 3;        // PD1_vert, PD2_vert
        int faceScalarSize = F_rows;           // PV1_face, PV2_face
        int vertScalarSize = V_rows;            // PV1_vert, PV2_vert
        
        int totalFloats = hdrSize 
                        + faceDataSize * 2      // PD1_face, PD2_face
                        + faceScalarSize * 2     // PV1_face, PV2_face
                        + vertDataSize * 2       // PD1_vert, PD2_vert
                        + vertScalarSize * 2;    // PV1_vert, PV2_vert
        
        result.buffer.resize(totalFloats);
        float* buf = result.buffer.data();
        
        // Fill header (as floats, but will be read as ints on JS side)
        // Use bit-cast through memcpy to be safe
        int* ibuf = reinterpret_cast<int*>(buf);
        int offset = 0;
        
        ibuf[offset++] = 0xDEADBEEF;   // magic
        ibuf[offset++] = V_rows;         // numVertices
        ibuf[offset++] = F_rows;         // numFaces
        ibuf[offset++] = hdrSize;        // PD1_face offset
        ibuf[offset++] = hdrSize + faceDataSize;   // PD2_face offset
        ibuf[offset++] = hdrSize + faceDataSize * 2; // PV1_face offset
        ibuf[offset++] = hdrSize + faceDataSize * 2 + faceScalarSize; // PV2_face offset
        ibuf[offset++] = hdrSize + faceDataSize * 2 + faceScalarSize * 2; // PD1_vert offset
        ibuf[offset++] = hdrSize + faceDataSize * 2 + faceScalarSize * 2 + vertDataSize; // PD2_vert offset
        ibuf[offset++] = hdrSize + faceDataSize * 2 + faceScalarSize * 2 + vertDataSize * 2; // PV1_vert offset
        ibuf[offset++] = hdrSize + faceDataSize * 2 + faceScalarSize * 2 + vertDataSize * 2 + vertScalarSize; // PV2_vert offset
        
        // Now fill float data starting after header
        int foff = hdrSize;
        
        // PD1_face (numFaces x 3)
        for (int i = 0; i < F_rows; i++)
            for (int j = 0; j < 3; j++)
                buf[foff++] = static_cast<float>(PD1_face(i, j));
        
        // PD2_face (numFaces x 3)
        for (int i = 0; i < F_rows; i++)
            for (int j = 0; j < 3; j++)
                buf[foff++] = static_cast<float>(PD2_face(i, j));
        
        // PV1_face (numFaces x 1)
        for (int i = 0; i < F_rows; i++)
            buf[foff++] = static_cast<float>(PV1_face(i));
        
        // PV2_face (numFaces x 1)
        for (int i = 0; i < F_rows; i++)
            buf[foff++] = static_cast<float>(PV2_face(i));
        
        // PD1_vert (numVert x 3)
        for (int i = 0; i < V_rows; i++)
            for (int j = 0; j < 3; j++)
                buf[foff++] = static_cast<float>(PD1_vert(i, j));
        
        // PD2_vert (numVert x 3)
        for (int i = 0; i < V_rows; i++)
            for (int j = 0; j < 3; j++)
                buf[foff++] = static_cast<float>(PD2_vert(i, j));
        
        // PV1_vert (numVert x 1)
        for (int i = 0; i < V_rows; i++)
            buf[foff++] = static_cast<float>(PV1_vert(i));
        
        // PV2_vert (numVert x 1)
        for (int i = 0; i < V_rows; i++)
            buf[foff++] = static_cast<float>(PV2_vert(i));
        
        // Return pointer to the buffer (it's in static storage, valid until next call)
        return reinterpret_cast<int*>(result.buffer.data());
        
    } catch (const std::exception& e) {
        std::cerr << "principal_curvature error: " << e.what() << std::endl;
        static int err_buf[2] = {0, -1};  // magic=0 means error
        return err_buf;
    }
}

/**
 * Get the total size of the result buffer in bytes.
 * Call this first to allocate memory if using a different approach.
 */
EMSCRIPTEN_KEEPALIVE
int get_result_buffer_size(
    float* V_ptr,
    int V_rows,
    int* F_ptr,
    int F_rows
) {
    int hdrSize = 11;
    int faceDataSize = F_rows * 3;
    int vertDataSize = V_rows * 3;
    int faceScalarSize = F_rows;
    int vertScalarSize = V_rows;
    
    int totalFloats = hdrSize 
                    + faceDataSize * 2 
                    + faceScalarSize * 2
                    + vertDataSize * 2 
                    + vertScalarSize * 2;
    
    return totalFloats * sizeof(float);  // total bytes
}

}  // extern "C"
