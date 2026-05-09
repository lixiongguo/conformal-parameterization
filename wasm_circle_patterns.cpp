/**
 * wasm_circle_patterns.cpp - WebAssembly 入口，Circle Patterns + 锥奇异点
 *
 * Exports:
 *   solve_cp(pos, posLen, faces, faceLen, optScheme,
 *            coneIdx, coneIdxLen, coneAngles, coneAnglesLen) → int
 *   get_cp_uv_result() → double*
 *   get_cp_uv_result_size() → int
 *   get_cp_last_time_ms() → double
 *   cp_dispose() → void
 */

#include <emscripten.h>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <cstdio>
#include <limits>

#include "Mesh.h"
#include "MeshIO.h"
#include "CirclePatternsWasm.h"

static Mesh* g_mesh = nullptr;
static double g_lastTimeMs = 0.0;
static std::vector<double> g_uv_result;

static bool loadMesh(const double* pos, int posLen, const int* faces, int faceLen) {
    if (g_mesh) { delete g_mesh; g_mesh = nullptr; }
    if (posLen < 9 || faceLen < 3) return false;
    g_mesh = new Mesh();
    size_t nV = posLen/3, nF = faceLen/3;
    std::stringstream ss;
    ss << "# CPSolver\n# " << nV << " v, " << nF << " f\n\n";
    for (size_t i=0;i<nV;i++) ss<<"v "<<pos[i*3]<<" "<<pos[i*3+1]<<" "<<pos[i*3+2]<<"\n";
    ss<<"\n";
    for (size_t i=0;i<nF;i++) ss<<"f "<<faces[i*3]+1<<" "<<faces[i*3+1]+1<<" "<<faces[i*3+2]+1<<"\n";
    std::string s=ss.str();
    FILE* fp=fopen("/tmp/input_cp.obj","wb");
    if(!fp)return false;
    fwrite(s.c_str(),1,s.size(),fp);fclose(fp);
    if(!g_mesh->read("/tmp/input_cp.obj")){delete g_mesh;g_mesh=nullptr;return false;}
    return true;
}

extern "C" {

EMSCRIPTEN_KEEPALIVE
int solve_cp(double* posPtr, int posLen, int* facePtr, int faceLen, int optScheme,
             int* coneIdx, int coneIdxLen, double* coneAngles, int coneAnglesLen) {
    if (!posPtr || !facePtr || posLen < 9 || faceLen < 3) return -1;
    if (!loadMesh(posPtr, posLen, facePtr, faceLen)) return -1;

    auto t0 = std::chrono::high_resolution_clock::now();
    if (optScheme != 0 && optScheme != 1 && optScheme != 3) optScheme = 1;

    CirclePatternsWasm cp(*g_mesh, optScheme);
    
    // Set cone singulars if provided
    if (coneIdx && coneIdxLen > 0 && coneAngles && coneAnglesLen > 0) {
        std::vector<int> ci(coneIdx, coneIdx + coneIdxLen);
        std::vector<double> ca(coneAngles, coneAngles + coneAnglesLen);
        cp.setConeSingulars(ci, ca);
    }
    
    cp.parameterize();

    auto t1 = std::chrono::high_resolution_clock::now();
    g_lastTimeMs = std::chrono::duration<double,std::milli>(t1-t0).count();

    g_uv_result.clear(); g_uv_result.reserve(g_mesh->vertices.size()*2);
    double min_u=std::numeric_limits<double>::infinity(),max_u=-min_u,min_v=min_u,max_v=-min_u;
    for(VertexCIter v=g_mesh->vertices.begin();v!=g_mesh->vertices.end();++v){
        double u=v->uv.x(),w=v->uv.y();
        if(u<min_u)min_u=u;if(u>max_u)max_u=u;
        if(w<min_v)min_v=w;if(w>max_v)max_v=w;
    }
    double ru=max_u-min_u,rv=max_v-min_v;
    if(ru<1e-10)ru=1.0;if(rv<1e-10)rv=1.0;
    for(VertexCIter v=g_mesh->vertices.begin();v!=g_mesh->vertices.end();++v){
        g_uv_result.push_back((v->uv.x()-min_u)/ru);
        g_uv_result.push_back((v->uv.y()-min_v)/rv);
    }
    return 0;
}

EMSCRIPTEN_KEEPALIVE double* get_cp_uv_result(){return g_uv_result.empty()?nullptr:g_uv_result.data();}
EMSCRIPTEN_KEEPALIVE int get_cp_uv_result_size(){return(int)g_uv_result.size();}
EMSCRIPTEN_KEEPALIVE double get_cp_last_time_ms(){return g_lastTimeMs;}
EMSCRIPTEN_KEEPALIVE void cp_dispose(){if(g_mesh){delete g_mesh;g_mesh=nullptr;}g_uv_result.clear();g_uv_result.shrink_to_fit();}

} // extern "C"
