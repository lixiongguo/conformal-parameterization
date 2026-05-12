/**
 * wasm_quadcover.cpp - QuadCover WASM 入口
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
#include "QuadCover.h"

static Mesh* g_mesh = nullptr;
static double g_lastTimeMs = 0.0;
static std::vector<double> g_uv_result;

static bool loadMesh(const double* pos, int posLen, const int* faces, int faceLen) {
    if(g_mesh){delete g_mesh;g_mesh=nullptr;}
    if(posLen<9||faceLen<3)return false;
    g_mesh=new Mesh();
    size_t nV=posLen/3,nF=faceLen/3;
    std::stringstream ss;
    ss<<"# QuadCover\n# "<<nV<<" v "<<nF<<" f\n\n";
    for(size_t i=0;i<nV;i++)ss<<"v "<<pos[i*3]<<" "<<pos[i*3+1]<<" "<<pos[i*3+2]<<"\n";
    ss<<"\n";
    for(size_t i=0;i<nF;i++)ss<<"f "<<faces[i*3]+1<<" "<<faces[i*3+1]+1<<" "<<faces[i*3+2]+1<<"\n";
    std::string s=ss.str();
    FILE* fp=fopen("/tmp/input_qc.obj","wb");
    if(!fp)return false;
    fwrite(s.c_str(),1,s.size(),fp);fclose(fp);
    if(!g_mesh->read("/tmp/input_qc.obj")){delete g_mesh;g_mesh=nullptr;return false;}
    return true;
}

extern "C" {
EMSCRIPTEN_KEEPALIVE
int solve_qc_with_field(double* posPtr, int posLen, int* facePtr, int faceLen, double* dirPtr, int dirLen) {
    if(!posPtr||!facePtr||posLen<9||faceLen<3)return -1;
    if(!loadMesh(posPtr,posLen,facePtr,faceLen))return -1;
    auto t0=std::chrono::high_resolution_clock::now();
    QuadCover qc(*g_mesh);
    if(dirPtr&&dirLen==faceLen){
        Eigen::MatrixXd dirs(faceLen/3,3);
        for(int i=0;i<dirLen/3;i++){
            dirs(i,0)=dirPtr[i*3];
            dirs(i,1)=dirPtr[i*3+1];
            dirs(i,2)=dirPtr[i*3+2];
        }
        qc.setFaceDirections(dirs);
    }
    qc.parameterize();
    auto t1=std::chrono::high_resolution_clock::now();
    g_lastTimeMs=std::chrono::duration<double,std::milli>(t1-t0).count();
    g_uv_result.clear();g_uv_result.reserve(g_mesh->vertices.size()*2);
    double mu=std::numeric_limits<double>::infinity(),Mu=-mu,mv=mu,Mv=-mu;
    for(VertexCIter v=g_mesh->vertices.begin();v!=g_mesh->vertices.end();++v){
        double u=v->uv.x(),w=v->uv.y();
        if(u<mu)mu=u;if(u>Mu)Mu=u;
        if(w<mv)mv=w;if(w>Mv)Mv=w;
    }
    double ru=Mu-mu,rv=Mv-mv;
    if(ru<1e-10)ru=1.0;if(rv<1e-10)rv=1.0;
    for(VertexCIter v=g_mesh->vertices.begin();v!=g_mesh->vertices.end();++v){
        g_uv_result.push_back((v->uv.x()-mu)/ru);
        g_uv_result.push_back((v->uv.y()-mv)/rv);
    }
    return 0;
}
EMSCRIPTEN_KEEPALIVE
int solve_qc(double* posPtr, int posLen, int* facePtr, int faceLen) {
    return solve_qc_with_field(posPtr,posLen,facePtr,faceLen,nullptr,0);
}
EMSCRIPTEN_KEEPALIVE double* get_qc_uv_result(){return g_uv_result.empty()?nullptr:g_uv_result.data();}
EMSCRIPTEN_KEEPALIVE int get_qc_uv_result_size(){return(int)g_uv_result.size();}
EMSCRIPTEN_KEEPALIVE double get_qc_last_time_ms(){return g_lastTimeMs;}
EMSCRIPTEN_KEEPALIVE void qc_dispose(){if(g_mesh){delete g_mesh;g_mesh=nullptr;}g_uv_result.clear();g_uv_result.shrink_to_fit();}
}
