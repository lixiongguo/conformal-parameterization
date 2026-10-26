#include <string>
#include <vector>

#include <emscripten.h>

#include "NaiveQuadMerge.h"

namespace {
NaiveQuadMerge::Result g_result;
std::vector<std::string> g_stage_names;
std::vector<double> g_stage_ms;

void publishStages(const std::vector<StageTime>& stages) {
    g_stage_names.clear();
    g_stage_ms.clear();
    for (const StageTime& stage : stages) {
        g_stage_names.push_back(stage.name);
        g_stage_ms.push_back(stage.ms);
    }
}
}

extern "C" {

EMSCRIPTEN_KEEPALIVE
int merge_tris_to_quads(double* pos, int nVerts, int* faces, int nFaces) {
    g_result = NaiveQuadMerge::merge(pos, nVerts, faces, nFaces);
    publishStages(g_result.stages);
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int get_quad_count() { return static_cast<int>(g_result.quads.size() / 4); }

EMSCRIPTEN_KEEPALIVE
int* get_quads() { return g_result.quads.empty() ? nullptr : g_result.quads.data(); }

EMSCRIPTEN_KEEPALIVE
int get_leftover_count() { return static_cast<int>(g_result.leftovers.size() / 3); }

EMSCRIPTEN_KEEPALIVE
int* get_leftovers() { return g_result.leftovers.empty() ? nullptr : g_result.leftovers.data(); }

EMSCRIPTEN_KEEPALIVE
int get_stage_count() { return static_cast<int>(g_stage_names.size()); }

EMSCRIPTEN_KEEPALIVE
const char* get_stage_name(int index) {
    if (index < 0 || index >= static_cast<int>(g_stage_names.size())) return "";
    return g_stage_names[static_cast<size_t>(index)].c_str();
}

EMSCRIPTEN_KEEPALIVE
double get_stage_ms(int index) {
    if (index < 0 || index >= static_cast<int>(g_stage_ms.size())) return 0;
    return g_stage_ms[static_cast<size_t>(index)];
}

EMSCRIPTEN_KEEPALIVE
int naive_stage_count() { return 2; }

EMSCRIPTEN_KEEPALIVE
const char* naive_stage_name(int index) {
    if (index == 0) return "建立邻接";
    if (index == 1) return "贪心配对";
    return "";
}

}
