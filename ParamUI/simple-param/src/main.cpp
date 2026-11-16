// Flat gallery of one simple parameterization on every topological disk.
// The algorithm is chosen from a dropdown; every disk is shown together.

#include "ARAP.h"
#include "AbfPlusPlus.h"
#include "BoundedDistortionMapping.hpp"
#include "TreeCotreeBasis.h"
#include "Lscm.h"
#include "LinAbf.h"
#include "Mesh.h"
#include "MeshIO.h"
#include "QcError.h"
#include "Scp.h"
#include "Tutte.h"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <OpenGL/gl3.h>

#include "imgui.h"
#include "imgui_chrome.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <condition_variable>
#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <set>
#include <tuple>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

enum class Algo { Lscm, Scp, Arap, TutteCircle, Abfpp, LinAbfAlgo, BdLscm };

struct DiskModel {
    std::string label;
    std::filesystem::path path;
    int faces = 0;
};

struct Cell {
    std::vector<float> uv;
    std::vector<float> pos;
    std::vector<float> nrm;
    std::vector<float> seamPos;
    std::vector<unsigned> flipTris;
    float cx = 0;
    float cy = 0;
    float cz = 0;
    float radius = 1;
    std::vector<unsigned> tris;
    int flipped = 0;
    double seconds = 0;
    bool ok = false;
    std::string note;
    std::string cutTag;
    bool puncture = false;
    std::vector<std::pair<int, int>> seams;
    std::vector<float> seamUv;
    std::vector<float> faceColor;
    std::vector<float> areaColor;
    float qcAvg = 1.f;
    float qcMax = 1.f;
    float areaAvg = 1.f;
    float areaMax = 1.f;
    GLuint distortVao = 0;
    GLuint distortVbo = 0;
    GLsizei distortCount = 0;
    GLuint surfVao = 0;
    GLuint surfVbo = 0;
    GLsizei surfCount = 0;
    GLuint qcVao = 0;
    GLuint qcVbo = 0;
    GLsizei qcCount = 0;
    GLuint areaVao = 0;
    GLuint areaVbo = 0;
    GLsizei areaCount = 0;
    GLuint flipVao = 0;
    GLuint flipVbo = 0;
    GLsizei flipCount = 0;
    GLuint flipUvVao = 0;
    GLuint flipUvVbo = 0;
    GLsizei flipUvCount = 0;
    int contentRev = 0;
    int uploadedRev = -1;
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint triEbo = 0;
    GLuint lineEbo = 0;
    GLsizei triCount = 0;
    GLsizei lineCount = 0;
    bool uploaded = false;
};

namespace {

const Algo kAlgos[] = {
    Algo::Lscm, Algo::Scp, Algo::Arap, Algo::TutteCircle, Algo::Abfpp, Algo::LinAbfAlgo, Algo::BdLscm,
};
const int kAlgoCount = 7;
const float kCell = 260.f;
const float kGap = 16.f;
const float kMargin = 16.f;
const float kHud = 52.f;
// 自适应适配时单格不缩到这个尺寸以下（点），宁可让画布略微溢出、靠拖动查看。
const float kMinCellPoints = 300.f;

std::filesystem::path executableDir() {
#ifdef __APPLE__
    char buf[4096];
    uint32_t size = sizeof(buf);
    if (_NSGetExecutablePath(buf, &size) == 0) return std::filesystem::weakly_canonical(buf).parent_path();
#endif
    return {};
}

void addDir(std::vector<std::filesystem::path>& dirs, const std::filesystem::path& path) {
    std::error_code error;
    if (!std::filesystem::is_directory(path, error)) return;
    const std::filesystem::path canonical = std::filesystem::weakly_canonical(path, error);
    if (error) return;
    if (std::find(dirs.begin(), dirs.end(), canonical) == dirs.end()) dirs.push_back(canonical);
}

void addNearby(std::vector<std::filesystem::path>& dirs, const std::filesystem::path& start) {
    const char* suffixes[] = {"assets/Models", "assets/Model1", "input", "assets/Models2"};
    std::error_code error;
    std::filesystem::path dir = std::filesystem::weakly_canonical(start, error);
    if (error || dir.empty()) return;
    for (int depth = 0; depth < 6 && !dir.empty(); ++depth) {
        for (const char* suffix : suffixes) addDir(dirs, dir / suffix);
        error.clear();
        for (const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator(dir, error)) {
            if (error || !entry.is_directory()) continue;
            for (const char* suffix : suffixes) addDir(dirs, entry.path() / suffix);
        }
        const std::filesystem::path parent = dir.parent_path();
        if (parent.empty() || parent == dir) break;
        dir = parent;
    }
}

bool isTopologicalDisk(const Mesh& mesh) {
    const int boundaryLoops = static_cast<int>(mesh.boundaries.size());
    const int euler = static_cast<int>(mesh.vertices.size()) - static_cast<int>(mesh.edges.size()) +
                      static_cast<int>(mesh.faces.size());
    return boundaryLoops == 1 && euler == 2;
}

int realFaceCount(const Mesh& mesh) {
    int count = 0;
    for (const Face& face : mesh.faces) {
        if (!face.isBoundary()) ++count;
    }
    return count;
}

std::vector<DiskModel> listDisks() {
    std::vector<std::filesystem::path> dirs;
    const std::filesystem::path exe = executableDir();
    const char* relatives[] = {
        "../../../../assets/Models", "../../../../assets/Model1", "assets/Models", "assets/Model1", "input",
    };
    for (const char* relative : relatives) addDir(dirs, relative);
    if (!exe.empty()) {
        addDir(dirs, exe / "../../../../../../assets/Models");
        addDir(dirs, exe / "../../../../../../assets/Model1");
    }
    std::error_code cwdError;
    const std::filesystem::path cwd = std::filesystem::current_path(cwdError);
    if (!cwdError) addNearby(dirs, cwd);
    if (!exe.empty()) addNearby(dirs, exe);

    std::map<std::string, DiskModel> found;
    for (const std::filesystem::path& dir : dirs) {
        for (const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator(dir)) {
            if (!entry.is_regular_file() || entry.path().extension() != ".obj") continue;
            Mesh mesh;
            if (!mesh.read(entry.path().string())) {
                std::cout << "跳过（读入失败）  " << entry.path().filename().string() << std::endl;
                continue;
            }
            if (!isTopologicalDisk(mesh)) {
                std::cout << "非圆盘，展开前切循环边  " << entry.path().filename().string()
                          << "  边界 " << mesh.boundaries.size() << std::endl;
            }
            DiskModel model;
            model.label = entry.path().filename().string();
            model.path = entry.path();
            model.faces = realFaceCount(mesh);
            found[model.label] = std::move(model);
        }
    }
    std::vector<DiskModel> models;
    for (auto& item : found) models.push_back(std::move(item.second));
    std::sort(models.begin(), models.end(), [](const DiskModel& a, const DiskModel& b) {
        return a.label < b.label;
    });
    return models;
}

bool uvHasSpread(const Mesh& mesh) {
    if (mesh.vertices.empty()) return false;
    double minU = mesh.vertices[0].uv.x(), maxU = minU;
    double minV = mesh.vertices[0].uv.y(), maxV = minV;
    for (const Vertex& vertex : mesh.vertices) {
        minU = std::min(minU, vertex.uv.x());
        maxU = std::max(maxU, vertex.uv.x());
        minV = std::min(minV, vertex.uv.y());
        maxV = std::max(maxV, vertex.uv.y());
    }
    return (maxU - minU) >= 1e-12 || (maxV - minV) >= 1e-12;
}

std::pair<int, int> farthestBoundaryPins(const Mesh& mesh) {
    std::vector<int> boundary;
    for (const Vertex& vertex : mesh.vertices) {
        if (vertex.isBoundary()) boundary.push_back(vertex.index);
    }
    std::pair<int, int> pins = {-1, -1};
    double best = -1.0;
    for (size_t i = 0; i < boundary.size(); ++i) {
        for (size_t j = i + 1; j < boundary.size(); ++j) {
            const double distance = (mesh.vertices[boundary[i]].position - mesh.vertices[boundary[j]].position).norm();
            if (distance > best) {
                best = distance;
                pins = {boundary[i], boundary[j]};
            }
        }
    }
    return pins;
}

const char* algoName(Algo algo) {
    switch (algo) {
        case Algo::Lscm: return "LSCM";
        case Algo::Scp: return "SCP";
        case Algo::Arap: return "ARAP";
        case Algo::TutteCircle: return "Tutte";
        case Algo::Abfpp: return "ABF++";
        case Algo::LinAbfAlgo: return "linABF";
        case Algo::BdLscm: return "BD-LSCM";
    }
    return "";
}

struct FormulaText {
    const char* title;
    const char* body;
};

FormulaText formulaText(Algo algo) {
    switch (algo) {
        case Algo::Lscm:
            return {"LSCM · 最小二乘共形映射",
                    "E(u) = ∫ |∂z̄u|² dA = ½ ∫ ‖∇u‖² dA − A(u)\n"
                    "∂z̄u = ½(∂xu + i·∂yu) 是 Cauchy–Riemann 算子，A(u) 是 UV 有向面积。能量对 u 是二次的，解一次线性方程就得到全局最优。"};
        case Algo::Scp:
            return {"SCP · 谱共形参数化",
                    "Q u = λ M u，λmin = min (uᵀ Q u) / (uᵀ M u)\n"
                    "Q 是共形能量矩阵（cot-Laplacian 减去 2 倍有向面积），M 是顶点面积质量矩阵。取最小非零广义特征值对应的特征向量。"};
        case Algo::Arap:
            return {"ARAP · 尽可能刚性",
                    "E(u) = min Σf Σ(i,j)∈f  wij ‖(ui−uj) − Rf (xi−xj)‖²，Rf ∈ SO(2)\n"
                    "local：逐面用 SVD 求最优旋转 Rf。global：固定 Rf 后解泊松方程。用 Tutte 做初值，两步交替迭代。"};
        case Algo::TutteCircle:
            return {"Tutte · 圆边界调和映射",
                    "内部顶点 Σj∈N(i) wij (uj − ui) = 0，边界固定到圆周，wij 为余切权重\n"
                    "这是把边界钉在圆上的调和映射，内部由稀疏线性方程唯一确定，也用作 ARAP 的初值。"};
        case Algo::Abfpp:
            return {"ABF++ · 基于角度的展开",
                    "min Σf Σk (αk − βk)² / βk²\n"
                    "约束：三角形内角和为 π；内部顶点角和为 2π，边界为 θv；内部顶点满足正弦定律。α 是待求角，β 是三维目标角。角度求出后再线性重建 UV。"};
        case Algo::LinAbfAlgo:
            return {"linABF · 线性化角度展开",
                    "在目标角附近把 ABF 的角度约束线性化，直接解稀疏线性系统\n"
                    "得到的角度再重建 UV。它保留 ABF 的角度误差，但不再做非线性迭代。"};
        case Algo::BdLscm:
            return {"BD-LSCM · 有界失真",
                    "min ½ ∫ ‖∇u‖² dA − A(u)，约束 K = σ1/σ2 ≤ σ 且 J(u) > 0\n"
                    "K 是共形失真上界，J 是雅可比，要求保向、不翻面。以 LSCM 为初值，用增广拉格朗日把失真压进上界。"};
    }
    return {"", ""};
}

std::string shortName(const std::string& filename) {
    std::string stem = filename;
    const auto dot = stem.rfind('.');
    if (dot != std::string::npos) stem.erase(dot);
    const auto mark = stem.rfind("_F");
    if (mark != std::string::npos) {
        bool digits = mark + 2 < stem.size();
        for (size_t i = mark + 2; i < stem.size(); ++i) digits = digits && std::isdigit(static_cast<unsigned char>(stem[i]));
        if (digits) stem.erase(mark);
    }
    return stem;
}

std::string runAlgorithm(Mesh& mesh, Algo algo) {
    if (algo == Algo::Lscm) {
        Lscm lscm(mesh);
        lscm.parameterize();
        if (!lscm.succeeded()) return "LSCM 失败";
        return "";
    }
    if (algo == Algo::Scp) {
        Scp(mesh).parameterize();
        return "";
    }
    if (algo == Algo::Arap) {
        ARAP(mesh, 5, TutteBoundary::CIRCLE, TutteWeight::COTAN).parameterize();
        return "";
    }
    if (algo == Algo::TutteCircle) {
        Tutte(mesh, TutteBoundary::CIRCLE, TutteWeight::COTAN).parameterize();
        return "";
    }
    if (algo == Algo::Abfpp) {
        AbfPlusPlus(mesh).parameterize();
        return "";
    }
    if (algo == Algo::LinAbfAlgo) {
        LinAbf(mesh).parameterize();
        return "";
    }

    Lscm lscm(mesh);
    lscm.parameterize();
    if (!lscm.succeeded() || !uvHasSpread(mesh)) return "LSCM 初值失败";
    const std::pair<int, int> pins = farthestBoundaryPins(mesh);
    if (pins.first < 0) return "没有边界锚点";
    const int nV = static_cast<int>(mesh.vertices.size());
    Eigen::MatrixXd positions(nV, 3);
    Eigen::MatrixXd initial(nV, 2);
    for (const Vertex& vertex : mesh.vertices) {
        positions.row(vertex.index) = vertex.position.transpose();
        initial.row(vertex.index) = vertex.uv.transpose();
    }
    std::vector<std::array<int, 3>> triangles;
    for (const Face& face : mesh.faces) {
        if (face.isBoundary()) continue;
        triangles.push_back({face.he->vertex->index, face.he->next->vertex->index, face.he->next->next->vertex->index});
    }
    Eigen::MatrixXi faces(static_cast<int>(triangles.size()), 3);
    for (int face = 0; face < faces.rows(); ++face) {
        faces(face, 0) = triangles[face][0];
        faces(face, 1) = triangles[face][1];
        faces(face, 2) = triangles[face][2];
    }
    bounded_distortion::Options options;
    options.distortion_bound = 1.5;
    options.outer_iterations = 4;
    options.inner_iterations = 8;
    options.lscm_weight = 1.0;
    options.distortion_penalty = 5000.0;
    options.positivity_penalty = 5000.0;
    options.initial_step = 1e-2;
    const std::vector<bounded_distortion::Anchor> anchors = {
        {pins.first, mesh.vertices[pins.first].uv},
        {pins.second, mesh.vertices[pins.second].uv},
    };
    const bounded_distortion::SolveResult result =
        bounded_distortion::solveBoundedDistortionLscm(positions, faces, initial, anchors, options);
    if (result.uv.rows() != nV) return "UV 行数不对";
    for (int vertex = 0; vertex < nV; ++vertex) {
        mesh.vertices[vertex].uv.x() = result.uv(vertex, 0);
        mesh.vertices[vertex].uv.y() = result.uv(vertex, 1);
    }
    return "";
}

void fillCell(const Mesh& mesh, Cell& cell) {
    const int nV = static_cast<int>(mesh.vertices.size());
    cell.uv.assign(nV * 2, 0.f);
    cell.pos.assign(static_cast<size_t>(nV) * 3, 0.f);
    cell.nrm.assign(static_cast<size_t>(nV) * 3, 0.f);
    float minX = 1e30f, minY = 1e30f, minZ = 1e30f;
    float maxX = -1e30f, maxY = -1e30f, maxZ = -1e30f;
    for (const Vertex& vertex : mesh.vertices) {
        if (vertex.index < 0 || vertex.index >= nV) continue;
        cell.uv[vertex.index * 2] = static_cast<float>(vertex.uv.x());
        cell.uv[vertex.index * 2 + 1] = static_cast<float>(vertex.uv.y());
        const float x = static_cast<float>(vertex.position.x());
        const float y = static_cast<float>(vertex.position.y());
        const float z = static_cast<float>(vertex.position.z());
        cell.pos[static_cast<size_t>(vertex.index) * 3] = x;
        cell.pos[static_cast<size_t>(vertex.index) * 3 + 1] = y;
        cell.pos[static_cast<size_t>(vertex.index) * 3 + 2] = z;
        minX = std::min(minX, x);
        minY = std::min(minY, y);
        minZ = std::min(minZ, z);
        maxX = std::max(maxX, x);
        maxY = std::max(maxY, y);
        maxZ = std::max(maxZ, z);
    }
    cell.cx = (minX + maxX) * 0.5f;
    cell.cy = (minY + maxY) * 0.5f;
    cell.cz = (minZ + maxZ) * 0.5f;
    const float dx = maxX - minX, dy = maxY - minY, dz = maxZ - minZ;
    cell.radius = std::max(0.05f, 0.5f * std::sqrt(dx * dx + dy * dy + dz * dz));
    for (const Face& face : mesh.faces) {
        if (face.isBoundary()) continue;
        cell.tris.push_back(static_cast<unsigned>(face.he->vertex->index));
        cell.tris.push_back(static_cast<unsigned>(face.he->next->vertex->index));
        cell.tris.push_back(static_cast<unsigned>(face.he->next->next->vertex->index));
        const Eigen::Vector3d e1 = face.he->next->vertex->position - face.he->vertex->position;
        const Eigen::Vector3d e2 = face.he->next->next->vertex->position - face.he->vertex->position;
        const Eigen::Vector3d normal = e1.cross(e2);
        const int ids[3] = {face.he->vertex->index, face.he->next->vertex->index, face.he->next->next->vertex->index};
        for (int id : ids) {
            if (id < 0 || id >= nV) continue;
            cell.nrm[static_cast<size_t>(id) * 3] += static_cast<float>(normal.x());
            cell.nrm[static_cast<size_t>(id) * 3 + 1] += static_cast<float>(normal.y());
            cell.nrm[static_cast<size_t>(id) * 3 + 2] += static_cast<float>(normal.z());
        }
    }
    for (int vertex = 0; vertex < nV; ++vertex) {
        const float x = cell.nrm[static_cast<size_t>(vertex) * 3];
        const float y = cell.nrm[static_cast<size_t>(vertex) * 3 + 1];
        const float z = cell.nrm[static_cast<size_t>(vertex) * 3 + 2];
        const float length = std::sqrt(x * x + y * y + z * z);
        if (length < 1e-8f) continue;
        cell.nrm[static_cast<size_t>(vertex) * 3] = x / length;
        cell.nrm[static_cast<size_t>(vertex) * 3 + 1] = y / length;
        cell.nrm[static_cast<size_t>(vertex) * 3 + 2] = z / length;
    }
    cell.flipTris.clear();
    const double eps = 1e-12;
    int flipped = 0;
    for (const Face& face : mesh.faces) {
        if (face.isBoundary()) continue;
        const int idx[3] = {face.he->vertex->index, face.he->next->vertex->index, face.he->next->next->vertex->index};
        const Eigen::Vector2d& a = face.he->vertex->uv;
        const Eigen::Vector2d& b = face.he->next->vertex->uv;
        const Eigen::Vector2d& c = face.he->next->next->vertex->uv;
        if ((b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y()) >= -eps) continue;
        ++flipped;
        for (int corner : idx) cell.flipTris.push_back(static_cast<unsigned>(corner));
    }
    cell.flipped = flipped;
    cell.faceColor.clear();
    cell.areaColor.clear();
    double qcSum = 0.0;
    double qcWeight = 0.0;
    double qcPeak = 1.0;
    std::vector<double> uvAreas;
    std::vector<double> meshAreas;
    for (const Face& face : mesh.faces) {
        if (face.isBoundary()) continue;
        std::vector<Eigen::Vector3d> position(3);
        std::vector<Eigen::Vector3d> uv(3);
        HalfEdgeCIter he = face.he;
        for (int corner = 0; corner < 3; ++corner) {
            position[static_cast<size_t>(corner)] = he->vertex->position;
            uv[static_cast<size_t>(corner)] = Eigen::Vector3d(he->vertex->uv.x(), he->vertex->uv.y(), 0.0);
            he = he->next;
        }
        double distortion = QuasiConformalError::compute(position, uv);
        if (!std::isfinite(distortion) || distortion < 1.0) distortion = 1.0;
        const Eigen::Vector3d color = QuasiConformalError::color(distortion);
        cell.faceColor.push_back(static_cast<float>(color.x()));
        cell.faceColor.push_back(static_cast<float>(color.y()));
        cell.faceColor.push_back(static_cast<float>(color.z()));
        const double meshArea = face.area();
        const double signedUv = 0.5 * ((uv[1].x() - uv[0].x()) * (uv[2].y() - uv[0].y()) -
                                        (uv[2].x() - uv[0].x()) * (uv[1].y() - uv[0].y()));
        uvAreas.push_back(std::abs(signedUv));
        meshAreas.push_back(std::isfinite(meshArea) && meshArea > 0.0 ? meshArea : 0.0);
        if (meshAreas.back() > 0.0) {
            qcSum += meshAreas.back() * distortion;
            qcWeight += meshAreas.back();
        }
        qcPeak = std::max(qcPeak, distortion);
    }
    cell.qcAvg = qcWeight > 0.0 ? static_cast<float>(qcSum / qcWeight) : 1.f;
    cell.qcMax = qcWeight > 0.0 ? static_cast<float>(qcPeak) : 1.f;
    double sumUv = 0.0;
    double sumMesh = 0.0;
    for (size_t face = 0; face < uvAreas.size(); ++face) {
        sumUv += uvAreas[face];
        sumMesh += meshAreas[face];
    }
    const double scale = sumUv > 1e-18 && sumMesh > 1e-18 ? sumMesh / sumUv : 1.0;
    double areaSum = 0.0;
    double areaWeight = 0.0;
    double areaPeak = 1.0;
    for (size_t face = 0; face < uvAreas.size(); ++face) {
        double stretch = 1.0;
        if (meshAreas[face] > 1e-18 && uvAreas[face] > 1e-18) {
            const double rho = uvAreas[face] * scale / meshAreas[face];
            stretch = std::max(rho, 1.0 / rho);
        } else if (meshAreas[face] > 1e-18) {
            stretch = 1.0e6;
        }
        if (!std::isfinite(stretch) || stretch < 1.0) stretch = 1.0;
        const Eigen::Vector3d color = QuasiConformalError::color(stretch);
        cell.areaColor.push_back(static_cast<float>(color.x()));
        cell.areaColor.push_back(static_cast<float>(color.y()));
        cell.areaColor.push_back(static_cast<float>(color.z()));
        if (meshAreas[face] > 0.0) {
            areaSum += meshAreas[face] * stretch;
            areaWeight += meshAreas[face];
        }
        areaPeak = std::max(areaPeak, stretch);
    }
    cell.areaAvg = areaWeight > 0.0 ? static_cast<float>(areaSum / areaWeight) : 1.f;
    cell.areaMax = areaWeight > 0.0 ? static_cast<float>(areaPeak) : 1.f;
}

std::pair<int, int> canonEdge(int a, int b) {
    return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
}

void appendSeamUv(const Edge& edge, Cell& cell) {
    cell.seamUv.push_back(static_cast<float>(edge.he->vertex->uv.x()));
    cell.seamUv.push_back(static_cast<float>(edge.he->vertex->uv.y()));
    cell.seamUv.push_back(static_cast<float>(edge.he->flip->vertex->uv.x()));
    cell.seamUv.push_back(static_cast<float>(edge.he->flip->vertex->uv.y()));
    cell.seamPos.push_back(static_cast<float>(edge.he->vertex->position.x()));
    cell.seamPos.push_back(static_cast<float>(edge.he->vertex->position.y()));
    cell.seamPos.push_back(static_cast<float>(edge.he->vertex->position.z()));
    cell.seamPos.push_back(static_cast<float>(edge.he->flip->vertex->position.x()));
    cell.seamPos.push_back(static_cast<float>(edge.he->flip->vertex->position.y()));
    cell.seamPos.push_back(static_cast<float>(edge.he->flip->vertex->position.z()));
}

void recordCutLines(const CutSeamMesh& cut, Cell& cell) {
    cell.seams = cut.seamVertexPairs();
    cell.puncture = false;
    cell.seamUv.clear();
    cell.seamPos.clear();
    cell.cutTag = "切缝 " + std::to_string(static_cast<int>(cell.seams.size()));
    for (const Edge& edge : cut.edges) {
        if (cut.isCutBoundaryEdge(edge.index)) appendSeamUv(edge, cell);
    }
}

void recordOpenedBoundary(const Mesh& mesh, Cell& cell) {
    cell.seamUv.clear();
    cell.seamPos.clear();
    for (const Edge& edge : mesh.edges) {
        if (edge.isBoundary()) appendSeamUv(edge, cell);
    }
}

std::vector<int> cotreeEdges(const Mesh& mesh, bool pruneLeaves) {
    const int nFaces = static_cast<int>(mesh.faces.size());
    const int nEdges = static_cast<int>(mesh.edges.size());
    const int nVerts = static_cast<int>(mesh.vertices.size());
    std::vector<std::vector<std::pair<int, int>>> adjacent(static_cast<size_t>(nFaces));
    std::vector<int> interior;
    for (const Edge& edge : mesh.edges) {
        if (edge.isBoundary()) continue;
        if (edge.he->face->isBoundary() || edge.he->flip->face->isBoundary()) continue;
        const int f0 = edge.he->face->index;
        const int f1 = edge.he->flip->face->index;
        adjacent[static_cast<size_t>(f0)].push_back({f1, edge.index});
        adjacent[static_cast<size_t>(f1)].push_back({f0, edge.index});
        interior.push_back(edge.index);
    }

    int root = -1;
    for (const Face& face : mesh.faces) {
        if (!face.isBoundary()) {
            root = face.index;
            break;
        }
    }
    std::vector<char> inTree(static_cast<size_t>(nEdges), 0);
    if (root >= 0) {
        std::vector<char> visited(static_cast<size_t>(nFaces), 0);
        std::queue<int> pending;
        pending.push(root);
        visited[static_cast<size_t>(root)] = 1;
        while (!pending.empty()) {
            const int face = pending.front();
            pending.pop();
            for (const std::pair<int, int>& neighbor : adjacent[static_cast<size_t>(face)]) {
                if (visited[static_cast<size_t>(neighbor.first)]) continue;
                visited[static_cast<size_t>(neighbor.first)] = 1;
                inTree[static_cast<size_t>(neighbor.second)] = 1;
                pending.push(neighbor.first);
            }
        }
    }

    std::vector<char> cut(static_cast<size_t>(nEdges), 0);
    std::vector<int> edges;
    for (int edge : interior) {
        if (inTree[static_cast<size_t>(edge)]) continue;
        cut[static_cast<size_t>(edge)] = 1;
        edges.push_back(edge);
    }
    if (!pruneLeaves) return edges;

    std::vector<int> degree(static_cast<size_t>(nVerts), 0);
    std::vector<char> onBoundary(static_cast<size_t>(nVerts), 0);
    for (const Vertex& vertex : mesh.vertices) {
        if (vertex.isBoundary()) onBoundary[static_cast<size_t>(vertex.index)] = 1;
    }
    auto endpoint = [&](int edge, int side) {
        const Edge& item = mesh.edges[static_cast<size_t>(edge)];
        return side == 0 ? item.he->vertex->index : item.he->flip->vertex->index;
    };
    for (int edge : edges) {
        degree[static_cast<size_t>(endpoint(edge, 0))]++;
        degree[static_cast<size_t>(endpoint(edge, 1))]++;
    }
    bool removed = true;
    while (removed) {
        removed = false;
        for (size_t i = 0; i < edges.size();) {
            const int edge = edges[i];
            const int a = endpoint(edge, 0);
            const int b = endpoint(edge, 1);
            const bool leafA = degree[static_cast<size_t>(a)] == 1 && !onBoundary[static_cast<size_t>(a)];
            const bool leafB = degree[static_cast<size_t>(b)] == 1 && !onBoundary[static_cast<size_t>(b)];
            if (!leafA && !leafB) {
                ++i;
                continue;
            }
            degree[static_cast<size_t>(a)]--;
            degree[static_cast<size_t>(b)]--;
            cut[static_cast<size_t>(edge)] = 0;
            edges.erase(edges.begin() + static_cast<std::ptrdiff_t>(i));
            removed = true;
        }
    }
    return edges;
}

bool punctureToDisk(const Mesh& source, Mesh& out, std::vector<std::pair<int, int>>& seam) {
    MeshData data;
    data.positions.resize(source.vertices.size());
    for (const Vertex& vertex : source.vertices) data.positions[static_cast<size_t>(vertex.index)] = vertex.position;
    bool skipped = false;
    for (const Face& face : source.faces) {
        if (face.isBoundary()) continue;
        if (!skipped) {
            skipped = true;
            int corners[3];
            int count = 0;
            HalfEdgeCIter he = face.he;
            do {
                if (count < 3) corners[count] = he->vertex->index;
                ++count;
                he = he->next;
            } while (he != face.he);
            if (count >= 3) {
                seam.push_back(canonEdge(corners[0], corners[1]));
                seam.push_back(canonEdge(corners[1], corners[2]));
                seam.push_back(canonEdge(corners[2], corners[0]));
            }
            continue;
        }
        std::vector<Index> corners;
        HalfEdgeCIter he = face.he;
        do {
            corners.emplace_back(he->vertex->index, -1, -1);
            he = he->next;
        } while (he != face.he);
        data.indices.push_back(std::move(corners));
    }
    return MeshIO::buildMesh(data, out) && isTopologicalDisk(out);
}

bool cutBoundaryToDisk(const Mesh& source, CutSeamMesh& cut, std::string& tag) {
    const std::vector<int> loops = cotreeEdges(source, true);
    if (loops.empty() || !cut.buildFromSeamEdges(source, loops) || !isTopologicalDisk(cut)) return false;
    tag = "切缝 " + std::to_string(static_cast<int>(loops.size()));
    return true;
}

Cell solveOne(const DiskModel& model, Algo algo, const std::vector<std::pair<int, int>>* customSeams) {
    Cell cell;
    const auto started = std::chrono::steady_clock::now();
    try {
        Mesh mesh;
        if (!mesh.read(model.path.string())) {
            cell.note = "无法读取";
            return cell;
        }
        CutSeamMesh homology;
        CutSeamMesh boundaryCut;
        CutSeamMesh customCut;
        Mesh punctured;
        Mesh* target = &mesh;
        if (customSeams) {
            cell.seams = *customSeams;
            if (customSeams->empty()) {
                if (!isTopologicalDisk(mesh)) {
                    cell.note = "切缝后不是圆盘";
                    return cell;
                }
            } else if (!customCut.buildFromSeamVertexPairs(mesh, *customSeams) || !isTopologicalDisk(customCut)) {
                cell.note = "切缝后不是圆盘";
                return cell;
            } else {
                target = &customCut;
            }
        } else if (!isTopologicalDisk(mesh)) {
            bool opened = false;
            if (mesh.boundaries.empty() &&
                topology::TreeCotreeBasis::buildCutSeamMesh(mesh, homology) &&
                isTopologicalDisk(homology)) {
                target = &homology;
                opened = true;
            } else if (!mesh.boundaries.empty() && cutBoundaryToDisk(mesh, boundaryCut, cell.cutTag)) {
                target = &boundaryCut;
                opened = true;
            } else if (mesh.boundaries.empty()) {
                std::vector<std::pair<int, int>> hole;
                if (punctureToDisk(mesh, punctured, hole)) {
                    target = &punctured;
                    cell.seams = std::move(hole);
                    cell.puncture = true;
                    cell.cutTag = "开边界";
                    opened = true;
                }
            }
            if (!opened) {
                cell.note = "切缝失败";
                return cell;
            }
        }
        cell.note = runAlgorithm(*target, algo);
        cell.ok = cell.note.empty() && uvHasSpread(*target);
        if (cell.note.empty() && !cell.ok) cell.note = "UV 退化";
        if (cell.ok) {
            fillCell(*target, cell);
            if (target == &homology || target == &boundaryCut || target == &customCut) {
                recordCutLines(static_cast<const CutSeamMesh&>(*target), cell);
            } else if (cell.puncture) {
                recordOpenedBoundary(*target, cell);
            }
        }
    } catch (const std::exception& error) {
        cell.ok = false;
        cell.note = error.what();
    }
    cell.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    return cell;
}

int gRows = 1;
std::vector<float> gYaw;
std::vector<float> gPitch;
int gDragCell = -1;
int gPressDetail = -1;

int galleryColumns(int count) {
    if (count <= 1) return 1;
    if (count <= 3) return count;
    if (count <= 8) return 4;
    return 5;
}

float gridWidth() {
    const int columns = galleryColumns(gRows);
    return kMargin * 2.f + columns * kCell + (columns - 1) * kGap;
}
float gridHeight(int count) {
    const int columns = galleryColumns(count);
    const int rows = (count + columns - 1) / columns;
    return kMargin * 2.f + rows * kCell + std::max(0, rows - 1) * kGap;
}
float cellX(int index) { return kMargin + (index % galleryColumns(gRows)) * (kCell + kGap); }
float cellY(int index) { return kMargin + (index / galleryColumns(gRows)) * (kCell + kGap); }

struct View {
    float originX = 0;
    float originY = 0;
    float scale = 1;
    bool drag = false;
    double lastX = 0;
    double lastY = 0;
    bool fit = true;
    int fbW = 1;
    int fbH = 1;
    float contentScale = 1;
};

View gView;

struct Batch {
    std::mutex mu;
    std::condition_variable cv;
    bool quit = false;
    bool busy = false;
    int algo = 0;
    int selected = 0;
    int next = 0;
    int count = 0;
    int jobs = 1;
    int focus = 0;
    int finished = 0;
    bool again = false;
};

Batch gBatch;

enum class QualityView { Checker, Conformal, Area };
QualityView gQuality = QualityView::Checker;
const std::vector<DiskModel>* gModels = nullptr;
std::vector<Cell>* gCells = nullptr;
std::atomic<unsigned char>* gReady = nullptr;

struct ModelCut {
    bool user = false;
    bool haveAuto = false;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::pair<int, int>> automatic;
};

std::mutex gCutMu;
std::vector<ModelCut> gCuts;
std::unique_ptr<std::atomic<int>[]> gEpoch;
std::unique_ptr<std::atomic<unsigned char>[]> gEditOwns;

struct EditJob {
    int model = -1;
    int algo = 0;
    int epoch = 0;
    bool custom = false;
    std::vector<std::pair<int, int>> edges;
};

std::mutex gEditMu;
std::condition_variable gEditCv;
bool gEditJobOn = false;
EditJob gEditJob;
std::thread gEditThread;

std::mutex gPendingMu;
bool gPendingOn = false;
int gPendingModel = -1;
int gPendingAlgo = -1;
int gPendingEpoch = -1;
Cell gPendingCell;

struct Vec3 {
    float x = 0, y = 0, z = 0;
};

struct Mat4 {
    float m[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
};

struct SourceMesh {
    int model = -1;
    std::vector<float> pos;
    std::vector<float> normal;
    std::vector<unsigned> tris;
    std::vector<std::vector<std::pair<int, int>>> edges;
    Vec3 center;
    float radius = 1;
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint nbo = 0;
    GLuint ebo = 0;
    GLsizei triCount = 0;
    bool uploaded = false;
};

struct Editor {
    bool open = false;
    int model = -1;
    bool pick = false;
    bool painting = false;
    int paintButton = -1;
    bool paintMoved = false;
    double paintX = 0;
    double paintY = 0;
    bool firstHit = false;
    bool firstWasSelected = false;
    std::pair<int, int> firstEdge = {-1, -1};
    std::set<std::pair<int, int>> seams;
    std::set<std::pair<int, int>> strokeBefore;
    std::vector<std::pair<int, int>> autoSeams;
    std::vector<std::pair<int, int>> lastCommitted;
    std::vector<std::vector<std::pair<int, int>>> history;
    bool orbit = false;
    float yaw = 0.7f;
    float pitch = 0.35f;
    float dist = 2.6f;
    std::string status;
    SourceMesh source;
};

Editor gEd;

struct Rect {
    float x = 0, y = 0, w = 0, h = 0;
    bool contains(float px, float py) const { return px >= x && py >= y && px < x + w && py < y + h; }
};

enum class SeamCmd { Back, Pick, Trunk, Clear };
void beginBatchLocked() {
    gBatch.algo = gBatch.selected;
    gBatch.next = 0;
    gBatch.finished = 0;
    gBatch.jobs = 1;
    gBatch.busy = true;
    gBatch.again = false;
    gBatch.cv.notify_all();
}

void finishBatchJobLocked() {
    ++gBatch.finished;
    if (gBatch.finished < gBatch.jobs) return;
    gBatch.busy = false;
    const int slot = gBatch.selected * gBatch.count + gBatch.focus;
    const bool have = gReady && gBatch.count > 0 && slot >= 0 &&
                      gReady[static_cast<size_t>(slot)].load(std::memory_order_acquire) != 0;
    if (!gBatch.quit && (gBatch.again || gBatch.selected != gBatch.algo) && !have) beginBatchLocked();
    else gBatch.again = false;
}

void selectAlgo(int index) {
    std::lock_guard<std::mutex> lock(gBatch.mu);
    if (index < 0 || index >= kAlgoCount) return;
    if (index == gBatch.selected && !gBatch.busy) return;
    gBatch.selected = index;
    const int slot = index * gBatch.count + gBatch.focus;
    const bool have = gReady && gBatch.count > 0 &&
                      gReady[static_cast<size_t>(slot)].load(std::memory_order_acquire) != 0;
    if (have) return;
    if (!gBatch.busy) beginBatchLocked();
    else gBatch.again = true;
}

Vec3 vecAdd(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 vecSub(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 vecScale(Vec3 a, float s) { return {a.x * s, a.y * s, a.z * s}; }
float vecDot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 vecCross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
Vec3 vecNorm(Vec3 a) {
    const float length = std::sqrt(vecDot(a, a));
    if (length < 1e-8f) return {0.f, 1.f, 0.f};
    return vecScale(a, 1.f / length);
}

Mat4 matMul(const Mat4& a, const Mat4& b) {
    Mat4 result{};
    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            float sum = 0.f;
            for (int k = 0; k < 4; ++k) sum += a.m[k * 4 + row] * b.m[col * 4 + k];
            result.m[col * 4 + row] = sum;
        }
    }
    return result;
}

Mat4 perspective(float fovy, float aspect, float nearPlane, float farPlane) {
    Mat4 result{};
    const float t = std::tan(fovy * 0.5f);
    result.m[0] = 1.f / (aspect * t);
    result.m[5] = 1.f / t;
    result.m[10] = (farPlane + nearPlane) / (nearPlane - farPlane);
    result.m[11] = -1.f;
    result.m[14] = (2.f * farPlane * nearPlane) / (nearPlane - farPlane);
    result.m[15] = 0.f;
    return result;
}

Mat4 lookAt(Vec3 eye, Vec3 target, Vec3 up) {
    const Vec3 forward = vecNorm(vecSub(target, eye));
    const Vec3 side = vecNorm(vecCross(forward, up));
    const Vec3 lifted = vecCross(side, forward);
    Mat4 result{};
    result.m[0] = side.x;
    result.m[1] = lifted.x;
    result.m[2] = -forward.x;
    result.m[4] = side.y;
    result.m[5] = lifted.y;
    result.m[6] = -forward.y;
    result.m[8] = side.z;
    result.m[9] = lifted.z;
    result.m[10] = -forward.z;
    result.m[12] = -vecDot(side, eye);
    result.m[13] = -vecDot(lifted, eye);
    result.m[14] = vecDot(forward, eye);
    result.m[15] = 1.f;
    return result;
}

void cameraBasis(Vec3& eye, Vec3& forward, Vec3& right, Vec3& up) {
    const float cp = std::cos(gEd.pitch);
    const float sp = std::sin(gEd.pitch);
    const float cy = std::cos(gEd.yaw);
    const float sy = std::sin(gEd.yaw);
    const Vec3 offset = vecScale({cp * sy, sp, cp * cy}, std::max(gEd.source.radius, 1e-4f) * gEd.dist);
    eye = vecAdd(gEd.source.center, offset);
    forward = vecNorm(vecSub(gEd.source.center, eye));
    right = vecNorm(vecCross(forward, {0.f, 1.f, 0.f}));
    up = vecCross(right, forward);
}

Mat4 meshMvp(float aspect) {
    Vec3 eye, forward, right, up;
    cameraBasis(eye, forward, right, up);
    const float radius = std::max(gEd.source.radius, 1e-4f);
    const Mat4 projection = perspective(0.9f, std::max(aspect, 1e-3f), radius * 0.02f, radius * gEd.dist * 12.f);
    return matMul(projection, lookAt(eye, gEd.source.center, {0.f, 1.f, 0.f}));
}

Mat4 previewMvp(Vec3 center, float radius, float yaw, float pitch, float dist, float aspect) {
    const float cp = std::cos(pitch);
    const float sp = std::sin(pitch);
    const float cy = std::cos(yaw);
    const float sy = std::sin(yaw);
    const float safeRadius = std::max(radius, 1e-4f);
    const Vec3 eye = vecAdd(center, vecScale({cp * sy, sp, cp * cy}, safeRadius * dist));
    const Mat4 projection = perspective(0.9f, std::max(aspect, 1e-3f), safeRadius * 0.02f, safeRadius * dist * 12.f);
    return matMul(projection, lookAt(eye, center, {0.f, 1.f, 0.f}));
}

float toolbarPoints() {
    const float font = 16.f;
    const float button = font + 12.f;
    const float pad = 12.f;
    const float gap = 8.f;
    return pad + button + gap + button + gap + font + pad;
}

Rect panel3d() {
    const float hud = toolbarPoints() * gView.contentScale;
    const float uvWidth = std::min(static_cast<float>(gView.fbW) * 0.42f, 680.f * gView.contentScale);
    return {0.f, hud, std::max(1.f, static_cast<float>(gView.fbW) - uvWidth),
            std::max(1.f, static_cast<float>(gView.fbH) - hud)};
}

Rect panelUv() {
    const Rect left = panel3d();
    return {left.x + left.w, left.y, std::max(1.f, static_cast<float>(gView.fbW) - left.w), left.h};
}

bool projectToPanel(const Mat4& mvp, const Rect& panel, Vec3 point, float& sx, float& sy) {
    const float x = mvp.m[0] * point.x + mvp.m[4] * point.y + mvp.m[8] * point.z + mvp.m[12];
    const float y = mvp.m[1] * point.x + mvp.m[5] * point.y + mvp.m[9] * point.z + mvp.m[13];
    const float w = mvp.m[3] * point.x + mvp.m[7] * point.y + mvp.m[11] * point.z + mvp.m[15];
    if (w <= 1e-5f) return false;
    const float ndcX = x / w;
    const float ndcY = y / w;
    sx = panel.x + (ndcX * 0.5f + 0.5f) * panel.w;
    sy = panel.y + (1.f - (ndcY * 0.5f + 0.5f)) * panel.h;
    return true;
}

bool rayTriangle(Vec3 origin, Vec3 direction, Vec3 v0, Vec3 v1, Vec3 v2, float& t) {
    const Vec3 edge1 = vecSub(v1, v0);
    const Vec3 edge2 = vecSub(v2, v0);
    const Vec3 p = vecCross(direction, edge2);
    const float det = vecDot(edge1, p);
    if (std::fabs(det) < 1e-7f) return false;
    const float inv = 1.f / det;
    const Vec3 s = vecSub(origin, v0);
    const float u = vecDot(s, p) * inv;
    if (u < 0.f || u > 1.f) return false;
    const Vec3 q = vecCross(s, edge1);
    const float v = vecDot(direction, q) * inv;
    if (v < 0.f || u + v > 1.f) return false;
    t = vecDot(edge2, q) * inv;
    return t > 1e-4f;
}

float pointSegmentDistance2(float px, float py, float ax, float ay, float bx, float by) {
    const float abx = bx - ax;
    const float aby = by - ay;
    const float apx = px - ax;
    const float apy = py - ay;
    const float len2 = abx * abx + aby * aby;
    float t = len2 > 1e-8f ? (apx * abx + apy * aby) / len2 : 0.f;
    t = std::max(0.f, std::min(1.f, t));
    const float dx = apx - t * abx;
    const float dy = apy - t * aby;
    return dx * dx + dy * dy;
}

bool rayFromPanel(float sx, float sy, const Rect& panel, Vec3& origin, Vec3& direction) {
    if (!panel.contains(sx, sy)) return false;
    const float ndcX = ((sx - panel.x) / panel.w) * 2.f - 1.f;
    const float ndcY = 1.f - ((sy - panel.y) / panel.h) * 2.f;
    Vec3 eye, forward, right, up;
    cameraBasis(eye, forward, right, up);
    const float t = std::tan(0.45f);
    const float aspect = panel.w / std::max(panel.h, 1.f);
    direction = vecNorm(vecAdd(forward, vecAdd(vecScale(right, ndcX * t * aspect), vecScale(up, ndcY * t))));
    origin = eye;
    return true;
}

Vec3 positionOf(int index) {
    return {gEd.source.pos[static_cast<size_t>(index) * 3], gEd.source.pos[static_cast<size_t>(index) * 3 + 1],
            gEd.source.pos[static_cast<size_t>(index) * 3 + 2]};
}

bool pickSeamEdge(float sx, float sy, const Rect& panel, const Mat4& mvp, std::pair<int, int>& edge) {
    Vec3 origin, direction;
    if (!rayFromPanel(sx, sy, panel, origin, direction)) return false;
    const auto& tris = gEd.source.tris;
    float bestT = 1e30f;
    int best = -1;
    for (size_t i = 0; i + 2 < tris.size(); i += 3) {
        float t = 0.f;
        if (!rayTriangle(origin, direction, positionOf(static_cast<int>(tris[i])), positionOf(static_cast<int>(tris[i + 1])),
                         positionOf(static_cast<int>(tris[i + 2])), t)) {
            continue;
        }
        if (t < bestT) {
            bestT = t;
            best = static_cast<int>(i);
        }
    }
    if (best < 0) return false;
    const int corners[3] = {static_cast<int>(tris[static_cast<size_t>(best)]),
                            static_cast<int>(tris[static_cast<size_t>(best) + 1]),
                            static_cast<int>(tris[static_cast<size_t>(best) + 2])};
    std::set<std::pair<int, int>> candidates;
    for (int corner : corners) {
        if (corner < 0 || corner >= static_cast<int>(gEd.source.edges.size())) continue;
        for (const auto& item : gEd.source.edges[static_cast<size_t>(corner)]) candidates.insert(item);
    }
    const float limit = 15.f * gView.contentScale;
    float bestDistance = limit * limit;
    bool found = false;
    for (const auto& candidate : candidates) {
        float ax = 0, ay = 0, bx = 0, by = 0;
        if (!projectToPanel(mvp, panel, positionOf(candidate.first), ax, ay)) continue;
        if (!projectToPanel(mvp, panel, positionOf(candidate.second), bx, by)) continue;
        const float distance = pointSegmentDistance2(sx, sy, ax, ay, bx, by);
        if (distance <= bestDistance) {
            bestDistance = distance;
            edge = candidate;
            found = true;
        }
    }
    return found;
}

std::vector<std::pair<int, int>> seamVector(const std::set<std::pair<int, int>>& seams) {
    return {seams.begin(), seams.end()};
}

void releaseCellGl(Cell& cell) {
    if (cell.vao) glDeleteVertexArrays(1, &cell.vao);
    if (cell.vbo) glDeleteBuffers(1, &cell.vbo);
    if (cell.triEbo) glDeleteBuffers(1, &cell.triEbo);
    if (cell.lineEbo) glDeleteBuffers(1, &cell.lineEbo);
    if (cell.distortVao) glDeleteVertexArrays(1, &cell.distortVao);
    if (cell.distortVbo) glDeleteBuffers(1, &cell.distortVbo);
    if (cell.surfVao) glDeleteVertexArrays(1, &cell.surfVao);
    if (cell.surfVbo) glDeleteBuffers(1, &cell.surfVbo);
    if (cell.qcVao) glDeleteVertexArrays(1, &cell.qcVao);
    if (cell.qcVbo) glDeleteBuffers(1, &cell.qcVbo);
    if (cell.areaVao) glDeleteVertexArrays(1, &cell.areaVao);
    if (cell.areaVbo) glDeleteBuffers(1, &cell.areaVbo);
    if (cell.flipVao) glDeleteVertexArrays(1, &cell.flipVao);
    if (cell.flipVbo) glDeleteBuffers(1, &cell.flipVbo);
    if (cell.flipUvVao) glDeleteVertexArrays(1, &cell.flipUvVao);
    if (cell.flipUvVbo) glDeleteBuffers(1, &cell.flipUvVbo);
    cell.vao = 0;
    cell.vbo = 0;
    cell.triEbo = 0;
    cell.lineEbo = 0;
    cell.distortVao = 0;
    cell.distortVbo = 0;
    cell.distortCount = 0;
    cell.surfVao = 0;
    cell.surfVbo = 0;
    cell.surfCount = 0;
    cell.qcVao = 0;
    cell.qcVbo = 0;
    cell.qcCount = 0;
    cell.areaVao = 0;
    cell.areaVbo = 0;
    cell.areaCount = 0;
    cell.flipVao = 0;
    cell.flipVbo = 0;
    cell.flipCount = 0;
    cell.flipUvVao = 0;
    cell.flipUvVbo = 0;
    cell.flipUvCount = 0;
    cell.uploaded = false;
}

void releaseSourceGl() {
    if (gEd.source.vao) glDeleteVertexArrays(1, &gEd.source.vao);
    if (gEd.source.vbo) glDeleteBuffers(1, &gEd.source.vbo);
    if (gEd.source.nbo) glDeleteBuffers(1, &gEd.source.nbo);
    if (gEd.source.ebo) glDeleteBuffers(1, &gEd.source.ebo);
    gEd.source.vao = 0;
    gEd.source.vbo = 0;
    gEd.source.nbo = 0;
    gEd.source.ebo = 0;
    gEd.source.uploaded = false;
}

bool loadSource(int model) {
    if (gEd.source.model == model && !gEd.source.pos.empty()) return true;
    Mesh mesh;
    if (!mesh.read((*gModels)[static_cast<size_t>(model)].path.string())) return false;
    SourceMesh source;
    source.model = model;
    const int count = static_cast<int>(mesh.vertices.size());
    source.pos.assign(static_cast<size_t>(count) * 3, 0.f);
    source.normal.assign(static_cast<size_t>(count) * 3, 0.f);
    source.edges.assign(static_cast<size_t>(count), {});
    Vec3 low{1e30f, 1e30f, 1e30f};
    Vec3 high{-1e30f, -1e30f, -1e30f};
    for (const Vertex& vertex : mesh.vertices) {
        if (vertex.index < 0 || vertex.index >= count) continue;
        const float x = static_cast<float>(vertex.position.x());
        const float y = static_cast<float>(vertex.position.y());
        const float z = static_cast<float>(vertex.position.z());
        source.pos[static_cast<size_t>(vertex.index) * 3] = x;
        source.pos[static_cast<size_t>(vertex.index) * 3 + 1] = y;
        source.pos[static_cast<size_t>(vertex.index) * 3 + 2] = z;
        low.x = std::min(low.x, x);
        low.y = std::min(low.y, y);
        low.z = std::min(low.z, z);
        high.x = std::max(high.x, x);
        high.y = std::max(high.y, y);
        high.z = std::max(high.z, z);
    }
    std::set<std::pair<int, int>> seen;
    for (const Face& face : mesh.faces) {
        if (face.isBoundary()) continue;
        int idx[3] = {-1, -1, -1};
        int corner = 0;
        HalfEdgeCIter he = face.he;
        do {
            if (corner < 3) idx[corner] = he->vertex->index;
            ++corner;
            he = he->next;
        } while (he != face.he);
        if (corner < 3 || idx[0] < 0 || idx[1] < 0 || idx[2] < 0) continue;
        source.tris.push_back(static_cast<unsigned>(idx[0]));
        source.tris.push_back(static_cast<unsigned>(idx[1]));
        source.tris.push_back(static_cast<unsigned>(idx[2]));
        const Vec3 a{source.pos[static_cast<size_t>(idx[0]) * 3], source.pos[static_cast<size_t>(idx[0]) * 3 + 1],
                     source.pos[static_cast<size_t>(idx[0]) * 3 + 2]};
        const Vec3 b{source.pos[static_cast<size_t>(idx[1]) * 3], source.pos[static_cast<size_t>(idx[1]) * 3 + 1],
                     source.pos[static_cast<size_t>(idx[1]) * 3 + 2]};
        const Vec3 c{source.pos[static_cast<size_t>(idx[2]) * 3], source.pos[static_cast<size_t>(idx[2]) * 3 + 1],
                     source.pos[static_cast<size_t>(idx[2]) * 3 + 2]};
        const Vec3 normal = vecNorm(vecCross(vecSub(b, a), vecSub(c, a)));
        for (int vertex : idx) {
            source.normal[static_cast<size_t>(vertex) * 3] += normal.x;
            source.normal[static_cast<size_t>(vertex) * 3 + 1] += normal.y;
            source.normal[static_cast<size_t>(vertex) * 3 + 2] += normal.z;
        }
        for (int edge = 0; edge < 3; ++edge) {
            const auto key = canonEdge(idx[edge], idx[(edge + 1) % 3]);
            if (!seen.insert(key).second) continue;
            source.edges[static_cast<size_t>(key.first)].push_back(key);
            source.edges[static_cast<size_t>(key.second)].push_back(key);
        }
    }
    for (int vertex = 0; vertex < count; ++vertex) {
        Vec3 normal{source.normal[static_cast<size_t>(vertex) * 3], source.normal[static_cast<size_t>(vertex) * 3 + 1],
                    source.normal[static_cast<size_t>(vertex) * 3 + 2]};
        normal = vecNorm(normal);
        source.normal[static_cast<size_t>(vertex) * 3] = normal.x;
        source.normal[static_cast<size_t>(vertex) * 3 + 1] = normal.y;
        source.normal[static_cast<size_t>(vertex) * 3 + 2] = normal.z;
    }
    source.center = vecScale(vecAdd(low, high), 0.5f);
    source.radius = std::max(0.05f, vecDot(vecSub(high, low), vecSub(high, low)));
    source.radius = std::sqrt(source.radius) * 0.5f;
    const GLuint vao = gEd.source.vao;
    const GLuint vbo = gEd.source.vbo;
    const GLuint nbo = gEd.source.nbo;
    const GLuint ebo = gEd.source.ebo;
    const bool uploaded = gEd.source.uploaded;
    gEd.source = std::move(source);
    gEd.source.vao = vao;
    gEd.source.vbo = vbo;
    gEd.source.nbo = nbo;
    gEd.source.ebo = ebo;
    gEd.source.uploaded = uploaded;
    return !gEd.source.tris.empty();
}

void uploadSource() {
    if (gEd.source.uploaded || gEd.source.pos.empty() || gEd.source.tris.empty()) return;
    if (!gEd.source.vao) glGenVertexArrays(1, &gEd.source.vao);
    if (!gEd.source.vbo) glGenBuffers(1, &gEd.source.vbo);
    if (!gEd.source.nbo) glGenBuffers(1, &gEd.source.nbo);
    if (!gEd.source.ebo) glGenBuffers(1, &gEd.source.ebo);
    glBindVertexArray(gEd.source.vao);
    glBindBuffer(GL_ARRAY_BUFFER, gEd.source.vbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(gEd.source.pos.size() * sizeof(float)), gEd.source.pos.data(),
                 GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glBindBuffer(GL_ARRAY_BUFFER, gEd.source.nbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(gEd.source.normal.size() * sizeof(float)), gEd.source.normal.data(),
                 GL_STATIC_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gEd.source.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(gEd.source.tris.size() * sizeof(unsigned)),
                 gEd.source.tris.data(), GL_STATIC_DRAW);
    gEd.source.triCount = static_cast<GLsizei>(gEd.source.tris.size());
    gEd.source.uploaded = true;
}

void requestSeamSolve(const std::vector<std::pair<int, int>>& edges, bool user) {
    if (!gEpoch || !gReady || !gCells || gEd.model < 0) return;
    int algo = 0;
    int count = 0;
    {
        std::lock_guard<std::mutex> lock(gBatch.mu);
        algo = gBatch.selected;
        count = gBatch.count;
    }
    int epoch = 0;
    {
        std::lock_guard<std::mutex> lock(gCutMu);
        gCuts[static_cast<size_t>(gEd.model)].user = user;
        gCuts[static_cast<size_t>(gEd.model)].edges = edges;
        epoch = gEpoch[static_cast<size_t>(gEd.model)].fetch_add(1) + 1;
    }
    gEditOwns[static_cast<size_t>(gEd.model)].store(1);
    for (int index = 0; index < kAlgoCount; ++index) {
        gReady[static_cast<size_t>(index * count + gEd.model)].store(0);
    }
    {
        std::lock_guard<std::mutex> lock(gEditMu);
        gEditJob.model = gEd.model;
        gEditJob.algo = algo;
        gEditJob.epoch = epoch;
        gEditJob.custom = user;
        gEditJob.edges = edges;
        gEditJobOn = true;
        gEditCv.notify_one();
    }
}

void commitSeams() {
    const std::vector<std::pair<int, int>> edges = seamVector(gEd.seams);
    if (edges == gEd.lastCommitted) return;
    gEd.lastCommitted = edges;
    requestSeamSolve(edges, edges != gEd.autoSeams);
}

void pushSeamHistory(const std::set<std::pair<int, int>>& before) {
    gEd.history.push_back(seamVector(before));
    if (gEd.history.size() > 100) gEd.history.erase(gEd.history.begin());
}

void extractTrunk() {
    if (gEd.seams.empty()) {
        gEd.status = "还没有切割边";
        return;
    }
    std::map<int, std::set<int>> adjacent;
    for (const auto& edge : gEd.seams) {
        adjacent[edge.first].insert(edge.second);
        adjacent[edge.second].insert(edge.first);
    }
    std::set<int> alive;
    std::queue<int> pending;
    for (const auto& item : adjacent) {
        alive.insert(item.first);
        if (item.second.size() <= 1) pending.push(item.first);
    }
    while (!pending.empty()) {
        const int vertex = pending.front();
        pending.pop();
        if (!alive.count(vertex) || adjacent[vertex].size() > 1) continue;
        alive.erase(vertex);
        const std::set<int> neighbors = adjacent[vertex];
        for (int neighbor : neighbors) {
            if (!alive.count(neighbor)) continue;
            adjacent[neighbor].erase(vertex);
            if (adjacent[neighbor].size() <= 1) pending.push(neighbor);
        }
        adjacent[vertex].clear();
    }
    std::set<std::pair<int, int>> keep;
    for (const auto& edge : gEd.seams) {
        if (alive.count(edge.first) && alive.count(edge.second)) keep.insert(edge);
    }
    if (keep.empty()) {
        gEd.status = "没有形成环路，无法抽取主干";
        return;
    }
    std::set<int> seen;
    std::set<int> best;
    int bestCount = -1;
    for (int start : alive) {
        if (!seen.insert(start).second) continue;
        std::set<int> component;
        std::vector<int> stack{start};
        while (!stack.empty()) {
            const int vertex = stack.back();
            stack.pop_back();
            component.insert(vertex);
            for (int neighbor : adjacent[vertex]) {
                if (!alive.count(neighbor) || !seen.insert(neighbor).second) continue;
                stack.push_back(neighbor);
            }
        }
        int count = 0;
        for (const auto& edge : keep) {
            if (component.count(edge.first) && component.count(edge.second)) ++count;
        }
        if (count > bestCount) {
            bestCount = count;
            best = std::move(component);
        }
    }
    if (bestCount >= 0 && bestCount < static_cast<int>(keep.size())) {
        std::set<std::pair<int, int>> filtered;
        for (const auto& edge : keep) {
            if (best.count(edge.first) && best.count(edge.second)) filtered.insert(edge);
        }
        keep.swap(filtered);
    }
    if (keep == gEd.seams) {
        gEd.status = "已是干净环路，共 " + std::to_string(keep.size()) + " 条";
        return;
    }
    const int removed = static_cast<int>(gEd.seams.size() - keep.size());
    pushSeamHistory(gEd.seams);
    gEd.seams = std::move(keep);
    gEd.status = "抽取主干，保留 " + std::to_string(gEd.seams.size()) + " 条，去掉 " + std::to_string(removed) + " 条";
    commitSeams();
}

void undoSeam() {
    if (gEd.history.empty()) {
        gEd.status = "没有可回退的选边";
        return;
    }
    const std::vector<std::pair<int, int>> previous = gEd.history.back();
    gEd.history.pop_back();
    gEd.seams.clear();
    gEd.seams.insert(previous.begin(), previous.end());
    gEd.status = "已回退，共 " + std::to_string(gEd.seams.size()) + " 条";
    commitSeams();
}

void openEditor(int model) {
    if (gEd.source.model != model) releaseSourceGl();
    if (!loadSource(model)) {
        gEd.status = "无法读取模型";
        return;
    }
    int selected = 0;
    int count = 0;
    {
        std::lock_guard<std::mutex> lock(gBatch.mu);
        selected = gBatch.selected;
        count = gBatch.count;
    }
    gEd.open = true;
    gEd.model = model;
    gEd.pick = false;
    gEd.painting = false;
    gEd.orbit = false;
    gEd.history.clear();
    gEd.seams.clear();
    const Cell& cell = (*gCells)[static_cast<size_t>(selected * count + model)];
    {
        std::lock_guard<std::mutex> lock(gCutMu);
        const ModelCut& cut = gCuts[static_cast<size_t>(model)];
        const std::vector<std::pair<int, int>>& shown = cut.user ? cut.edges : cell.seams;
        for (const auto& edge : shown) gEd.seams.insert(canonEdge(edge.first, edge.second));
        gEd.autoSeams = cut.haveAuto ? cut.automatic : seamVector(gEd.seams);
    }
    gEd.lastCommitted = seamVector(gEd.seams);
    gEd.yaw = model < static_cast<int>(gYaw.size()) ? gYaw[static_cast<size_t>(model)] : 0.7f;
    gEd.pitch = model < static_cast<int>(gPitch.size()) ? gPitch[static_cast<size_t>(model)] : 0.35f;
    gEd.dist = 2.6f;
    gEd.status = "拖拽旋转，滚轮缩放，C 选切割边。粉红线是当前切缝";
}

void finishPaint() {
    if (!gEd.paintMoved && gEd.paintButton == GLFW_MOUSE_BUTTON_LEFT && gEd.firstHit && gEd.firstWasSelected) {
        gEd.seams.erase(gEd.firstEdge);
    }
    if (!gEd.paintMoved && gEd.paintButton == GLFW_MOUSE_BUTTON_RIGHT && gEd.firstHit && !gEd.firstWasSelected) {
        gEd.status = "该边未被选中，右键用于擦除";
        gEd.painting = false;
        return;
    }
    gEd.painting = false;
    if (gEd.seams != gEd.strokeBefore) {
        pushSeamHistory(gEd.strokeBefore);
        gEd.status = "已选 " + std::to_string(gEd.seams.size()) + " 条，正在重新展开";
        commitSeams();
    } else {
        gEd.status = "已选 " + std::to_string(gEd.seams.size()) + " 条";
    }
}

void showModel(int index) {
    if (!gModels || gModels->empty() || !gReady) return;
    const int count = static_cast<int>(gModels->size());
    index = (index % count + count) % count;
    {
        std::lock_guard<std::mutex> lock(gBatch.mu);
        gBatch.focus = index;
        const bool have = gReady[static_cast<size_t>(gBatch.selected * gBatch.count + index)].load(std::memory_order_acquire) != 0;
        if (!have) {
            if (!gBatch.busy) beginBatchLocked();
            else if (gBatch.next > 0) gBatch.again = true;
        }
    }
    openEditor(index);
}

void handleSeamCmd(SeamCmd cmd) {
    if (cmd == SeamCmd::Back) {
        gEd.pick = false;
        gEd.painting = false;
        gEd.orbit = false;
        gEd.status = "拖拽旋转，滚轮缩放，C 选切割边。粉红线是当前切缝";
        return;
    }
    if (cmd == SeamCmd::Pick) {
        gEd.pick = !gEd.pick;
        gEd.painting = false;
        gEd.orbit = false;
        gEd.status = gEd.pick ? "左键绘制，右键擦除，S 回退，再按 C 可旋转" : "拖拽旋转，C 选切割边";
        return;
    }
    if (cmd == SeamCmd::Clear) {
        if (gEd.seams.empty()) {
            gEd.status = "当前没有切割线";
            return;
        }
        pushSeamHistory(gEd.seams);
        gEd.seams.clear();
        gEd.status = "已清空切割线";
        commitSeams();
        return;
    }
    extractTrunk();
}

void editLoop() {
    while (true) {
        EditJob job;
        {
            std::unique_lock<std::mutex> lock(gEditMu);
            gEditCv.wait(lock, [] { return gBatch.quit || gEditJobOn; });
            if (gBatch.quit) return;
            job = gEditJob;
            gEditJobOn = false;
        }
        Cell solved = solveOne((*gModels)[static_cast<size_t>(job.model)], kAlgos[job.algo], job.custom ? &job.edges : nullptr);
        if (!job.custom) {
            std::lock_guard<std::mutex> cutLock(gCutMu);
            if (!gCuts[static_cast<size_t>(job.model)].haveAuto) {
                gCuts[static_cast<size_t>(job.model)].automatic = solved.seams;
                std::sort(gCuts[static_cast<size_t>(job.model)].automatic.begin(),
                          gCuts[static_cast<size_t>(job.model)].automatic.end());
                gCuts[static_cast<size_t>(job.model)].haveAuto = true;
            }
        }
        solved.contentRev = job.epoch + 1;
        std::lock_guard<std::mutex> lock(gPendingMu);
        if (!gEpoch || gEpoch[static_cast<size_t>(job.model)].load() != job.epoch) continue;
        gPendingCell = std::move(solved);
        gPendingModel = job.model;
        gPendingAlgo = job.algo;
        gPendingEpoch = job.epoch;
        gPendingOn = true;
    }
}

void applyPendingEdit(std::vector<Cell>& cells, int modelCount) {
    Cell incoming;
    int model = -1;
    int algo = -1;
    {
        std::lock_guard<std::mutex> lock(gPendingMu);
        if (!gPendingOn) return;
        if (!gEpoch || gEpoch[static_cast<size_t>(gPendingModel)].load() != gPendingEpoch) {
            gPendingOn = false;
            return;
        }
        incoming = std::move(gPendingCell);
        model = gPendingModel;
        algo = gPendingAlgo;
        gPendingOn = false;
    }
    const int slot = algo * modelCount + model;
    releaseCellGl(cells[static_cast<size_t>(slot)]);
    cells[static_cast<size_t>(slot)] = std::move(incoming);
    gReady[static_cast<size_t>(slot)].store(1, std::memory_order_release);
    gEditOwns[static_cast<size_t>(model)].store(0, std::memory_order_release);
}

struct UvFit {
    float minU = 0, minV = 0, maxU = 1, maxV = 1, fit = 1, padX = 12, padY = 36;
};

UvFit fitUv(const std::vector<float>& uv, float cellW, float cellH) {
    UvFit fit;
    fit.minU = uv[0];
    fit.maxU = uv[0];
    fit.minV = uv[1];
    fit.maxV = uv[1];
    for (size_t i = 0; i + 1 < uv.size(); i += 2) {
        fit.minU = std::min(fit.minU, uv[i]);
        fit.maxU = std::max(fit.maxU, uv[i]);
        fit.minV = std::min(fit.minV, uv[i + 1]);
        fit.maxV = std::max(fit.maxV, uv[i + 1]);
    }
    const float spanX = std::max(fit.maxU - fit.minU, 1e-6f);
    const float spanY = std::max(fit.maxV - fit.minV, 1e-6f);
    const float innerX = cellW - 24.f;
    const float innerY = cellH - 72.f;
    fit.fit = std::min(innerX / spanX, innerY / spanY);
    fit.padX = 12.f + (innerX - spanX * fit.fit) * 0.5f;
    fit.padY = 36.f + (innerY - spanY * fit.fit) * 0.5f;
    return fit;
}

void uvToLocal(const UvFit& fit, float u, float v, float& x, float& y) {
    x = fit.padX + (u - fit.minU) * fit.fit;
    y = fit.padY + (fit.maxV - v) * fit.fit;
}

void cursorFramebuffer(GLFWwindow* window, float& sx, float& sy) {
    double x = 0, y = 0;
    glfwGetCursorPos(window, &x, &y);
    int winW = 1, winH = 1;
    glfwGetWindowSize(window, &winW, &winH);
    sx = static_cast<float>(x * gView.fbW / std::max(winW, 1));
    sy = static_cast<float>(y * gView.fbH / std::max(winH, 1));
}

void fitView() {
    const float availW = static_cast<float>(gView.fbW);
    const float availH = static_cast<float>(gView.fbH) - kHud * gView.contentScale;
    const float gw = gridWidth();
    const float gh = gridHeight(gRows);
    gView.scale = std::max(0.05f, std::min(availW / gw, availH / gh) * 0.98f);
    // 结果优先：缩到看不清时改为默认放大，允许画布溢出，按 F 可再回到全览。
    const float readable = kMinCellPoints * gView.contentScale / kCell;
    gView.scale = std::min(8.f, std::max(gView.scale, readable));
    const float cx = gw * 0.5f;
    const float cy = gh * 0.5f;
    const float sx = availW * 0.5f;
    const float sy = kHud * gView.contentScale + availH * 0.5f;
    gView.originX = cx - sx / gView.scale;
    gView.originY = cy - sy / gView.scale;
}

bool gPressDown = false;
int gPressButton = -1;
double gPressX = 0;
double gPressY = 0;
bool gPressDragged = false;

Rect detailButton(int model) {
    const float scale = gView.contentScale;
    const float sx = (cellX(model) - gView.originX) * gView.scale;
    const float sy = (cellY(model) - gView.originY) * gView.scale;
    const float sw = kCell * gView.scale;
    const float bw = 58.f * scale;
    const float bh = 24.f * scale;
    const float inset = 8.f * scale;
    return {sx + sw - bw - inset, sy + inset, bw, bh};
}

int frontCell(const std::function<bool(int)>& hit) {
    for (int model = gRows - 1; model >= 0; --model) {
        if (hit(model)) return model;
    }
    return -1;
}

int detailAt(float sx, float sy) {
    if (sy < kHud * gView.contentScale) return -1;
    return frontCell([&](int model) { return detailButton(model).contains(sx, sy); });
}

int cellAt(float sx, float sy) {
    if (sy < kHud * gView.contentScale) return -1;
    const float worldX = gView.originX + sx / std::max(gView.scale, 1e-4f);
    const float worldY = gView.originY + sy / std::max(gView.scale, 1e-4f);
    return frontCell([&](int model) {
        const Rect cell{cellX(model), cellY(model), kCell, kCell};
        return cell.contains(worldX, worldY);
    });
}

void scrollCallback(GLFWwindow* window, double xoff, double dy) {
    imguiScroll(window, xoff, dy);
    if (imguiCapturesMouse()) return;
    if (gEd.open) {
        gEd.dist = std::max(0.35f, std::min(18.f, gEd.dist * static_cast<float>(std::pow(0.9, dy))));
        return;
    }
    double x = 0, y = 0;
    glfwGetCursorPos(window, &x, &y);
    int winW = 1, winH = 1;
    glfwGetWindowSize(window, &winW, &winH);
    const float sx = static_cast<float>(x * gView.fbW / std::max(winW, 1));
    const float sy = static_cast<float>(y * gView.fbH / std::max(winH, 1));
    const float old = gView.scale;
    gView.scale = std::max(0.05f, std::min(8.f, gView.scale * static_cast<float>(std::pow(1.12, dy))));
    const float worldX = gView.originX + sx / old;
    const float worldY = gView.originY + sy / old;
    gView.originX = worldX - sx / gView.scale;
    gView.originY = worldY - sy / gView.scale;
    gView.fit = false;
}
void mouseButton(GLFWwindow* window, int button, int action, int mods) {
    imguiMouseButton(window, button, action, mods);
    float sx = 0, sy = 0;
    cursorFramebuffer(window, sx, sy);
    if (imguiCapturesMouse()) {
        if (action == GLFW_RELEASE) {
            gPressDown = false;
            gEd.orbit = false;
            gView.drag = false;
            if (gEd.painting) finishPaint();
        }
        return;
    }
    if (action == GLFW_PRESS) {
        gPressDown = true;
        gPressButton = button;
        gPressDragged = false;
        glfwGetCursorPos(window, &gPressX, &gPressY);
        gView.lastX = gPressX;
        gView.lastY = gPressY;
        if (gEd.open) {
            if (gEd.pick && (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_RIGHT)) {
                const Rect panel = panel3d();
                const Mat4 mvp = meshMvp(panel.w / std::max(panel.h, 1.f));
                gEd.painting = true;
                gEd.paintButton = button;
                gEd.paintMoved = false;
                gEd.paintX = sx;
                gEd.paintY = sy;
                gEd.strokeBefore = gEd.seams;
                gEd.firstHit = false;
                std::pair<int, int> edge;
                if (pickSeamEdge(sx, sy, panel, mvp, edge)) {
                    gEd.firstHit = true;
                    gEd.firstEdge = edge;
                    gEd.firstWasSelected = gEd.seams.count(edge) > 0;
                    if (button == GLFW_MOUSE_BUTTON_LEFT) gEd.seams.insert(edge);
                    else gEd.seams.erase(edge);
                }
                return;
            }
            if (button == GLFW_MOUSE_BUTTON_LEFT) gEd.orbit = true;
            return;
        }
        gPressDetail = detailAt(sx, sy);
        gDragCell = -1;
        if (gPressDetail >= 0) return;
        const int grabbed = cellAt(sx, sy);
        if (grabbed >= 0) gDragCell = grabbed;
        return;
    }
    if (!gPressDown || button != gPressButton) return;
    gPressDown = false;
    gEd.orbit = false;
    gView.drag = false;
    if (gEd.open) {
        if (gEd.painting) finishPaint();
        return;
    }
    if (button == GLFW_MOUSE_BUTTON_LEFT && !gPressDragged && gPressDetail >= 0 && detailAt(sx, sy) == gPressDetail) {
        int selected = 0;
        int count = 0;
        {
            std::lock_guard<std::mutex> lock(gBatch.mu);
            selected = gBatch.selected;
            count = gBatch.count;
        }
        if (gReady && count > 0 && gReady[static_cast<size_t>(selected * count + gPressDetail)].load()) openEditor(gPressDetail);
    } else if (button == GLFW_MOUSE_BUTTON_LEFT && !gPressDragged && gDragCell < 0 && gPressDetail < 0) {
        static double lastClick = 0;
        const double now = glfwGetTime();
        if (now - lastClick < 0.28 && cellAt(sx, sy) < 0) gView.fit = true;
        lastClick = now;
    }
    gDragCell = -1;
    gPressDetail = -1;
}

void cursorPos(GLFWwindow* window, double x, double y) {
    if (!gPressDown) return;
    const float dx = static_cast<float>(x - gView.lastX);
    const float dy = static_cast<float>(y - gView.lastY);
    if (!gPressDragged && (x - gPressX) * (x - gPressX) + (y - gPressY) * (y - gPressY) > 16.0) gPressDragged = true;
    if (gEd.open && gEd.painting) {
        float sx = 0, sy = 0;
        cursorFramebuffer(window, sx, sy);
        const float moved = (sx - static_cast<float>(gEd.paintX)) * (sx - static_cast<float>(gEd.paintX)) +
                            (sy - static_cast<float>(gEd.paintY)) * (sy - static_cast<float>(gEd.paintY));
        if (!gEd.paintMoved && moved > 16.f) gEd.paintMoved = true;
        if (gEd.paintMoved) {
            const Rect panel = panel3d();
            const Mat4 mvp = meshMvp(panel.w / std::max(panel.h, 1.f));
            std::pair<int, int> edge;
            if (pickSeamEdge(sx, sy, panel, mvp, edge)) {
                if (gEd.paintButton == GLFW_MOUSE_BUTTON_LEFT) gEd.seams.insert(edge);
                else gEd.seams.erase(edge);
            }
        }
        return;
    }
    if (gEd.open && gEd.orbit && gPressDragged) {
        gEd.yaw += dx * 0.01f;
        gEd.pitch = std::max(-1.15f, std::min(1.15f, gEd.pitch + dy * 0.01f));
        gView.lastX = x;
        gView.lastY = y;
        return;
    }
    if (!gEd.open && gDragCell >= 0 && gPressDragged && gPressButton == GLFW_MOUSE_BUTTON_LEFT) {
        if (gDragCell < static_cast<int>(gYaw.size())) {
            gYaw[static_cast<size_t>(gDragCell)] += dx * 0.01f;
            gPitch[static_cast<size_t>(gDragCell)] =
                std::max(-1.15f, std::min(1.15f, gPitch[static_cast<size_t>(gDragCell)] + dy * 0.01f));
        }
        gView.lastX = x;
        gView.lastY = y;
        return;
    }
    if (!gEd.open && gDragCell < 0 && gPressDetail < 0 && gPressDragged && gPressButton == GLFW_MOUSE_BUTTON_LEFT) {
        int winW = 1, winH = 1;
        glfwGetWindowSize(window, &winW, &winH);
        const float sdx = static_cast<float>((x - gView.lastX) * gView.fbW / std::max(winW, 1));
        const float sdy = static_cast<float>((y - gView.lastY) * gView.fbH / std::max(winH, 1));
        gView.originX -= sdx / gView.scale;
        gView.originY -= sdy / gView.scale;
        gView.lastX = x;
        gView.lastY = y;
        gView.fit = false;
        gView.drag = true;
    }
}
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    imguiKey(window, key, scancode, action, mods);
    if (imguiCapturesKeyboard()) return;
    if (action != GLFW_PRESS) return;
    if (gEd.open && key == GLFW_KEY_ESCAPE) {
        gEd.pick = false;
        gEd.painting = false;
        gEd.orbit = false;
        gEd.status = "拖拽旋转，滚轮缩放，C 选切割边。粉红线是当前切缝";
        return;
    }
    if (gEd.open && key == GLFW_KEY_C) {
        handleSeamCmd(SeamCmd::Pick);
        return;
    }
    if (gEd.open && key == GLFW_KEY_S) {
        undoSeam();
        return;
    }
    if (gEd.open && key == GLFW_KEY_T) {
        extractTrunk();
        return;
    }
    if (key == GLFW_KEY_F || key == GLFW_KEY_HOME) gView.fit = true;
}

GLuint compileShader(GLenum type, const char* source) {
    const GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);
    GLint ok = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[512];
        glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
        std::cerr << log << std::endl;
    }
    return shader;
}
void uploadFaceColors(Cell& cell, const std::vector<float>& colors, GLuint& vao, GLuint& vbo, GLsizei& count) {
    if (colors.size() != cell.tris.size() || cell.pos.size() < 3 || cell.tris.empty()) return;
    const size_t vertices = cell.pos.size() / 3;
    std::vector<float> packed(cell.tris.size() * 9);
    const size_t faceCount = cell.tris.size() / 3;
    for (size_t face = 0; face < faceCount; ++face) {
        unsigned index[3];
        float corner[3][3];
        for (int c = 0; c < 3; ++c) {
            index[c] = cell.tris[face * 3 + static_cast<size_t>(c)];
            if (index[c] >= vertices) return;
            corner[c][0] = cell.pos[static_cast<size_t>(index[c]) * 3];
            corner[c][1] = cell.pos[static_cast<size_t>(index[c]) * 3 + 1];
            corner[c][2] = cell.pos[static_cast<size_t>(index[c]) * 3 + 2];
        }
        const float e1x = corner[1][0] - corner[0][0];
        const float e1y = corner[1][1] - corner[0][1];
        const float e1z = corner[1][2] - corner[0][2];
        const float e2x = corner[2][0] - corner[0][0];
        const float e2y = corner[2][1] - corner[0][1];
        const float e2z = corner[2][2] - corner[0][2];
        float nx = e1y * e2z - e1z * e2y;
        float ny = e1z * e2x - e1x * e2z;
        float nz = e1x * e2y - e1y * e2x;
        const float length = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (length > 1e-8f) {
            nx /= length;
            ny /= length;
            nz /= length;
        }
        for (int c = 0; c < 3; ++c) {
            const size_t dst = (face * 3 + static_cast<size_t>(c)) * 9;
            packed[dst] = corner[c][0];
            packed[dst + 1] = corner[c][1];
            packed[dst + 2] = corner[c][2];
            packed[dst + 3] = nx;
            packed[dst + 4] = ny;
            packed[dst + 5] = nz;
            packed[dst + 6] = colors[face * 3];
            packed[dst + 7] = colors[face * 3 + 1];
            packed[dst + 8] = colors[face * 3 + 2];
        }
    }
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(packed.size() * sizeof(float)), packed.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), reinterpret_cast<void*>(6 * sizeof(float)));
    count = static_cast<GLsizei>(cell.tris.size());
}

void uploadCell(Cell& cell) {
    if (cell.uploaded || cell.uv.empty() || cell.tris.empty()) return;
    std::vector<unsigned> lines;
    lines.reserve(cell.tris.size() * 2);
    for (size_t i = 0; i < cell.tris.size(); i += 3) {
        const unsigned a = cell.tris[i], b = cell.tris[i + 1], c = cell.tris[i + 2];
        lines.insert(lines.end(), {a, b, b, c, c, a});
    }
    glGenVertexArrays(1, &cell.vao);
    glBindVertexArray(cell.vao);
    glGenBuffers(1, &cell.vbo);
    glBindBuffer(GL_ARRAY_BUFFER, cell.vbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(cell.uv.size() * sizeof(float)), cell.uv.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), nullptr);
    glGenBuffers(1, &cell.triEbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell.triEbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(cell.tris.size() * sizeof(unsigned)), cell.tris.data(), GL_STATIC_DRAW);
    cell.triCount = static_cast<GLsizei>(cell.tris.size());
    glGenBuffers(1, &cell.lineEbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell.lineEbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(lines.size() * sizeof(unsigned)), lines.data(), GL_STATIC_DRAW);
    cell.lineCount = static_cast<GLsizei>(lines.size());
    if (cell.faceColor.size() == cell.tris.size() && !cell.tris.empty()) {
        std::vector<float> packed(cell.tris.size() * 5);
        const size_t faceCount = cell.tris.size() / 3;
        for (size_t face = 0; face < faceCount; ++face) {
            const float red = cell.faceColor[face * 3];
            const float green = cell.faceColor[face * 3 + 1];
            const float blue = cell.faceColor[face * 3 + 2];
            for (int corner = 0; corner < 3; ++corner) {
                const unsigned vertex = cell.tris[face * 3 + static_cast<size_t>(corner)];
                const size_t dst = (face * 3 + static_cast<size_t>(corner)) * 5;
                packed[dst] = cell.uv[static_cast<size_t>(vertex) * 2];
                packed[dst + 1] = cell.uv[static_cast<size_t>(vertex) * 2 + 1];
                packed[dst + 2] = red;
                packed[dst + 3] = green;
                packed[dst + 4] = blue;
            }
        }
        glGenVertexArrays(1, &cell.distortVao);
        glBindVertexArray(cell.distortVao);
        glGenBuffers(1, &cell.distortVbo);
        glBindBuffer(GL_ARRAY_BUFFER, cell.distortVbo);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(packed.size() * sizeof(float)), packed.data(), GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), nullptr);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), reinterpret_cast<void*>(2 * sizeof(float)));
        cell.distortCount = static_cast<GLsizei>(cell.tris.size());
    }
    if (cell.pos.size() == cell.nrm.size() && cell.pos.size() == cell.uv.size() / 2 * 3 && !cell.tris.empty()) {
        std::vector<float> packed(cell.pos.size() / 3 * 8);
        const size_t vertices = cell.pos.size() / 3;
        for (size_t vertex = 0; vertex < vertices; ++vertex) {
            const size_t dst = vertex * 8;
            packed[dst] = cell.pos[vertex * 3];
            packed[dst + 1] = cell.pos[vertex * 3 + 1];
            packed[dst + 2] = cell.pos[vertex * 3 + 2];
            packed[dst + 3] = cell.nrm[vertex * 3];
            packed[dst + 4] = cell.nrm[vertex * 3 + 1];
            packed[dst + 5] = cell.nrm[vertex * 3 + 2];
            packed[dst + 6] = cell.uv[vertex * 2];
            packed[dst + 7] = cell.uv[vertex * 2 + 1];
        }
        glGenVertexArrays(1, &cell.surfVao);
        glBindVertexArray(cell.surfVao);
        glGenBuffers(1, &cell.surfVbo);
        glBindBuffer(GL_ARRAY_BUFFER, cell.surfVbo);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(packed.size() * sizeof(float)), packed.data(), GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), nullptr);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), reinterpret_cast<void*>(6 * sizeof(float)));
        cell.surfCount = static_cast<GLsizei>(cell.tris.size());
    }
    uploadFaceColors(cell, cell.faceColor, cell.qcVao, cell.qcVbo, cell.qcCount);
    uploadFaceColors(cell, cell.areaColor, cell.areaVao, cell.areaVbo, cell.areaCount);
    if (!cell.flipTris.empty() && cell.pos.size() >= 3 && cell.uv.size() >= 2) {
        const size_t vertices = cell.pos.size() / 3;
        std::vector<float> meshPacked;
        std::vector<float> uvPacked;
        meshPacked.reserve(cell.flipTris.size() * 6);
        uvPacked.reserve(cell.flipTris.size() * 6);
        bool complete = true;
        for (size_t face = 0; face + 2 < cell.flipTris.size(); face += 3) {
            unsigned index[3];
            float corner[3][3];
            float uvCorner[3][2];
            for (int c = 0; c < 3; ++c) {
                index[c] = cell.flipTris[face + static_cast<size_t>(c)];
                if (index[c] >= vertices || static_cast<size_t>(index[c]) * 2 + 1 >= cell.uv.size()) {
                    complete = false;
                    break;
                }
                corner[c][0] = cell.pos[static_cast<size_t>(index[c]) * 3];
                corner[c][1] = cell.pos[static_cast<size_t>(index[c]) * 3 + 1];
                corner[c][2] = cell.pos[static_cast<size_t>(index[c]) * 3 + 2];
                uvCorner[c][0] = cell.uv[static_cast<size_t>(index[c]) * 2];
                uvCorner[c][1] = cell.uv[static_cast<size_t>(index[c]) * 2 + 1];
            }
            if (!complete) break;
            const float e1x = corner[1][0] - corner[0][0];
            const float e1y = corner[1][1] - corner[0][1];
            const float e1z = corner[1][2] - corner[0][2];
            const float e2x = corner[2][0] - corner[0][0];
            const float e2y = corner[2][1] - corner[0][1];
            const float e2z = corner[2][2] - corner[0][2];
            float nx = e1y * e2z - e1z * e2y;
            float ny = e1z * e2x - e1x * e2z;
            float nz = e1x * e2y - e1y * e2x;
            const float length = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (length > 1e-8f) {
                nx /= length;
                ny /= length;
                nz /= length;
            }
            for (int c = 0; c < 3; ++c) {
                meshPacked.insert(meshPacked.end(), {corner[c][0], corner[c][1], corner[c][2], nx, ny, nz});
                uvPacked.insert(uvPacked.end(), {uvCorner[c][0], uvCorner[c][1], 1.f, 0.86f, 0.05f, 1.f});
            }
        }
        if (complete && !meshPacked.empty()) {
            glGenVertexArrays(1, &cell.flipVao);
            glBindVertexArray(cell.flipVao);
            glGenBuffers(1, &cell.flipVbo);
            glBindBuffer(GL_ARRAY_BUFFER, cell.flipVbo);
            glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(meshPacked.size() * sizeof(float)), meshPacked.data(), GL_STATIC_DRAW);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
            cell.flipCount = static_cast<GLsizei>(cell.flipTris.size());
            glGenVertexArrays(1, &cell.flipUvVao);
            glBindVertexArray(cell.flipUvVao);
            glGenBuffers(1, &cell.flipUvVbo);
            glBindBuffer(GL_ARRAY_BUFFER, cell.flipUvVbo);
            glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(uvPacked.size() * sizeof(float)), uvPacked.data(), GL_STATIC_DRAW);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);
            glEnableVertexAttribArray(2);
            glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(2 * sizeof(float)));
            cell.flipUvCount = static_cast<GLsizei>(cell.flipTris.size());
        }
    }
    cell.uploaded = true;
}

struct ScreenVert {
    float x, y, u, v, r, g, b, a;
};

void appendRect(std::vector<ScreenVert>& verts, float x, float y, float w, float h, float r, float g, float b, float a) {
    const ScreenVert corners[4] = {
        {x, y, 0, 0, r, g, b, a},
        {x + w, y, 1, 0, r, g, b, a},
        {x + w, y + h, 1, 1, r, g, b, a},
        {x, y + h, 0, 1, r, g, b, a},
    };
    const int order[6] = {0, 1, 2, 0, 2, 3};
    for (int index : order) verts.push_back(corners[index]);
}

void appendScreenLine(std::vector<ScreenVert>& verts, float x0, float y0, float x1, float y1, float width, float r,
                      float g, float b) {
    const float dx = x1 - x0;
    const float dy = y1 - y0;
    const float length = std::sqrt(dx * dx + dy * dy);
    if (length < 0.5f) return;
    const float nx = -dy / length * width;
    const float ny = dx / length * width;
    const ScreenVert corners[4] = {
        {x0 + nx, y0 + ny, 0, 0, r, g, b, 1},
        {x1 + nx, y1 + ny, 0, 0, r, g, b, 1},
        {x1 - nx, y1 - ny, 0, 0, r, g, b, 1},
        {x0 - nx, y0 - ny, 0, 0, r, g, b, 1},
    };
    const int order[6] = {0, 1, 2, 0, 2, 3};
    for (int index : order) verts.push_back(corners[index]);
}

void appendProjectedLines(std::vector<ScreenVert>& verts, const Mat4& mvp, float px, float py, float pw, float ph,
                          const std::vector<float>& segments, float half, float red, float green, float blue) {
    const Rect panel{px, py, pw, ph};
    for (size_t i = 0; i + 5 < segments.size(); i += 6) {
        float ax = 0, ay = 0, bx = 0, by = 0;
        const Vec3 a{segments[i], segments[i + 1], segments[i + 2]};
        const Vec3 bpt{segments[i + 3], segments[i + 4], segments[i + 5]};
        if (!projectToPanel(mvp, panel, a, ax, ay)) continue;
        if (!projectToPanel(mvp, panel, bpt, bx, by)) continue;
        appendScreenLine(verts, ax, ay, bx, by, half, red, green, blue);
    }
}

void workerLoop() {
    while (true) {
        int model = -1;
        int algo = 0;
        {
            std::unique_lock<std::mutex> lock(gBatch.mu);
            gBatch.cv.wait(lock, [] { return gBatch.quit || (gBatch.busy && gBatch.next < gBatch.jobs); });
            if (gBatch.quit) return;
            model = gBatch.focus;
            algo = gBatch.algo;
            ++gBatch.next;
        }
        const int slot = algo * gBatch.count + model;
        if (gReady[static_cast<size_t>(slot)].load(std::memory_order_acquire) == 0) {
            std::vector<std::pair<int, int>> custom;
            const std::vector<std::pair<int, int>>* customPtr = nullptr;
            int epoch = 0;
            {
                std::lock_guard<std::mutex> cutLock(gCutMu);
                epoch = gEpoch[static_cast<size_t>(model)].load(std::memory_order_acquire);
                if (gCuts[static_cast<size_t>(model)].user) {
                    custom = gCuts[static_cast<size_t>(model)].edges;
                    customPtr = &custom;
                }
            }
            if (gEditOwns[static_cast<size_t>(model)].load(std::memory_order_acquire)) {
                std::lock_guard<std::mutex> lock(gBatch.mu);
                finishBatchJobLocked();
                continue;
            }
            Cell solved = solveOne((*gModels)[static_cast<size_t>(model)], kAlgos[algo], customPtr);
            bool stale = gEditOwns[static_cast<size_t>(model)].load(std::memory_order_acquire) != 0;
            {
                std::lock_guard<std::mutex> cutLock(gCutMu);
                if (gEpoch[static_cast<size_t>(model)].load(std::memory_order_acquire) != epoch) stale = true;
                if (!stale && !gCuts[static_cast<size_t>(model)].user && !gCuts[static_cast<size_t>(model)].haveAuto) {
                    gCuts[static_cast<size_t>(model)].automatic = solved.seams;
                    std::sort(gCuts[static_cast<size_t>(model)].automatic.begin(),
                              gCuts[static_cast<size_t>(model)].automatic.end());
                    gCuts[static_cast<size_t>(model)].haveAuto = true;
                }
            }
            if (!stale) {
                Cell& destination = (*gCells)[static_cast<size_t>(slot)];
                destination.uv = std::move(solved.uv);
                destination.pos = std::move(solved.pos);
                destination.nrm = std::move(solved.nrm);
                destination.cx = solved.cx;
                destination.cy = solved.cy;
                destination.cz = solved.cz;
                destination.radius = solved.radius;
                destination.tris = std::move(solved.tris);
                destination.flipped = solved.flipped;
                destination.seconds = solved.seconds;
                destination.ok = solved.ok;
                destination.note = std::move(solved.note);
                destination.cutTag = std::move(solved.cutTag);
                destination.puncture = solved.puncture;
                destination.seams = std::move(solved.seams);
                destination.seamUv = std::move(solved.seamUv);
                destination.seamPos = std::move(solved.seamPos);
                destination.flipTris = std::move(solved.flipTris);
                destination.faceColor = std::move(solved.faceColor);
                destination.areaColor = std::move(solved.areaColor);
                destination.qcAvg = solved.qcAvg;
                destination.qcMax = solved.qcMax;
                destination.areaAvg = solved.areaAvg;
                destination.areaMax = solved.areaMax;
                destination.contentRev += 1;
                gReady[static_cast<size_t>(slot)].store(1, std::memory_order_release);
                std::cout << (*gModels)[static_cast<size_t>(model)].label << "  " << algoName(kAlgos[algo])
                          << "  " << (destination.ok ? "ok" : destination.note) << std::endl;
            }
        }
        std::lock_guard<std::mutex> lock(gBatch.mu);
        finishBatchJobLocked();
    }
}

} // namespace

int main(int argc, char** argv) {
    const std::vector<DiskModel> models = listDisks();
    const bool listOnly = argc > 1 && std::string(argv[1]) == "--list";
    if (listOnly) {
        for (const DiskModel& model : models) std::cout << model.label << "  " << model.faces << " 面" << std::endl;
        std::cout << "模型 " << models.size() << " 个" << std::endl;
        return models.empty() ? 1 : 0;
    }
    if (argc > 1 && std::string(argv[1]) == "--check") {
        if (models.empty()) return 1;
        int smallest = 0;
        for (int i = 1; i < static_cast<int>(models.size()); ++i) {
            if (models[i].faces < models[smallest].faces) smallest = i;
        }
        int failed = 0;
        for (Algo algo : kAlgos) {
            const Cell cell = solveOne(models[smallest], algo, nullptr);
            std::cout << models[smallest].label << "  " << algoName(algo) << "  "
                      << (cell.ok ? "ok" : cell.note) << "  翻转 " << cell.flipped
                      << "  K " << cell.qcAvg << "/" << cell.qcMax
                      << "  A " << cell.areaAvg << "/" << cell.areaMax
                      << "  " << std::llround(cell.seconds * 1000.0) << " ms" << std::endl;
            if (!cell.ok) ++failed;
        }
        return failed ? 1 : 0;
    }
    if (argc > 1 && std::string(argv[1]) == "--seams") {
        int failed = 0;
        for (const DiskModel& model : models) {
            const Cell cell = solveOne(model, Algo::Lscm, nullptr);
            std::cout << model.label << "  " << (cell.cutTag.empty() ? "圆盘" : cell.cutTag) << "  "
                      << (cell.ok ? "ok" : cell.note) << "  翻转 " << cell.flipped << "  缝线 "
                      << (cell.seamUv.size() / 4) << "  K " << cell.qcAvg << "/" << cell.qcMax
                      << "  A " << cell.areaAvg << "/" << cell.areaMax
                      << "  " << std::llround(cell.seconds * 1000.0) << " ms" << std::endl;
            if (!cell.ok) ++failed;
        }
        return failed ? 1 : 0;
    }
    if (models.empty()) {
        std::cerr << "assets/Models 里没有可读模型" << std::endl;
        return 1;
    }

    const int modelCount = static_cast<int>(models.size());
    std::vector<Cell> cells(static_cast<size_t>(modelCount * kAlgoCount));
    auto ready = std::make_unique<std::atomic<unsigned char>[]>(static_cast<size_t>(modelCount * kAlgoCount));
    for (int i = 0; i < modelCount * kAlgoCount; ++i) ready[static_cast<size_t>(i)].store(0);
    gModels = &models;
    gCells = &cells;
    gReady = ready.get();
    gCuts.assign(static_cast<size_t>(modelCount), {});
    gEpoch = std::make_unique<std::atomic<int>[]>(static_cast<size_t>(modelCount));
    gEditOwns = std::make_unique<std::atomic<unsigned char>[]>(static_cast<size_t>(modelCount));
    for (int index = 0; index < modelCount; ++index) {
        gEpoch[static_cast<size_t>(index)].store(0);
        gEditOwns[static_cast<size_t>(index)].store(0);
    }
    {
        std::lock_guard<std::mutex> lock(gBatch.mu);
        gBatch.count = modelCount;
        gBatch.selected = 0;
        beginBatchLocked();
    }
    gEditThread = std::thread(editLoop);
    const unsigned threads = std::max(2u, std::min(4u, std::thread::hardware_concurrency()));
    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) workers.emplace_back(workerLoop);
    auto stopWorkers = [&]() {
        {
            std::lock_guard<std::mutex> lock(gBatch.mu);
            gBatch.quit = true;
            gBatch.cv.notify_all();
        }
        gEditCv.notify_all();
        if (gEditThread.joinable()) gEditThread.join();
        for (std::thread& worker : workers) {
            if (worker.joinable()) worker.join();
        }
    };

    if (!glfwInit()) {
        std::cerr << "无法打开窗口" << std::endl;
        stopWorkers();
        return 1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_DEPTH_BITS, 24);
    GLFWwindow* window = glfwCreateWindow(1440, 900, "简单参数化", nullptr, nullptr);
    if (!window) {
        std::cerr << "无法打开窗口" << std::endl;
        glfwTerminate();
        stopWorkers();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetMouseButtonCallback(window, mouseButton);
    glfwSetCursorPosCallback(window, cursorPos);
    glfwSetCharCallback(window, imguiChar);
    imguiInit(window);
    glfwSetKeyCallback(window, keyCallback);
    gRows = 1;
    gYaw.assign(models.size(), 0.7f);
    gPitch.assign(models.size(), 0.35f);
    openEditor(0);

    const char* vertexSource = R"GLSL(
      #version 330 core
      layout(location = 0) in vec2 aPos;
      layout(location = 1) in vec2 aUv;
      layout(location = 2) in vec4 aColor;
      uniform int uMode;
      uniform vec2 uUvMin;
      uniform vec2 uUvMax;
      uniform vec2 uCellOrigin;
      uniform vec2 uCellSize;
      uniform vec2 uView;
      uniform float uScale;
      uniform vec2 uScreen;
      out vec2 vUv;
      out vec4 vColor;
      flat out int vMode;
      void main() {
        vMode = uMode;
        vColor = aColor;
        vec2 screen = aPos;
        vUv = aUv;
        if (uMode == 0 || uMode == 3 || uMode == 4) {
          vec2 span = max(uUvMax - uUvMin, vec2(1e-6));
          vec2 inner = vec2(uCellSize.x - 24.0, uCellSize.y - 72.0);
          float fit = min(inner.x / span.x, inner.y / span.y);
          vec2 draw = span * fit;
          vec2 pad = vec2(12.0, 36.0) + (inner - draw) * 0.5;
          vec2 local = pad + vec2((aPos.x - uUvMin.x) * fit, (uUvMax.y - aPos.y) * fit);
          vec2 world = uCellOrigin + local;
          screen = (world - uView) * uScale;
          vUv = aPos;
        }
        vec2 ndc = vec2(screen.x / uScreen.x * 2.0 - 1.0, 1.0 - screen.y / uScreen.y * 2.0);
        gl_Position = vec4(ndc, 0.0, 1.0);
      }
    )GLSL";
    const char* fragmentSource = R"GLSL(
      #version 330 core
      in vec2 vUv;
      in vec4 vColor;
      flat in int vMode;
      uniform sampler2D uTex;
      out vec4 frag;
      void main() {
        if (vMode == 0) {
          frag = vec4(0.90, 0.93, 0.96, 1.0);
        } else if (vMode == 3) {
          frag = vec4(0.22, 0.28, 0.36, 0.55);
        } else if (vMode == 2) {
          frag = texture(uTex, vUv);
        } else {
          frag = vColor;
        }
      }
    )GLSL";
    const GLuint program = glCreateProgram();
    const GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertexSource);
    const GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentSource);
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    GLuint screenVao = 0, screenVbo = 0;
    glGenVertexArrays(1, &screenVao);
    glGenBuffers(1, &screenVbo);
    glBindVertexArray(screenVao);
    glBindBuffer(GL_ARRAY_BUFFER, screenVbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(ScreenVert), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(ScreenVert), reinterpret_cast<void*>(2 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(ScreenVert), reinterpret_cast<void*>(4 * sizeof(float)));

    const GLint modeLoc = glGetUniformLocation(program, "uMode");
    const GLint uvMinLoc = glGetUniformLocation(program, "uUvMin");
    const GLint uvMaxLoc = glGetUniformLocation(program, "uUvMax");
    const GLint cellLoc = glGetUniformLocation(program, "uCellOrigin");
    const GLint cellSizeLoc = glGetUniformLocation(program, "uCellSize");
    const GLint viewLoc = glGetUniformLocation(program, "uView");
    const GLint scaleLoc = glGetUniformLocation(program, "uScale");
    const GLint screenLoc = glGetUniformLocation(program, "uScreen");
    const char* meshVertex = R"GLSL(
      #version 330 core
      layout(location = 0) in vec3 aPos;
      layout(location = 1) in vec3 aNrm;
      layout(location = 2) in vec2 aUv;
      layout(location = 3) in vec3 aColor;
      uniform mat4 uMvp;
      out vec3 vN;
      out vec2 vUv;
      flat out vec3 vColor;
      void main() {
        vN = aNrm;
        vUv = aUv;
        vColor = aColor;
        gl_Position = uMvp * vec4(aPos, 1.0);
      }
    )GLSL";
    const char* meshFragment = R"GLSL(
      #version 330 core
      in vec3 vN;
      in vec2 vUv;
      flat in vec3 vColor;
      uniform int uShade;
      out vec4 frag;
      void main() {
        if (uShade == 3) {
          frag = vec4(1.0, 0.86, 0.05, 1.0);
          return;
        }
        vec3 n = normalize(vN);
        vec3 light = normalize(vec3(0.25, 0.85, 0.45));
        float shade = abs(dot(n, light));
        float lambert = uShade == 2 ? 0.82 + 0.18 * shade : 0.42 + 0.58 * shade;
        float checker = mod(floor(vUv.x * 8.0) + floor(vUv.y * 8.0), 2.0);
        vec3 albedo = uShade == 0 ? mix(vec3(0.96, 0.97, 0.98), vec3(0.42, 0.60, 0.82), checker)
                     : uShade == 2 ? vColor
                                   : vec3(0.63, 0.72, 0.82);
        frag = vec4(albedo * lambert, 1.0);
      }
    )GLSL";
    const GLuint meshProgram = glCreateProgram();
    const GLuint meshVs = compileShader(GL_VERTEX_SHADER, meshVertex);
    const GLuint meshFs = compileShader(GL_FRAGMENT_SHADER, meshFragment);
    glAttachShader(meshProgram, meshVs);
    glAttachShader(meshProgram, meshFs);
    glLinkProgram(meshProgram);
    glDeleteShader(meshVs);
    glDeleteShader(meshFs);
    const GLint mvpLoc = glGetUniformLocation(meshProgram, "uMvp");
    const GLint shadeLoc = glGetUniformLocation(meshProgram, "uShade");
    GLFWcursor* crossCursor = glfwCreateStandardCursor(GLFW_CROSSHAIR_CURSOR);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    auto drawScreen = [&](const std::vector<ScreenVert>& verts, int mode, GLuint tex) {
        if (verts.empty()) return;
        glUniform1i(modeLoc, mode);
        if (mode == 2) {
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, tex);
        }
        glBindVertexArray(screenVao);
        glBindBuffer(GL_ARRAY_BUFFER, screenVbo);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(verts.size() * sizeof(ScreenVert)), verts.data(), GL_STREAM_DRAW);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(verts.size()));
    };

    while (!glfwWindowShouldClose(window)) {
        int winW = 1, winH = 1;
        glfwGetFramebufferSize(window, &gView.fbW, &gView.fbH);
        glfwGetWindowSize(window, &winW, &winH);
        float sx = 1.f, sy = 1.f;
        glfwGetWindowContentScale(window, &sx, &sy);
        gView.contentScale = std::max(sx, 1.f);
        if (gView.fbW < 1) gView.fbW = 1;
        if (gView.fbH < 1) gView.fbH = 1;
        applyPendingEdit(cells, modelCount);
        if (gView.fit && !gEd.open) fitView();
        glfwSetCursor(window, gEd.open && gEd.pick ? crossCursor : nullptr);

        int selected = 0;
        {
            std::lock_guard<std::mutex> lock(gBatch.mu);
            selected = gBatch.selected;
        }
        const int shown = (gEd.model >= 0 && gEd.model < modelCount) ? gEd.model : 0;
        const bool shownReady = ready[static_cast<size_t>(selected * modelCount + shown)].load(std::memory_order_acquire) != 0;
        std::string title = models[static_cast<size_t>(shown)].label + "  " + algoName(kAlgos[selected]);
        if (!shownReady) title += "  计算中";
        glfwSetWindowTitle(window, title.c_str());
        if (shownReady && gEd.open && gEd.model == shown && !gEd.painting && gEd.history.empty()) {
            bool userCut = false;
            {
                std::lock_guard<std::mutex> lock(gCutMu);
                userCut = gCuts[static_cast<size_t>(shown)].user;
            }
            const Cell& arrived = cells[static_cast<size_t>(selected * modelCount + shown)];
            if (!userCut && arrived.ok) {
                std::set<std::pair<int, int>> incoming;
                for (const auto& edge : arrived.seams) incoming.insert(canonEdge(edge.first, edge.second));
                if (incoming != gEd.seams) {
                    gEd.seams = std::move(incoming);
                    gEd.autoSeams = seamVector(gEd.seams);
                    gEd.lastCommitted = gEd.autoSeams;
                }
            }
        }

        glViewport(0, 0, gView.fbW, gView.fbH);
        glDisable(GL_DEPTH_TEST);
        glClearColor(0.90f, 0.90f, 0.888f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(program);
        glUniform2f(viewLoc, gView.originX, gView.originY);
        glUniform1f(scaleLoc, gView.scale);
        glUniform2f(screenLoc, static_cast<float>(gView.fbW), static_cast<float>(gView.fbH));

        auto worldRect = [&](float x, float y, float w, float h, float r, float g, float b) {
            const float s0x = (x - gView.originX) * gView.scale;
            const float s0y = (y - gView.originY) * gView.scale;
            std::vector<ScreenVert> verts;
            appendRect(verts, s0x, s0y, w * gView.scale, h * gView.scale, r, g, b, 1.f);
            drawScreen(verts, 1, 0);
        };
        auto ensureUploaded = [&](Cell& cell) {
            if (cell.contentRev == cell.uploadedRev) return;
            releaseCellGl(cell);
            if (cell.ok && !cell.uv.empty() && !cell.tris.empty()) uploadCell(cell);
            cell.uploadedRev = cell.contentRev;
        };
        auto drawUvLines = [&](const Cell& cell, const std::vector<float>& lines, float originX, float originY, float cellW,
                               float cellH, bool screenSpace, float red, float green, float blue, float half) {
            if (lines.empty() || cell.uv.empty()) return;
            const UvFit fit = fitUv(cell.uv, cellW, cellH);
            std::vector<ScreenVert> verts;
            for (size_t i = 0; i + 3 < lines.size(); i += 4) {
                float lx0 = 0, ly0 = 0, lx1 = 0, ly1 = 0;
                uvToLocal(fit, lines[i], lines[i + 1], lx0, ly0);
                uvToLocal(fit, lines[i + 2], lines[i + 3], lx1, ly1);
                float x0 = 0, y0 = 0, x1 = 0, y1 = 0;
                if (screenSpace) {
                    x0 = originX + lx0;
                    y0 = originY + ly0;
                    x1 = originX + lx1;
                    y1 = originY + ly1;
                } else {
                    x0 = (originX + lx0 - gView.originX) * gView.scale;
                    y0 = (originY + ly0 - gView.originY) * gView.scale;
                    x1 = (originX + lx1 - gView.originX) * gView.scale;
                    y1 = (originY + ly1 - gView.originY) * gView.scale;
                }
                appendScreenLine(verts, x0, y0, x1, y1, half, red, green, blue);
            }
            drawScreen(verts, 1, 0);
        };

        if (!gEd.open) for (int model = 0; model < gRows; ++model) {
            worldRect(cellX(model), cellY(model), kCell, kCell, 1.f, 1.f, 1.f);
            const int slot = selected * modelCount + model;
            if (!ready[static_cast<size_t>(slot)].load(std::memory_order_acquire)) continue;
            Cell& cell = cells[static_cast<size_t>(slot)];
            if (!cell.ok) {
                worldRect(cellX(model) + 8, cellY(model) + 8, kCell - 16, kCell - 16, 0.98f, 0.90f, 0.90f);
                continue;
            }
            ensureUploaded(cell);
            if (!cell.uploaded) continue;
            const bool quality = gQuality != QualityView::Checker;
            const GLsizei colorCount = gQuality == QualityView::Area ? cell.areaCount : cell.qcCount;
            if (cell.surfCount > 0 && (!quality || colorCount > 0)) {
                const float padTop = 36.f;
                const float padBot = 36.f;
                const float sx = (cellX(model) - gView.originX) * gView.scale;
                const float sy = (cellY(model) + padTop - gView.originY) * gView.scale;
                const float sw = kCell * gView.scale;
                const float sh = (kCell - padTop - padBot) * gView.scale;
                const int x0 = static_cast<int>(std::floor(sx));
                const int y0 = static_cast<int>(std::floor(sy));
                const int x1 = static_cast<int>(std::ceil(sx + sw));
                const int y1 = static_cast<int>(std::ceil(sy + sh));
                const int vx = std::max(0, x0);
                const int vy = std::max(0, static_cast<int>(gView.fbH) - y1);
                const int vw = std::min(gView.fbW, x1) - vx;
                const int vh = std::min(gView.fbH, static_cast<int>(gView.fbH) - y0) - vy;
                if (vw > 1 && vh > 1) {
                    const float yaw = model < static_cast<int>(gYaw.size()) ? gYaw[static_cast<size_t>(model)] : 0.7f;
                    const float pitch = model < static_cast<int>(gPitch.size()) ? gPitch[static_cast<size_t>(model)] : 0.35f;
                    const Mat4 mvp = previewMvp({cell.cx, cell.cy, cell.cz}, cell.radius, yaw, pitch, 2.6f, sw / sh);
                    glUseProgram(meshProgram);
                    glUniformMatrix4fv(mvpLoc, 1, GL_FALSE, mvp.m);
                    glUniform1i(shadeLoc, quality ? 2 : 0);
                    glViewport(vx, vy, vw, vh);
                    glEnable(GL_SCISSOR_TEST);
                    glScissor(vx, vy, vw, vh);
                    glEnable(GL_DEPTH_TEST);
                    glClear(GL_DEPTH_BUFFER_BIT);
                    if (quality) {
                        glBindVertexArray(gQuality == QualityView::Area ? cell.areaVao : cell.qcVao);
                        glDrawArrays(GL_TRIANGLES, 0, colorCount);
                    } else {
                        glBindVertexArray(cell.surfVao);
                        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell.triEbo);
                        glDrawElements(GL_TRIANGLES, cell.surfCount, GL_UNSIGNED_INT, nullptr);
                    }
                    if (cell.flipCount > 0) {
                        glEnable(GL_POLYGON_OFFSET_FILL);
                        glPolygonOffset(-1.f, -1.f);
                        glUniform1i(shadeLoc, 3);
                        glBindVertexArray(cell.flipVao);
                        glDrawArrays(GL_TRIANGLES, 0, cell.flipCount);
                        glDisable(GL_POLYGON_OFFSET_FILL);
                    }
                    glDisable(GL_DEPTH_TEST);
                    glViewport(0, 0, gView.fbW, gView.fbH);
                    glUseProgram(program);
                    glUniform2f(viewLoc, gView.originX, gView.originY);
                    glUniform1f(scaleLoc, gView.scale);
                    glUniform2f(screenLoc, static_cast<float>(gView.fbW), static_cast<float>(gView.fbH));
                    std::vector<ScreenVert> overlays;
                    const float seamHalf = std::max(2.4f, 2.2f * gView.contentScale);
                    appendProjectedLines(overlays, mvp, sx, sy, sw, sh, cell.seamPos, seamHalf, 1.f, 0.231f, 0.42f);
                    drawScreen(overlays, 1, 0);
                    glDisable(GL_SCISSOR_TEST);
                }
            } else {
                float minU = cell.uv[0], maxU = cell.uv[0], minV = cell.uv[1], maxV = cell.uv[1];
                for (size_t i = 0; i < cell.uv.size(); i += 2) {
                    minU = std::min(minU, cell.uv[i]);
                    maxU = std::max(maxU, cell.uv[i]);
                    minV = std::min(minV, cell.uv[i + 1]);
                    maxV = std::max(maxV, cell.uv[i + 1]);
                }
                glUniform2f(uvMinLoc, minU, minV);
                glUniform2f(uvMaxLoc, maxU, maxV);
                glUniform2f(cellLoc, cellX(model), cellY(model));
                glUniform2f(cellSizeLoc, kCell, kCell);
                if (gQuality == QualityView::Conformal && cell.distortCount > 0) {
                    glUniform1i(modeLoc, 4);
                    glBindVertexArray(cell.distortVao);
                    glDrawArrays(GL_TRIANGLES, 0, cell.distortCount);
                } else {
                    glUniform1i(modeLoc, 0);
                    glBindVertexArray(cell.vao);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell.triEbo);
                    glDrawElements(GL_TRIANGLES, cell.triCount, GL_UNSIGNED_INT, nullptr);
                }
                if (cell.flipUvCount > 0) {
                    glUniform1i(modeLoc, 4);
                    glBindVertexArray(cell.flipUvVao);
                    glDrawArrays(GL_TRIANGLES, 0, cell.flipUvCount);
                }
                glUniform1i(modeLoc, 3);
                glBindVertexArray(cell.vao);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell.lineEbo);
                glDrawElements(GL_LINES, cell.lineCount, GL_UNSIGNED_INT, nullptr);
                const float seamHalf = std::max(2.4f, 2.2f * gView.contentScale);
                drawUvLines(cell, cell.seamUv, cellX(model), cellY(model), kCell, kCell, false, 1.f, 0.231f, 0.42f, seamHalf);
            }
        }
        const float cellScreen = kCell * gView.scale;
        const float pointScale = std::max(gView.contentScale, 1.f);
        std::string uvNote;
        const FormulaText editorFormula = gEd.open ? formulaText(kAlgos[selected]) : FormulaText{"", ""};
        struct CellLabel {
            ImVec2 pos;
            ImVec2 clip0;
            ImVec2 clip1;
            std::string text;
        };
        std::vector<CellLabel> labels;
        auto addLabel = [&](const std::string& text, float worldX, float worldY, float clipWorldX, float clipWorldY,
                            float clipWorldW, float clipWorldH) {
            const float sx0 = (worldX - gView.originX) * gView.scale / pointScale;
            const float sy0 = (worldY - gView.originY) * gView.scale / pointScale;
            const float cx0 = (clipWorldX - gView.originX) * gView.scale / pointScale;
            const float cy0 = (clipWorldY - gView.originY) * gView.scale / pointScale;
            const float cw = clipWorldW * gView.scale / pointScale;
            const float ch = clipWorldH * gView.scale / pointScale;
            labels.push_back(CellLabel{ImVec2(sx0, sy0), ImVec2(cx0, cy0), ImVec2(cx0 + cw, cy0 + ch), text});
        };
        if (!gEd.open) for (int model = 0; model < gRows; ++model) {
            const auto detail = detailButton(model);
            std::vector<ScreenVert> detailVerts;
            appendRect(detailVerts, detail.x, detail.y, detail.w, detail.h, 0.86f, 0.91f, 0.98f, 1.f);
            drawScreen(detailVerts, 1, 0);
            const float nameRight = cellX(model) + kCell - (detail.w + 16.f) / std::max(gView.scale, 1e-4f);
            addLabel(shortName(models[static_cast<size_t>(model)].label), cellX(model) + 12.f, cellY(model) + 8.f,
                     cellX(model) + 8.f, cellY(model) + 4.f, std::max(8.f, nameRight - (cellX(model) + 8.f)), 28.f);
            labels.push_back(CellLabel{ImVec2((detail.x + 8.f) / pointScale, (detail.y + 2.f) / pointScale),
                                       ImVec2(detail.x / pointScale, detail.y / pointScale),
                                       ImVec2((detail.x + detail.w) / pointScale, (detail.y + detail.h) / pointScale),
                                       "详情"});
            if (cellScreen < 70.f) continue;
            const int slot = selected * modelCount + model;
            const float captionTop = cellY(model) + kCell - 42.f;
            if (!ready[static_cast<size_t>(slot)].load(std::memory_order_acquire)) {
                addLabel("计算中", cellX(model) + 12.f, captionTop, cellX(model) + 8.f, captionTop, kCell - 16.f, 36.f);
                continue;
            }
            const Cell& cell = cells[static_cast<size_t>(slot)];
            if (!cell.ok) {
                addLabel(cell.note.empty() ? "失败" : cell.note, cellX(model) + 12.f, captionTop, cellX(model) + 8.f,
                         captionTop, kCell - 16.f, 36.f);
                continue;
            }
            std::string caption = std::to_string(models[static_cast<size_t>(model)].faces) + "面  翻" +
                                  std::to_string(cell.flipped);
            if (!cell.cutTag.empty()) caption = cell.cutTag + "  " + caption;
            addLabel(caption, cellX(model) + 12.f, captionTop, cellX(model) + 8.f, captionTop, kCell - 16.f, 18.f);
            char second[80];
            if (gQuality == QualityView::Area) {
                std::snprintf(second, sizeof(second), "%lld ms  A %.2f/%.2f", std::llround(cell.seconds * 1000.0),
                              cell.areaAvg, cell.areaMax);
            } else if (gQuality == QualityView::Conformal) {
                std::snprintf(second, sizeof(second), "%lld ms  K %.2f/%.2f", std::llround(cell.seconds * 1000.0),
                              cell.qcAvg, cell.qcMax);
            } else {
                std::snprintf(second, sizeof(second), "%lld ms", std::llround(cell.seconds * 1000.0));
            }
            addLabel(second, cellX(model) + 12.f, captionTop + 16.f, cellX(model) + 8.f, captionTop + 16.f, kCell - 16.f,
                     18.f);
        } else if (gEd.model >= 0 && gEd.model < modelCount) {
            const auto view3d = panel3d();
            const auto viewUv = panelUv();
            uploadSource();
            if (gEd.source.uploaded) {
                Cell* paintedCell = nullptr;
                const int editSlot = selected * modelCount + gEd.model;
                if (ready[static_cast<size_t>(editSlot)].load(std::memory_order_acquire)) {
                    Cell& candidate = cells[static_cast<size_t>(editSlot)];
                    if (candidate.ok) {
                        ensureUploaded(candidate);
                        if (candidate.uploaded && candidate.surfCount > 0) paintedCell = &candidate;
                    }
                }
                std::vector<ScreenVert> meshBack;
                appendRect(meshBack, view3d.x, view3d.y, view3d.w, view3d.h, 0.86f, 0.89f, 0.92f, 1.f);
                drawScreen(meshBack, 1, 0);
                glUseProgram(meshProgram);
                const Mat4 mvp = meshMvp(view3d.w / std::max(view3d.h, 1.f));
                glUniformMatrix4fv(mvpLoc, 1, GL_FALSE, mvp.m);
                glUniform1i(shadeLoc, paintedCell ? 0 : 1);
                glViewport(static_cast<GLint>(view3d.x), static_cast<GLint>(gView.fbH - (view3d.y + view3d.h)),
                           static_cast<GLsizei>(view3d.w), static_cast<GLsizei>(view3d.h));
                glEnable(GL_SCISSOR_TEST);
                glScissor(static_cast<GLint>(view3d.x), static_cast<GLint>(gView.fbH - (view3d.y + view3d.h)),
                          static_cast<GLsizei>(view3d.w), static_cast<GLsizei>(view3d.h));
                glEnable(GL_DEPTH_TEST);
                glClear(GL_DEPTH_BUFFER_BIT);
                if (paintedCell) {
                    glBindVertexArray(paintedCell->surfVao);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, paintedCell->triEbo);
                    glDrawElements(GL_TRIANGLES, paintedCell->surfCount, GL_UNSIGNED_INT, nullptr);
                } else {
                    glBindVertexArray(gEd.source.vao);
                    glDrawElements(GL_TRIANGLES, gEd.source.triCount, GL_UNSIGNED_INT, nullptr);
                }
                if (paintedCell && paintedCell->flipCount > 0) {
                    glEnable(GL_POLYGON_OFFSET_FILL);
                    glPolygonOffset(-1.f, -1.f);
                    glUniform1i(shadeLoc, 3);
                    glBindVertexArray(paintedCell->flipVao);
                    glDrawArrays(GL_TRIANGLES, 0, paintedCell->flipCount);
                    glDisable(GL_POLYGON_OFFSET_FILL);
                }
                glDisable(GL_DEPTH_TEST);
                glDisable(GL_SCISSOR_TEST);
                glViewport(0, 0, gView.fbW, gView.fbH);
                glUseProgram(program);
                glUniform2f(viewLoc, gView.originX, gView.originY);
                glUniform1f(scaleLoc, gView.scale);
                glUniform2f(screenLoc, static_cast<float>(gView.fbW), static_cast<float>(gView.fbH));
                std::vector<ScreenVert> seamVerts;
                const float half = std::max(1.8f, 2.1f * gView.contentScale);
                for (const auto& edge : gEd.seams) {
                    float ax = 0, ay = 0, bx = 0, by = 0;
                    if (!projectToPanel(mvp, view3d, positionOf(edge.first), ax, ay)) continue;
                    if (!projectToPanel(mvp, view3d, positionOf(edge.second), bx, by)) continue;
                    appendScreenLine(seamVerts, ax, ay, bx, by, half, 1.f, 0.231f, 0.42f);
                }
                glEnable(GL_SCISSOR_TEST);
                glScissor(static_cast<GLint>(view3d.x), static_cast<GLint>(gView.fbH - (view3d.y + view3d.h)),
                          static_cast<GLsizei>(view3d.w), static_cast<GLsizei>(view3d.h));
                drawScreen(seamVerts, 1, 0);
                glDisable(GL_SCISSOR_TEST);
            }
            std::vector<ScreenVert> uvBack;
            appendRect(uvBack, viewUv.x, viewUv.y, viewUv.w, viewUv.h, 0.96f, 0.96f, 0.94f, 1.f);
            drawScreen(uvBack, 1, 0);
            const float cardH = 156.f * gView.contentScale;
            const float chartY = viewUv.y + cardH;
            const float chartH = std::max(1.f, viewUv.h - cardH);
            const int slot = selected * modelCount + gEd.model;
            if (ready[static_cast<size_t>(slot)].load(std::memory_order_acquire)) {
                Cell& cell = cells[static_cast<size_t>(slot)];
                if (cell.ok) {
                    ensureUploaded(cell);
                    if (cell.uploaded) {
                        const UvFit fit = fitUv(cell.uv, viewUv.w, chartH);
                        glUniform2f(uvMinLoc, fit.minU, fit.minV);
                        glUniform2f(uvMaxLoc, fit.maxU, fit.maxV);
                        glUniform2f(viewLoc, 0.f, 0.f);
                        glUniform1f(scaleLoc, 1.f);
                        glUniform2f(cellLoc, viewUv.x, chartY);
                        glUniform2f(cellSizeLoc, viewUv.w, chartH);
                        glUniform1i(modeLoc, 0);
                        glBindVertexArray(cell.vao);
                        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell.triEbo);
                        glDrawElements(GL_TRIANGLES, cell.triCount, GL_UNSIGNED_INT, nullptr);
                        if (cell.flipUvCount > 0) {
                            glUniform1i(modeLoc, 4);
                            glBindVertexArray(cell.flipUvVao);
                            glDrawArrays(GL_TRIANGLES, 0, cell.flipUvCount);
                        }
                        glUniform1i(modeLoc, 3);
                        glBindVertexArray(cell.vao);
                        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell.lineEbo);
                        glDrawElements(GL_LINES, cell.lineCount, GL_UNSIGNED_INT, nullptr);
                        glUniform2f(viewLoc, gView.originX, gView.originY);
                        glUniform1f(scaleLoc, gView.scale);
                        const float seamHalf = std::max(2.4f, 2.2f * gView.contentScale);
                        drawUvLines(cell, cell.seamUv, viewUv.x, chartY, viewUv.w, chartH, true, 1.f, 0.231f, 0.42f, seamHalf);
                    }
                }
                uvNote = cell.ok ? (cell.cutTag.empty() ? "展开" : cell.cutTag) : (cell.note.empty() ? "失败" : cell.note);
            } else {
                uvNote = "正在按切缝重新展开";
            }
        }
        const Cell& shownCell = cells[static_cast<size_t>(selected * modelCount + shown)];
        std::string hudText = "正在展开当前模型";
        if (shownReady && gEd.pick) hudText = gEd.status;
        else if (shownReady && !shownCell.ok) hudText = shownCell.note.empty() ? "展开失败" : shownCell.note;
        else if (shownReady) {
            hudText = std::string(algoName(kAlgos[selected])) + " · " +
                      std::to_string(std::llround(shownCell.seconds * 1000.0)) + " ms · 翻转 " +
                      std::to_string(shownCell.flipped);
            if (!shownCell.cutTag.empty()) hudText += " · " + shownCell.cutTag;
            if (!gEd.status.empty() && gEd.status.rfind("拖拽", 0) != 0) hudText = gEd.status;
        }
        std::string modelLabel = models[static_cast<size_t>(shown)].label;
        const auto dot = modelLabel.rfind(".obj");
        if (dot != std::string::npos) modelLabel.erase(dot, 4);
        imguiBegin();
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(16.f, 12.f));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(12.f, 6.f));
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 7.f);
        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.f);
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(8.f, 8.f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.f);
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.975f, 0.972f, 0.964f, 1.f));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.84f, 0.835f, 0.816f, 1.f));
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(1.f, 1.f, 1.f, 1.f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.955f, 0.953f, 0.945f, 1.f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.91f, 0.905f, 0.89f, 1.f));
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.14f, 0.15f, 0.17f, 1.f));
        ImGui::SetNextWindowPos(ImVec2(0.f, 0.f));
        ImGui::SetNextWindowSize(ImVec2(static_cast<float>(gView.fbW) / pointScale, toolbarPoints()));
        ImGui::Begin("工具", nullptr,
                     ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings |
                         ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
        ImGui::SetWindowFontScale(1.f / pointScale);
        auto darkButton = [](const char* label) {
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.15f, 0.18f, 0.22f, 1.f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.24f, 0.28f, 0.33f, 1.f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.09f, 0.11f, 0.14f, 1.f));
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.98f, 0.98f, 0.96f, 1.f));
            const bool pressed = ImGui::Button(label);
            ImGui::PopStyleColor(4);
            return pressed;
        };
        auto toggleButton = [&](const char* label, bool on) { return on ? darkButton(label) : ImGui::Button(label); };
        if (ImGui::Button("<##prev")) showModel(shown - 1);
        ImGui::SameLine();
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted(modelLabel.c_str());
        ImGui::SameLine();
        if (ImGui::Button(">##next")) showModel(shown + 1);
        ImGui::SameLine(0.f, 20.f);
        for (int index = 0; index < kAlgoCount; ++index) {
            if (index > 0) ImGui::SameLine();
            if (toggleButton(algoName(kAlgos[index]), index == selected)) selectAlgo(index);
        }
        if (toggleButton("棋盘格", gQuality == QualityView::Checker)) gQuality = QualityView::Checker;
        ImGui::SameLine();
        if (toggleButton("共形畸变", gQuality == QualityView::Conformal)) gQuality = QualityView::Conformal;
        ImGui::SameLine();
        if (toggleButton("面积畸变", gQuality == QualityView::Area)) gQuality = QualityView::Area;
        ImGui::SameLine(0.f, 20.f);
        if (toggleButton(gEd.pick ? "结束选边" : "选切割边", gEd.pick)) handleSeamCmd(SeamCmd::Pick);
        ImGui::SameLine();
        if (ImGui::Button("抽取主干")) handleSeamCmd(SeamCmd::Trunk);
        ImGui::SameLine();
        if (ImGui::Button("清空")) handleSeamCmd(SeamCmd::Clear);
        const float legendReserve = 72.f;
        const float statusWidth = std::max(40.f, ImGui::GetContentRegionAvail().x - legendReserve);
        std::string status = hudText;
        if (ImGui::CalcTextSize(status.c_str()).x > statusWidth) {
            const std::string ellipsis = "…";
            while (!status.empty() && ImGui::CalcTextSize((status + ellipsis).c_str()).x > statusWidth) {
                status.pop_back();
                while (!status.empty() && (static_cast<unsigned char>(status.back()) & 0xC0) == 0x80) status.pop_back();
            }
            status += ellipsis;
        }
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.40f, 0.42f, 0.45f, 1.f));
        ImGui::TextUnformatted(status.c_str());
        ImGui::PopStyleColor();
        ImDrawList* barDraw = ImGui::GetWindowDrawList();
        const ImVec2 legend(ImGui::GetWindowPos().x + ImGui::GetWindowWidth() - legendReserve, ImGui::GetItemRectMin().y);
        barDraw->AddCircleFilled(ImVec2(legend.x + 4.f, legend.y + ImGui::GetFontSize() * 0.5f), 4.f, IM_COL32(255, 59, 107, 255), 16);
        barDraw->AddText(ImVec2(legend.x + 14.f, legend.y), IM_COL32(92, 84, 88, 255), "切缝");
        const ImVec2 bar = ImGui::GetWindowPos();
        barDraw->AddLine(ImVec2(bar.x, bar.y + ImGui::GetWindowHeight() - 1.f),
                         ImVec2(bar.x + ImGui::GetWindowWidth(), bar.y + ImGui::GetWindowHeight() - 1.f),
                         IM_COL32(214, 212, 206, 255), 1.f);
        ImGui::End();
        ImGui::PopStyleColor(6);
        ImGui::PopStyleVar(6);
        ImDrawList* overlay = ImGui::GetForegroundDrawList();
        for (const CellLabel& label : labels) {
            overlay->PushClipRect(label.clip0, label.clip1, true);
            overlay->AddText(label.pos, IM_COL32(30, 34, 40, 255), label.text.c_str());
            overlay->PopClipRect();
        }
        if (gEd.open && gEd.model >= 0) {
            const Rect viewUv = panelUv();
            ImGui::SetNextWindowPos(ImVec2(viewUv.x / pointScale, viewUv.y / pointScale));
            ImGui::SetNextWindowSize(ImVec2(viewUv.w / pointScale, 156.f));
            ImGui::Begin("公式", nullptr,
                         ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
            ImGui::SetWindowFontScale(1.f / pointScale);
            ImGui::TextUnformatted(editorFormula.title);
            ImGui::TextWrapped("%s", editorFormula.body);
            ImGui::End();
            if (!uvNote.empty()) {
                overlay->AddText(ImVec2((viewUv.x + 16.f) / pointScale, (viewUv.y + 156.f * pointScale + 8.f) / pointScale),
                                 IM_COL32(30, 34, 40, 255), uvNote.c_str());
            }
        }
        imguiEnd();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    imguiShutdown();
    stopWorkers();
    if (crossCursor) glfwDestroyCursor(crossCursor);
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
