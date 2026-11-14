#define GL_SILENCE_DEPRECATION
#define GLFW_INCLUDE_NONE
#include <OpenGL/gl3.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_chrome.h"

#include "Cetm.h"
#include "CirclePatterns.h"
#include "IncrementalFlattening.h"
#include "Mesh.h"
#include "QcError.h"
#include "RicciFlow.h"
#include "models.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr float kPi = 3.14159265358979323846f;
constexpr double kDefaultCone = 4.0 * M_PI;

enum class Method { Circle, Cetm, Ricci, Incremental };
enum class Texture { Checker, Distortion };

struct DrawVertex {
    float px, py, pz, nx, ny, nz, u, v, r, g, b;
};

struct Mat4 {
    float m[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
};

struct Vec3 {
    float x = 0, y = 0, z = 0;
};

Mesh gMesh;
std::vector<ModelFile> gModels;
int gModel = 0;
Method gMethod = Method::Circle;
int gScheme = NEWTON;
Texture gTexture = Texture::Checker;
bool gPickCones = false;
bool gShowExtra = false;
bool gAutoRun = true;
bool gHasResult = false;
bool gDrag = false;
int gDragPanel = -1;
double gLastX = 0, gLastY = 0;
float gYaw = 0.7f, gPitch = 0.35f, gRadius = 3.2f;
float gUvZoom = 1.f, gUvPanX = 0.f, gUvPanY = 0.f;
Vec3 gCenter;
std::unordered_map<int, double> gCones;
std::string gStatus = "选择模型后点「展开」。";
double gSeconds = 0;
int gFlips = 0;
int gWidth = 1440, gHeight = 860;
float gScale = 1.f;
std::vector<DrawVertex> gTriangles, gUvTriangles, gFlipTris, gFlips3d, gConeMarks;
bool gSolving = false;
bool gUnfoldNow = false;
GLFWwindow* gWindow = nullptr;

void drawFrame();
void reportCircleProgress(int completed, int total);
void reportIncrementalProgress(int phase, int completed, int total);

GLuint gProgram = 0, gVao = 0, gBuffer = 0;

Mat4 multiply(const Mat4& a, const Mat4& b) {
    Mat4 out;
    for (int column = 0; column < 4; ++column) {
        for (int row = 0; row < 4; ++row) {
            out.m[column * 4 + row] = 0;
            for (int k = 0; k < 4; ++k) out.m[column * 4 + row] += a.m[k * 4 + row] * b.m[column * 4 + k];
        }
    }
    return out;
}

Mat4 perspective(float fov, float aspect, float nearPlane, float farPlane) {
    const float f = 1.f / std::tan(fov * 0.5f);
    Mat4 out{};
    out.m[0] = f / std::max(aspect, 1e-4f);
    out.m[5] = f;
    out.m[10] = (farPlane + nearPlane) / (nearPlane - farPlane);
    out.m[11] = -1.f;
    out.m[14] = (2.f * farPlane * nearPlane) / (nearPlane - farPlane);
    out.m[15] = 0.f;
    return out;
}

Mat4 lookAt(Vec3 eye, Vec3 center) {
    Vec3 forward{center.x - eye.x, center.y - eye.y, center.z - eye.z};
    const float fl = std::sqrt(forward.x * forward.x + forward.y * forward.y + forward.z * forward.z);
    forward.x /= fl; forward.y /= fl; forward.z /= fl;
    Vec3 up{0, 1, 0};
    Vec3 side{forward.y * up.z - forward.z * up.y, forward.z * up.x - forward.x * up.z, forward.x * up.y - forward.y * up.x};
    const float sl = std::sqrt(side.x * side.x + side.y * side.y + side.z * side.z);
    side.x /= sl; side.y /= sl; side.z /= sl;
    up = {side.y * forward.z - side.z * forward.y, side.z * forward.x - side.x * forward.z, side.x * forward.y - side.y * forward.x};
    Mat4 out;
    out.m[0] = side.x; out.m[4] = side.y; out.m[8] = side.z; out.m[12] = -(side.x * eye.x + side.y * eye.y + side.z * eye.z);
    out.m[1] = up.x; out.m[5] = up.y; out.m[9] = up.z; out.m[13] = -(up.x * eye.x + up.y * eye.y + up.z * eye.z);
    out.m[2] = -forward.x; out.m[6] = -forward.y; out.m[10] = -forward.z;
    out.m[14] = forward.x * eye.x + forward.y * eye.y + forward.z * eye.z;
    out.m[15] = 1;
    return out;
}

Mat4 ortho(float left, float right, float bottom, float top) {
    Mat4 out{};
    out.m[0] = 2.f / (right - left);
    out.m[5] = 2.f / (top - bottom);
    out.m[10] = -1.f;
    out.m[12] = -(right + left) / (right - left);
    out.m[13] = -(top + bottom) / (top - bottom);
    out.m[15] = 1.f;
    return out;
}

Vec3 eyePosition() {
    const float cp = std::cos(gPitch), sp = std::sin(gPitch);
    return {gCenter.x + gRadius * cp * std::sin(gYaw), gCenter.y + gRadius * sp,
            gCenter.z + gRadius * cp * std::cos(gYaw)};
}

void fitCamera() {
    if (gMesh.vertices.empty()) return;
    Eigen::Vector3d lo = gMesh.vertices[0].position, hi = lo;
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    for (VertexCIter vertex = gMesh.vertices.begin(); vertex != gMesh.vertices.end(); ++vertex) {
        lo = lo.cwiseMin(vertex->position);
        hi = hi.cwiseMax(vertex->position);
        sum += vertex->position;
    }
    const Eigen::Vector3d mid = sum / std::max(1, static_cast<int>(gMesh.vertices.size()));
    gCenter = {static_cast<float>(mid.x()), static_cast<float>(mid.y()), static_cast<float>(mid.z())};
    gRadius = static_cast<float>((hi - lo).norm()) * 1.8f;
    if (gRadius < 0.2f) gRadius = 2.f;
}

const char* methodName(Method method) {
    if (method == Method::Cetm) return "CETM";
    if (method == Method::Ricci) return "Ricci 流";
    if (method == Method::Incremental) return "Incremental Flattening";
    return "Circle Pattern";
}

const char* schemeName(int scheme) {
    if (scheme == GRAD_DESCENT) return "梯度下降";
    if (scheme == LBFGS) return "LBFGS";
    return "牛顿法";
}

void rebuildGeometry();

void loadModel(int index) {
    if (gModels.empty()) return;
    gModel = (index % static_cast<int>(gModels.size()) + static_cast<int>(gModels.size())) % static_cast<int>(gModels.size());
    gMesh.read(gModels[gModel].path);
    gCones.clear();
    gHasResult = false;
    gSeconds = 0;
    gFlips = 0;
    gPickCones = false;
    gUvZoom = 1.f;
    gUvPanX = gUvPanY = 0.f;
    fitCamera();
    gStatus = gModels[gModel].label + " · " + std::to_string(gMesh.faces.size()) + " 面。点「展开」才计算。";
    rebuildGeometry();
}

double cornerAngle(const Eigen::Vector3d& at, const Eigen::Vector3d& b, const Eigen::Vector3d& c) {
    const Eigen::Vector3d u = b - at, v = c - at;
    const double lu = u.norm(), lv = v.norm();
    if (lu < 1e-12 || lv < 1e-12) return 0;
    return std::acos(std::clamp(u.dot(v) / (lu * lv), -1.0, 1.0));
}

void oneRing(int index, std::vector<int>& neighbors) {
    neighbors.clear();
    const Vertex& vertex = gMesh.vertices[index];
    if (vertex.isIsolated()) return;
    HalfEdgeCIter halfEdge = vertex.he;
    const HalfEdgeCIter first = halfEdge;
    do {
        neighbors.push_back(halfEdge->next->vertex->index);
        halfEdge = halfEdge->flip->next;
    } while (halfEdge != first);
}

void placeBenChenCones() {
    if (gMesh.vertices.empty()) return;
    const int count = static_cast<int>(gMesh.vertices.size());
    std::vector<double> angleSum(count, 0.0);
    std::vector<char> boundary(count, 0);
    for (VertexCIter vertex = gMesh.vertices.begin(); vertex != gMesh.vertices.end(); ++vertex) {
        boundary[vertex->index] = vertex->isBoundary() ? 1 : 0;
    }
    for (FaceCIter face = gMesh.faces.begin(); face != gMesh.faces.end(); ++face) {
        if (face->isBoundary()) continue;
        HalfEdgeCIter he = face->he;
        const int ids[3] = {he->vertex->index, he->next->vertex->index, he->next->next->vertex->index};
        const Eigen::Vector3d p[3] = {he->vertex->position, he->next->vertex->position, he->next->next->vertex->position};
        for (int corner = 0; corner < 3; ++corner) {
            angleSum[ids[corner]] += cornerAngle(p[corner], p[(corner + 1) % 3], p[(corner + 2) % 3]);
        }
    }
    std::vector<double> curvature(count, 0.0);
    double curvatureSum = 0.0;
    int interiorCount = 0;
    for (int index = 0; index < count; ++index) {
        if (boundary[index]) continue;
        curvature[index] = 2.0 * M_PI - angleSum[index];
        curvatureSum += curvature[index];
        ++interiorCount;
    }

    struct Candidate { int index; double score; Eigen::Vector3d position; };
    std::vector<Candidate> ranked;
    ranked.reserve(interiorCount);
    Eigen::Vector3d lo = gMesh.vertices[0].position, hi = lo;
    for (int index = 0; index < count; ++index) {
        lo = lo.cwiseMin(gMesh.vertices[index].position);
        hi = hi.cwiseMax(gMesh.vertices[index].position);
        if (boundary[index]) continue;
        ranked.push_back({index, std::abs(curvature[index]), gMesh.vertices[index].position});
    }
    std::sort(ranked.begin(), ranked.end(), [](const Candidate& a, const Candidate& b) {
        return a.score > b.score;
    });
    const int limit = std::clamp(static_cast<int>(std::ceil(std::sqrt(static_cast<double>(count)) / 8.0)), 4, 24);
    const double diagonal = std::max((hi - lo).norm(), 1e-6);
    std::vector<int> chosen;
    auto pickSeparated = [&](double separation) {
        chosen.clear();
        for (const Candidate& candidate : ranked) {
            if (static_cast<int>(chosen.size()) >= limit) break;
            bool crowded = false;
            for (int index : chosen) {
                if ((gMesh.vertices[index].position - candidate.position).norm() < separation) crowded = true;
            }
            if (!crowded) chosen.push_back(candidate.index);
        }
    };
    double separation = 0.08 * diagonal;
    const int want = std::min(limit, interiorCount);
    pickSeparated(separation);
    while (static_cast<int>(chosen.size()) < want && separation > diagonal * 0.01) {
        separation *= 0.5;
        pickSeparated(separation);
    }

    std::vector<char> isCone(count, 0);
    for (int index : chosen) isCone[index] = 1;
    std::vector<char> seen(count, 0);
    std::vector<int> stack;
    std::vector<int> ring;
    for (int index = 0; index < count; ++index) {
        if (boundary[index] || seen[index]) continue;
        stack.clear();
        stack.push_back(index);
        seen[index] = 1;
        bool hasCone = false;
        int strongest = index;
        for (size_t cursor = 0; cursor < stack.size(); ++cursor) {
            const int vertex = stack[cursor];
            if (isCone[vertex]) hasCone = true;
            if (std::abs(curvature[vertex]) > std::abs(curvature[strongest])) strongest = vertex;
            oneRing(vertex, ring);
            for (int neighbor : ring) {
                if (boundary[neighbor] || seen[neighbor]) continue;
                seen[neighbor] = 1;
                stack.push_back(neighbor);
            }
        }
        if (!hasCone) {
            chosen.push_back(strongest);
            isCone[strongest] = 1;
        }
    }

    std::vector<int> transientOf(count, -1);
    std::vector<int> transient;
    transient.reserve(interiorCount);
    for (int index = 0; index < count; ++index) {
        if (boundary[index] || isCone[index]) continue;
        transientOf[index] = static_cast<int>(transient.size());
        transient.push_back(index);
    }

    bool solved = false;
    std::vector<double> target = curvature;
    if (!transient.empty() && !chosen.empty()) {
        std::vector<Eigen::Triplet<double>> entries;
        entries.reserve(transient.size() * 7);
        Eigen::VectorXd rhs(transient.size());
        std::vector<int> neighbors;
        for (int row = 0; row < static_cast<int>(transient.size()); ++row) {
            const int vertex = transient[row];
            oneRing(vertex, neighbors);
            rhs[row] = curvature[vertex];
            entries.emplace_back(row, row, 1.0);
            if (neighbors.empty()) continue;
            const double step = 1.0 / static_cast<double>(neighbors.size());
            for (int neighbor : neighbors) {
                const int column = transientOf[neighbor];
                if (column >= 0) entries.emplace_back(column, row, -step);
            }
        }
        Eigen::SparseMatrix<double> walk(transient.size(), transient.size());
        walk.setFromTriplets(entries.begin(), entries.end());
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(walk);
        if (solver.info() == Eigen::Success) {
            const Eigen::VectorXd potential = solver.solve(rhs);
            if (solver.info() == Eigen::Success) {
                std::unordered_map<int, double> absorbed;
                for (int row = 0; row < static_cast<int>(transient.size()); ++row) {
                    oneRing(transient[row], neighbors);
                    if (neighbors.empty()) continue;
                    const double share = potential[row] / static_cast<double>(neighbors.size());
                    for (int neighbor : neighbors) {
                        if (isCone[neighbor]) absorbed[neighbor] += share;
                    }
                }
                for (int index : chosen) target[index] = curvature[index] + absorbed[index];
                solved = true;
            }
        }
    } else if (!chosen.empty()) {
        solved = true;
    }

    gCones.clear();
    double prescribed = 0.0;
    for (int index : chosen) {
        const double angle = solved ? 2.0 * M_PI - target[index] : std::clamp(2.0 * M_PI - curvature[index], M_PI * 0.5, M_PI * 2.5);
        gCones[index] = angle;
        prescribed += 2.0 * M_PI - angle;
    }
    std::ostringstream text;
    text.setf(std::ios::fixed);
    text.precision(3);
    text << "Ben-Chen 2008：" << gCones.size() << " 个锥点，曲率和 " << prescribed
         << "，网格曲率和 " << curvatureSum << "。点「展开」才计算。";
    if (!solved) text << " 曲率重分配没有解出，锥角暂用局部亏缺。";
    gStatus = text.str();
}

bool negativeUv(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c) {
    return (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x()) < 0;
}

void rebuildGeometry() {
    gTriangles.clear();
    gUvTriangles.clear();
    gFlipTris.clear();
    gFlips3d.clear();
    gConeMarks.clear();
    gFlips = 0;
    std::vector<Eigen::Vector3d> distortion(gMesh.vertices.size(), Eigen::Vector3d(0.55, 0.62, 0.72));
    if (gHasResult && gTexture == Texture::Distortion) {
        std::vector<double> weight(gMesh.vertices.size(), 0);
        std::vector<Eigen::Vector3d> color(gMesh.vertices.size(), Eigen::Vector3d::Zero());
        for (FaceCIter face = gMesh.faces.begin(); face != gMesh.faces.end(); ++face) {
            if (face->isBoundary()) continue;
            HalfEdgeCIter he = face->he;
            VertexCIter verts[3] = {he->vertex, he->next->vertex, he->next->next->vertex};
            std::vector<Eigen::Vector3d> p = {verts[0]->position, verts[1]->position, verts[2]->position};
            std::vector<Eigen::Vector3d> q;
            for (int corner = 0; corner < 3; ++corner) {
                q.emplace_back(verts[corner]->uv.x(), verts[corner]->uv.y(), 0);
            }
            const Eigen::Vector3d faceColor = QuasiConformalError::color(QuasiConformalError::compute(p, q));
            for (int corner = 0; corner < 3; ++corner) {
                color[verts[corner]->index] += faceColor;
                weight[verts[corner]->index] += 1;
            }
        }
        for (size_t index = 0; index < color.size(); ++index) {
            if (weight[index] > 0) distortion[index] = color[index] / weight[index];
        }
    }
    auto push = [](std::vector<DrawVertex>& target, const Eigen::Vector3d& p, const Eigen::Vector3d& n,
                   const Eigen::Vector2d& uv, const Eigen::Vector3d& color) {
        target.push_back(DrawVertex{static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()),
                                 static_cast<float>(n.x()), static_cast<float>(n.y()), static_cast<float>(n.z()),
                                 static_cast<float>(uv.x()), static_cast<float>(uv.y()),
                                 static_cast<float>(color.x()), static_cast<float>(color.y()), static_cast<float>(color.z())});
    };
    for (FaceCIter face = gMesh.faces.begin(); face != gMesh.faces.end(); ++face) {
        if (face->isBoundary()) continue;
        HalfEdgeCIter he = face->he;
        VertexCIter verts[3] = {he->vertex, he->next->vertex, he->next->next->vertex};
        Eigen::Vector3d positions[3] = {verts[0]->position, verts[1]->position, verts[2]->position};
        Eigen::Vector3d normal = (positions[1] - positions[0]).cross(positions[2] - positions[0]);
        if (normal.norm() > 1e-12) normal.normalize();
        else normal = Eigen::Vector3d(0, 1, 0);
        const bool flipped = gHasResult && negativeUv(verts[0]->uv, verts[1]->uv, verts[2]->uv);
        if (flipped) ++gFlips;
        for (int corner = 0; corner < 3; ++corner) {
            const Eigen::Vector3d color = gTexture == Texture::Distortion && gHasResult
                ? distortion[verts[corner]->index] : Eigen::Vector3d(0.72, 0.78, 0.84);
            push(gTriangles, positions[corner], normal, gHasResult ? verts[corner]->uv : Eigen::Vector2d::Zero(), color);
            if (gHasResult) {
                push(gUvTriangles, Eigen::Vector3d(verts[corner]->uv.x(), verts[corner]->uv.y(), 0),
                     Eigen::Vector3d(0, 0, 1), verts[corner]->uv, color);
                if (flipped) {
                    push(gFlipTris, Eigen::Vector3d(verts[corner]->uv.x(), verts[corner]->uv.y(), 0.001),
                         Eigen::Vector3d(0, 0, 1), verts[corner]->uv, Eigen::Vector3d(0.95, 0.78, 0.15));
                    push(gFlips3d, positions[corner] + normal * 0.002, normal, verts[corner]->uv,
                         Eigen::Vector3d(0.95, 0.78, 0.15));
                }
            }
        }
    }
    for (const auto& cone : gCones) {
        if (cone.first < 0 || cone.first >= static_cast<int>(gMesh.vertices.size())) continue;
        const Eigen::Vector3d p = gMesh.vertices[cone.first].position;
        push(gConeMarks, p, Eigen::Vector3d(0, 1, 0), Eigen::Vector2d::Zero(), Eigen::Vector3d(0.85, 0.2, 0.15));
    }
}

void unfold() {
    if (gMesh.vertices.empty() || gSolving) return;
    gSolving = true;
    const auto started = std::chrono::steady_clock::now();
    std::string error;
    try {
        if (gMethod == Method::Incremental) {
            IncrementalFlattening solver(gMesh);
            if (!gCones.empty()) {
                std::vector<int> indices;
                std::vector<double> angles;
                indices.reserve(gCones.size());
                angles.reserve(gCones.size());
                for (const auto& cone : gCones) {
                    indices.push_back(cone.first);
                    angles.push_back(cone.second);
                }
                solver.setPrescribedCones(indices, angles);
            }
            solver.setProgressCallback(reportIncrementalProgress);
            gStatus = "Incremental Flattening · 摊平 0/24";
            reportIncrementalProgress(0, 0, 24);
            solver.parameterize();
            gCones.clear();
            const std::vector<int>& cones = solver.coneVertices();
            const Eigen::VectorXd& curvature = solver.coneCurvatures();
            for (size_t index = 0; index < cones.size(); ++index) {
                const double angle = index < static_cast<size_t>(curvature.size()) ? 2.0 * M_PI - curvature[static_cast<Eigen::Index>(index)] : kDefaultCone;
                gCones[cones[index]] = angle;
            }
        } else {
            std::vector<int> indices;
            std::vector<double> angles;
            indices.reserve(gCones.size());
            angles.reserve(gCones.size());
            for (const auto& cone : gCones) {
                indices.push_back(cone.first);
                angles.push_back(cone.second);
            }
            if (gMethod == Method::Circle) {
                CirclePatterns solver(gMesh, gScheme);
                solver.setProgressCallback(reportCircleProgress);
                gStatus = "Circle Pattern · 角度二次规划 0/8000";
                reportCircleProgress(0, 8000);
                solver.setConeSingulars(indices, angles);
                solver.parameterize();
            } else if (gMethod == Method::Cetm) {
                Cetm solver(gMesh, gScheme);
                solver.setConeSingulars(indices, angles);
                solver.parameterize();
            } else {
                RicciFlow solver(gMesh, gScheme);
                solver.setConeSingulars(indices, angles);
                solver.parameterize();
            }
        }
        int placed = 0;
        double minU = 1e9, maxU = -1e9, minV = 1e9, maxV = -1e9;
        for (VertexCIter vertex = gMesh.vertices.begin(); vertex != gMesh.vertices.end(); ++vertex) {
            if (!std::isfinite(vertex->uv.x()) || !std::isfinite(vertex->uv.y())) continue;
            ++placed;
            minU = std::min(minU, vertex->uv.x());
            maxU = std::max(maxU, vertex->uv.x());
            minV = std::min(minV, vertex->uv.y());
            maxV = std::max(maxV, vertex->uv.y());
        }
        gHasResult = placed > 2 && (maxU - minU) > 1e-6 && (maxV - minV) > 1e-6;
        if (!gHasResult) {
            error = gCones.empty()
                ? "展开没有得到坐标。先放锥点，让目标角和与欧拉示性数一致后再试。"
                : "展开没有得到坐标。";
        }
    } catch (const std::exception& exception) {
        error = exception.what();
        gHasResult = false;
    }
    gSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    rebuildGeometry();
    gStatus = error.empty()
        ? std::string(methodName(gMethod)) + " · " + schemeName(gScheme) + " · " +
          std::to_string(gSeconds).substr(0, 5) + " 秒 · 翻转 " + std::to_string(gFlips) +
          " · 锥点 " + std::to_string(gCones.size())
        : error;
    gSolving = false;
}

struct Panel {
    float x, y, w, h;
};

float toolbarHeight() { return gShowExtra ? 118.f : 86.f; }

Panel leftPanel() {
    const float top = toolbarHeight();
    const float gap = 8.f;
    // 结果优先：UV 面板比 3D 预览更宽。
    const float width = (static_cast<float>(gWidth) - 3.f * gap) * 0.44f;
    return {gap, top, width, static_cast<float>(gHeight) - top - gap};
}

Panel rightPanel() {
    const Panel left = leftPanel();
    const float width = static_cast<float>(gWidth) - left.x - left.w - 16.f;
    return {left.x + left.w + 8.f, left.y, width, left.h};
}

Mat4 modelViewProjection(const Panel& panel) {
    const Vec3 eye = eyePosition();
    const float aspect = std::max(panel.w, 1.f) / std::max(panel.h, 1.f);
    const float radius = std::max(gRadius, 0.2f);
    return multiply(perspective(0.9f, aspect, radius * 0.02f, radius * 8.f), lookAt(eye, gCenter));
}

void project(const Eigen::Vector3d& point, const Panel& panel, float& sx, float& sy) {
    const Mat4 mvp = modelViewProjection(panel);
    const float x = static_cast<float>(point.x()), y = static_cast<float>(point.y()), z = static_cast<float>(point.z());
    const float cx = mvp.m[0] * x + mvp.m[4] * y + mvp.m[8] * z + mvp.m[12];
    const float cy = mvp.m[1] * x + mvp.m[5] * y + mvp.m[9] * z + mvp.m[13];
    const float cw = mvp.m[3] * x + mvp.m[7] * y + mvp.m[11] * z + mvp.m[15];
    sx = panel.x + (cx / cw * 0.5f + 0.5f) * panel.w;
    sy = panel.y + (0.5f - cy / cw * 0.5f) * panel.h;
}

bool rayHit(double mouseX, double mouseY, int& vertexIndex) {
    const Panel panel = leftPanel();
    if (mouseX < panel.x || mouseY < panel.y || mouseX > panel.x + panel.w || mouseY > panel.y + panel.h) return false;
    const float nx = ((static_cast<float>(mouseX) - panel.x) / panel.w) * 2.f - 1.f;
    const float ny = 1.f - ((static_cast<float>(mouseY) - panel.y) / panel.h) * 2.f;
    const Mat4 mvp = modelViewProjection(panel);
    Eigen::Matrix4f matrix;
    for (int column = 0; column < 4; ++column) {
        for (int row = 0; row < 4; ++row) matrix(row, column) = mvp.m[column * 4 + row];
    }
    const Eigen::Matrix4f inverse = matrix.inverse();
    auto unproject = [&](float clipZ) {
        const Eigen::Vector4f world = inverse * Eigen::Vector4f(nx, ny, clipZ, 1.f);
        return Eigen::Vector3d(world.x() / world.w(), world.y() / world.w(), world.z() / world.w());
    };
    const Eigen::Vector3d origin = unproject(-1.f);
    Eigen::Vector3d direction = unproject(1.f) - origin;
    if (direction.norm() < 1e-8) return false;
    direction.normalize();
    double best = 1e18;
    for (FaceCIter face = gMesh.faces.begin(); face != gMesh.faces.end(); ++face) {
        if (face->isBoundary()) continue;
        HalfEdgeCIter he = face->he;
        const Eigen::Vector3d a = he->vertex->position, b = he->next->vertex->position, c = he->next->next->vertex->position;
        const Eigen::Vector3d edge1 = b - a, edge2 = c - a;
        const Eigen::Vector3d pvec = direction.cross(edge2);
        const double det = edge1.dot(pvec);
        if (std::abs(det) < 1e-10) continue;
        const double invDet = 1.0 / det;
        const Eigen::Vector3d tvec = origin - a;
        const double u = tvec.dot(pvec) * invDet;
        if (u < 0 || u > 1) continue;
        const Eigen::Vector3d qvec = tvec.cross(edge1);
        const double v = direction.dot(qvec) * invDet;
        if (v < 0 || u + v > 1) continue;
        const double t = edge2.dot(qvec) * invDet;
        if (t > 1e-5 && t < best) {
            best = t;
            const int ids[3] = {he->vertex->index, he->next->vertex->index, he->next->next->vertex->index};
            const Eigen::Vector3d corners[3] = {a, b, c};
            float bestScreen = 1e9f;
            vertexIndex = ids[0];
            for (int corner = 0; corner < 3; ++corner) {
                float sx = 0, sy = 0;
                project(corners[corner], panel, sx, sy);
                const float distance = std::hypot(sx - static_cast<float>(mouseX), sy - static_cast<float>(mouseY));
                if (distance < bestScreen) {
                    bestScreen = distance;
                    vertexIndex = ids[corner];
                }
            }
        }
    }
    return best < 1e17;
}

void toggleCone(int index) {
    const auto found = gCones.find(index);
    if (found == gCones.end()) gCones[index] = kDefaultCone;
    else gCones.erase(found);
    gStatus = "锥点 " + std::to_string(gCones.size()) + " 个，默认锥角 4π。点「展开」才计算。";
    rebuildGeometry();
}

void clickButton(int id) {
    if (gSolving) return;
    if (id == 1) loadModel(gModel - 1);
    else if (id == 2) loadModel(gModel + 1);
    else if (id == 10) gMethod = Method::Circle;
    else if (id == 11) gMethod = Method::Cetm;
    else if (id == 12) gMethod = Method::Ricci;
    else if (id == 13) gMethod = Method::Incremental;
    else if (id == 20) {
        if (gScheme == NEWTON) gScheme = GRAD_DESCENT;
        else if (gScheme == GRAD_DESCENT) gScheme = LBFGS;
        else gScheme = NEWTON;
    } else if (id == 30) {
        gPickCones = !gPickCones;
    } else if (id == 31) {
        placeBenChenCones();
    } else if (id == 40) {
        gUnfoldNow = true;
    } else if (id == 60) {
        gShowExtra = !gShowExtra;
    } else if (id == 41) {
        gHasResult = false;
        gSeconds = 0;
        gFlips = 0;
        gStatus = "已清除展开结果。锥点还在，点「展开」再算。";
    } else if (id == 50) gTexture = Texture::Checker;
    else if (id == 51) gTexture = Texture::Distortion;
    if (id != 40) rebuildGeometry();
}

GLuint compile(GLenum type, const char* source) {
    const GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);
    return shader;
}

void setupGL() {
    const char* vertexSource = R"(#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;
layout(location=2) in vec2 aUv;
layout(location=3) in vec3 aCol;
uniform mat4 uMvp;
out vec3 vNrm; out vec2 vUv; out vec3 vCol;
void main(){ gl_Position = uMvp * vec4(aPos,1.0); vNrm=aNrm; vUv=aUv; vCol=aCol; })";
    const char* fragmentSource = R"(#version 330 core
in vec3 vNrm; in vec2 vUv; in vec3 vCol; out vec4 frag;
uniform int uShade; uniform vec3 uLight;
void main(){
 vec3 n = normalize(vNrm);
 float light = 0.42 + 0.58*max(dot(n, normalize(uLight)), 0.0);
 vec3 color = vCol;
 if(uShade==0){
   float cell = step(0.5, fract(vUv.x*8.0)) == step(0.5, fract(vUv.y*8.0)) ? 0.92 : 0.28;
   color = vec3(cell);
 }
 if(uShade==3) color = vec3(0.95,0.78,0.15);
 frag = vec4(color * light, 1.0);
})";
    gProgram = glCreateProgram();
    glAttachShader(gProgram, compile(GL_VERTEX_SHADER, vertexSource));
    glAttachShader(gProgram, compile(GL_FRAGMENT_SHADER, fragmentSource));
    glLinkProgram(gProgram);
    glGenVertexArrays(1, &gVao);
    glGenBuffers(1, &gBuffer);
}

void drawMesh(const std::vector<DrawVertex>& vertices, const Mat4& mvp, int shade, float pointSize) {
    if (vertices.empty()) return;
    glUseProgram(gProgram);
    glUniformMatrix4fv(glGetUniformLocation(gProgram, "uMvp"), 1, GL_FALSE, mvp.m);
    glUniform1i(glGetUniformLocation(gProgram, "uShade"), shade);
    const Vec3 eye = eyePosition();
    Vec3 light{eye.x - gCenter.x, eye.y - gCenter.y, eye.z - gCenter.z};
    glUniform3f(glGetUniformLocation(gProgram, "uLight"), light.x, light.y, light.z);
    glBindVertexArray(gVao);
    glBindBuffer(GL_ARRAY_BUFFER, gBuffer);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(DrawVertex), vertices.data(), GL_DYNAMIC_DRAW);
    const GLsizei stride = sizeof(DrawVertex);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(12));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(24));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(32));
    if (pointSize > 0) {
        glPointSize(pointSize);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(vertices.size()));
    } else {
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices.size()));
    }
}

Mat4 uvMatrix(const Panel& panel) {
    float minU = 0, minV = 0, maxU = 1, maxV = 1;
    if (!gUvTriangles.empty()) {
        minU = maxU = gUvTriangles[0].px;
        minV = maxV = gUvTriangles[0].py;
        for (const DrawVertex& vertex : gUvTriangles) {
            minU = std::min(minU, vertex.px); maxU = std::max(maxU, vertex.px);
            minV = std::min(minV, vertex.py); maxV = std::max(maxV, vertex.py);
        }
    }
    const float cx = (minU + maxU) * 0.5f + gUvPanX;
    const float cy = (minV + maxV) * 0.5f + gUvPanY;
    const float span = std::max(maxU - minU, maxV - minV) * 0.62f / std::max(gUvZoom, 0.05f);
    const float aspect = panel.w / std::max(panel.h, 1.f);
    return ortho(cx - span * aspect, cx + span * aspect, cy - span, cy + span);
}

bool toolButton(const char* label, bool on) {
    if (on) ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.70f, 0.82f, 0.95f, 1.f));
    const bool pressed = ImGui::Button(label);
    if (on) ImGui::PopStyleColor();
    return pressed;
}

void drawInterface() {
    imguiBegin();
    ImGui::SetNextWindowPos(ImVec2(0.f, 0.f));
    ImGui::SetNextWindowSize(ImVec2(static_cast<float>(gWidth), toolbarHeight()));
    ImGui::Begin("工具", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
    if (ImGui::Button("<")) clickButton(1);
    ImGui::SameLine();
    ImGui::TextUnformatted(gModels.empty() ? "没有模型" : gModels[gModel].label.c_str());
    ImGui::SameLine();
    if (ImGui::Button(">")) clickButton(2);
    ImGui::SameLine();
    if (toolButton("Circle Pattern", gMethod == Method::Circle)) clickButton(10);
    ImGui::SameLine();
    if (toolButton("CETM", gMethod == Method::Cetm)) clickButton(11);
    ImGui::SameLine();
    if (toolButton("Ricci 流", gMethod == Method::Ricci)) clickButton(12);
    ImGui::SameLine();
    if (toolButton("Incremental Flattening", gMethod == Method::Incremental)) clickButton(13);
    ImGui::SameLine();
    if (ImGui::Button("展开")) clickButton(40);
    ImGui::SameLine();
    if (ImGui::Button(gShowExtra ? "收起" : "更多")) clickButton(60);
    if (gShowExtra) {
        if (toolButton(gPickCones ? "锥点点选中" : "锥奇异点", gPickCones)) clickButton(30);
        ImGui::SameLine();
        if (ImGui::Button("Ben-Chen 2008")) clickButton(31);
        ImGui::SameLine();
        if (ImGui::Button("清除")) clickButton(41);
        ImGui::SameLine();
        if (toolButton("棋盘格", gTexture == Texture::Checker)) clickButton(50);
        ImGui::SameLine();
        if (toolButton("共形畸变", gTexture == Texture::Distortion)) clickButton(51);
        ImGui::SameLine();
        if (ImGui::Button(schemeName(gScheme))) clickButton(20);
    }
    ImGui::TextWrapped("%s", gStatus.c_str());
    ImGui::End();
    const Panel left = leftPanel();
    const Panel right = rightPanel();
    ImDrawList* overlay = ImGui::GetForegroundDrawList();
    overlay->AddText(ImVec2(left.x + 10.f, left.y + 8.f), IM_COL32(30, 34, 40, 255), "3D");
    overlay->AddText(ImVec2(right.x + 10.f, right.y + 8.f), IM_COL32(30, 34, 40, 255), gHasResult ? "UV" : "UV · 点展开");
    imguiEnd();
}

void drawFrame() {
    glViewport(0, 0, static_cast<int>(gWidth * gScale), static_cast<int>(gHeight * gScale));
    glClearColor(0.93f, 0.94f, 0.95f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    const Panel left = leftPanel();
    const Panel right = rightPanel();
    auto viewportOf = [](const Panel& panel) {
        glViewport(static_cast<int>(panel.x * gScale), static_cast<int>((gHeight - panel.y - panel.h) * gScale),
                   static_cast<int>(panel.w * gScale), static_cast<int>(panel.h * gScale));
        glScissor(static_cast<int>(panel.x * gScale), static_cast<int>((gHeight - panel.y - panel.h) * gScale),
                  static_cast<int>(panel.w * gScale), static_cast<int>(panel.h * gScale));
    };
    glEnable(GL_SCISSOR_TEST);
    viewportOf(left);
    glClearColor(0.86f, 0.89f, 0.92f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    const int shade = gHasResult && gTexture == Texture::Checker ? 0 : 2;
    drawMesh(gTriangles, modelViewProjection(left), shade, 0);
    drawMesh(gFlips3d, modelViewProjection(left), 3, 0);
    drawMesh(gConeMarks, modelViewProjection(left), 2, 11.f * gScale);

    viewportOf(right);
    glClearColor(0.96f, 0.96f, 0.94f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if (gHasResult) {
        const Mat4 uv = uvMatrix(right);
        drawMesh(gUvTriangles, uv, shade, 0);
        glDisable(GL_DEPTH_TEST);
        drawMesh(gFlipTris, uv, 3, 0);
        glEnable(GL_DEPTH_TEST);
    }
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_DEPTH_TEST);
    glViewport(0, 0, static_cast<int>(gWidth * gScale), static_cast<int>(gHeight * gScale));
    drawInterface();
}

void reportCircleProgress(int completed, int total) {
    gStatus = "Circle Pattern · 角度二次规划 " + std::to_string(completed) + "/" + std::to_string(total);
    if (!gWindow) return;
    drawFrame();
    glfwSwapBuffers(gWindow);
    glfwPollEvents();
}

void reportIncrementalProgress(int phase, int completed, int total) {
    const char* name = "计算";
    if (phase == 0) name = "摊平";
    else if (phase == 1) name = "曲率取整";
    else if (phase == 2) name = "割缝";
    else if (phase == 3) name = "交叉场";
    else if (phase == 4) name = "ARAP";
    gStatus = std::string("Incremental Flattening · ") + name + " " + std::to_string(completed) + "/" + std::to_string(total);
    if (!gWindow) return;
    drawFrame();
    glfwSwapBuffers(gWindow);
    glfwPollEvents();
}

void mouseButton(GLFWwindow* window, int button, int action, int mods) {
    imguiMouseButton(window, button, action, mods);
    if (button != GLFW_MOUSE_BUTTON_LEFT) return;
    if (action == GLFW_RELEASE) { gDrag = false; return; }
    if (imguiCapturesMouse()) return;
    glfwGetCursorPos(window, &gLastX, &gLastY);
    int vertex = -1;
    if (gPickCones && rayHit(gLastX, gLastY, vertex)) {
        toggleCone(vertex);
        return;
    }
    const Panel left = leftPanel();
    const Panel right = rightPanel();
    if (gLastX >= left.x && gLastX <= left.x + left.w && gLastY >= left.y && gLastY <= left.y + left.h) {
        gDrag = true; gDragPanel = 0;
    } else if (gLastX >= right.x && gLastX <= right.x + right.w && gLastY >= right.y && gLastY <= right.y + right.h) {
        gDrag = true; gDragPanel = 1;
    }
}

void cursor(GLFWwindow*, double x, double y) {
    if (!gDrag) return;
    const double dx = x - gLastX, dy = y - gLastY;
    gLastX = x; gLastY = y;
    if (gDragPanel == 0) {
        gYaw += static_cast<float>(dx) * 0.008f;
        gPitch = std::clamp(gPitch + static_cast<float>(dy) * 0.008f, -1.2f, 1.2f);
    } else {
        gUvPanX -= static_cast<float>(dx) * 0.002f / std::max(gUvZoom, 0.05f);
        gUvPanY += static_cast<float>(dy) * 0.002f / std::max(gUvZoom, 0.05f);
    }
}

void scroll(GLFWwindow* window, double, double yoff) {
    imguiScroll(window, 0, yoff);
    if (imguiCapturesMouse()) return;
    double x = 0, y = 0;
    glfwGetCursorPos(glfwGetCurrentContext(), &x, &y);
    const Panel right = rightPanel();
    if (x >= right.x && x <= right.x + right.w && y >= right.y && y <= right.y + right.h) {
        gUvZoom = std::clamp(gUvZoom * (yoff > 0 ? 1.1f : 0.9f), 0.2f, 40.f);
    } else {
        gRadius = std::clamp(gRadius * (yoff > 0 ? 0.9f : 1.1f), 0.2f, 400.f);
    }
}

} // namespace

int main() {
    gModels = listModels();
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
    GLFWwindow* window = glfwCreateWindow(gWidth, gHeight, "锥奇异点展开", nullptr, nullptr);
    if (!window) return 1;
    gWindow = window;
    glfwMakeContextCurrent(window);
    glfwSetMouseButtonCallback(window, mouseButton);
    glfwSetCursorPosCallback(window, cursor);
    glfwSetScrollCallback(window, scroll);
    glfwSetKeyCallback(window, imguiKey);
    glfwSetCharCallback(window, imguiChar);
    glfwSwapInterval(1);
    setupGL();
    imguiInit(window);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    if (!gModels.empty()) loadModel(0);
    else gStatus = "没有找到模型。";
    while (!glfwWindowShouldClose(window)) {
        int framebufferWidth = 0, framebufferHeight = 0;
        glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);
        glfwGetWindowSize(window, &gWidth, &gHeight);
        gScale = gWidth > 0 ? static_cast<float>(framebufferWidth) / static_cast<float>(gWidth) : 1.f;
        drawFrame();
        glfwSwapBuffers(window);
        glfwPollEvents();
        if (gAutoRun) {
            gAutoRun = false;
            if (!gModels.empty() && !gSolving) {
                placeBenChenCones();
                unfold();
            }
        }
        if (gUnfoldNow && !gSolving) {
            gUnfoldNow = false;
            unfold();
        }
    }
    imguiShutdown();
    glfwDestroyWindow(window);
    gWindow = nullptr;
    glfwTerminate();
    return 0;
}
