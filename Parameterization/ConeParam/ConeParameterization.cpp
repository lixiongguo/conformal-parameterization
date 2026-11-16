#include "ConeParameterization.h"
#include <algorithm>
#include <cmath>
#include <stack>

ConeParameterization::ConeParameterization(Mesh& mesh0)
    : Parameterization(mesh0)
{
}

void ConeParameterization::setConeSingulars(const std::vector<int>& coneIdx,
                                            const std::vector<double>& coneAngles)
{
    coneSingulars.clear();
    for (size_t i = 0; i < coneIdx.size() && i < coneAngles.size(); i++) {
        coneSingulars[coneIdx[i]] = coneAngles[i];
    }
}

bool ConeParameterization::isConeVertex(int vertexIndex) const
{
    return coneSingulars.find(vertexIndex) != coneSingulars.end();
}

double ConeParameterization::getConeTargetAngle(int vertexIndex) const
{
    auto it = coneSingulars.find(vertexIndex);
    return (it != coneSingulars.end()) ? it->second : 0.0;
}

const std::unordered_map<int, double>& ConeParameterization::getConeSingulars() const
{
    return coneSingulars;
}

void ConeParameterization::clearConeSingulars()
{
    coneSingulars.clear();
}

void ConeParameterization::layoutIntrinsicMetric(const Eigen::VectorXd& edgeLengths)
{
    const int nV = static_cast<int>(mesh.vertices.size());
    const int nF = static_cast<int>(mesh.faces.size());
    if (nV == 0 || nF == 0) return;

    std::vector<char> placed(static_cast<size_t>(nV), 0);
    std::vector<char> faceDone(static_cast<size_t>(nF), 0);
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        v->uv.setZero();
    }

    auto lenAt = [&](int edgeIndex) -> double {
        if (edgeIndex < 0 || edgeIndex >= edgeLengths.size()) return 1e-8;
        const double length = edgeLengths[edgeIndex];
        if (!std::isfinite(length) || length < 1e-12) return 1e-12;
        return length;
    };

    double cursor = 0.0;
    for (FaceCIter face = mesh.faces.begin(); face != mesh.faces.end(); ++face) {
        if (face->isBoundary()) continue;
        const int faceIndex = face->index;
        if (faceIndex < 0 || faceIndex >= nF || faceDone[static_cast<size_t>(faceIndex)]) continue;

        HalfEdgeCIter seed = face->he;
        const int ia = seed->vertex->index;
        const int ib = seed->next->vertex->index;
        if (ia < 0 || ib < 0 || ia >= nV || ib >= nV) continue;
        if (placed[static_cast<size_t>(ia)] || placed[static_cast<size_t>(ib)]) continue;

        const double seedLen = lenAt(seed->edge->index);
        seed->vertex->uv = Eigen::Vector2d::Zero();
        seed->next->vertex->uv = Eigen::Vector2d(seedLen, 0.0);
        placed[static_cast<size_t>(ia)] = 1;
        placed[static_cast<size_t>(ib)] = 1;

        std::vector<int> component;
        component.push_back(ia);
        component.push_back(ib);

        std::stack<HalfEdgeCIter> stack;
        stack.push(seed);
        stack.push(seed->flip);
        while (!stack.empty()) {
            const HalfEdgeCIter he = stack.top();
            stack.pop();
            if (he->onBoundary) continue;

            const int fIdx = he->face->index;
            if (fIdx < 0 || fIdx >= nF || faceDone[static_cast<size_t>(fIdx)]) continue;

            const HalfEdgeCIter next = he->next;
            const HalfEdgeCIter prev = next->next;
            const int a = he->vertex->index;
            const int b = next->vertex->index;
            const int c = prev->vertex->index;
            if (a < 0 || b < 0 || c < 0 || a >= nV || b >= nV || c >= nV) continue;
            if (!placed[static_cast<size_t>(a)] || !placed[static_cast<size_t>(b)]) continue;

            faceDone[static_cast<size_t>(fIdx)] = 1;
            if (!placed[static_cast<size_t>(c)]) {
                const Eigen::Vector2d A = he->vertex->uv;
                const Eigen::Vector2d B = next->vertex->uv;
                const Eigen::Vector2d AB = B - A;
                const double lab = AB.norm();
                if (lab < 1e-14) continue;
                const double lac = lenAt(prev->edge->index);
                const double lbc = lenAt(next->edge->index);
                const double x = (lab * lab + lac * lac - lbc * lbc) / (2.0 * lab);
                const double h2 = lac * lac - x * x;
                const double h = h2 > 0.0 ? std::sqrt(h2) : 0.0;
                const Eigen::Vector2d dir = AB / lab;
                const Eigen::Vector2d left(-dir.y(), dir.x());
                prev->vertex->uv = A + dir * x + left * h;
                placed[static_cast<size_t>(c)] = 1;
                component.push_back(c);
            }

            stack.push(next->flip);
            stack.push(prev->flip);
        }

        double minX = 1e300;
        double maxX = -1e300;
        for (int id : component) {
            const double x = mesh.vertices[id].uv.x();
            minX = std::min(minX, x);
            maxX = std::max(maxX, x);
        }
        if (component.empty() || !(minX < 1e299)) continue;
        const double shift = cursor - minX;
        for (int id : component) mesh.vertices[id].uv.x() += shift;
        cursor += std::max(maxX - minX, seedLen) * 1.1;
    }

    Eigen::Vector2d anchor = Eigen::Vector2d::Zero();
    int placedCount = 0;
    for (int i = 0; i < nV; ++i) {
        if (!placed[static_cast<size_t>(i)]) continue;
        anchor += mesh.vertices[i].uv;
        ++placedCount;
    }
    if (placedCount > 0) anchor /= static_cast<double>(placedCount);
    for (int i = 0; i < nV; ++i) {
        if (!placed[static_cast<size_t>(i)]) mesh.vertices[i].uv = anchor;
    }
}
