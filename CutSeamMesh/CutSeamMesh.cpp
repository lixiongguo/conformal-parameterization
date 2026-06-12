#include "CutSeamMesh.h"
#include "../BaseMesh/MeshIO.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <unordered_map>

CutSeamMesh::CutSeamMesh(const Mesh& source)
{
    buildFromSeamEdges(source, std::vector<int>());
}

bool CutSeamMesh::buildFromSeamEdges(const Mesh& source,
                                     const std::vector<int>& seamEdgeIndices)
{
    std::vector<std::pair<int, int>> seamPairs;
    seamPairs.reserve(seamEdgeIndices.size());

    for (int edgeIndex : seamEdgeIndices) {
        if (edgeIndex < 0 || edgeIndex >= static_cast<int>(source.edges.size())) {
            std::cerr << "Error: seam edge index " << edgeIndex
                      << " is outside the source mesh edge range." << std::endl;
            return false;
        }

        EdgeCIter edge = source.edges.begin() + edgeIndex;
        int a = edge->he->vertex->index;
        int b = edge->he->flip->vertex->index;
        seamPairs.push_back(canonicalPair(a, b));
    }

    return buildFromSeamPairSet(source, seamPairs);
}

bool CutSeamMesh::buildFromSeamVertexPairs(
    const Mesh& source,
    const std::vector<std::pair<int, int>>& seamVertexPairs)
{
    std::vector<std::pair<int, int>> seamPairs;
    seamPairs.reserve(seamVertexPairs.size());

    for (const std::pair<int, int>& pair : seamVertexPairs) {
        if (pair.first < 0 || pair.second < 0 ||
            pair.first >= static_cast<int>(source.vertices.size()) ||
            pair.second >= static_cast<int>(source.vertices.size())) {
            std::cerr << "Error: seam vertex pair (" << pair.first << ", "
                      << pair.second << ") is outside the source mesh vertex range."
                      << std::endl;
            return false;
        }

        seamPairs.push_back(canonicalPair(pair.first, pair.second));
    }

    return buildFromSeamPairSet(source, seamPairs);
}

bool CutSeamMesh::buildFromSeamPairSet(
    const Mesh& source,
    const std::vector<std::pair<int, int>>& seamPairs)
{
    clearMetadata();

    std::set<std::pair<int, int>> seamSet;
    for (const std::pair<int, int>& pair : seamPairs) {
        seamSet.insert(canonicalPair(pair.first, pair.second));
    }
    seamVertexPairs_.assign(seamSet.begin(), seamSet.end());

    std::set<std::pair<int, int>> sourceEdges;
    for (EdgeCIter edge = source.edges.begin(); edge != source.edges.end(); ++edge) {
        int a = edge->he->vertex->index;
        int b = edge->he->flip->vertex->index;
        sourceEdges.insert(canonicalPair(a, b));
    }

    for (const std::pair<int, int>& pair : seamVertexPairs_) {
        if (sourceEdges.find(pair) == sourceEdges.end()) {
            std::cerr << "Warning: seam vertex pair (" << pair.first << ", "
                      << pair.second << ") is not an edge in the source mesh."
                      << std::endl;
        }
    }

    std::vector<std::vector<int>> cornersByVertex(source.vertices.size());
    std::vector<std::vector<int>> faceCorners;
    std::vector<int> sourceFaceIndices;

    for (FaceCIter face = source.faces.begin(); face != source.faces.end(); ++face) {
        if (face->isBoundary()) {
            continue;
        }

        std::vector<int> corners;
        HalfEdgeCIter he = face->he;
        do {
            corners.push_back(he->index);
            cornersByVertex[he->vertex->index].push_back(he->index);
            he = he->next;
        } while (he != face->he);

        faceCorners.push_back(corners);
        sourceFaceIndices.push_back(face->index);
    }

    MeshData cutData;
    cutData.positions.reserve(source.vertices.size());

    std::vector<Eigen::Vector2d> duplicateUvs;
    std::vector<int> duplicateForHalfEdge(source.halfEdges.size(), -1);

    for (size_t vertexIndex = 0; vertexIndex < cornersByVertex.size(); ++vertexIndex) {
        const std::vector<int>& corners = cornersByVertex[vertexIndex];
        if (corners.empty()) {
            continue;
        }

        std::unordered_map<int, int> localIndex;
        localIndex.reserve(corners.size());
        for (size_t i = 0; i < corners.size(); ++i) {
            localIndex[corners[i]] = static_cast<int>(i);
        }

        std::vector<std::vector<int>> adjacency(corners.size());
        for (int cornerIndex : corners) {
            HalfEdgeCIter he = source.halfEdges.begin() + cornerIndex;
            const std::pair<int, int> edgePair =
                canonicalPair(he->vertex->index, he->flip->vertex->index);

            if (he->edge->isBoundary() || seamSet.find(edgePair) != seamSet.end()) {
                continue;
            }

            HalfEdgeCIter oppositeCorner = he->flip->next;
            if (oppositeCorner->face->isBoundary()) {
                continue;
            }

            std::unordered_map<int, int>::const_iterator neighborIt =
                localIndex.find(oppositeCorner->index);
            if (neighborIt == localIndex.end()) {
                continue;
            }

            int a = localIndex[cornerIndex];
            int b = neighborIt->second;
            adjacency[a].push_back(b);
            adjacency[b].push_back(a);
        }

        std::vector<bool> visited(corners.size(), false);
        for (size_t seed = 0; seed < corners.size(); ++seed) {
            if (visited[seed]) {
                continue;
            }

            int duplicateIndex = static_cast<int>(cutData.positions.size());
            VertexCIter sourceVertex = source.vertices.begin() + vertexIndex;
            cutData.positions.push_back(sourceVertex->position);
            duplicateUvs.push_back(sourceVertex->uv);
            originalVertexForVertex_.push_back(static_cast<int>(vertexIndex));

            std::vector<int> stack;
            stack.push_back(static_cast<int>(seed));
            visited[seed] = true;

            while (!stack.empty()) {
                int current = stack.back();
                stack.pop_back();

                duplicateForHalfEdge[corners[current]] = duplicateIndex;

                for (int neighbor : adjacency[current]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        stack.push_back(neighbor);
                    }
                }
            }
        }
    }

    for (const std::vector<int>& corners : faceCorners) {
        std::vector<Index> faceIndices;
        faceIndices.reserve(corners.size());

        for (int cornerIndex : corners) {
            int duplicateIndex = duplicateForHalfEdge[cornerIndex];
            if (duplicateIndex < 0) {
                std::cerr << "Error: failed to assign a cut vertex for halfedge "
                          << cornerIndex << "." << std::endl;
                clearMetadata();
                return false;
            }

            faceIndices.push_back(Index(duplicateIndex, duplicateIndex, duplicateIndex));
        }

        cutData.indices.push_back(faceIndices);
    }

    if (!MeshIO::buildMesh(cutData, *this)) {
        clearMetadata();
        return false;
    }

    for (size_t i = 0; i < duplicateUvs.size() && i < vertices.size(); ++i) {
        vertices[i].uv = duplicateUvs[i];
    }

    originalFaceForFace_.assign(faces.size(), -1);
    for (size_t i = 0; i < sourceFaceIndices.size() && i < originalFaceForFace_.size(); ++i) {
        originalFaceForFace_[i] = sourceFaceIndices[i];
    }

    markCutBoundaryEdges(seamVertexPairs_);
    return true;
}

int CutSeamMesh::originalVertexIndex(int cutVertexIndex) const
{
    if (cutVertexIndex < 0 ||
        cutVertexIndex >= static_cast<int>(originalVertexForVertex_.size())) {
        return -1;
    }

    return originalVertexForVertex_[cutVertexIndex];
}

int CutSeamMesh::originalFaceIndex(int cutFaceIndex) const
{
    if (cutFaceIndex < 0 ||
        cutFaceIndex >= static_cast<int>(originalFaceForFace_.size())) {
        return -1;
    }

    return originalFaceForFace_[cutFaceIndex];
}

bool CutSeamMesh::isCutBoundaryEdge(int cutEdgeIndex) const
{
    if (cutEdgeIndex < 0 ||
        cutEdgeIndex >= static_cast<int>(cutBoundaryEdges_.size())) {
        return false;
    }

    return cutBoundaryEdges_[cutEdgeIndex];
}

bool CutSeamMesh::isCutBoundaryEdge(EdgeCIter edge) const
{
    return isCutBoundaryEdge(edge->index);
}

const std::vector<int>& CutSeamMesh::originalVertexForVertex() const
{
    return originalVertexForVertex_;
}

const std::vector<int>& CutSeamMesh::originalFaceForFace() const
{
    return originalFaceForFace_;
}

const std::vector<bool>& CutSeamMesh::cutBoundaryEdges() const
{
    return cutBoundaryEdges_;
}

const std::vector<std::pair<int, int>>& CutSeamMesh::seamVertexPairs() const
{
    return seamVertexPairs_;
}

void CutSeamMesh::clearMetadata()
{
    halfEdges.clear();
    vertices.clear();
    edges.clear();
    faces.clear();
    boundaries.clear();

    originalVertexForVertex_.clear();
    originalFaceForFace_.clear();
    cutBoundaryEdges_.clear();
    seamVertexPairs_.clear();
}

void CutSeamMesh::markCutBoundaryEdges(
    const std::vector<std::pair<int, int>>& seamPairs)
{
    std::set<std::pair<int, int>> seamSet;
    for (const std::pair<int, int>& pair : seamPairs) {
        seamSet.insert(canonicalPair(pair.first, pair.second));
    }

    cutBoundaryEdges_.assign(edges.size(), false);
    for (EdgeCIter edge = edges.begin(); edge != edges.end(); ++edge) {
        if (!edge->isBoundary()) {
            continue;
        }

        int a = originalVertexIndex(edge->he->vertex->index);
        int b = originalVertexIndex(edge->he->flip->vertex->index);
        if (seamSet.find(canonicalPair(a, b)) != seamSet.end()) {
            cutBoundaryEdges_[edge->index] = true;
        }
    }
}

std::pair<int, int> CutSeamMesh::canonicalPair(int a, int b)
{
    if (a > b) {
        std::swap(a, b);
    }

    return std::pair<int, int>(a, b);
}
