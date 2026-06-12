#include "CutSeamParameterization.h"
#include "TreeCotreeBasis.h"

CutSeamParameterization::CutSeamParameterization(Mesh& mesh0)
    : Parameterization(mesh0)
    , cutMesh_(mesh0)
    , hasSeams_(false)
{
}

void CutSeamParameterization::setSeamEdges(const std::vector<int>& seamEdgeIndices)
{
    hasSeams_ = cutMesh_.buildFromSeamEdges(mesh, seamEdgeIndices);
}

void CutSeamParameterization::setSeamVertexPairs(const std::vector<std::pair<int, int>>& seamVertexPairs)
{
    hasSeams_ = cutMesh_.buildFromSeamVertexPairs(mesh, seamVertexPairs);
}

void CutSeamParameterization::clearSeams()
{
    cutMesh_ = CutSeamMesh(mesh);
    hasSeams_ = false;
}

bool CutSeamParameterization::prepareCutMesh()
{
    if (hasSeams_) {
        return true;
    }

    if (mesh.boundaries.empty()) {
        hasSeams_ = topology::TreeCotreeBasis::buildCutSeamMesh(mesh, cutMesh_);
        return hasSeams_;
    }

    cutMesh_ = CutSeamMesh(mesh);
    return !cutMesh_.vertices.empty();
}

void CutSeamParameterization::copyCutMeshUvsToOriginal()
{
    std::vector<Eigen::Vector2d> uvSums(mesh.vertices.size(), Eigen::Vector2d::Zero());
    std::vector<int> uvCounts(mesh.vertices.size(), 0);

    for (VertexCIter cutVertex = cutMesh_.vertices.begin();
         cutVertex != cutMesh_.vertices.end();
         ++cutVertex) {
        const int originalIndex = cutMesh_.originalVertexIndex(cutVertex->index);
        if (originalIndex < 0 || originalIndex >= static_cast<int>(mesh.vertices.size())) {
            continue;
        }

        uvSums[originalIndex] += cutVertex->uv;
        uvCounts[originalIndex] += 1;
    }

    for (VertexIter vertex = mesh.vertices.begin(); vertex != mesh.vertices.end(); ++vertex) {
        if (vertex->index < 0 || vertex->index >= static_cast<int>(uvCounts.size())) {
            continue;
        }

        if (uvCounts[vertex->index] > 0) {
            vertex->uv = uvSums[vertex->index] / static_cast<double>(uvCounts[vertex->index]);
        }
    }
}
