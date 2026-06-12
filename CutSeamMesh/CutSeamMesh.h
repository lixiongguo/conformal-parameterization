#ifndef CUT_SEAM_MESH_H
#define CUT_SEAM_MESH_H

#include "../BaseMesh/Mesh.h"
#include <utility>
#include <vector>

class CutSeamMesh : public Mesh {
public:
    CutSeamMesh() = default;
    explicit CutSeamMesh(const Mesh& source);

    bool buildFromSeamEdges(const Mesh& source,
                            const std::vector<int>& seamEdgeIndices);
    bool buildFromSeamVertexPairs(const Mesh& source,
                                  const std::vector<std::pair<int, int>>& seamVertexPairs);

    int originalVertexIndex(int cutVertexIndex) const;
    int originalFaceIndex(int cutFaceIndex) const;
    bool isCutBoundaryEdge(int cutEdgeIndex) const;
    bool isCutBoundaryEdge(EdgeCIter edge) const;

    const std::vector<int>& originalVertexForVertex() const;
    const std::vector<int>& originalFaceForFace() const;
    const std::vector<bool>& cutBoundaryEdges() const;
    const std::vector<std::pair<int, int>>& seamVertexPairs() const;

private:
    bool buildFromSeamPairSet(const Mesh& source,
                              const std::vector<std::pair<int, int>>& seamPairs);
    void clearMetadata();
    void markCutBoundaryEdges(const std::vector<std::pair<int, int>>& seamPairs);

    static std::pair<int, int> canonicalPair(int a, int b);

    std::vector<int> originalVertexForVertex_;
    std::vector<int> originalFaceForFace_;
    std::vector<bool> cutBoundaryEdges_;
    std::vector<std::pair<int, int>> seamVertexPairs_;
};

#endif
