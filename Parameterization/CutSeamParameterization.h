#ifndef CUT_SEAM_PARAMETERIZATION_H
#define CUT_SEAM_PARAMETERIZATION_H

#include "Parameterization.h"
#include "CutSeamMesh.h"
#include <utility>
#include <vector>

/**
 * Extended base class for parameterization algorithms that operate on
 * cut-seam meshes (meshes with artificial boundaries created by cutting
 * along homology / seam edges).
 *
 * The CutSeamMesh duplicates vertices along cut edges, turning closed
 * meshes into open meshes with boundary. The original mesh is available
 * via the base class Parameterization::mesh.
 *
 * Typical pipeline:
 *   1. Construct the subclass with the original mesh.
 *   2. Call setSeamEdges() or setSeamVertexPairs() to define the cut.
 *   3. Call parameterize() — subclass uses cutMesh_ internally.
 *   4. Results are written back to mesh.vertices[i].uv.
 *
 * Algorithms supported:
 *   - AbelJacobiParameterization  — genus-aware global conformal
 *   - HolomorphicOneForm           — Hodge / cotangent harmonic 1-form
 *   - BFF                          — Boundary First Flattening
 */
class CutSeamParameterization : public Parameterization {
public:
    // constructor
    CutSeamParameterization(Mesh& mesh0);

    // destructor
    virtual ~CutSeamParameterization() {}

    /**
     * Define cut seams by edge indices in the original mesh.
     * @param seamEdgeIndices  indices of edges to cut along
     */
    virtual void setSeamEdges(const std::vector<int>& seamEdgeIndices);

    /**
     * Define cut seams by vertex pair (each pair defines a seam edge).
     * The vertices must be connected by an edge in the original mesh.
     * @param seamVertexPairs  pairs (v_i, v_j) of vertex indices
     */
    virtual void setSeamVertexPairs(const std::vector<std::pair<int, int>>& seamVertexPairs);

    /** Access the cut mesh for reading or modification. */
    CutSeamMesh& cutMesh() { return cutMesh_; }

    /** Const access to the cut mesh. */
    const CutSeamMesh& cutMesh() const { return cutMesh_; }

    /** Whether cut seams have been defined. */
    bool hasSeams() const { return hasSeams_; }

    /** Clear seams and reset the cut mesh to a copy of the original mesh. */
    void clearSeams();

protected:
    /** Ensure cutMesh_ is ready. Closed meshes are cut with tree-cotree seams. */
    bool prepareCutMesh();

    /** Copy/average UVs from duplicated cut vertices back to original mesh vertices. */
    void copyCutMeshUvsToOriginal();

    // cut mesh: a derived Mesh with duplicated vertices along seams
    CutSeamMesh cutMesh_;

    // true if setSeamEdges() or setSeamVertexPairs() was called with valid data
    bool hasSeams_;
};

#endif
