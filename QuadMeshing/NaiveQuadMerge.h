#ifndef NAIVE_QUAD_MERGE_H
#define NAIVE_QUAD_MERGE_H

#include <vector>

#include "StageLog.h"

namespace NaiveQuadMerge {

/**
 * Greedy pairing of a triangle mesh into quads.
 *
 * Each triangle is used at most once. For an unused triangle, the partner is the
 * unused face across its longest shared edge. That edge is the deleted diagonal,
 * which is what puts a diagonal-split quad back together. Vertices are not moved.
 *
 * A triangle (a → b → c) whose shared edge is a → b, with the partner's opposite
 * vertex d, becomes the boundary a → d → b → c.
 *
 * quads holds 4 indices per face. leftovers holds triangles that had no free partner.
 */
struct Result {
    std::vector<int> quads;
    std::vector<int> leftovers;
    std::vector<StageTime> stages;
};

Result merge(const double* positions, int nVerts, const int* faces, int nFaces);

} // namespace NaiveQuadMerge

#endif
