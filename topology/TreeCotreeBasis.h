#ifndef TREE_COTREE_BASIS_H
#define TREE_COTREE_BASIS_H

#include <vector>
#include <utility>

namespace topology {

/**
 * Build a tree–cotree homology basis on a closed, connected triangle mesh.
 *
 * Input mesh is given by:
 *  - nV vertices, nF faces, nE edges
 *  - each edge e has endpoints (v0,v1) with v0 < v1
 *  - each edge e has adjacent faces (f0,f1) with f0,f1 in [0,nF),
 *    and for a closed mesh every edge must have two incident faces.
 *
 * Output:
 *  - a set of edge cycles. Each cycle is a list of (edgeIndex, sign) pairs,
 *    representing an oriented 1-cycle in the primal graph:
 *      sum_i sign_i * omega[edgeIndex_i] = 0
 *
 * For a closed genus g surface, the number of returned cycles should be 2g.
 */
class TreeCotreeBasis {
public:
    struct Edge {
        int v0 = 0;
        int v1 = 0;
        int f0 = -1;
        int f1 = -1;
    };

    using SignedEdge = std::pair<int, int>; // (edgeIdx, sign in {+1,-1})
    using Cycle = std::vector<SignedEdge>;

    /**
     * Build a tree–cotree homology basis.
     *
     * @param nV number of vertices
     * @param nF number of faces (interior faces)
     * @param edges edge list size nE, each edge has v0<v1 and f0,f1>=0
     * @param genus genus g (expects g>=0)
     * @param findEdge callback returning edge index for (a,b), or -1 if missing
     * @param edgeSign callback returning +1 if (from->to) matches stored (v0->v1),
     *        -1 if opposite, 0 if invalid
     */
    template <typename FindEdgeFn, typename EdgeSignFn>
    static std::vector<Cycle> buildClosedMeshBasis(
        int nV,
        int nF,
        const std::vector<Edge>& edges,
        int genus,
        FindEdgeFn findEdge,
        EdgeSignFn edgeSign);
};

} // namespace topology

// Implementation (header-only template)
namespace topology {

namespace detail {

inline void buildVertexAdj(int nV, const std::vector<TreeCotreeBasis::Edge>& edges,
                           std::vector<std::vector<int>>& vEdges) {
    vEdges.assign(nV, {});
    for (int ei = 0; ei < (int)edges.size(); ++ei) {
        const auto& e = edges[ei];
        if (e.v0 >= 0 && e.v0 < nV) vEdges[e.v0].push_back(ei);
        if (e.v1 >= 0 && e.v1 < nV) vEdges[e.v1].push_back(ei);
    }
}

inline void buildFaceAdj(int nF, const std::vector<TreeCotreeBasis::Edge>& edges,
                         std::vector<std::vector<std::pair<int,int>>>& fAdj) {
    fAdj.assign(nF, {});
    for (int ei = 0; ei < (int)edges.size(); ++ei) {
        const auto& e = edges[ei];
        if (e.f0 < 0 || e.f1 < 0) continue;
        if (e.f0 >= nF || e.f1 >= nF) continue;
        fAdj[e.f0].push_back({e.f1, ei});
        fAdj[e.f1].push_back({e.f0, ei});
    }
}

inline bool pathToRoot(int v, const std::vector<int>& parent, std::vector<int>& out) {
    out.clear();
    int cur = v;
    std::vector<char> seen(parent.size(), 0);
    while (cur >= 0 && cur < (int)parent.size() && !seen[cur]) {
        seen[cur] = 1;
        out.push_back(cur);
        if (parent[cur] == cur) break;
        cur = parent[cur];
    }
    return !out.empty();
}

} // namespace detail

template <typename FindEdgeFn, typename EdgeSignFn>
std::vector<typename TreeCotreeBasis::Cycle> TreeCotreeBasis::buildClosedMeshBasis(
    int nV,
    int nF,
    const std::vector<Edge>& edges,
    int genus,
    FindEdgeFn findEdge,
    EdgeSignFn edgeSign) {

    std::vector<Cycle> cycles;
    if (genus <= 0) return cycles;
    if (nV <= 0 || nF <= 0 || edges.empty()) return cycles;

    const int nE = (int)edges.size();

    // 1) Build primal spanning tree T on vertex graph
    std::vector<std::vector<int>> vEdges;
    detail::buildVertexAdj(nV, edges, vEdges);

    std::vector<int> parentV(nV, -1);
    std::vector<int> parentEdge(nV, -1);
    std::vector<int> depth(nV, 0);
    std::vector<char> visitedV(nV, 0);
    std::vector<char> inPrimalTree(nE, 0);

    std::vector<int> q;
    q.reserve(nV);
    parentV[0] = 0;
    visitedV[0] = 1;
    q.push_back(0);

    for (size_t qi = 0; qi < q.size(); ++qi) {
        const int u = q[qi];
        for (int ei : vEdges[u]) {
            const auto& e = edges[ei];
            const int v = (e.v0 == u) ? e.v1 : e.v0;
            if (v < 0 || v >= nV) continue;
            if (visitedV[v]) continue;
            visitedV[v] = 1;
            parentV[v] = u;
            parentEdge[v] = ei;
            depth[v] = depth[u] + 1;
            inPrimalTree[ei] = 1;
            q.push_back(v);
        }
    }

    // 2) Build dual spanning tree T* on faces, forbidding dual edges whose primal is in T
    std::vector<std::vector<std::pair<int,int>>> fAdj;
    detail::buildFaceAdj(nF, edges, fAdj);

    std::vector<int> parentF(nF, -1);
    std::vector<int> parentDualEdge(nF, -1); // primal edge index used to connect in dual tree
    std::vector<char> visitedF(nF, 0);
    std::vector<char> inDualTree(nE, 0);

    std::vector<int> fq;
    fq.reserve(nF);
    parentF[0] = 0;
    visitedF[0] = 1;
    fq.push_back(0);

    for (size_t qi = 0; qi < fq.size(); ++qi) {
        const int f = fq[qi];
        for (const auto& nb : fAdj[f]) {
            const int g = nb.first;
            const int ei = nb.second;
            if (g < 0 || g >= nF) continue;
            if (visitedF[g]) continue;
            if (ei < 0 || ei >= nE) continue;
            if (inPrimalTree[ei]) continue; // forbidden
            visitedF[g] = 1;
            parentF[g] = f;
            parentDualEdge[g] = ei;
            inDualTree[ei] = 1;
            fq.push_back(g);
        }
    }

    // 3) Select remaining edges R = E \ (T ∪ T*)
    std::vector<int> remainingEdges;
    remainingEdges.reserve(2 * genus);
    for (int ei = 0; ei < nE; ++ei) {
        const auto& e = edges[ei];
        if (e.f0 < 0 || e.f1 < 0) continue; // boundary (shouldn't happen for closed)
        if (inPrimalTree[ei]) continue;
        if (inDualTree[ei]) continue;
        remainingEdges.push_back(ei);
    }

    // Expect size==2g for a closed connected mesh, but keep robust.
    const int target = 2 * genus;
    if ((int)remainingEdges.size() > target) remainingEdges.resize(target);

    // Helper: LCA by climbing with depth (tree is rooted at 0)
    auto lca = [&](int a, int b) {
        int x = a, y = b;
        while (x != y) {
            if (x < 0 || y < 0) return 0;
            if (depth[x] > depth[y]) x = parentV[x];
            else if (depth[y] > depth[x]) y = parentV[y];
            else { x = parentV[x]; y = parentV[y]; }
        }
        return x;
    };

    // 4) For each remaining edge, build fundamental cycle: e + path_T(u,v)
    for (int ei : remainingEdges) {
        const auto& e = edges[ei];
        const int u = e.v0;
        const int v = e.v1;
        if (u < 0 || u >= nV || v < 0 || v >= nV) continue;

        const int w = lca(u, v);
        // Path u -> w
        std::vector<int> pathU;
        int cur = u;
        while (cur != w && cur >= 0 && cur < nV && parentV[cur] >= 0) {
            pathU.push_back(cur);
            cur = parentV[cur];
        }
        pathU.push_back(w);

        // Path v -> w
        std::vector<int> pathV;
        cur = v;
        while (cur != w && cur >= 0 && cur < nV && parentV[cur] >= 0) {
            pathV.push_back(cur);
            cur = parentV[cur];
        }
        pathV.push_back(w);

        // Build ordered vertex loop: u -> ... -> w -> ... -> v -> u
        std::vector<int> loopVerts;
        loopVerts.reserve(pathU.size() + pathV.size());
        for (int vv : pathU) loopVerts.push_back(vv);
        for (int i = (int)pathV.size() - 2; i >= 0; --i) loopVerts.push_back(pathV[i]); // skip w duplicate
        loopVerts.push_back(u); // close using e (v->u eventually handled in edges)

        Cycle cyc;
        // Convert consecutive vertex pairs into signed edges.
        // Note: the closing edge in loopVerts is u at end; the last segment is from v to u.
        for (size_t i = 0; i + 1 < loopVerts.size(); ++i) {
            const int a = loopVerts[i];
            const int b = loopVerts[i + 1];
            const int eidx = findEdge(a, b);
            const int s = edgeSign(a, b);
            if (eidx >= 0 && s != 0) {
                cyc.push_back({eidx, s});
            }
        }

        if (!cyc.empty()) cycles.push_back(std::move(cyc));
        if ((int)cycles.size() >= target) break;
    }

    return cycles;
}

} // namespace topology

#endif

