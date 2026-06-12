#include "GaussianCurvature.h"

#include "../Mesh.h"
#include "../Vertex.h"
#include "../Face.h"
#include "../HalfEdge.h"

#include <cmath>
#include <vector>

using Eigen::Vector3d;
using Eigen::VectorXd;

namespace {

double angleAtVertexInFace(const Vertex& v, const HalfEdge& heAtV) {
    // heAtV is a halfedge whose tail vertex is v, inside a (non-boundary) face.
    const Vector3d& p = v.position;
    const Vector3d& p1 = heAtV.next->vertex->position;
    const Vector3d& p2 = heAtV.next->next->vertex->position;

    Vector3d e1 = p1 - p;
    Vector3d e2 = p2 - p;
    const double n1 = e1.norm();
    const double n2 = e2.norm();
    if (n1 < 1e-14 || n2 < 1e-14) return 0.0;
    e1 /= n1;
    e2 /= n2;
    double d = e1.dot(e2);
    d = std::max(-1.0, std::min(1.0, d));
    return std::acos(d);
}

} // namespace

namespace geometry {

VectorXd gaussianCurvatureAngleDeficit(const Mesh& mesh) {
    const int nV = (int)mesh.vertices.size();
    VectorXd K = VectorXd::Zero(nV);

    for (VertexCIter vIt = mesh.vertices.begin(); vIt != mesh.vertices.end(); ++vIt) {
        if (vIt->isIsolated()) continue;

        const bool isBnd = vIt->isBoundary();
        const double target = isBnd ? M_PI : (2.0 * M_PI);

        double sumAngles = 0.0;
        HalfEdgeCIter h = vIt->he;
        // walk around one-ring (may visit boundary face halfedges too; skip boundary faces)
        do {
            if (!h->face->isBoundary()) {
                sumAngles += angleAtVertexInFace(*vIt, *h);
            }
            h = h->flip->next;
        } while (h != vIt->he);

        K(vIt->index) = target - sumAngles;
    }

    return K;
}

VectorXd gaussianCurvaturePerArea(const Mesh& mesh, double eps) {
    const int nV = (int)mesh.vertices.size();
    VectorXd K = gaussianCurvatureAngleDeficit(mesh);
    VectorXd A = VectorXd::Zero(nV);

    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        HalfEdgeCIter he = f->he;
        const Vector3d p0 = he->vertex->position;
        const Vector3d p1 = he->next->vertex->position;
        const Vector3d p2 = he->next->next->vertex->position;
        const double area = 0.5 * (p1 - p0).cross(p2 - p0).norm();
        if (!std::isfinite(area) || area <= 0.0) continue;
        const double oneThird = area / 3.0;
        A(he->vertex->index) += oneThird;
        A(he->next->vertex->index) += oneThird;
        A(he->next->next->vertex->index) += oneThird;
    }

    for (int i = 0; i < nV; ++i) {
        const double denom = std::max(A(i), eps);
        K(i) = K(i) / denom;
    }
    return K;
}

} // namespace geometry

