#include "MIQQuad.h"
#include <cmath>
#include <map>

using namespace Eigen;

namespace {
const double kQuarterTurn = std::acos(-1.0) / 2.0;
} // namespace

MIQQuad::MIQQuad(Mesh& mesh0)
    : GlobalFieldsParameterization(mesh0)
{
}

void MIQQuad::buildLocalFrame(const Vector3d& n, Vector3d& t1, Vector3d& t2)
{
    Vector3d ref(1, 0, 0);
    if (std::fabs(n.dot(ref)) > 0.9) ref = Vector3d(0, 1, 0);
    t1 = (ref - ref.dot(n) * n).normalized();
    t2 = n.cross(t1).normalized();
}

double MIQQuad::cotan(const Vector3d& a, const Vector3d& b, const Vector3d& c)
{
    const Vector3d u = a - b;
    const Vector3d v = c - b;
    const double d = u.dot(v);
    const double cr = u.cross(v).norm();
    if (cr < 1e-12) return 0.0;
    return d / cr;
}

bool MIQQuad::initMeshData()
{
    nV = static_cast<int>(mesh.vertices.size());
    nF = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (!f->isBoundary()) ++nF;
    }
    if (nV < 3 || nF < 1) return false;

    vertPos.resize(nV, 3);
    faces.resize(nF * 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        vertPos.row(v->index) = v->position;
    }

    int fi = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        faces[fi * 3] = f->he->vertex->index;
        faces[fi * 3 + 1] = f->he->next->vertex->index;
        faces[fi * 3 + 2] = f->he->next->next->vertex->index;
        ++fi;
    }

    edgeList.clear();
    std::map<std::pair<int, int>, int> em;
    for (int f = 0; f < nF; ++f) {
        for (int k = 0; k < 3; ++k) {
            int v1 = faces[f * 3 + k];
            int v2 = faces[f * 3 + (k + 1) % 3];
            if (v1 > v2) std::swap(v1, v2);
            const auto key = std::make_pair(v1, v2);
            auto it = em.find(key);
            if (it == em.end()) {
                Edge e;
                e.v1 = v1;
                e.v2 = v2;
                e.f1 = f;
                e.f2 = -1;
                e.idx = static_cast<int>(edgeList.size());
                em[key] = e.idx;
                edgeList.push_back(e);
            } else {
                edgeList[it->second].f2 = f;
            }
        }
    }
    nE = static_cast<int>(edgeList.size());
    jump.resize(nE);
    jump.setZero();

    mixedIntegerProgram_.reset();
    return true;
}

void MIQQuad::initCrossField()
{
    faceN.resize(nF, 3);
    theta.resize(nF);
    theta.setZero();

    for (int fi = 0; fi < nF; ++fi) {
        const int v0 = faces[fi * 3];
        const int v1 = faces[fi * 3 + 1];
        const int v2 = faces[fi * 3 + 2];
        const Vector3d p0 = vertPos.row(v0);
        const Vector3d p1 = vertPos.row(v1);
        const Vector3d p2 = vertPos.row(v2);
        const Vector3d fn = (p1 - p0).cross(p2 - p0).normalized();
        faceN.row(fi) = fn;

        Vector3d t1, t2;
        buildLocalFrame(fn, t1, t2);
        const Vector3d e0 = p1 - p0;
        const Vector3d e1 = p2 - p1;
        const Vector3d e2 = p0 - p2;
        const double l0 = e0.squaredNorm();
        const double l1 = e1.squaredNorm();
        const double l2 = e2.squaredNorm();
        Vector3d d;
        if (l0 >= l1 && l0 >= l2) d = e0;
        else if (l1 >= l0 && l1 >= l2) d = e1;
        else d = e2;
        const double dx = d.dot(t1);
        const double dy = d.dot(t2);
        theta(fi) = std::atan2(dy, dx);
    }
}

std::vector<MixedIntegerProgram::Constraint> MIQQuad::buildFaceAdjacencyConstraints() const
{
    std::vector<MixedIntegerProgram::Constraint> constraints;
    constraints.reserve(edgeList.size());
    for (const auto& e : edgeList) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        MixedIntegerProgram::Constraint c;
        c.i = e.f1;
        c.j = e.f2;
        c.idx = e.idx;
        constraints.push_back(c);
    }
    return constraints;
}

double MIQQuad::crossFieldEnergy() const
{
    if (!mixedIntegerProgram_) return 0.0;
    return mixedIntegerProgram_->energy(theta, jump);
}

bool MIQQuad::solveCrossFieldIP()
{
    auto constraints = buildFaceAdjacencyConstraints();
    mixedIntegerProgram_ = std::make_unique<MixedIntegerProgram>(
        nF,
        std::move(constraints),
        kQuarterTurn);
    mixedIntegerProgram_->setIntegerBounds(jumpLo_, jumpHi_);
    mixedIntegerProgram_->setAlternatingIterations(crossIters_);
    mixedIntegerProgram_->setRefinePasses(jumpRefinePasses_);
    if (!mixedIntegerProgram_->solve(theta, jump)) {
        return false;
    }

    wrapCrossFieldAngles();
    return true;
}

void MIQQuad::wrapCrossFieldAngles()
{
    for (int fi = 0; fi < theta.size(); ++fi) {
        theta(fi) = std::fmod(theta(fi), kQuarterTurn);
        if (theta(fi) < 0.0) {
            theta(fi) += kQuarterTurn;
        }
    }
}

void MIQQuad::buildTargetDirs()
{
    faceD1.resize(nF, 3);
    faceD2.resize(nF, 3);

    for (int fi = 0; fi < nF; ++fi) {
        const Vector3d n = faceN.row(fi);
        Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);
        const double c = std::cos(theta(fi));
        const double s = std::sin(theta(fi));
        faceD1.row(fi) = (c * t1 + s * t2).normalized();
        faceD2.row(fi) = (-s * t1 + c * t2).normalized();
    }
}

void MIQQuad::buildCotLaplacian(SparseMatrix<double>& L)
{
    L.resize(nV, nV);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nV);

    for (const auto& e : edgeList) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        const int vi = e.v1;
        const int vj = e.v2;

        int o1 = -1;
        int o2 = -1;
        for (int k = 0; k < 3; ++k) {
            const int v = faces[e.f1 * 3 + k];
            if (v != vi && v != vj) {
                o1 = v;
                break;
            }
        }
        for (int k = 0; k < 3; ++k) {
            const int v = faces[e.f2 * 3 + k];
            if (v != vi && v != vj) {
                o2 = v;
                break;
            }
        }

        double w = 0.5 * (cotan(vertPos.row(o1), vertPos.row(vi), vertPos.row(vj)) +
                          cotan(vertPos.row(o2), vertPos.row(vi), vertPos.row(vj)));
        if (w < 0) w = 0.0;

        diag(vi) += w;
        diag(vj) += w;
        trips.emplace_back(vi, vj, -w);
        trips.emplace_back(vj, vi, -w);
    }
    for (int i = 0; i < nV; ++i) {
        trips.emplace_back(i, i, diag(i) + 1e-8);
    }
    L.setFromTriplets(trips.begin(), trips.end());
}

void MIQQuad::solvePoisson()
{
    MatrixXd Vd1(nV, 3);
    MatrixXd Vd2(nV, 3);
    Vd1.setZero();
    Vd2.setZero();
    VectorXd vA = VectorXd::Zero(nV);

    for (int fi = 0; fi < nF; ++fi) {
        const int v0 = faces[fi * 3];
        const int v1 = faces[fi * 3 + 1];
        const int v2 = faces[fi * 3 + 2];
        const Vector3d a = vertPos.row(v0);
        const Vector3d b = vertPos.row(v1);
        const Vector3d c = vertPos.row(v2);
        const double area = 0.5 * (b - a).cross(c - a).norm();
        for (int k = 0; k < 3; ++k) {
            const int v = faces[fi * 3 + k];
            vA(v) += area;
            Vd1.row(v) += area * faceD1.row(fi);
            Vd2.row(v) += area * faceD2.row(fi);
        }
    }
    for (int i = 0; i < nV; ++i) {
        if (vA(i) > 1e-12) {
            Vd1.row(i) /= vA(i);
            Vd2.row(i) /= vA(i);
        }
    }

    VectorXd div1 = VectorXd::Zero(nV);
    VectorXd div2 = VectorXd::Zero(nV);
    for (int fi = 0; fi < nF; ++fi) {
        const int v0 = faces[fi * 3];
        const int v1 = faces[fi * 3 + 1];
        const int v2 = faces[fi * 3 + 2];
        const Vector3d p0 = vertPos.row(v0);
        const Vector3d p1 = vertPos.row(v1);
        const Vector3d p2 = vertPos.row(v2);
        const Vector3d fn = (p1 - p0).cross(p2 - p0);
        const Vector3d e01 = p1 - p0;
        const Vector3d e12 = p2 - p1;
        const Vector3d e20 = p0 - p2;
        const Vector3d en01 = fn.cross(e01).normalized();
        const Vector3d en12 = fn.cross(e12).normalized();
        const Vector3d en20 = fn.cross(e20).normalized();
        const double cot0 = cotan(p2, p0, p1);
        const double cot1 = cotan(p0, p1, p2);
        const double cot2 = cotan(p1, p2, p0);
        const Vector3d fd1 = faceD1.row(fi);
        const Vector3d fd2 = faceD2.row(fi);
        const double d1 = 0.5 * (fd1.dot(en01) * cot2 + fd1.dot(en12) * cot0 + fd1.dot(en20) * cot1);
        const double d2 = 0.5 * (fd2.dot(en01) * cot2 + fd2.dot(en12) * cot0 + fd2.dot(en20) * cot1);
        div1(v0) += d1;
        div1(v1) += d1;
        div1(v2) += d1;
        div2(v0) += d2;
        div2(v1) += d2;
        div2(v2) += d2;
    }

    SparseMatrix<double> L;
    buildCotLaplacian(L);

    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        if (v->isBoundary()) {
            L.coeffRef(v->index, v->index) += 1e6;
            div1(v->index) = 0.0;
            div2(v->index) = 0.0;
        }
    }
    bool hasBnd = false;
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        if (v->isBoundary()) {
            hasBnd = true;
            break;
        }
    }
    if (!hasBnd) {
        L.coeffRef(0, 0) += 1e6;
        div1(0) = 0.0;
        div2(0) = 0.0;
    }

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Success) {
        UV = MatrixXd::Zero(nV, 2);
        return;
    }
    const VectorXd u = solver.solve(div1);
    const VectorXd v = solver.solve(div2);

    const double uMin = u.minCoeff();
    const double uMax = u.maxCoeff();
    const double vMin = v.minCoeff();
    const double vMax = v.maxCoeff();
    double uR = uMax - uMin;
    double vR = vMax - vMin;
    if (uR < 1e-10) uR = 1.0;
    if (vR < 1e-10) vR = 1.0;

    UV.resize(nV, 2);
    for (int i = 0; i < nV; ++i) {
        UV(i, 0) = (u(i) - uMin) / uR;
        UV(i, 1) = (v(i) - vMin) / vR;
    }
}

void MIQQuad::parameterize()
{
    if (!initMeshData()) return;

    initCrossField();
    if (!solveCrossFieldIP()) return;
    buildTargetDirs();
    solvePoisson();

    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        v->uv = Vector2d(UV(v->index, 0), UV(v->index, 1));
    }

    normalize();
}
