#include "PGP.h"
#include <cmath>
#include <map>
#include <algorithm>

using namespace Eigen;

namespace {

// Rotate a 2D vector by k·π/2, k ∈ {0,1,2,3}
Vector2d rotateK(int k, const Vector2d& v)
{
    switch (k & 3) {
        case 0: return v;
        case 1: return Vector2d(-v(1),  v(0));
        case 2: return Vector2d(-v(0), -v(1));
        case 3: return Vector2d( v(1), -v(0));
    }
    return v;
}

double cotan(const Vector3d& a, const Vector3d& b, const Vector3d& c)
{
    Vector3d u = a - b;
    Vector3d v = c - b;
    double d  = u.dot(v);
    double cr = u.cross(v).norm();
    return (cr < 1e-12) ? 0.0 : d / cr;
}

} // anonymous namespace

PGP::PGP(Mesh& mesh0)
    : Parameterization(mesh0),
      nV(0), nF(0),
      hasExternalDirs(false),
      maxIters(3),
      doLocalSearch(true),
      solverReady(false)
{
}

void PGP::setCrossField(const MatrixXd& dirs)
{
    faceDirs    = dirs;
    hasExternalDirs = (dirs.cols() == 3 && dirs.rows() > 0);
}

// ============================================================
//  buildLocalFrame
// ============================================================
void PGP::buildLocalFrame(const Vector3d& n, Vector3d& t1, Vector3d& t2)
{
    Vector3d ref(1.0, 0.0, 0.0);
    if (std::fabs(n.dot(ref)) > 0.9) ref = Vector3d(0.0, 1.0, 0.0);
    t1 = (ref - ref.dot(n) * n).normalized();
    t2 = n.cross(t1).normalized();
}

// ============================================================
//  initMeshData
// ============================================================
bool PGP::initMeshData()
{
    nV = static_cast<int>(mesh.vertices.size());
    nF = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f)
        if (!f->isBoundary()) ++nF;
    if (nV < 3 || nF < 1) return false;

    vertPos.resize(nV, 3);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v)
        vertPos.row(v->index) = v->position;

    faceIndices.resize(nF, 3);
    int fi = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        faceIndices(fi, 0) = f->he->vertex->index;
        faceIndices(fi, 1) = f->he->next->vertex->index;
        faceIndices(fi, 2) = f->he->next->next->vertex->index;
        ++fi;
    }

    // internal edges
    edges.clear();
    std::map<std::pair<int,int>, int> emap;
    for (int f = 0; f < nF; ++f) {
        for (int k = 0; k < 3; ++k) {
            int a = faceIndices(f, k);
            int b = faceIndices(f, (k + 1) % 3);
            if (a > b) std::swap(a, b);
            auto key = std::make_pair(a, b);
            auto it  = emap.find(key);
            if (it == emap.end()) {
                InternalEdge e;
                e.v1 = a; e.v2 = b; e.f1 = f; e.f2 = -1;
                e.idx = static_cast<int>(edges.size());
                emap[key] = e.idx;
                edges.push_back(e);
            } else {
                edges[it->second].f2 = f;
            }
        }
    }

    // per-face data
    faceNormals.resize(nF, 3);
    faceCenters.resize(nF, 3);
    faceAreas.resize(nF);

    for (int f = 0; f < nF; ++f) {
        Vector3d p0 = vertPos.row(faceIndices(f, 0));
        Vector3d p1 = vertPos.row(faceIndices(f, 1));
        Vector3d p2 = vertPos.row(faceIndices(f, 2));
        Vector3d n  = (p1 - p0).cross(p2 - p0);
        double area = 0.5 * n.norm();
        faceAreas[f]   = std::max(area, 1e-16);
        faceNormals.row(f) = n.normalized();
        faceCenters.row(f) = (p0 + p1 + p2) / 3.0;
    }

    // integer variables
    kVar.resize(nF);   kVar.setZero();
    pVar.resize(nF, 2); pVar.setZero();

    solverReady = false;
    return true;
}

// ============================================================
//  buildLocalFrames + buildGradientOperators
// ============================================================
void PGP::buildLocalFrames()
{
    faceGrad.resize(nF);
    faceD1.resize(nF);
    faceD2.resize(nF);

    for (int f = 0; f < nF; ++f) {
        Vector3d p0 = vertPos.row(faceIndices(f, 0));
        Vector3d p1 = vertPos.row(faceIndices(f, 1));
        Vector3d p2 = vertPos.row(faceIndices(f, 2));
        Vector3d n  = faceNormals.row(f);

        Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);

        // project to 2D  (q0 at origin)
        Vector3d d1 = p1 - p0, d2 = p2 - p0;
        Vector2d q1(d1.dot(t1), d1.dot(t2));
        Vector2d q2(d2.dot(t1), d2.dot(t2));

        double A2 = 0.5 * std::fabs(q1(0) * q2(1) - q1(1) * q2(0));
        if (A2 < 1e-16) A2 = 1e-16;

        // gradient operator  G_T ∈ R^{2×3}
        // ∇u = G_T · [u0,u1,u2]^T
        // row 0 (∂/∂x):  (y1-y2,  y2-y0,  y0-y1) / (2A)
        // row 1 (∂/∂y):  (x2-x1,  x0-x2,  x1-x0) / (2A)
        Matrix2d G;
        G(0,0) =  q1(1) - q2(1);
        G(0,1) =  q2(1);
        G(0,2) = -q1(1);
        G(1,0) =  q2(0) - q1(0);
        G(1,1) = -q2(0);
        G(1,2) =  q1(0);
        G /= (2.0 * A2);

        faceGrad[f] = G;
    }
}

// ============================================================
//  computeCrossField
// ============================================================
void PGP::computeCrossField()
{
    faceDirs.resize(nF, 3);

    if (hasExternalDirs && faceDirs.rows() == nF) {
        for (int f = 0; f < nF; ++f) {
            Vector3d n = faceNormals.row(f);
            Vector3d d = faceDirs.row(f);
            d -= d.dot(n) * n;
            double len = d.norm();
            faceDirs.row(f) = (len > 1e-12) ? (d / len) : Vector3d::UnitX();
        }
    } else {
        // heuristic: longest edge direction per face
        for (int f = 0; f < nF; ++f) {
            Vector3d p0 = vertPos.row(faceIndices(f, 0));
            Vector3d p1 = vertPos.row(faceIndices(f, 1));
            Vector3d p2 = vertPos.row(faceIndices(f, 2));
            Vector3d e01 = p1 - p0, e12 = p2 - p1, e20 = p0 - p2;
            double l01 = e01.squaredNorm();
            double l12 = e12.squaredNorm();
            double l20 = e20.squaredNorm();

            Vector3d best = e01;
            if      (l12 > l01 && l12 > l20) best = e12;
            else if (l20 > l01 && l20 > l12) best = e20;

            Vector3d n = faceNormals.row(f);
            best -= best.dot(n) * n;
            double len = best.norm();
            faceDirs.row(f) = (len > 1e-12) ? (best / len) : Vector3d::UnitX();
        }
    }

    // build faceD1, faceD2 in 2D local frames
    for (int f = 0; f < nF; ++f) {
        Vector3d n = faceNormals.row(f);
        Vector3d t1, t2;
        buildLocalFrame(n, t1, t2);

        Vector3d d = faceDirs.row(f);
        d -= d.dot(n) * n;
        if (d.norm() < 1e-12) d = t1;
        d.normalize();

        faceD1[f] = Vector2d(d.dot(t1), d.dot(t2));
        faceD2[f] = Vector2d(-faceD1[f](1), faceD1[f](0));
    }
}

// ============================================================
//  assembleLaplacian
// ============================================================
void PGP::assembleLaplacian(SparseMatrix<double>& L)
{
    L.resize(nV, nV);
    std::vector<Triplet<double>> trips;
    trips.reserve(nF * 9);

    for (int f = 0; f < nF; ++f) {
        double A = faceAreas[f];
        const Matrix2d& G = faceGrad[f];          // 2×3
        Matrix3d Mf = 2.0 * A * G.transpose() * G; // 3×3

        int vi[3] = { faceIndices(f,0), faceIndices(f,1), faceIndices(f,2) };
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (std::fabs(Mf(i,j)) > 1e-16)
                    trips.emplace_back(vi[i], vi[j], Mf(i,j));
    }
    L.setFromTriplets(trips.begin(), trips.end());
}

// ============================================================
//  assembleRHS
//  bu = Σ_f A_f · G_T^T · (R(k·π/2)·d1_f + p_T)
//  bv = Σ_f A_f · G_T^T · (R(k·π/2)·d2_f + p_T)
// ============================================================
void PGP::assembleRHS(const VectorXi& kT, const MatrixXi& pT,
                      VectorXd& bu, VectorXd& bv)
{
    bu.resize(nV); bu.setZero();
    bv.resize(nV); bv.setZero();

    for (int f = 0; f < nF; ++f) {
        double A    = faceAreas[f];
        int    k    = kT(f);
        Vector2d pu(static_cast<double>(pT(f,0)),
                     static_cast<double>(pT(f,1)));

        Vector2d d1k = rotateK(k, faceD1[f]);
        Vector2d d2k = rotateK(k, faceD2[f]);

        Vector2d tu = d1k + pu;  // target gradient for u
        Vector2d tv = d2k + pu;  // target gradient for v

        const Matrix2d& G = faceGrad[f];          // 2×3
        Vector3d rhsu = A * G.transpose() * tu;    // 3×1
        Vector3d rhsv = A * G.transpose() * tv;    // 3×1

        int vi[3] = { faceIndices(f,0), faceIndices(f,1), faceIndices(f,2) };
        for (int i = 0; i < 3; ++i) {
            bu(vi[i]) += rhsu(i);
            bv(vi[i]) += rhsv(i);
        }
    }
}

// ============================================================
//  solveForUV
// ============================================================
void PGP::solveForUV(const VectorXi& kT, const MatrixXi& pT,
                     VectorXd& u, VectorXd& v)
{
    if (!solverReady) {
        SparseMatrix<double> L;
        assembleLaplacian(L);
        // pin vertex 0 to (0,0)
        for (SparseMatrix<double>::InnerIterator it(L, 0); it; ++it) {
            if (it.row() != it.col())
                L.coeffRef(it.row(), it.col()) = 0.0;
        }
        L.coeffRef(0, 0) = 1.0;
        solver.compute(L);
        solverReady = true;
    }
    if (solver.info() != Success) return;

    VectorXd bu, bv;
    assembleRHS(kT, pT, bu, bv);
    bu(0) = 0.0; bv(0) = 0.0;   // pin

    u = solver.solve(bu);
    v = solver.solve(bv);

    if (solver.info() != Success) {
        // fallback: CG
        SparseMatrix<double> L;
        assembleLaplacian(L);
        ConjugateGradient<SparseMatrix<double>> cg;
        cg.compute(L);
        u = cg.solveWithGuess(bu, u);
        v = cg.solveWithGuess(bv, v);
    }
}

// ============================================================
//  computeGradient
// ============================================================
Vector2d PGP::computeGradient(const VectorXd& u, int f) const
{
    const Matrix2d& G = faceGrad[f];
    Vector3d uf{ u(faceIndices(f,0)), u(faceIndices(f,1)), u(faceIndices(f,2)) };
    return G * uf;
}

// ============================================================
//  faceEnergy
//  E_f = A_f [||∇u - (R(k·π/2)·d1 + p)||² + ||∇v - (R(k·π/2)·d2 + p)||²]
// ============================================================
double PGP::faceEnergy(int f, const Vector2d& gu, const Vector2d& gv,
                        int k, const Vector2i& p) const
{
    Vector2d d1k = rotateK(k, faceD1[f]);
    Vector2d d2k = rotateK(k, faceD2[f]);
    Vector2d pv(static_cast<double>(p(0)), static_cast<double>(p(1)));

    double eu = (gu - d1k - pv).squaredNorm();
    double ev = (gv - d2k - pv).squaredNorm();
    return faceAreas[f] * (eu + ev);
}

// ============================================================
//  optimizePerFace
// ============================================================
void PGP::optimizePerFace(const VectorXd& u, const VectorXd& v,
                           VectorXi& kT, MatrixXi& pT)
{
    for (int f = 0; f < nF; ++f) {
        Vector2d gu = computeGradient(u, f);
        Vector2d gv = computeGradient(v, f);

        double bestE = 1e30;
        int    bestK = 0;
        Vector2i bestP(0, 0);

        for (int k = 0; k < 4; ++k) {
            Vector2d d1k = rotateK(k, faceD1[f]);
            Vector2d d2k = rotateK(k, faceD2[f]);

            // estimate continuous optimal p (average of u and v targets)
            Vector2d pc = (gu - d1k + gv - d2k) * 0.5;

            // try rounding neighborhood
            for (int dp0 = -1; dp0 <= 1; ++dp0) {
                for (int dp1 = -1; dp1 <= 1; ++dp1) {
                    Vector2i pi;
                    pi(0) = static_cast<int>(std::round(pc(0))) + dp0;
                    pi(1) = static_cast<int>(std::round(pc(1))) + dp1;
                    double e = faceEnergy(f, gu, gv, k, pi);
                    if (e < bestE) { bestE = e; bestK = k; bestP = pi; }
                }
            }
        }
        kT(f)       = bestK;
        pT(f, 0)    = bestP(0);
        pT(f, 1)    = bestP(1);
    }
}

// ============================================================
//  localSearchP
// ============================================================
void PGP::localSearchP(const VectorXd& u, const VectorXd& v,
                        VectorXi& kT, MatrixXi& pT)
{
    static const int delta[8][2] = {{1,0},{-1,0},{0,1},{0,-1},
                                     {1,1},{1,-1},{-1,1},{-1,-1}};
    for (int f = 0; f < nF; ++f) {
        Vector2d gu = computeGradient(u, f);
        Vector2d gv = computeGradient(v, f);
        int    k    = kT(f);
        Vector2i p(pT(f,0), pT(f,1));

        double curE = faceEnergy(f, gu, gv, k, p);
        bool improved = true;
        while (improved) {
            improved = false;
            for (int d = 0; d < 8; ++d) {
                Vector2i np(p(0)+delta[d][0], p(1)+delta[d][1]);
                double ne = faceEnergy(f, gu, gv, k, np);
                if (ne < curE) {
                    p = np; curE = ne; improved = true; break;
                }
            }
        }
        pT(f,0) = p(0); pT(f,1) = p(1);
    }
}

// ============================================================
//  writeUV
// ============================================================
void PGP::writeUV(const VectorXd& u, const VectorXd& v)
{
    for (VertexIter vit = mesh.vertices.begin(); vit != mesh.vertices.end(); ++vit)
        vit->uv = Vector2d(u(vit->index), v(vit->index));
    normalize();
}

// ============================================================
//  parameterize  (main entry)
// ============================================================
void PGP::parameterize()
{
    if (!initMeshData()) return;

    buildLocalFrames();
    computeCrossField();

    // initial solve with all zeros
    VectorXd u(nV), v(nV);
    u.setZero(); v.setZero();
    solveForUV(kVar, pVar, u, v);

    for (int iter = 0; iter < maxIters; ++iter) {
        optimizePerFace(u, v, kVar, pVar);
        solveForUV(kVar, pVar, u, v);

        if (doLocalSearch) {
            localSearchP(u, v, kVar, pVar);
            solveForUV(kVar, pVar, u, v);
        }
    }

    writeUV(u, v);
}
