#include "AbfPlusPlus.h"
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_set>
#include <vector>

namespace {

constexpr double kPi     = 3.14159265358979323846;
constexpr double kTwoPi  = 6.28318530717958647692;

// --------------- helpers ---------------

double angle3D(const Eigen::Vector3d& a,
               const Eigen::Vector3d& b,
               const Eigen::Vector3d& c)
{
    Eigen::Vector3d u = b - a;
    Eigen::Vector3d v = c - a;
    double nu = u.norm();
    double nv = v.norm();
    if (nu < 1e-14 || nv < 1e-14) return kPi / 3.0;
    u /= nu;
    v /= nv;
    double d = std::max(-1.0, std::min(1.0, u.dot(v)));
    return std::acos(d);
}

void collectBoundaryLoops(const Mesh& mesh,
                          std::vector<std::vector<int>>& loops)
{
    loops.clear();
    loops.reserve(mesh.boundaries.size());
    for (size_t i = 0; i < mesh.boundaries.size(); i++) {
        HalfEdgeCIter start = mesh.boundaries[i];
        HalfEdgeCIter he    = start;
        std::vector<int> loop;
        do {
            loop.push_back(he->vertex->index);
            he = he->next;
        } while (he != start);
        if (loop.size() >= 3) loops.push_back(loop);
    }
}

double loopPerimeter(const Mesh& mesh, const std::vector<int>& loop)
{
    double p = 0.0;
    int n = static_cast<int>(loop.size());
    for (int i = 0; i < n; i++) {
        int a = loop[i];
        int b = loop[(i + 1) % n];
        p += (mesh.vertices[b].position - mesh.vertices[a].position).norm();
    }
    return p;
}

void pinLoopOnCircle(Mesh& mesh,
                     const std::vector<int>& loop,
                     const Eigen::Vector2d& center,
                     double radius)
{
    int n = static_cast<int>(loop.size());
    if (n < 3) return;
    for (int i = 0; i < n; i++) {
        double t     = static_cast<double>(i) / static_cast<double>(n);
        double theta = kTwoPi * t;
        mesh.vertices[loop[i]].uv = Eigen::Vector2d(
            center.x() + radius * std::cos(theta),
            center.y() + radius * std::sin(theta));
    }
}

inline double safeCot(double angle)
{
    double s = std::sin(angle);
    if (std::abs(s) < 1e-12)
        return (s >= 0.0) ? 1e12 : -1e12;
    return std::cos(angle) / s;
}

} // anonymous namespace


// ================================================================
//  AbfPlusPlus — public
// ================================================================

AbfPlusPlus::AbfPlusPlus(Mesh& mesh0, int maxIterations)
    : Parameterization(mesh0)
    , m_maxIterations(std::max(1, maxIterations))
{
}

void AbfPlusPlus::parameterize()
{
    if (mesh.boundaries.empty()) {
        std::cerr << "[AbfPlusPlus] closed mesh is not supported.\n";
        return;
    }

    // Phase 1 — angle variable mapping
    buildAngleMapping();
    if (nAngles_ == 0) {
        normalize();
        return;
    }

    // Phase 2 — build constraint system
    buildConstraints();

    // Phase 3 — configure Augmented Lagrangian solver
    // ABF++ constraints are linear (angle sums = constants),
    // so ∇²c_j ≡ 0  →  constraint Hessian callback returns empty.
    using namespace std::placeholders;
    auto conFunc  = std::bind(&AbfPlusPlus::evalConstraint,        this, _1);
    auto jacFunc  = std::bind(&AbfPlusPlus::evalJacobian,          this, _1);
    auto hessFunc = std::bind(&AbfPlusPlus::evalConstraintHessian, this, _1, _2, _3);

    AugmentedLagrangian::Config cfg;
    cfg.maxOuter = std::min(m_maxIterations, 40);
    al_ = std::make_unique<AugmentedLagrangian>(
        nAngles_, nConstraints_,
        conFunc, jacFunc, hessFunc, cfg);

    // Phase 4 — solve in angle space
    Eigen::VectorXd alpha = beta_;               // initialise at target angles
    al_->solve(alpha, beta_);

    // Phase 5 — reconstruct UV from optimal angles
    angleToUv(alpha);

    normalize();
}


// ================================================================
//  Phase 1 : build angle variable mapping
// ================================================================

void AbfPlusPlus::buildAngleMapping()
{
    nInteriorFaces_ = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++)
        if (!f->isBoundary()) nInteriorFaces_++;

    nAngles_ = 3 * nInteriorFaces_;
    faceAngles_.assign(mesh.faces.size(), {-1, -1, -1});
    angleToVertex_.resize(nAngles_);
    beta_.resize(nAngles_);

    int idx = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;

        int fi = f->index;

        HalfEdgeCIter he = f->he;
        int v0 = he->vertex->index;
        const Eigen::Vector3d& p0 = he->vertex->position;
        he = he->next;
        int v1 = he->vertex->index;
        const Eigen::Vector3d& p1 = he->vertex->position;
        he = he->next;
        int v2 = he->vertex->index;
        const Eigen::Vector3d& p2 = he->vertex->position;

        faceAngles_[fi] = {{idx, idx + 1, idx + 2}};
        angleToVertex_[idx]     = v0;
        angleToVertex_[idx + 1] = v1;
        angleToVertex_[idx + 2] = v2;

        // 3‑D target angles (angle at v0 is ∠(v1,v0,v2), etc.)
        beta_(idx)     = angle3D(p0, p1, p2);
        beta_(idx + 1) = angle3D(p1, p2, p0);
        beta_(idx + 2) = angle3D(p2, p0, p1);

        idx += 3;
    }
}


// ================================================================
//  Phase 2 : build constraint system
// ================================================================

void AbfPlusPlus::buildConstraints()
{
    int nVerts = static_cast<int>(mesh.vertices.size());

    // which vertices participate? (at least one incident interior face)
    std::vector<bool> participates(nVerts, false);
    for (int a = 0; a < nAngles_; a++)
        participates[angleToVertex_[a]] = true;

    // count
    int nTriConstr  = nInteriorFaces_;
    int nVertConstr = 0;
    for (int v = 0; v < nVerts; v++)
        if (participates[v]) nVertConstr++;

    nConstraints_ = nTriConstr + nVertConstr;
    constraintRhs_.resize(nConstraints_);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * nTriConstr + 8 * nVertConstr);

    // ---- triangle constraints : α₁+α₂+α₃ = π ----
    int row = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        auto& fa = faceAngles_[f->index];
        triplets.emplace_back(row, fa[0], 1.0);
        triplets.emplace_back(row, fa[1], 1.0);
        triplets.emplace_back(row, fa[2], 1.0);
        constraintRhs_(row) = kPi;
        row++;
    }

    // ---- vertex constraint RHS ----
    std::vector<double> vertRhs(nVerts, kTwoPi);

    // boundary vertices: 3‑D curvature angle instead of 2π
    std::vector<std::vector<int>> loops;
    collectBoundaryLoops(mesh, loops);
    for (const auto& loop : loops) {
        int n = static_cast<int>(loop.size());
        for (int i = 0; i < n; i++) {
            int v    = loop[i];
            int prev = loop[(i - 1 + n) % n];
            int next = loop[(i + 1) % n];
            vertRhs[v] = angle3D(mesh.vertices[prev].position,
                                 mesh.vertices[v].position,
                                 mesh.vertices[next].position);
        }
    }

    // map vertex → constraint row
    std::vector<int> vertRow(nVerts, -1);
    for (int v = 0; v < nVerts; v++) {
        if (!participates[v]) continue;
        vertRow[v] = row;
        constraintRhs_(row) = vertRhs[v];
        row++;
    }

    // fill J for vertex constraints
    for (int a = 0; a < nAngles_; a++) {
        int v = angleToVertex_[a];
        int r = vertRow[v];
        if (r >= 0)
            triplets.emplace_back(r, a, 1.0);
    }

    constraintJacobian_.resize(nConstraints_, nAngles_);
    constraintJacobian_.setFromTriplets(triplets.begin(), triplets.end());
}


// ================================================================
//  Augmented Lagrangian callbacks
//    (ABF++ constraints are linear:  c(α) = J·α − rhs)
//    ∇²c_j ≡ 0 for all j  →  constraint Hessian returns empty.
// ================================================================

Eigen::VectorXd AbfPlusPlus::evalConstraint(const Eigen::VectorXd& alpha) const
{
    return constraintJacobian_ * alpha - constraintRhs_;
}

Eigen::SparseMatrix<double> AbfPlusPlus::evalJacobian(const Eigen::VectorXd& /*alpha*/) const
{
    // Jacobian is constant: each row sums a subset of α variables.
    return constraintJacobian_;
}

Eigen::SparseMatrix<double> AbfPlusPlus::evalConstraintHessian(
    const Eigen::VectorXd& /*alpha*/,
    const Eigen::VectorXd& /*lambda*/,
    double /*rho*/) const
{
    // All constraints are linear in α  →  ∇²c_j = 0 for every j.
    // Return an empty matrix — the solver will add nothing to the Hessian.
    return Eigen::SparseMatrix<double>(nAngles_, nAngles_);
}


// ================================================================
//  Helper  (kept for API completeness)
// ================================================================

double AbfPlusPlus::boundaryCurvature(int v0, int v1, int v2) const
{
    return angle3D(mesh.vertices[v0].position,
                   mesh.vertices[v1].position,
                   mesh.vertices[v2].position);
}


// ================================================================
//  Phase 5 : angle → UV  (cotangent‑Laplacian reconstruction)
// ================================================================

void AbfPlusPlus::angleToUv(const Eigen::VectorXd& alpha)
{
    int nVerts = static_cast<int>(mesh.vertices.size());

    // ---- accumulate cot‑weights from optimised angles ----
    std::map<std::pair<int,int>, double> edgeWeights;

    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        auto& fa = faceAngles_[f->index];
        if (fa[0] < 0) continue;

        int    v0 = angleToVertex_[fa[0]];
        int    v1 = angleToVertex_[fa[1]];
        int    v2 = angleToVertex_[fa[2]];
        double a0 = alpha(fa[0]);
        double a1 = alpha(fa[1]);
        double a2 = alpha(fa[2]);

        auto p01 = std::minmax(v0, v1);
        auto p12 = std::minmax(v1, v2);
        auto p20 = std::minmax(v2, v0);
        edgeWeights[p01] += safeCot(a2);
        edgeWeights[p12] += safeCot(a0);
        edgeWeights[p20] += safeCot(a1);
    }

    // ---- boundary circle placement (largest loop) ----
    std::vector<std::vector<int>> loops;
    collectBoundaryLoops(mesh, loops);
    if (loops.empty()) return;

    std::sort(loops.begin(), loops.end(),
              [&](const std::vector<int>& a, const std::vector<int>& b) {
                  return loopPerimeter(mesh, a) > loopPerimeter(mesh, b);
              });

    pinLoopOnCircle(mesh, loops[0], Eigen::Vector2d(0.0, 0.0), 1.0);

    int holes = static_cast<int>(loops.size()) - 1;
    if (holes > 0) {
        double holeRadius = std::max(0.04, std::min(0.20, 0.35 / static_cast<double>(holes)));
        double centerRing = (holes == 1) ? 0.0 : 0.55;
        for (int i = 0; i < holes; i++) {
            double angle       = kTwoPi * static_cast<double>(i) / static_cast<double>(holes);
            Eigen::Vector2d c  = Eigen::Vector2d(centerRing * std::cos(angle),
                                                  centerRing * std::sin(angle));
            pinLoopOnCircle(mesh, loops[i + 1], c, holeRadius);
        }
    }

    // mark boundary vertices as fixed
    std::unordered_set<int> boundarySet;
    for (const auto& loop : loops)
        for (int v : loop) boundarySet.insert(v);

    // ---- map: interior vertex → system index ----
    std::vector<int> interiorIdx(nVerts, -1);
    int nInterior = 0;
    for (int v = 0; v < nVerts; v++)
        if (!boundarySet.count(v)) interiorIdx[v] = nInterior++;

    if (nInterior == 0) return;

    // ---- build sparse Laplacian  L·u = 0, L·v = 0  (fixed boundary) ----
    std::vector<Eigen::Triplet<double>> lapTriplets;
    lapTriplets.reserve(edgeWeights.size() * 4);

    Eigen::VectorXd bu = Eigen::VectorXd::Zero(nInterior);
    Eigen::VectorXd bv = Eigen::VectorXd::Zero(nInterior);

    for (const auto& kv : edgeWeights) {
        int i = kv.first.first;
        int j = kv.first.second;
        double w = kv.second;

        int ii = interiorIdx[i];
        int jj = interiorIdx[j];

        if (ii >= 0 && jj >= 0) {
            lapTriplets.emplace_back(ii, ii,  w);
            lapTriplets.emplace_back(ii, jj, -w);
            lapTriplets.emplace_back(jj, ii, -w);
            lapTriplets.emplace_back(jj, jj,  w);
        } else if (ii >= 0) {
            lapTriplets.emplace_back(ii, ii, w);
            bu(ii) += w * mesh.vertices[j].uv.x();
            bv(ii) += w * mesh.vertices[j].uv.y();
        } else if (jj >= 0) {
            lapTriplets.emplace_back(jj, jj, w);
            bu(jj) += w * mesh.vertices[i].uv.x();
            bv(jj) += w * mesh.vertices[i].uv.y();
        }
    }

    Eigen::SparseMatrix<double> L(nInterior, nInterior);
    L.setFromTriplets(lapTriplets.begin(), lapTriplets.end());

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Eigen::Success) return;

    Eigen::VectorXd uSol = solver.solve(bu);
    Eigen::VectorXd vSol = solver.solve(bv);
    if (solver.info() != Eigen::Success) return;

    for (int v = 0; v < nVerts; v++) {
        if (interiorIdx[v] >= 0) {
            mesh.vertices[v].uv = Eigen::Vector2d(uSol(interiorIdx[v]),
                                                   vSol(interiorIdx[v]));
        }
    }
}
