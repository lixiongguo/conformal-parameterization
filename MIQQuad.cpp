#include "MIQQuad.h"
#include <map>
#include <algorithm>
#include <Eigen/SparseCholesky>
using namespace Eigen;

MIQQuad::MIQQuad(Mesh& mesh0): Parameterization(mesh0) {}

int MIQQuad::roundToInt(double x) { return (int)floor(x + 0.5); }

// ============================================================
// Helpers
// ============================================================
void MIQQuad::buildLocalFrame(const Vector3d& n, Vector3d& t1, Vector3d& t2) {
    Vector3d ref(1,0,0);
    if (fabs(n.dot(ref)) > 0.9) ref = Vector3d(0,1,0);
    t1 = (ref - ref.dot(n)*n).normalized();
    t2 = n.cross(t1).normalized();
}
double MIQQuad::cotan(const Vector3d& a, const Vector3d& b, const Vector3d& c) {
    Vector3d u = a-b, v = c-b;
    double d = u.dot(v), cr = u.cross(v).norm();
    if (cr < 1e-12) return 0.0;
    return d / cr;
}

// ============================================================
// Phase 1a: Initialize cross field from geometry
// ============================================================
void MIQQuad::initCrossField() {
    faceN.resize(nF, 3);
    theta.resize(nF); theta.setZero();

    for (int fi = 0; fi < nF; fi++) {
        int v0=faces[fi*3], v1=faces[fi*3+1], v2=faces[fi*3+2];
        Vector3d p0=vertPos.row(v0), p1=vertPos.row(v1), p2=vertPos.row(v2);
        Vector3d fn = (p1-p0).cross(p2-p0).normalized();
        faceN.row(fi) = fn;

        Vector3d t1,t2; buildLocalFrame(fn,t1,t2);
        // Use longest edge as reference direction
        Vector3d e0=p1-p0, e1=p2-p1, e2=p0-p2;
        double l0=e0.squaredNorm(), l1=e1.squaredNorm(), l2=e2.squaredNorm();
        Vector3d d;
        if (l0>=l1&&l0>=l2) d=e0; else if(l1>=l0&&l1>=l2) d=e1; else d=e2;
        double dx=d.dot(t1), dy=d.dot(t2);
        theta(fi) = atan2(dy, dx);
    }
}

// ============================================================
// Phase 1b: Face graph Laplacian
// ============================================================
void MIQQuad::buildFaceLaplacian(SparseMatrix<double>& L) {
    L.resize(nF, nF);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nF);

    for (auto& e : edgeList) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        trips.push_back(Triplet<double>(e.f1, e.f2, -1.0));
        trips.push_back(Triplet<double>(e.f2, e.f1, -1.0));
        diag(e.f1)++; diag(e.f2)++;
    }
    for (int i = 0; i < nF; i++)
        trips.push_back(Triplet<double>(i, i, diag(i) + 1e-8));
    L.setFromTriplets(trips.begin(), trips.end());
}

// ============================================================
// Phase 1c: Alternating MIQP solver
// ============================================================
bool MIQQuad::optimizeCrossField(int maxIter) {
    const double K = M_PI / 2.0; // quarter-turn

    // Build face Laplacian once
    SparseMatrix<double> Lf;
    buildFaceLaplacian(Lf);

    // Precompute LDLT
    SimplicialCholesky<SparseMatrix<double>> solver;
    solver.compute(Lf);
    if (solver.info() != Success) return false;

    for (int iter = 0; iter < maxIter; iter++) {
        // Step A: Fix integer jumps p, solve for continuous θ
        VectorXd b = VectorXd::Zero(nF);
        for (auto& e : edgeList) {
            if (e.f1 < 0 || e.f2 < 0) continue;
            // Contribution: -K·p_ij from each adjacent face
            // b_i += K·p_ij for edge (i,j) where θ_i - θ_j ≈ -K·p_ij
            // L·θ = b where b_i = Σ_j K·p_ij (with orientation)
            // Actually: θ_i - θ_j = -K·p_ij  →  b_i = Σ_j (-K·p_ij)
            // But careful: L·θ = b resolves θ_i with fixed differences
            // b_i = Σ_j K·p_ij where p_ij = p for edge (i→j), -p for (j→i)
            int f1 = e.f1, f2 = e.f2;
            // Convention: p_ij is jump from f_i to f_j, so θ_j = θ_i + K·p_ij
            // → θ_i·deg(i) - Σ_j θ_j = -K Σ_j p_ij  → b_i = -K·Σ_j p_ij
            // Actually let me reconsider:
            // Energy: (θ_i - θ_j + K·p_ij)²
            // ∂/∂θ_i: 2(θ_i - θ_j + K·p_ij) = 0  → θ_i - θ_j = -K·p_ij
            // For fixed p, this is: L·θ = b where b_i = -K Σ_j p_ij
            // Note: p_ij from i→j. From j→i, p_ji = -p_ij so the signs work out.
            
            double val = K * jump[&e - &edgeList[0]];
            b(f1) += val;
            b(f2) -= val; // opposite direction
        }

        VectorXd newTheta = solver.solve(b);

        // Step B: Fix θ, round integer jumps
        for (auto& e : edgeList) {
            if (e.f1 < 0 || e.f2 < 0) continue;
            int idx = (int)(&e - &edgeList[0]);
            double diff = newTheta(e.f2) - newTheta(e.f1);
            jump[idx] = roundToInt(diff / K);
        }

        // Check convergence
        double change = (newTheta - theta).norm();
        theta = newTheta;
        if (change < 1e-6) break;
    }

    // Normalize θ to [0, π/2)
    for (int i = 0; i < nF; i++) {
        theta(i) = fmod(theta(i), K);
        if (theta(i) < 0) theta(i) += K;
    }

    return true;
}

// ============================================================
// Phase 1d: Build target direction vectors from θ
// ============================================================
void MIQQuad::buildTargetDirs() {
    faceD1.resize(nF, 3); faceD2.resize(nF, 3);

    for (int fi = 0; fi < nF; fi++) {
        Vector3d n = faceN.row(fi);
        Vector3d t1, t2; buildLocalFrame(n, t1, t2);

        double c = cos(theta(fi)), s = sin(theta(fi));
        faceD1.row(fi) = ( c*t1 + s*t2).normalized();
        faceD2.row(fi) = (-s*t1 + c*t2).normalized();
    }
}

// ============================================================
// Phase 2: Cotangent Laplacian (vertex-based)
// ============================================================
void MIQQuad::buildCotLaplacian(SparseMatrix<double>& L) {
    L.resize(nV, nV);
    std::vector<Triplet<double>> trips;
    VectorXd diag = VectorXd::Zero(nV);

    for (auto& e : edgeList) {
        if (e.f1 < 0 || e.f2 < 0) continue;
        int vi = e.v1, vj = e.v2;

        // Find opposite vertices
        int o1 = -1, o2 = -1;
        for (int k = 0; k < 3; k++) {
            int v = faces[e.f1*3+k];
            if (v != vi && v != vj) { o1 = v; break; }
        }
        for (int k = 0; k < 3; k++) {
            int v = faces[e.f2*3+k];
            if (v != vi && v != vj) { o2 = v; break; }
        }

        double w = 0.5 * (cotan(vertPos.row(o1), vertPos.row(vi), vertPos.row(vj)) +
                          cotan(vertPos.row(o2), vertPos.row(vi), vertPos.row(vj)));
        if (w < 0) w = 0;

        diag(vi) += w; diag(vj) += w;
        trips.push_back(Triplet<double>(vi, vj, -w));
        trips.push_back(Triplet<double>(vj, vi, -w));
    }
    for (int i = 0; i < nV; i++)
        trips.push_back(Triplet<double>(i, i, diag(i) + 1e-8));
    L.setFromTriplets(trips.begin(), trips.end());
}

// ============================================================
// Phase 2b: Poisson solve for UV
// ============================================================
void MIQQuad::solvePoisson() {
    // 1. Interpolate face directions to vertices (area-weighted)
    MatrixXd Vd1(nV,3), Vd2(nV,3); Vd1.setZero(); Vd2.setZero();
    VectorXd vA = VectorXd::Zero(nV);

    for (int fi = 0; fi < nF; fi++) {
        int v0=faces[fi*3], v1=faces[fi*3+1], v2=faces[fi*3+2];
        Vector3d a=vertPos.row(v0),b=vertPos.row(v1),c=vertPos.row(v2);
        double area=0.5*(b-a).cross(c-a).norm();
        for (int k=0;k<3;k++) {
            int v=faces[fi*3+k]; vA(v)+=area;
            Vd1.row(v)+=area*faceD1.row(fi); Vd2.row(v)+=area*faceD2.row(fi);
        }
    }
    for (int i=0;i<nV;i++) if(vA(i)>1e-12) { Vd1.row(i)/=vA(i); Vd2.row(i)/=vA(i); }

    // 2. Discrete divergence
    VectorXd div1=VectorXd::Zero(nV), div2=VectorXd::Zero(nV);
    for (int fi=0;fi<nF;fi++){
        int v0=faces[fi*3],v1=faces[fi*3+1],v2=faces[fi*3+2];
        Vector3d p0=vertPos.row(v0),p1=vertPos.row(v1),p2=vertPos.row(v2);
        Vector3d fn=(p1-p0).cross(p2-p0);
        Vector3d e01=p1-p0,e12=p2-p1,e20=p0-p2;
        Vector3d en01=fn.cross(e01).normalized(),en12=fn.cross(e12).normalized(),en20=fn.cross(e20).normalized();
        double cot0=cotan(p2,p0,p1),cot1=cotan(p0,p1,p2),cot2=cotan(p1,p2,p0);
        Vector3d fd1=faceD1.row(fi),fd2=faceD2.row(fi);
        double d1=0.5*(fd1.dot(en01)*cot2+fd1.dot(en12)*cot0+fd1.dot(en20)*cot1);
        double d2=0.5*(fd2.dot(en01)*cot2+fd2.dot(en12)*cot0+fd2.dot(en20)*cot1);
        div1(v0)+=d1;div1(v1)+=d1;div1(v2)+=d1;
        div2(v0)+=d2;div2(v1)+=d2;div2(v2)+=d2;
    }

    // 3. Build Laplacian + boundary conditions
    SparseMatrix<double> L;
    buildCotLaplacian(L);

    // Fix boundary vertices via penalty
    for (VertexCIter v=mesh.vertices.begin();v!=mesh.vertices.end();v++)
        if(v->isBoundary()){ L.coeffRef(v->index,v->index)+=1e6; div1(v->index)=0; div2(v->index)=0; }
    // If closed mesh, fix vertex 0
    bool hasBnd=false;
    for(VertexCIter v=mesh.vertices.begin();v!=mesh.vertices.end();v++) if(v->isBoundary()){hasBnd=true;break;}
    if(!hasBnd){L.coeffRef(0,0)+=1e6;div1(0)=0;div2(0)=0;}

    // 4. Solve
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(L);
    if(solver.info()!=Success){UV=MatrixXd::Zero(nV,2);return;}
    VectorXd u=solver.solve(div1), v=solver.solve(div2);

    // 5. Normalize
    double uMin=u.minCoeff(),uMax=u.maxCoeff(),vMin=v.minCoeff(),vMax=v.maxCoeff();
    double uR=uMax-uMin,vR=vMax-vMin;
    if(uR<1e-10)uR=1.0;if(vR<1e-10)vR=1.0;
    UV.resize(nV,2);
    for(int i=0;i<nV;i++){UV(i,0)=(u(i)-uMin)/uR;UV(i,1)=(v(i)-vMin)/vR;}
}

// ============================================================
// Main entry point
// ============================================================
void MIQQuad::parameterize() {
    // Extract mesh data
    nV=(int)mesh.vertices.size();nF=0;
    for(FaceCIter f=mesh.faces.begin();f!=mesh.faces.end();f++)
        if(!f->isBoundary())nF++;
    vertPos.resize(nV,3);faces.resize(nF*3);
    for(VertexCIter v=mesh.vertices.begin();v!=mesh.vertices.end();v++)
        vertPos.row(v->index)=v->position;
    int fi=0;
    for(FaceCIter f=mesh.faces.begin();f!=mesh.faces.end();f++){
        if(!f->isBoundary()){
            faces[fi*3]=f->he->vertex->index;
            faces[fi*3+1]=f->he->next->vertex->index;
            faces[fi*3+2]=f->he->next->next->vertex->index;
            fi++;
        }
    }

    // Build edge list
    std::map<std::pair<int,int>,int> em;
    for(int fi2=0;fi2<nF;fi2++){
        for(int k=0;k<3;k++){
            int v1=faces[fi2*3+k],v2=faces[fi2*3+(k+1)%3];
            if(v1>v2)std::swap(v1,v2);
            auto key=std::make_pair(v1,v2);
            auto it=em.find(key);
            if(it==em.end()){em[key]=(int)edgeList.size();Edge e;e.v1=v1;e.v2=v2;e.f1=fi2;e.f2=-1;edgeList.push_back(e);}
            else edgeList[it->second].f2=fi2;
        }
    }
    nE=(int)edgeList.size();
    jump.resize(nE);jump.setZero();

    // Phase 1: Cross field optimization (MIQP)
    initCrossField();
    optimizeCrossField(8);
    buildTargetDirs();

    // Phase 2: Poisson UV parameterization
    solvePoisson();

    // Copy UV back
    for(VertexIter v=mesh.vertices.begin();v!=mesh.vertices.end();v++)
        v->uv=Vector2d(UV(v->index,0),UV(v->index,1));
}
