#include "HolomorphicOneForm.h"
#include <map>
#include <queue>
using namespace Eigen;
using namespace std;

HolomorphicOneForm::HolomorphicOneForm(Mesh& mesh0): Parameterization(mesh0) {}

double HolomorphicOneForm::cotan(const Vector3d& a, const Vector3d& b, const Vector3d& c) {
    Vector3d u=a-b, v=c-b; double d=u.dot(v), cr=u.cross(v).norm();
    if(cr<1e-12)return 0.0; return d/cr;
}
void HolomorphicOneForm::localFrame(const Vector3d& n, Vector3d& t1, Vector3d& t2) {
    Vector3d ref(1,0,0); if(fabs(n.dot(ref))>0.9)ref=Vector3d(0,1,0);
    t1=(ref-ref.dot(n)*n).normalized(); t2=n.cross(t1).normalized();
}

// ============================================================
// Phase 1: Build |E|×|E| system: [closedness; harmonicity; cut]
// ============================================================
void HolomorphicOneForm::buildClosedness(SparseMatrix<double>& A, int& off) {
    for(auto& e:edgeList){ if(e.f1<0)continue; // skip boundary
        for(int k=0;k<3;k++) {
            int v0=F[e.f1*3+k], v1=F[e.f1*3+(k+1)%3];
            int eIdx=-1;
            for(int ei=0;ei<nE;ei++){
                if(edgeList[ei].f1>=0&&((edgeList[ei].v1==v0&&edgeList[ei].v2==v1)||(edgeList[ei].v1==v1&&edgeList[ei].v2==v0))){eIdx=ei;break;}
                if(edgeList[ei].f2>=0&&((edgeList[ei].v1==v0&&edgeList[ei].v2==v1)||(edgeList[ei].v1==v1&&edgeList[ei].v2==v0))){eIdx=ei;break;}
            }
            int sign=(edgeList[eIdx].v1==v0&&edgeList[eIdx].v2==v1)?1:-1;
            if(v0>v1)sign=-sign;
            A.coeffRef(off,eIdx)+=sign;
        }
        off++;
    }
}
void HolomorphicOneForm::buildHarmonicity(SparseMatrix<double>& A, int& off) {
    for(int vi=0;vi<nV-1;vi++){ // last vertex redundant
        // Sum of cot-weighted oriented edges around vi
        for(auto& e:edgeList){
            if(e.f1<0)continue;
            int sign=0;
            // Check if e is incident to vi with correct orientation
            if(e.v1==vi||e.v2==vi){
                // edge orientation from v1→v2
                // For cot weights, need to find the two opposite vertices
                int o1=-1,o2=-1;
                for(int k=0;k<3;k++){int vv=F[e.f1*3+k];if(vv!=e.v1&&vv!=e.v2){o1=vv;break;}}
                if(e.f2>=0) for(int k=0;k<3;k++){int vv=F[e.f2*3+k];if(vv!=e.v1&&vv!=e.v2){o2=vv;break;}}
                
                double w=0;
                if(o1>=0) w+=cotan(V.row(o1),V.row(e.v1),V.row(e.v2));
                if(o2>=0) w+=cotan(V.row(o2),V.row(e.v1),V.row(e.v2));
                
                int ei=(int)(&e-&edgeList[0]);
                if(e.v1==vi) A.coeffRef(off,ei)+=w; // outflow
                else A.coeffRef(off,ei)-=w; // inflow
            }
        }
        off++;
    }
}
void HolomorphicOneForm::fixCutConstraint(SparseMatrix<double>& A) {
    // Fix the integral along cut edge to 1.0
    A.coeffRef(nF-1+nV-1, cutEdge) = 1.0;
}

// ============================================================
// Phase 2: Integrate 1-form → UV coordinates (BFS)
// ============================================================
void HolomorphicOneForm::integrateUV(const VectorXd& real, const VectorXd& imag) {
    UV.resize(nV,2); UV.setZero();
    vector<bool> visited(nV,false);
    queue<int> q;
    q.push(0); visited[0]=true; UV(0,0)=0; UV(0,1)=0;

    while(!q.empty()){
        int vi=q.front();q.pop();
        for(auto& e:edgeList){
            if(e.f1<0)continue;
            int ei=(int)(&e-&edgeList[0]);
            int other=-1;
            if(e.v1==vi)other=e.v2; else if(e.v2==vi)other=e.v1; else continue;
            if(visited[other])continue;
            visited[other]=true;
            // Orientation: edge length × form value
            int sign=(e.v1==vi&&e.v2==other)?1:-1;
            UV(other,0)=UV(vi,0)+sign*real(ei);
            UV(other,1)=UV(vi,1)+sign*imag(ei);
            q.push(other);
        }
    }

    // Normalize
    double uMin=UV.col(0).minCoeff(), uMax=UV.col(0).maxCoeff();
    double vMin=UV.col(1).minCoeff(), vMax=UV.col(1).maxCoeff();
    double uR=uMax-uMin, vR=vMax-vMin;
    if(uR<1e-10)uR=1.0;if(vR<1e-10)vR=1.0;
    for(int i=0;i<nV;i++){ UV(i,0)=(UV(i,0)-uMin)/uR; UV(i,1)=(UV(i,1)-vMin)/vR; }
}

// ============================================================
// Main entry
// ============================================================
void HolomorphicOneForm::parameterize() {
    // Extract mesh
    nV=(int)mesh.vertices.size(); nF=0;
    for(FaceCIter f=mesh.faces.begin();f!=mesh.faces.end();f++) if(!f->isBoundary())nF++;
    V.resize(nV,3); F.resize(nF*3);
    for(VertexCIter v=mesh.vertices.begin();v!=mesh.vertices.end();v++) V.row(v->index)=v->position;
    int fi=0;
    for(FaceCIter f=mesh.faces.begin();f!=mesh.faces.end();f++){
        if(!f->isBoundary()){F[fi*3]=f->he->vertex->index;F[fi*3+1]=f->he->next->vertex->index;F[fi*3+2]=f->he->next->next->vertex->index;fi++;}
    }

    // Build edge list
    map<pair<int,int>,int> em;
    for(int fi2=0;fi2<nF;fi2++) for(int k=0;k<3;k++){
        int v1=F[fi2*3+k],v2=F[fi2*3+(k+1)%3];if(v1>v2)swap(v1,v2);
        auto it=em.find({v1,v2});
        if(it==em.end()){em[{v1,v2}]=(int)edgeList.size();Edge e;e.v1=v1;e.v2=v2;e.f1=fi2;e.f2=-1;e.len=(V.row(v1)-V.row(v2)).norm();edgeList.push_back(e);}
        else edgeList[it->second].f2=fi2;
    }
    nE=(int)edgeList.size();

    // Pick cut edge: longest edge
    cutEdge=0; double best=0;
    for(int i=0;i<nE;i++) if(edgeList[i].len>best){best=edgeList[i].len;cutEdge=i;}

    // Build system: A·ω = b
    SparseMatrix<double> A(nE,nE);
    vector<Triplet<double>> trips;
    int off=0;
    buildClosedness(A,off);         // |F|-1 rows
    buildHarmonicity(A,off);        // |V|-1 rows
    fixCutConstraint(A);            // 1 row: ω[cutEdge] = 1

    // RHS
    VectorXd b=VectorXd::Zero(nE);
    b(nF-1+nV-1)=1.0; // cut constraint

    // Solve real part
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(A);
    if(solver.info()!=Success){ UV.setZero(nV,2); return; }
    VectorXd real=solver.solve(b);

    // Compute conjugate via Hodge star: rotate edge values 90° in tangent space
    VectorXd imag=VectorXd::Zero(nE);
    for(int ei=0;ei<nE;ei++){
        auto& e=edgeList[ei];
        if(e.f1<0)continue;
        // Get face normal and compute 90° rotation
        int v0=F[e.f1*3],v1=F[e.f1*3+1],v2=F[e.f1*3+2];
        Vector3d p0=V.row(v0),p1=V.row(v1),p2=V.row(v2);
        Vector3d fn=(p1-p0).cross(p2-p0).normalized();
        Vector3d ed=V.row(e.v2)-V.row(e.v1);
        Vector3d rot=fn.cross(ed);
        // Project rotated vector back to edge direction to get scalar scaling
        double s=rot.dot(ed)/ed.squaredNorm();
        if(e.v1>e.v2)s=-s; // orientation
        imag(ei)=s*real(ei);
    }

    // Integrate for UV
    integrateUV(real,imag);

    // Copy UV to mesh
    for(VertexIter v=mesh.vertices.begin();v!=mesh.vertices.end();v++)
        v->uv=Vector2d(UV(v->index,0),UV(v->index,1));
}
