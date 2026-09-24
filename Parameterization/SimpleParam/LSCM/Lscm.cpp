#include "Lscm.h"
#include <Eigen/SparseCholesky>

Lscm::Lscm(Mesh& mesh0):
Parameterization(mesh0),
lastRunOk(false)
{
    
}

// 由调用方(wasm 绑定 / 前端「改pin点」)指定 pin 顶点。
// 内部统一以 2*vertexIndex 作为无约束列的下标基准(见 isPinnedVertex / buildMassMatrix),
// 所以这里同样存 2*v。pin 的目标 UV 与自动选点时保持一致:(-0.5,0) 与 (0.5,0)。
void Lscm::setPins(int v0, int v1)
{
    if (v0 < 0 || v1 < 0 || v0 == v1) return;   // 非法输入:保持未设置,交给 parameterize() 自动选点

    pinnedPositions = {Eigen::Vector2d(-0.5, 0), Eigen::Vector2d(0.5, 0)};
    pinnedVertices.assign(2, -1);
    pinnedVertices[0] = 2*v0;
    pinnedVertices[1] = 2*v1;
}

void Lscm::pinVertices()
{
    pinnedVertices.assign(2, -1);
    pinnedPositions = {Eigen::Vector2d(-0.5, 0), Eigen::Vector2d(0.5, 0)};

    // 注意:比较阈值初值必须是负数,不能是 0.0。
    //
    // 原先用 0.0 作初值:一旦边界上任何一对顶点的距离平方都不 > 0(边界顶点重合、
    // 边界环退化、或 boundaries 里只有单个顶点),pinnedVertices 会保持 {0, 0}。
    // 此时 isPinnedVertex() 对几乎所有顶点都会累加 shift(变成 4),于是
    //     vIndex = indices[i] - shift
    // 得到负数,再交给 M.setFromTriplets();或者 fIndex 超出 b 的长度。
    // 这类越界写不会当场崩溃(写进的仍是已映射的线性内存),但会破坏 dlmalloc
    // 的堆元数据,之后任意一次 free()/malloc() 才抛出
    //     RuntimeError: memory access out of bounds
    // 所以必须在这里就杜绝"pin 未赋值"的可能。
    double max = -1.0;

    // pin diameter vertices on longest boundary loop
    for (std::vector<HalfEdgeIter>::const_iterator it1 = mesh.boundaries.begin();
                                                   it1 != mesh.boundaries.end();
                                                   it1++) {
        HalfEdgeCIter he1 = *it1;
        do {
            const int& vIdx1(he1->vertex->index);
            const Eigen::Vector3d& p1(he1->vertex->position);

            for (std::vector<HalfEdgeIter>::const_iterator it2 = mesh.boundaries.begin();
                                                           it2 != mesh.boundaries.end();
                                                           it2++) {
                HalfEdgeCIter he2 = *it2;
                do {
                    const int& vIdx2(he2->vertex->index);
                    const Eigen::Vector3d& p2(he2->vertex->position);
                    double l = (p2-p1).squaredNorm();

                    if (l > max) {
                        max = l;
                        pinnedVertices[0] = 2*vIdx1;
                        pinnedVertices[1] = 2*vIdx2;
                    }

                    he2 = he2->next;
                } while (he2 != *it2);
            }

            he1 = he1->next;
        } while (he1 != *it1);
    }

    // 没有边界,或边界搜索没能给出可用 pin(退化边界):退化为"顶点 0 + 离它最远的顶点"
    if (pinnedVertices[0] < 0) {
        pinnedVertices[0] = 0;
        if (!mesh.vertices.empty()) {
            const Eigen::Vector3d& p0(mesh.vertices[0].position);
            max = -1.0;
            for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
                double l = (v->position - p0).squaredNorm();
                if (l > max) {
                    max = l;
                    pinnedVertices[1] = 2*v->index;
                }
            }
        }
    }

    // 最后兜底:两个 pin 必须存在、下标合法且互不相同
    const int nV2 = 2*static_cast<int>(mesh.vertices.size());
    if (pinnedVertices[0] < 0 || pinnedVertices[0] >= nV2) pinnedVertices[0] = 0;
    if (pinnedVertices[1] < 0 || pinnedVertices[1] >= nV2 ||
        pinnedVertices[1] == pinnedVertices[0]) {
        pinnedVertices[1] = (pinnedVertices[0] == 0 && nV2 > 2) ? 2 : 0;
    }
}

void computeLocalBasisCoordinates(std::vector<Eigen::Vector2d>& coords,
                                  const std::vector<Eigen::Vector3d>& positions)
{
    Eigen::Vector3d x = positions[1] - positions[0];
    Eigen::Vector3d y = positions[2] - positions[0];
    
    // compute orthonormal basis
    Eigen::Vector3d xhat = x; xhat.normalize();
    Eigen::Vector3d zhat = xhat.cross(y); zhat.normalize();
    Eigen::Vector3d yhat = zhat.cross(xhat); yhat.normalize();
    
    // compute coordinates in local basis
    coords.push_back(Eigen::Vector2d(0, 0));
    coords.push_back(Eigen::Vector2d(x.norm(), 0));
    coords.push_back(Eigen::Vector2d(y.dot(xhat), y.dot(yhat)));
}

bool Lscm::isPinnedVertex(const int& vIndex, int& shift, Eigen::Vector2d& pinnedPosition) const
{
    shift = 0;
    for (size_t i = 0; i < pinnedVertices.size(); i++) {
        if (pinnedVertices[i] == vIndex) {
            pinnedPosition = pinnedPositions[i];
            return true;
        }
        
        if (pinnedVertices[i] < vIndex) shift += 2;
    }
    
    return false;
}

bool Lscm::buildMassMatrix(Eigen::SparseMatrix<double>& M, Eigen::VectorXd& b) const
{
    std::vector<Eigen::Triplet<double>> MTriplets;

    const int nRows = static_cast<int>(M.rows());
    const int nCols = static_cast<int>(M.cols());
    
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fIndex = 2*f->index;

            // 行索引越界防护:setFromTriplets / b() 越界会静默写坏堆,见 pinVertices() 注释
            if (fIndex < 0 || fIndex + 1 >= nRows) return false;
            
            // get indices and positions of adjacent vertices
            std::vector<int> indices;
            std::vector<Eigen::Vector3d> positions;
            HalfEdgeCIter h = f->he;
            do {
                indices.push_back(2*h->vertex->index);
                positions.push_back(h->vertex->position);
                h = h->next;
                
            } while (h != f->he);

            // 少于 3 个角无法构造局部基
            if (indices.size() < 3 || positions.size() < 3) return false;
            
            // compute weights
            std::vector<Eigen::Vector2d> coords;
            computeLocalBasisCoordinates(coords, positions);
            std::vector<Eigen::Vector2d> ws = {coords[2] - coords[1],
                                               coords[0] - coords[2],
                                               coords[1] - coords[0]};
            
            // set M entries
            for (int i = 0; i < 3; i++) {
                int shift;
                Eigen::Vector2d pinnedPosition;
                const Eigen::Vector2d& w(ws[i]);
                
                if (isPinnedVertex(indices[i], shift, pinnedPosition)) {
                    b(fIndex) -= (w.x()*pinnedPosition.x() - w.y()*pinnedPosition.y());
                    b(fIndex+1) -= (w.y()*pinnedPosition.x() + w.x()*pinnedPosition.y());
                    
                } else {
                    int vIndex = indices[i] - shift;

                    // 列索引越界防护(负数下标是历史崩溃的元凶)
                    if (vIndex < 0 || vIndex + 1 >= nCols) return false;
                    
                    // set real components
                    MTriplets.push_back(Eigen::Triplet<double>(fIndex, vIndex, w.x()));
                    MTriplets.push_back(Eigen::Triplet<double>(fIndex+1, vIndex+1, w.x()));
                    
                    if (i != 2) {
                        // set imaginary components
                        MTriplets.push_back(Eigen::Triplet<double>(fIndex, vIndex+1, -w.y()));
                        MTriplets.push_back(Eigen::Triplet<double>(fIndex+1, vIndex, w.y()));
                    }
                }
            }
        }
    }
    
    M.setFromTriplets(MTriplets.begin(), MTriplets.end());
    return true;
}

void solveLeastSquares(const Eigen::SparseMatrix<double>& A,
                       const Eigen::VectorXd& b,
                       Eigen::VectorXd& x)
{
    // initialize
    Eigen::SparseMatrix<double> At = A.transpose();
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(At * A);
    Eigen::VectorXd Atb = At * b;
    
    // backsolve
    x = solver.solve(Atb);
}

void Lscm::setUvs(const Eigen::VectorXd& x)
{
    const int n = static_cast<int>(x.size());

    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int vIndex = 2*v->index;
        
        int shift;
        Eigen::Vector2d pinnedPosition;
        if (isPinnedVertex(vIndex, shift, pinnedPosition)) {
            v->uv(0) = pinnedPosition.x();
            v->uv(1) = pinnedPosition.y();
            
        } else {
            vIndex -= shift;
            if (vIndex < 0 || vIndex + 1 >= n) {   // 防护:绝不越界读
                v->uv.setZero();
                continue;
            }
            v->uv(0) = x(vIndex);
            v->uv(1) = x(vIndex+1);
        }
    }
    
    normalize();
}

void Lscm::parameterize()
{
    lastRunOk = false;

    // 退化网格直接放弃,不要进入下面的 Eigen 写入路径
    if (mesh.vertices.size() < 3 || mesh.faces.size() <= mesh.boundaries.size()) {
        return;
    }

    // pin vertices:调用方已用 setPins() 指定时,跳过自动选点
    if (pinnedVertices.size() != 2) {
        pinVertices();
    }

    // 校验 pin:必须是两个互异且下标合法的顶点,否则 isPinnedVertex 的 shift 语义不成立
    const int nV2 = 2*static_cast<int>(mesh.vertices.size());
    if (pinnedVertices.size() != 2 ||
        pinnedVertices[0] == pinnedVertices[1] ||
        pinnedVertices[0] < 0 || pinnedVertices[0] >= nV2 ||
        pinnedVertices[1] < 0 || pinnedVertices[1] >= nV2) {
        return;
    }

    const int nFaces = static_cast<int>(mesh.faces.size());
    const int nBnd   = static_cast<int>(mesh.boundaries.size());
    int f = 2*(nFaces - nBnd);
    int v = 2*((int)mesh.vertices.size() - 2);
    if (f <= 0 || v <= 0) {
        return;
    }
    
    // build mass matrix and handle boundary conditions
    Eigen::SparseMatrix<double> M(f, v);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(f);
    if (!buildMassMatrix(M, b)) {
        return;
    }
    
    // solve least squares
    Eigen::VectorXd x(v);
    solveLeastSquares(M, b, x);
    
    // set uv coords
    setUvs(x);

    lastRunOk = true;
}
