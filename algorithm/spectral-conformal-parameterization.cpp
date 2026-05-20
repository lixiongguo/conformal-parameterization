// Implement member functions for SpectralConformalParameterization class.
#include "spectral-conformal-parameterization.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
SpectralConformalParameterization::SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh,
                                                                     VertexPositionGeometry* inputGeo) {

    this->mesh = inputMesh;
    this->geometry = inputGeo;
}

/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {
    
    // TODO
  SparseMatrix<std::complex<double>> Ec = 0.5 * geometry->complexLaplaceMatrix();

  std::vector<Eigen::Triplet<std::complex<double>>> A_entries;
  
  for(BoundaryLoop bl : mesh->boundaryLoops())
    for (Halfedge he : bl.adjacentHalfedges())
    {
      A_entries.push_back(
        Eigen::Triplet<std::complex<double>>(
          he.tailVertex().getIndex(),
          he.tipVertex().getIndex(), 
          std::complex<double>(0, 0.25)
          )
      );
      A_entries.push_back(
        Eigen::Triplet<std::complex<double>>(
          he.tipVertex().getIndex(),
          he.tailVertex().getIndex(),
          std::complex<double>(0, -0.25)
          )
      );
    }

  Eigen::SparseMatrix<std::complex<double>> A(Ec.rows(), Ec.cols());
  A.setFromTriplets(A_entries.begin(), A_entries.end());

  Ec -= A;

  return Ec; // placeholder
}


/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {

    // TODO
    SparseMatrix<std::complex<double>> Ec = buildConformalEnergy();

    // Solve the eigenvalue problem with the smallest eigenvalue using solveInversePowerMethod
    Vector<std::complex<double>> eig_vec;
    eig_vec = solveInversePowerMethod(Ec);

    VertexData<Vector2> result(*mesh);

    for (Vertex v: mesh->vertices())
    {
      std::complex<double> result_cmplx = eig_vec[v.getIndex()];
      result[v].x = result_cmplx.real();
      result[v].y = result_cmplx.imag();
    }

    return result; // placeholder
}

// 使用边界循环的相邻半边来获取有序的边界顶点
std::vector<Vertex>SpectralConformalParameterization:: getBoundaryVertices()const {
    std::vector<Vertex> orderedBoundaryVertices;
    for (BoundaryLoop bl : mesh->boundaryLoops()) {
        std::vector<Vertex> loopVertices;
        for (Halfedge he : bl.adjacentHalfedges()) {
            Vertex v = he.tailVertex(); 
            loopVertices.push_back(v);
        }
        orderedBoundaryVertices.insert(orderedBoundaryVertices.end(), loopVertices.begin(), loopVertices.end());
        break;//只需要第一个边界循环，可以在这里break
    }
    
    return orderedBoundaryVertices;
}

VertexData<Vector2> SpectralConformalParameterization::tutte_flatten() const {

    VertexData<Vector2> result(*mesh);
    //直接拍平Flatten
    //for (Vertex v : mesh->vertices())
    //{
    //    Vector3 pos = geometry->inputVertexPositions[v];
    //    result[v].x = pos.x;
    //    result[v].y = pos.y;
    //}
    // 
    // 
    //for (Vertex v : mesh->vertices()) {
    //    Vector3 pos = geometry->inputVertexPositions[v];
    //    bool bisBoundary = false;
    //   /* for (Edge e : v.adjacentEdges()) {
    //        if (e.isBoundary()) {
    //            isBoundary = true;
    //            break;
    //        }
    //    }*/
    //    bisBoundary = isBoundary[v];
    //    if (!bisBoundary) {
    //        result[v].x = pos.x;
    //        result[v].y = pos.y;
    //    }
    //}

    //1.确定边界
    std::vector<Vertex> boundaryVertices = getBoundaryVertices();
    VertexData<bool> isBoundary(*mesh, false);
    VertexData<Vector2> boundaryPos(*mesh);
    size_t boundaryCount = boundaryVertices.size();

    for (size_t i = 0; i < boundaryCount; i++) {
        Vertex v = boundaryVertices[i];
        double angle = 2.0 * M_PI * i / boundaryCount;
        boundaryPos[v].x = 10 * cos(angle);
        boundaryPos[v].y = 10 * sin(angle);
        isBoundary[v] = true;
    }
    //2.建立方程组Ax = b
    size_t nVerteices = mesh->nVertices();
    SparseMatrix<double> A(nVerteices,nVerteices);
    std::vector<Eigen::Triplet<double>> A_entries;
    Vector<double> bx(nVerteices);
    Vector<double> by(nVerteices);
    for (Vertex v : mesh->vertices())
    {
        int vIdx = v.getIndex();
        bool bisBoundary = isBoundary[v];
        if (bisBoundary)
        {
            A_entries.push_back(
                Eigen::Triplet<double>(vIdx, vIdx, 1.0)
            );
            bx[vIdx] =  boundaryPos[v].x;
            by[vIdx] =  boundaryPos[v].y;
        }
        else
        {
           size_t Nv = 0;
           for (auto w : v.adjacentVertices())
           {
              int wIdx = w.getIndex();
              A_entries.push_back(
                  Eigen::Triplet<double>(vIdx, wIdx,1.0)
              );
              Nv++;
           }
           A_entries.push_back(Eigen::Triplet<double>(vIdx,vIdx,-static_cast<double>(Nv)));
           bx[vIdx] = 0;
           by[vIdx] = 0;
        }
       
    }
    A.setFromTriplets(A_entries.begin(), A_entries.end());
    //3. Solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "LU 分解失败" << std::endl;
        return result;
    }
    Eigen::VectorXd x = solver.solve(bx);
    Eigen::VectorXd y = solver.solve(by);
    for (Vertex v : mesh->vertices())
    {
        result[v].x = x[v.getIndex()];
        result[v].y = y[v.getIndex()];
    }

    return result; 
}