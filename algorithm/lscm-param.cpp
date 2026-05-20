#include "spectral-conformal-parameterization.h"




// 计算每个面的局部坐标系
// 对于每个三角面，建立一个局部的2D坐标系，其中一个顶点在原点，另一个顶点在x轴上
std::map<Face, std::array<Vector2, 3>> computeLocalCoordinates(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry) {
    std::map<Face, std::array<Vector2, 3>> localCoords;
    
    for (Face f : mesh->faces()) {
        // 获取面的三个顶点
        std::array<Vertex, 3> vertices;
        size_t i = 0;
        for (Vertex v : f.adjacentVertices()) {
            vertices[i++] = v;
            if (i >= 3) break;
        }
        
        // 获取三个顶点的3D坐标
        Vector3 p0 = geometry->inputVertexPositions[vertices[0]];
        Vector3 p1 = geometry->inputVertexPositions[vertices[1]];
        Vector3 p2 = geometry->inputVertexPositions[vertices[2]];
        
        // 计算边向量
        Vector3 e0 = p1 - p0;
        Vector3 e1 = p2 - p0;
        
        // 在局部坐标系中，第一个顶点在原点(0,0)
        Vector2 local0{0.0, 0.0};
        // 第二个顶点在x轴上，距离为e0的长度
        double len0 = e0.norm();
        Vector2 local1{len0, 0.0};
        
        // 第三个顶点的位置需要通过几何计算得到
        // 使用点积计算角度
        double x_val = dot(e0, e1);
        double y_val = cross(e0, e1).norm();
        double angle = std::atan2(y_val, x_val);
        
        // 第三个顶点的坐标
        double len1 = e1.norm();
        Vector2 local2{ len1 * std::cos(angle), len1 * std::sin(angle) };
        
        // 存储局部坐标
        localCoords[f] = {local0, local1, local2};
    }
    
    return localCoords;
}


// 找到边界上距离最远的两个点
std::pair<Vertex, Vertex> findFarthestBoundaryPoints(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry) {
    // 获取边界顶点
    std::vector<Vertex> boundaryVertices;
    for (BoundaryLoop bl : mesh->boundaryLoops()) {
        for (Halfedge he : bl.adjacentHalfedges()) {
            Vertex v = he.tailVertex(); 
            boundaryVertices.push_back(v);
        }
        break; // 只处理第一个边界环
    }
    
    if (boundaryVertices.size() < 2) {
        // 如果边界点少于2个，返回前两个顶点（这种情况不应该发生）
        auto vertices = mesh->vertices();
        auto it = vertices.begin();
        Vertex v1 = *it;
        ++it;
        Vertex v2 = *it;
        return std::make_pair(v1, v2);
    }
    
    // 找到距离最远的两个点
    Vertex farthestV1 = boundaryVertices[0];
    Vertex farthestV2 = boundaryVertices[1];
    double maxDistance = 0.0;
    
    for (size_t i = 0; i < boundaryVertices.size(); i++) {
        for (size_t j = i + 1; j < boundaryVertices.size(); j++) {
            Vertex v1 = boundaryVertices[i];
            Vertex v2 = boundaryVertices[j];
            
            Vector3 pos1 = geometry->inputVertexPositions[v1];
            Vector3 pos2 = geometry->inputVertexPositions[v2];
            
            double distance = (pos1 - pos2).norm();
            if (distance > maxDistance) {
                maxDistance = distance;
                farthestV1 = v1;
                farthestV2 = v2;
            }
        }
    }
    
    return std::make_pair(farthestV1, farthestV2);
}



VertexData<Vector2> SpectralConformalParameterization::arap_flatten() const 
{
    VertexData<Vector2> result(*mesh);
    std::map<Face, std::array<Vector2, 3>> localCoords = computeLocalCoordinates(mesh, geometry);
    std::pair<Vertex, Vertex> pinPoints = findFarthestBoundaryPoints(mesh, geometry);
    Vertex pin1 = pinPoints.first;
    Vertex pin2 = pinPoints.second;
    std::map<Vertex, size_t> vertexIndex;
    size_t index = 0;
    for (Vertex v : mesh->vertices()) {
        if (v != pin1 && v != pin2) {
            vertexIndex[v] = index++;
        }
    }
    size_t nVertices = mesh->nVertices();
    size_t nFreeVertices = nVertices - 2; 
    // 初始化结果，先用局部坐标作为初始猜测
    for (Vertex v : mesh->vertices()) {
        Vector3 pos = geometry->inputVertexPositions[v];
        result[v].x = pos.x;
        result[v].y = pos.y;
    }

   

    //// 设置pin点的坐标：将两个最远点分别放在(0,0)和(d,0)，其中d是它们在3D空间中的距离
    //Vector3 pos1 = geometry->inputVertexPositions[pin1];
    //Vector3 pos2 = geometry->inputVertexPositions[pin2];
    //double distance = (pos1 - pos2).norm();
    //
    //result[pin1].x = 0.0;
    //result[pin1].y = 0.0;
    //result[pin2].x = distance;
    //result[pin2].y = 0.0;

    //const int maxIterations = 10;
    //for (int iter = 0; iter < maxIterations; iter++) {

    //  VertexData<Vector2> prevResult = result;
    //  std::map<Face, Eigen::Matrix2d> rotations;

    //  //local step : SVD
    //  for (Face f : mesh->faces()) {
    //      //获取三角面上的三个顶点
    //      std::array<Vertex, 3> vertices;
    //      size_t i = 0;
    //      for (Vertex v : f.adjacentVertices()) {
    //          vertices[i++] = v;
    //          if (i >= 3)break;
    //      }

    //      std::array<Vector3, 3> pos3D;
    //      std::array<Vector2, 3> pos2D;
    //      for (int j = 0; j < 3; j++) {
    //          pos3D[j] = geometry->inputVertexPositions[vertices[j]];
    //          pos2D[j] = result[vertices[j]];
    //      }
    //      // 计算3D和2D的边向量
    //        Vector3 e0_3D = pos3D[1] - pos3D[0];
    //        Vector3 e1_3D = pos3D[2] - pos3D[0];
    //        Vector2 e0_2D = pos2D[1] - pos2D[0];
    //        Vector2 e1_2D = pos2D[2] - pos2D[0];
    //        // 构建3D和2D的协方差矩阵
    //        Eigen::MatrixXd S_3D(3,2);
    //        S_3D << e0_3D.x, e0_3D.y,
    //                e1_3D.x, e1_3D.y,
    //                e0_3D.z, e1_3D.z;
    //        Eigen::Matrix2d S_2D;
    //        S_2D << e0_2D.x, e0_2D.y,
    //                e1_2D.x, e1_2D.y;
    //        // 计算协方差矩阵
    //        Eigen::MatrixXd C = S_3D * S_2D;
    //        
    //        Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
    //        Eigen::Matrix2d U = svd.matrixU();
    //        Eigen::Matrix3d V = svd.matrixV();
    //        Eigen::Matrix2d R = U * V.transpose();
    //        if (R.determinant() < 0) {
    //            U.col(1) *= -1;
    //            R = U * V.transpose();
    //        }
    //        rotations[f] = R;
    //    }
    
      ////global step  Ax = b
      //{
      //    SparseMatrix<double> A(2 * nFreeVertices, 2 * nFreeVertices);
      //    Vector<double> b(2 * nFreeVertices);
      //    std::vector<Eigen::Triplet<double>> A_entries;

      //    for (Face f : mesh->faces()) {
      //        //获取三角面上的三个顶点
      //        std::array<Vertex, 3> vertices;
      //        size_t i = 0;
      //        for (Vertex v : f.adjacentVertices()) {
      //            vertices[i++] = v;
      //            if (i >= 3)break;
      //        }
      //        // 获取局部坐标
      //        std::array<Vector2, 3> local = localCoords[f];
      //        Eigen::Matrix2d R = rotations[f];

      //        // 为每个边建立方程
      //        for (int j = 0; j < 3; j++)
      //        {
      //            int next_j = (j + 1) % 3;
      //            // 计算局部坐标的边向量
      //            Vector2 e_local = local[next_j] - local[j];
      //            // 应用旋转矩阵
      //            Vector2 e_rotated{ R(0,0) * e_local.x + R(0,1) * e_local.y,
      //                              R(1,0) * e_local.x + R(1,1) * e_local.y };

      //            // 计算cotangent权重
      //            double cot_weight = 0.5; // 简化处理，实际应该计算真实的cotangent权重
      //            // 为顶点j建立方程
      //            int k = next_j;
      //            if (vertices[j] != pin1 && vertices[j] != pin2) {
      //                size_t idx = vertexIndex[vertices[j]];
      //                A_entries.push_back(Eigen::Triplet<double>(2 * idx, 2 * idx, cot_weight));
      //                A_entries.push_back(Eigen::Triplet<double>(2 * idx + 1, 2 * idx + 1, cot_weight));
      //            }
      //            else {
      //                // pin点贡献到右边
      //                Vector2 pinPos = (vertices[j] == pin1) ? Vector2{ 0.0, 0.0 } : Vector2{ distance, 0.0 };
      //                b[2 * vertexIndex[vertices[k]]] += cot_weight * (result[vertices[k]].x - e_rotated.x);
      //                b[2 * vertexIndex[vertices[k]] + 1] += cot_weight * (result[vertices[k]].y - e_rotated.y);
      //            }
      //        }
      //    }
      //    A.setFromTriplets(A_entries.begin(), A_entries.end());
      //    Eigen::SimplicialLDLT<SparseMatrix<double>>solver(A);
      //    if (solver.info() == Eigen::Success) {
      //        Vector<double> solution = solver.solve(b);
      //        // 更新结果
      //        for (const auto& pair : vertexIndex) {
      //            Vertex v = pair.first;
      //            size_t idx = pair.second;
      //            result[v].x = solution[2 * idx];
      //            result[v].y = solution[2 * idx + 1];
      //        }
      //    }
      //}
    //}

    // 设置pin点的坐标
   /* result[pin1].x = 0.0;
    result[pin1].y = 0.0;
    result[pin2].x = distance;
    result[pin2].y = 0.0;*/
    return result;
}


VertexData<Vector2> SpectralConformalParameterization::lscm_flatten() const {

    VertexData<Vector2> result(*mesh);
    //0.建立局部坐标系
    //1.选择两个pin点，边界上距离最远的两个点
    
    // 计算每个面的局部坐标系
    std::map<Face, std::array<Vector2, 3>> localCoords = computeLocalCoordinates(mesh, geometry);
    
    // 找到边界上距离最远的两个点
    std::pair<Vertex, Vertex> pinPoints = findFarthestBoundaryPoints(mesh, geometry);
    Vertex pin1 = pinPoints.first;
    Vertex pin2 = pinPoints.second;
    

    // 创建顶点索引映射（排除pin点）
    std::map<Vertex, size_t> vertexIndex;
    size_t index = 0;
    for (Vertex v : mesh->vertices()) {
        if (v != pin1 && v != pin2) {
            vertexIndex[v] = index++;
        }
    }

    size_t nFaces = mesh->nFaces();
    size_t nVertices = mesh->nVertices();
    size_t nFreeVertices = nVertices - 2;
    SparseMatrix<double> A(2*nFaces,2* nFreeVertices);
    Vector<double> bRow(2*nFaces);
    std::vector<Eigen::Triplet<double>> A_entries;
    size_t faceIndex = 0;
    for (Face f: mesh->faces())
    {
        // 获取面的三个顶点
        std::array<Vertex, 3> vertices;
        size_t i = 0;
       for (Vertex v : f.adjacentVertices()) {
            vertices[i++] = v;
            if (i >= 3) break;
        }
        // 获取局部坐标
        std::array<Vector2, 3> local = localCoords[f];

        // 计算复数边向量
        std::complex<double> z01(local[1].x - local[0].x, local[1].y - local[0].y);
        std::complex<double> z02(local[2].x - local[0].x, local[2].y - local[0].y);
        
        // 计算雅可比矩阵 J = [a -b; b a]，其中 a 和 b 是复数 z01/z02 的实部和虚部
        std::complex<double> jacobian = z01 / z02;
        double a = jacobian.real();
        double b = jacobian.imag();
        
        int row_real = 2 * faceIndex;
        int row_imag = 2 * faceIndex + 1;
        
        // 为每个顶点建立方程
        for (int j = 0; j < 3; j++) {
            Vector2 localPos = local[j];
            if (vertices[j] != pin1 && vertices[j] != pin2)
            {
                size_t idx2 = vertexIndex[vertices[j]];
                A_entries.push_back(Eigen::Triplet<double>(row_real, 2*idx2, a));
                A_entries.push_back(Eigen::Triplet<double>(row_real, 2*idx2+1, -b));
                A_entries.push_back(Eigen::Triplet<double>(row_imag, 2*idx2, b));
                A_entries.push_back(Eigen::Triplet<double>(row_imag, 2*idx2+1, a));
            }
            else
            {
                 // 如果是pin点，将其贡献移到右边
                Vector2 pinPos = Vector2{0.0, 0.0};
                if (vertices[2] == pin1) {
                    pinPos = Vector2{0.0, 0.0};
                } else if (vertices[2] == pin2) {
                    Vector3 pos1 = geometry->inputVertexPositions[pin1];
                    Vector3 pos2 = geometry->inputVertexPositions[pin2];
                    double distance = (pos1 - pos2).norm();
                    pinPos = Vector2{distance, 0.0};
                }
                bRow[row_real] -= a * pinPos.x - b * pinPos.y;
                bRow[row_imag] -= b * pinPos.x + a * pinPos.y;
            }
        }
      faceIndex++;
    }
    A.setFromTriplets(A_entries.begin(), A_entries.end());
    //3. 求解最小二乘问题: min ||Ax - b||^2
    // 正规方程: A^T A x = A^T b  // 是否可以憨直地求逆解？ x = (ATA)_-1Ab
    SparseMatrix<double> AtA = A.transpose() * A;
    Eigen::SimplicialLDLT<SparseMatrix<double>>solver(AtA);
    Vector<double> Atb = A.transpose() * bRow;
    Vector<double> solution = solver.solve(Atb);
    for (const auto& pair : vertexIndex) {
        Vertex v = pair.first;
        size_t idx = pair.second;
        result[v].x = solution[2 * idx];
        result[v].y = solution[2 * idx + 1];
    }
    // 设置pin点的坐标：将两个最远点分别放在(0,0)和(d,0)，其中d是它们在3D空间中的距离
    Vector3 pos1 = geometry->inputVertexPositions[pin1];
    Vector3 pos2 = geometry->inputVertexPositions[pin2];
    double distance = (pos1 - pos2).norm();
    result[pin1].x = 0.0;
    result[pin1].y = 0.0;
    result[pin2].x = distance;
    result[pin2].y = 0.0;

    return result;
}