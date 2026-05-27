#include "gcp.h"
//#include <geometrycentral/surface/surface_mesh.h>
//#include <geometrycentral/surface/surface_mesh_io.h>
//#include <geometrycentral/surface/vertex_position_geometry.h>
//#include <geometrycentral/surface/edge_length_geometry.h>
//#include <geometrycentral/surface/vertex_position_geometry.h>
//#include <geometrycentral/surface/face_data.h>
//#include <geometrycentral/surface/vertex_data.h>
//#include <geometrycentral/surface/edge_data.h>
//#include <geometrycentral/surface/face_data.h>
//#include <geometrycentral/surface/vertex_position_geometry.h>
//#include <geometrycentral/surface/mesh_io.h>
//#include <geometrycentral/surface/surface_mesh.h>
//#include <geometrycentral/surface/vertex_position_geometry.h>
//#include <geometrycentral/surface/edge_length_geometry.h>
//#include <geometrycentral/surface/vertex_data.h>
//#include <geometrycentral/surface/edge_data.h>
//#include <geometrycentral/surface/face_data.h>
//#include <geometrycentral/surface/parameterize.h>
//
//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <Eigen/Dense>
//
//using namespace geometrycentral;
//using namespace geometrycentral::surface;
//
//// 1. 计算 cotangent 权重
//void computeCotangentWeights(const VertexPositionGeometry& geometry, 
//                            FaceData<double>& cotWeights) {
//    for (Face f : geometry.mesh().faces()) {
//        double w = 0.0;
//        for (Vertex v : f.adjacentVertices()) {
//            // 获取对边的长度
//            Edge e = f.edge(v);
//            double l = geometry.edgeLength(e);
//            // 获取对边对面的顶点
//            Vertex v2 = e.otherVertex(v);
//            // 获取对边对面的顶点
//            Vertex v3 = e.otherVertex(v2);
//            // 计算对边对面的角
//            double cosA = (geometry.vertexPosition(v2) - geometry.vertexPosition(v)).dot(
//                            geometry.vertexPosition(v3) - geometry.vertexPosition(v)) / 
//                          (geometry.vertexPosition(v2).norm() * geometry.vertexPosition(v3).norm());
//            // cotangent = cosA / sinA
//            double sinA = sqrt(1 - cosA*cosA);
//            if (sinA > 1e-6) {
//                w += cosA / sinA;
//            }
//        }
//        cotWeights[f] = w;
//    }
//}
//
//// 2. 构建拉普拉斯矩阵
//void buildLaplacian(const VertexPositionGeometry& geometry, 
//                   const FaceData<double>& cotWeights,
//                   SparseMatrix<double>& L) {
//    int n = geometry.mesh().nVertices();
//    L = SparseMatrix<double>(n, n);
//    vector<Triplet<double>> triplets;
//
//    for (Vertex v : geometry.mesh().vertices()) {
//        double sum = 0.0;
//        for (Vertex w : v.adjacentVertices()) {
//            // 获取边
//            Edge e = v.edge(w);
//            double weight = cotWeights[e];
//            triplets.emplace_back(v.index(), w.index(), -weight);
//            sum += weight;
//        }
//        triplets.emplace_back(v.index(), v.index(), sum);
//    }
//    L.setFromTriplets(triplets.begin(), triplets.end());
//}
//
//// 3. 求解调和 1-形式
//VectorXd solveHarmonic1Form(const VertexPositionGeometry& geometry, 
//                           const FaceData<double>& cotWeights,
//                           const vector<vector<int>>& basis) {
//    int n = geometry.mesh().nVertices();
//    int m = basis.size();
//    
//    // 构建约束矩阵 A
//    MatrixXd A(m, n);
//    for (int i = 0; i < m; ++i) {
//        for (int j = 0; j < n; ++j) {
//            A(i, j) = (basis[i][j] ? 1.0 : 0.0);
//        }
//    }
//
//    // 构建质量矩阵 M (对角，cotangent 权重)
//    VectorXd M_diag(n);
//    for (int i = 0; i < n; ++i) {
//        M_diag[i] = 1.0; // 实际中用 cotangent 权重
//    }
//
//    // 求解 min omega^T M omega s.t. A omega = b
//    // 使用拉格朗日乘子法
//    SparseMatrix<double> K(n + m, n + m);
//    vector<Triplet<double>> triplets;
//    
//    // 2M 块
//    for (int i = 0; i < n; ++i) {
//        triplets.emplace_back(i, i, 2 * M_diag[i]);
//    }
//    // A^T 块
//    for (int i = 0; i < m; ++i) {
//        for (int j = 0; j < n; ++j) {
//            if (A(i, j) != 0) {
//                triplets.emplace_back(j, n + i, A(i, j));
//                triplets.emplace_back(n + i, j, A(i, j));
//            }
//        }
//    }
//    
//    K.setFromTriplets(triplets.begin(), triplets.end());
//    
//    // 构建右边
//    VectorXd rhs(n + m);
//    rhs.setZero();
//    // 仅第一个约束设为1
//    rhs(n) = 1.0;
//
//    // 求解
//    SimplicialLDLT<SparseMatrix<double>> solver(K);
//    VectorXd sol = solver.solve(rhs);
//    return sol.head(n);
//}
//
//// 4. 积分 1-形式得到坐标
//VectorXd integrate1Form(const VertexPositionGeometry& geometry, 
//                      const VectorXd& omega) {
//    int n = geometry.mesh().nVertices();
//    VectorXd u = VectorXd::Zero(n);
//    u[0] = 0.0; // 固定原点
//    
//    // 使用 BFS 从顶点0开始积分
//    vector<bool> visited(n, false);
//    queue<int> q;
//    q.push(0);
//    visited[0] = true;
//    
//    while (!q.empty()) {
//        int v = q.front(); q.pop();
//        for (Vertex w : geometry.mesh().vertex(v).adjacentVertices()) {
//            if (!visited[w.index()]) {
//                visited[w.index()] = true;
//                // 计算边上的值
//                double val = 0.0;
//                for (Edge e : geometry.mesh().vertex(v).adjacentEdges()) {
//                    if (e.otherVertex(v) == w) {
//                        // 使用边权重
//                        val = omega[e.index()];
//                        break;
//                    }
//                }
//                u[w.index()] = u[v] + val;
//                q.push(w.index());
//            }
//        }
//    }
//    return u;
//}
//
//int main(int argc, char* argv[]) {
//    if (argc < 2) {
//        std::cerr << "Usage: " << argv[0] << " <input.obj>" << std::endl;
//        return 1;
//    }
//
//    // 1. 读取网格
//    std::string input = argv[1];
//    std::unique_ptr<SurfaceMesh> mesh = loadMeshFromFile(input);
//    if (!mesh) {
//        std::cerr << "Failed to load mesh from " << input << std::endl;
//        return 1;
//    }
//
//    // 2. 创建几何对象
//    VertexPositionGeometry geometry(*mesh);
//    geometry.requireEdgeLengths();
//    geometry.requireFaceAreas();
//    
//    // 3. 计算 cotangent 权重
//    FaceData<double> cotWeights(*mesh);
//    computeCotangentWeights(geometry, cotWeights);
//    
//    // 4. 计算上同调基（简化版，仅用于教学）
//    // 实际中需要 Tree-Cotree 算法
//    int genus = 0; // 估算亏格
//    vector<vector<int>> basis;
//    if (genus > 0) {
//        // 为每个基创建一个 1-形式
//        for (int i = 0; i < 2 * genus; ++i) {
//            basis.push_back(vector<int>(mesh->nEdges(), 0));
//            basis[i][i] = 1; // 简化：使用第 i 条边作为基
//        }
//    } else {
//        // 球面（亏格0）：使用球面坐标
//        // 这里跳过，直接输出球面坐标
//        std::cout << "Genus 0: Using spherical parameterization" << std::endl;
//        // 实现球面坐标（简化）
//        for (Vertex v : geometry.mesh().vertices()) {
//            Vector3d pos = geometry.vertexPosition(v);
//            double r = pos.norm();
//            double theta = atan2(pos.y(), pos.x());
//            double phi = acos(pos.z() / r);
//            // 保存到 UV
//            // ... (这里省略，球面坐标不是 GCP)
//        }
//        return 0;
//    }
//
//    // 5. 求解调和 1-形式并积分
//    MatrixXd UV(mesh->nVertices(), 2);
//    
//    // 为每个基求解
//    for (int i = 0; i < basis.size(); ++i) {
//        VectorXd omega = solveHarmonic1Form(geometry, cotWeights, basis);
//        VectorXd u = integrate1Form(geometry, omega);
//        UV.col(i) = u;
//    }
//
//    // 6. 保存 UV 坐标
//    std::ofstream uvFile("output_uv.txt");
//    for (int i = 0; i < mesh->nVertices(); ++i) {
//        uvFile << UV(i, 0) << " " << UV(i, 1) << std::endl;
//    }
//    uvFile.close();
//
//    std::cout << "GCP completed. UV saved to output_uv.txt" << std::endl;
//    std::cout << "Number of vertices: " << mesh->nVertices() << std::endl;
//    std::cout << "Number of faces: " << mesh->nFaces() << std::endl;
//    std::cout << "Genus: " << genus << std::endl;
//
//    return 0;
//}