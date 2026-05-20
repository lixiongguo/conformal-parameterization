
VertexData<Vector2> SpectralConformalParameterization::getUVFromConformalFactorWithBFS(const VertexData<double>& conformalFactor) const {

   VertexData<Vector2> uvCoordinates(*mesh);
    VertexData<bool> visited(*mesh, false);
    std::queue<Vertex> bfsQueue;
    
    // 1. 计算边的缩放系数和参数化后的边长
    EdgeData<double> edgeScaleFactors(*mesh);
    EdgeData<double> parameterizedLengths(*mesh);
    
    geometry->requireEdgeLengths(); // 确保边长度已计算
    
    for (Edge e : mesh->edges()) {
        Vertex v1 = e.firstVertex();
        Vertex v2 = e.secondVertex();
        
        // 计算边的缩放系数：exp((φ(v1) + φ(v2))/2)
        double scaleFactor = exp((conformalFactor[v1] + conformalFactor[v2]) / 2.0);
        edgeScaleFactors[e] = scaleFactor;
        
        // 计算参数化后的边长
        double originalLength = geometry->edgeLengths[e];
        parameterizedLengths[e] = originalLength * scaleFactor;
    }
    
    // 2. 设置起始点和边界条件
    // 找到第一个边界顶点作为起始点
    Vertex startVertex;
    bool foundStart = false;
    
    for (Vertex v : mesh->vertices()) {
        if (v.isBoundary()) {
            startVertex = v;
            foundStart = true;
            break;
        }
    }
    
    // 如果没有边界顶点，选择任意内部顶点作为起始点
    if (!foundStart) {
        startVertex = *mesh->vertices().begin();
    }
    
    // 设置起始点的UV坐标
    uvCoordinates[startVertex] = Vector2(0.0, 0.0);
    visited[startVertex] = true;
    bfsQueue.push(startVertex);
    
    // 3. 使用BFS构建UV坐标
    while (!bfsQueue.empty()) {
        Vertex current = bfsQueue.front();
        bfsQueue.pop();
        
        // 遍历当前顶点的所有相邻顶点
        for (Halfedge he : current.outgoingHalfedges()) {
            Vertex neighbor = he.tipVertex();
            
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                bfsQueue.push(neighbor);
                
                // 获取当前边
                Edge e = he.edge();
                
                // 计算相邻顶点的UV坐标
                // 基本思想是沿着边的方向，根据参数化后的边长设置UV坐标
                Vector2 currentUV = uvCoordinates[current];
                Vector2 neighborUV;
                
                // 确定边的方向
                if (e.firstVertex() == current && e.secondVertex() == neighbor) {
                    // 寻找一个已访问的相邻顶点来确定方向
                    bool foundDirection = false;
                    
                    for (Halfedge heAdj : current.outgoingHalfedges()) {
                        Vertex adjNeighbor = heAdj.tipVertex();
                        if (adjNeighbor != neighbor && visited[adjNeighbor]) {
                            // 计算边向量
                            Vector2 dirVec = uvCoordinates[adjNeighbor] - currentUV;
                            if (dirVec.norm() > 1e-6) {
                                // 计算垂直方向
                                Vector2 perpDir(-dirVec.y, dirVec.x);
                                perpDir = perpDir.normalize();
                                
                                // 设置邻居的UV坐标
                                neighborUV = currentUV + perpDir * parameterizedLengths[e];
                                foundDirection = true;
                                break;
                            }
                        }
                    }
                    
                    // 如果没有找到方向参考，使用默认方向（x轴正方向）
                    if (!foundDirection) {
                        neighborUV = currentUV + Vector2(parameterizedLengths[e], 0.0);
                    }
                } else {
                    // 处理反向边
                    neighborUV = currentUV - Vector2(parameterizedLengths[e], 0.0);
                }
                
                uvCoordinates[neighbor] = neighborUV;
            }
        }
    }
}
VertexData<Vector2> SpectralConformalParameterization::benchen_flatten() const {
    //yamabe 方程

    VertexData<Vector2> result(*mesh);
}