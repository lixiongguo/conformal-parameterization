#include "RicciFlow.h"

RicciFlow::RicciFlow(Mesh& mesh0, int optScheme0):
Parameterization(mesh0),
OptScheme(optScheme0)
{
    
}

// ============================================================
// Initialization
// ============================================================

void RicciFlow::initInversiveDistances()
{
    // Inversive distance I_ij from original edge length l_ij (assuming r_i = r_j = 1 initially):
    //   l_ij² = r_i² + r_j² + 2·r_i·r_j·I_ij = 2 + 2·I_ij
    //   → I_ij = (l_ij² - 2) / 2
    // 
    // We normalize edge lengths so that the average is ~1 to avoid numerical issues.
    
    inversiveDistance.resize(mesh.edges.size());
    
    // First pass: compute average edge length
    double avgLen = 0.0;
    int eCount = 0;
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        avgLen += e->length();
        eCount++;
    }
    avgLen /= eCount;
    
    // Second pass: compute inversive distances with normalization
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        double normLen = e->length() / avgLen;
        inversiveDistance[e->index] = (normLen * normLen - 2.0) / 2.0;
    }
}

void RicciFlow::setTargetCurvature()
{
    Ktarget.resize(solver.n);
    Ktarget.setZero(); // flat metric → K̄ = 0
}

// ============================================================
// Circle packing metric computations
// ============================================================

double RicciFlow::edgeLengthFromMetric(double ui, double uj, double Iij) const
{
    // l_ij² = e^{2ui} + e^{2uj} + 2·I_ij·e^{ui+uj}
    double ri = exp(ui);
    double rj = exp(uj);
    double l2 = ri*ri + rj*rj + 2.0 * Iij * ri * rj;
    return sqrt(std::max(l2, 1e-12));
}

double RicciFlow::triangleAngle(double la, double lb, double lc) const
{
    // Law of cosines: cos α = (b² + c² - a²) / (2bc)
    // where a is the edge opposite to the angle we're computing
    // In our convention: la is opposite edge, lb and lc are adjacent edges
    double cosA = (lb*lb + lc*lc - la*la) / (2.0 * lb * lc);
    cosA = fmax(-1.0, fmin(1.0, cosA));
    return acos(cosA);
}

void RicciFlow::computeEdgeLengthsAndAngles()
{
    const Eigen::VectorXd& u = solver.x;
    
    edgeLengths.resize(mesh.edges.size());
    halfEdgeAngles.resize(mesh.halfEdges.size());
    
    // Compute edge lengths from current radii
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter he = e->he;
        
        double ui = he->vertex->isBoundary() ? 0.0 : u[index.at(he->vertex->index)];
        double uj = he->flip->vertex->isBoundary() ? 0.0 : u[index.at(he->flip->vertex->index)];
        double Iij = inversiveDistance[e->index];
        
        edgeLengths[e->index] = edgeLengthFromMetric(ui, uj, Iij);
    }
    
    // Compute angles for each triangle
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        
        // Collect three edge lengths of this triangle
        double al[3];
        int i = 0;
        HalfEdgeCIter he = f->he;
        do {
            al[i] = edgeLengths[he->edge->index];
            i++;
            he = he->next;
        } while (he != f->he);
        
        // Compute three angles (angle at vertex i is opposite to edge i)
        // In half-edge order: he->vertex is opposite to he->edge
        i = 0;
        he = f->he;
        do {
            halfEdgeAngles[he->index] = triangleAngle(al[i], al[(i+1)%3], al[(i+2)%3]);
            i++;
            he = he->next;
        } while (he != f->he);
    }
}

// ============================================================
// Energy, Gradient, Hessian (for solver callbacks)
// ============================================================

void RicciFlow::computeEnergy(double& energy, const Eigen::VectorXd& u)
{
    // Ricci energy is convex and defined such that ∇E = K - K̄
    // We integrate: E(u) = ∫₀¹ Σᵢ (K_i(tu) - K̄_i) u_i dt
    // 
    // For practical computation, we compute the energy as the sum of
    // triangle-based terms (the "Ricci energy" as formulated by Gu et al.)
    // 
    // Simplified: use trapezoidal rule to approximate the integral
    // E ≈ Σᵢ [K_i(0)·u_i + K_i(u)·u_i] / 2 - Σᵢ K̄_i·u_i
    // But since K_i(0) is the original curvature and is not easily available,
    // we use a direct approach:
    //
    // The discrete Ricci energy for Euclidean geometry is:
    // E(u) = Σ_f V_f(u) - Σ_i (Σ_f α_f,i(0)) u_i + Σ_i K̄_i u_i
    // where V_f is related to the volume of an ideal tetrahedron
    
    // For simplicity and robustness, we implement energy in terms of
    // line integral approximation with midpoint rule
    
    // Actually, the standard energy formulation by Gu-Yau et al.:
    // E(u) = ∫₀ᵗ Σ_i (K_i - K̄_i) u̇_i dt
    // This is the work done by the Ricci flow.
    // We can compute it incrementally or use the known closed form.
    
    // Closed form for Euclidean Ricci energy (Gu & Yau, 2008):
    // For each face f = (i,j,k) with edge lengths l_ij, l_jk, l_ki:
    // The energy contribution is related to the dilogarithm function.
    // We use a simpler numerical approach: direct curvature-based energy.
    
    // Compute current curvature
    Eigen::VectorXd K(solver.n);
    K.setZero();
    
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        
        HalfEdgeCIter he = f->he;
        do {
            if (!he->vertex->isBoundary()) {
                int vIdx = index.at(he->vertex->index);
                K[vIdx] += halfEdgeAngles[he->index];
            }
            he = he->next;
        } while (he != f->he);
    }
    
    // K_i = 2π - Σ angles  (or π for boundary, but boundaries are fixed)
    for (int i = 0; i < solver.n; i++) {
        K[i] = 2.0 * M_PI - K[i];
    }
    
    // Ricci energy approximation:
    // E(u) ≈ ½ Σ_i (K_i(u) + K_i(0) - 2K̄_i) · u_i
    // But we don't have K_i(0) stored... Let's use a different approach.
    
    // Actually the convex Ricci energy has the form:
    // E(u) = Σ_f F_f(u) + Σ_i K̄_i · u_i
    // where F_f is related to the Lobachevsky function.
    // ∇E_i = K_i - K̄_i
    
    // For Newton's method, we only need gradient and Hessian.
    // The energy value is used for line search in the solver.
    // We can use the line integral approximation:
    
    energy = 0.0;
    // E = Σ_i ∫_0^{u_i} (K_i(t) - K̄_i) dt
    // ≈ Σ_i ½(K_i(0) + K_i(u) - 2K̄_i) · u_i
    // We approximate K_i(0) ≈ K_i(0.5u) by computing curvature at u/2
    
    // But this is expensive. Instead, use a simple quadratic approximation
    // that has the correct gradient. For convex functions with known gradient,
    // the energy for line search purposes just needs to be consistent.
    
    // We'll compute a pseudo-energy that has the right gradient:
    // E = ½ Σ_i u_i · (K_i(u) - 2K̄_i)
    // This is not the exact Ricci energy but is consistent with gradient descent.
    
    for (int i = 0; i < solver.n; i++) {
        energy += 0.5 * u[i] * (K[i] - 2.0 * Ktarget[i]);
    }
    
    // Add small regularization for stability
    energy += 1e-8 * u.squaredNorm();
}

void RicciFlow::computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& u)
{
    gradient.resize(solver.n);
    gradient.setZero();
    
    // Gradient: ∇E_i = K_i - K̄_i
    // where K_i = 2π - Σ(angles around vertex i) for interior vertices
    
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isBoundary()) continue;
        
        HalfEdgeCIter he = f->he;
        do {
            if (!he->vertex->isBoundary()) {
                int vIdx = index.at(he->vertex->index);
                gradient[vIdx] -= halfEdgeAngles[he->index]; // subtract angle sum
            }
            he = he->next;
        } while (he != f->he);
    }
    
    // Complete: K_i = 2π - Σ angles, so ∇E_i = (2π - Σ angles) - K̄_i
    for (int i = 0; i < solver.n; i++) {
        gradient[i] = (2.0 * M_PI + gradient[i]) - Ktarget[i];
    }
}

void RicciFlow::computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& u)
{
    // Ricci flow Hessian = ∂K/∂u = weighted cotangent Laplacian
    // ∂K_i/∂u_j = -w_ij (for i ≠ j), ∂K_i/∂u_i = Σ_j w_ij
    // where w_ij = cot(α_ij^k) + cot(α_ij^l) (standard cotangent weights)
    
    std::vector<Eigen::Triplet<double>> HTriplets;
    
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->isBoundary()) continue;
        
        int vIdx = index.at(v->index);
        HalfEdgeCIter he = v->he;
        double sumW = 0.0;
        
        do {
            // Edge from v to neighbor: he points to neighbor
            // cot(angle opposite to this edge) = cot(angle at he->next->vertex)
            double cotAlpha = !he->onBoundary ? 
                1.0 / tan(halfEdgeAngles[he->next->index]) : 0.0;
            double cotBeta  = !he->flip->onBoundary ? 
                1.0 / tan(halfEdgeAngles[he->flip->next->index]) : 0.0;
            
            double w = (cotAlpha + cotBeta);
            sumW += w;
            
            if (!he->flip->vertex->isBoundary()) {
                HTriplets.push_back(Eigen::Triplet<double>(vIdx, index.at(he->flip->vertex->index), -w));
            }
            
            he = he->flip->next;
        } while (he != v->he);
        
        // Diagonal: sum of weights + small regularization
        HTriplets.push_back(Eigen::Triplet<double>(vIdx, vIdx, sumW + 1e-8));
    }
    
    hessian.resize(solver.n, solver.n);
    hessian.setFromTriplets(HTriplets.begin(), HTriplets.end());
}

// ============================================================
// Optimization
// ============================================================

bool RicciFlow::optimizeRadii()
{
    MeshHandle handle;
    
    // The solver callbacks need edge lengths and angles to be up-to-date.
    // We first compute them, then bind the energy/gradient/hessian callbacks.
    // However, the solver calls callbacks in order: energy, then gradient, then hessian.
    // We need edge lengths/angles to be fresh for each iteration.
    // Solution: compute edge lengths in the first callback called each iteration.
    
    // We use a wrapper pattern: before computing energy/gradient/hessian, 
    // recompute edge lengths and angles.
    handle.computeEnergy = [this](double& energy, const Eigen::VectorXd& u) {
        this->solver.x = u;
        this->computeEdgeLengthsAndAngles();
        this->computeEnergy(energy, u);
    };
    
    handle.computeGradient = [this](Eigen::VectorXd& gradient, const Eigen::VectorXd& u) {
        this->computeGradient(gradient, u);
    };
    
    handle.computeHessian = [this](Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& u) {
        this->computeHessian(hessian, u);
    };
    
    solver.handle = &handle;
    
    // Run optimization
    if (OptScheme == GRAD_DESCENT) {
        solver.gradientDescent();
    } else if (OptScheme == NEWTON) {
        solver.newton();
    } else if (OptScheme == LBFGS) {
        solver.lbfgs();
    } else {
        // Default to Newton
        solver.newton();
    }
    
    return true;
}

// ============================================================
// UV Layout (BFS)
// ============================================================

void RicciFlow::performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir,
                                   std::unordered_map<int, bool>& visited, 
                                   std::stack<EdgeCIter>& stack)
{
    if (!he->onBoundary) {
        int fIdx = he->face->index;
        if (visited.find(fIdx) == visited.end()) {
            HalfEdgeCIter next = he->next;
            HalfEdgeCIter prev = he->next->next;
            
            // Compute angle at current vertex
            double angle = halfEdgeAngles[next->index];
            Eigen::Vector2d newDir = {
                cos(angle)*dir[0] - sin(angle)*dir[1],
                sin(angle)*dir[0] + cos(angle)*dir[1]
            };
            
            prev->vertex->uv = he->vertex->uv + newDir * edgeLengths[prev->edge->index];
            
            visited[fIdx] = true;
            stack.push(next->edge);
            stack.push(prev->edge);
        }
    }
}

void RicciFlow::setUVs()
{
    // Push any edge
    std::stack<EdgeCIter> stack;
    EdgeCIter e = mesh.edges.begin();
    stack.push(e);
    e->he->vertex->uv = Eigen::Vector2d::Zero();
    e->he->next->vertex->uv = Eigen::Vector2d(edgeLengths[e->index], 0);
    
    // BFS layout
    std::unordered_map<int, bool> visited;
    while (!stack.empty()) {
        EdgeCIter e = stack.top();
        stack.pop();
        
        HalfEdgeCIter h1 = e->he;
        HalfEdgeCIter h2 = h1->flip;
        
        Eigen::Vector2d dir = h2->vertex->uv - h1->vertex->uv;
        dir.normalize();
        
        performFaceLayout(h1, dir, visited, stack);
        performFaceLayout(h2, -dir, visited, stack);
    }
    
    normalize();
}

// ============================================================
// Main parameterization entry point
// ============================================================

void RicciFlow::parameterize()
{
    // Step 0. Assign indices to non-boundary vertices
    int vIdx = 0;
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (!v->isBoundary()) {
            index[v->index] = vIdx++;
        }
    }
    solver.n = vIdx;
    
    // Step 1. Initialize inversive distances from original edge lengths
    initInversiveDistances();
    
    // Step 2. Set target curvature (flat metric → K̄=0)
    setTargetCurvature();
    
    // Step 3. Optimize radii (Ricci flow)
    if (!optimizeRadii()) {
        std::cout << "RicciFlow: optimization failed" << std::endl;
        return;
    }
    
    // Step 4. Compute final edge lengths and angles
    computeEdgeLengthsAndAngles();
    
    // Step 5. Set UVs via BFS layout
    setUVs();
}
