#include "CirclePatternsWasm.h"

CirclePatternsWasm::CirclePatternsWasm(Mesh& mesh0, int optScheme0):
Parameterization(mesh0),
angles(mesh.halfEdges.size()),
thetas(mesh.edges.size()),
OptScheme(optScheme0)
{
    radii.resize(mesh.faces.size() - 1);
    solver.n = (int)mesh.faces.size() - 1;
}

void CirclePatternsWasm::setConeSingulars(const std::vector<int>& coneIdx,
                                           const std::vector<double>& coneAngle)
{
    coneSingulars.clear();
    for (size_t i = 0; i < coneIdx.size() && i < coneAngle.size(); i++) {
        coneSingulars[coneIdx[i]] = coneAngle[i];
    }
}

// ============================================================
// Phase 1: Direct angle computation (Mosek-free)
// ============================================================

void CirclePatternsWasm::computeAnglesDirect()
{
    // Copy original half-edge angles with positivity clamping
    for (HalfEdgeCIter he = mesh.halfEdges.begin(); he != mesh.halfEdges.end(); he++) {
        if (!he->onBoundary) {
            angles[he->index] = std::max(he->angle(), EPSILON);
        } else {
            angles[he->index] = 0.0;
        }
    }
    
    // Compute edge weights τ_e from clamped angles
    // For interior edge e_ij: τ_e = π - a_ij^k - a_ij^l
    // For boundary edge e_ij: τ_e = π - a_ij^k
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter he = e->he;
        double a1 = angles[he->index];
        
        if (e->isBoundary()) {
            thetas[e->index] = M_PI - a1;
        } else {
            double a2 = angles[he->flip->index];
            double sum = a1 + a2;
            
            // Local Delaunay: a_ij^k + a_ij^l < π
            if (sum >= M_PI) {
                double scale = (M_PI - EPSILON) / sum;
                a1 *= scale; a2 *= scale;
                angles[he->index] = a1;
                angles[he->flip->index] = a2;
                sum = M_PI - EPSILON;
            }
            thetas[e->index] = std::max(M_PI - sum, EPSILON);
        }
    }
}

// ============================================================
// Phase 2: Radius optimization (Bobenko-Springborn variational principle)
// ============================================================

double ImLi2SumWasm(double dp, double theta)
{
    double tStar = M_PI - theta;
    double x = 2.0 * atan(tanh(0.5 * dp) * tan(0.5 * tStar));
    return x * dp + Cl2(x + tStar) + Cl2(-x + tStar) - Cl2(2.0 * tStar);
}

double feWasm(double dp, double theta)
{
    return atan2(sin(theta), exp(dp) - cos(theta));
}

void CirclePatternsWasm::computeEnergy(double& energy, const Eigen::VectorXd& rho)
{
    energy = 0.0;
    
    // Sum over edges: Bobenko-Springborn variational energy
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        int fk = e->he->face->index;
        
        if (e->isBoundary()) {
            energy -= 2.0 * (M_PI - thetas[e->index]) * rho[fk];
        } else {
            int fl = e->he->flip->face->index;
            energy += ImLi2SumWasm(rho[fk] - rho[fl], thetas[e->index])
                    - (M_PI - thetas[e->index]) * (rho[fk] + rho[fl]);
        }
    }
    
    // Sum over faces: 2π · ρ_f for each face
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) energy += 2.0 * M_PI * rho[f->index];
    }
}

void CirclePatternsWasm::computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& rho)
{
    // Gradient: ∇E_f = 2π - Σ_{e∈∂f} 2·φ_e(rho)
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fk = f->index;
            gradient[fk] = 2.0 * M_PI;
            
            HalfEdgeCIter he = f->he;
            do {
                EdgeCIter e = he->edge;
                if (e->isBoundary()) {
                    gradient[fk] -= 2.0 * (M_PI - thetas[e->index]);
                } else {
                    HalfEdgeCIter h = e->he;
                    int fl = (fk == (int)h->face->index) ? h->flip->face->index : h->face->index;
                    gradient[fk] -= 2.0 * feWasm(rho[fk] - rho[fl], thetas[e->index]);
                }
                he = he->next;
            } while (he != f->he);
        }
    }
}

void CirclePatternsWasm::computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& rho)
{
    std::vector<Eigen::Triplet<double>> HTriplets;
    
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!e->isBoundary()) {
            int fk = e->he->face->index;
            int fl = e->he->flip->face->index;
            
            double h = sin(thetas[e->index]) / (cosh(rho(fk) - rho(fl)) - cos(thetas[e->index]));
            HTriplets.push_back(Eigen::Triplet<double>(fk, fk, h + 1e-8));
            HTriplets.push_back(Eigen::Triplet<double>(fl, fl, h + 1e-8));
            HTriplets.push_back(Eigen::Triplet<double>(fk, fl, -h));
            HTriplets.push_back(Eigen::Triplet<double>(fl, fk, -h));
        }
    }
    
    hessian.resize(solver.n, solver.n);
    hessian.setFromTriplets(HTriplets.begin(), HTriplets.end());
}

void CirclePatternsWasm::setRadii()
{
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) radii[f->index] = exp(solver.x[f->index]);
    }
}

bool CirclePatternsWasm::computeRadii()
{
    MeshHandle handle;
    handle.computeEnergy   = std::bind(&CirclePatternsWasm::computeEnergy,   this, _1, _2);
    handle.computeGradient = std::bind(&CirclePatternsWasm::computeGradient, this, _1, _2);
    handle.computeHessian  = std::bind(&CirclePatternsWasm::computeHessian,  this, _1, _2);
    
    solver.handle = &handle;
    if (OptScheme == GRAD_DESCENT) solver.gradientDescent();
    else if (OptScheme == NEWTON) solver.newton();
    else if (OptScheme == LBFGS) solver.lbfgs();
    else solver.newton();
    
    setRadii();
    return true;
}

// ============================================================
// Phase 3: UV layout
// ============================================================

void CirclePatternsWasm::computeAnglesAndEdgeLengths(Eigen::VectorXd& lengths)
{
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter h1 = e->he;
        if (e->isBoundary() && h1->onBoundary) {
            h1 = h1->flip;
        }
        
        if (e->isBoundary()) {
            angles[h1->index] = M_PI - thetas[e->index];
        } else {
            HalfEdgeCIter h2 = h1->flip;
            double dp = log(radii[h1->face->index]) - log(radii[h2->face->index]);
            angles[h1->index] = feWasm(dp, thetas[e->index]);
            angles[h2->index] = feWasm(-dp, thetas[e->index]);
        }
        lengths[e->index] = 2.0 * radii[h1->face->index] * sin(angles[h1->index]);
    }
}

void CirclePatternsWasm::performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir,
                                            Eigen::VectorXd& lengths,
                                            std::unordered_map<int, bool>& visited,
                                            std::stack<EdgeCIter>& stack)
{
    if (!he->onBoundary) {
        int fIdx = he->face->index;
        if (visited.find(fIdx) == visited.end()) {
            HalfEdgeCIter next = he->next;
            HalfEdgeCIter prev = he->next->next;
            
            double a = angles[next->index];
            Eigen::Vector2d newDir = {
                cos(a)*dir[0] - sin(a)*dir[1],
                sin(a)*dir[0] + cos(a)*dir[1]
            };
            prev->vertex->uv = he->vertex->uv + newDir * lengths[prev->edge->index];
            
            visited[fIdx] = true;
            stack.push(next->edge);
            stack.push(prev->edge);
        }
    }
}

void CirclePatternsWasm::setUVs()
{
    Eigen::VectorXd lengths(mesh.edges.size());
    computeAnglesAndEdgeLengths(lengths);
    
    std::stack<EdgeCIter> stack;
    EdgeCIter e = mesh.edges.begin();
    stack.push(e);
    e->he->vertex->uv = Eigen::Vector2d::Zero();
    e->he->next->vertex->uv = Eigen::Vector2d(lengths[e->index], 0);
    
    std::unordered_map<int, bool> visited;
    while (!stack.empty()) {
        EdgeCIter e = stack.top(); stack.pop();
        HalfEdgeCIter h1 = e->he, h2 = h1->flip;
        Eigen::Vector2d dir = h2->vertex->uv - h1->vertex->uv;
        dir.normalize();
        performFaceLayout(h1, dir, lengths, visited, stack);
        performFaceLayout(h2, -dir, lengths, visited, stack);
    }
    
    normalize();
}

// ============================================================
// Main entry point
// ============================================================

void CirclePatternsWasm::parameterize()
{
    // Phase 1: Compute edge weights τ_e (Mosek-free)
    computeAnglesDirect();
    
    // Phase 2: Optimize face radii (convex)
    if (!computeRadii()) {
        std::cout << "CirclePatternsWasm: radius optimization failed" << std::endl;
        return;
    }
    
    // Phase 3: UV layout from radii
    setUVs();
}
