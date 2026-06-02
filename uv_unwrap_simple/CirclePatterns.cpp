#include "CirclePatterns.h"
#include <algorithm>
#include "QPsolver/QPSolver.h"

CirclePatterns::CirclePatterns(Mesh& mesh0, int optScheme0):
Parameterization(mesh0),
angles(mesh.halfEdges.size()),
thetas(mesh.edges.size()),
radii(std::max(1, (int)mesh.faces.size() - (int)mesh.boundaries.size())),
eIntIndices(mesh.edges.size()),
imaginaryHe(0),
OptScheme(optScheme0)
{
  // One radius per interior (non-imaginary-boundary) face; indices use Face::index.
  solver.n = std::max(1, (int)mesh.faces.size() - (int)mesh.boundaries.size());
}

void CirclePatterns::setupAngleOptProblem()
{
    // set triangle sum constraint
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fIdx = f->index;
            mosekSolver.bkc[fIdx] = MSK_BK_FX;
            mosekSolver.blc[fIdx] = M_PI;
            mosekSolver.buc[fIdx] = M_PI;
        }
    }
    
    // set vertex sum constraint
    int shift = (int)(mesh.faces.size() - mesh.boundaries.size());
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int vIdx = v->index + shift;
        if (v->isBoundary()) {
            mosekSolver.bkc[vIdx] = MSK_BK_FX;
            mosekSolver.blc[vIdx] = M_PI;
            mosekSolver.buc[vIdx] = M_PI;
        } else {
            mosekSolver.bkc[vIdx] = MSK_BK_FX;
            mosekSolver.blc[vIdx] = 2*M_PI;
            mosekSolver.buc[vIdx] = 2*M_PI;
        }
    }
    
    // set local delaunay constraint
    shift += (int)mesh.vertices.size();
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!e->isBoundary()) {
            int eIdx = eIntIndices[e->index] + shift;
            mosekSolver.bkc[eIdx] = MSK_BK_RA;
            mosekSolver.blc[eIdx] = EPSILON;
            mosekSolver.buc[eIdx] = M_PI - EPSILON;
        }
    }
    
    // set angle matrices
    int a = 0;
    int numanz = 0;
    for (HalfEdgeCIter he = mesh.halfEdges.begin(); he != mesh.halfEdges.end(); he++) {
        if (!he->onBoundary) {
            double angle = he->angle();
            
            // subtract angle from constraints and set linear constraint matrix
            mosekSolver.ptrb[a] = numanz;
            
            int fIdx = he->face->index;
            mosekSolver.blc[fIdx] -= angle;
            mosekSolver.buc[fIdx] -= angle;
            mosekSolver.asub[numanz] = fIdx;
            mosekSolver.aval[numanz] = 1.0;
            numanz++;
            
            int vIdx = he->next->next->vertex->index + (int)(mesh.faces.size() - mesh.boundaries.size());
            mosekSolver.blc[vIdx] -= angle;
            mosekSolver.buc[vIdx] -= angle;
            mosekSolver.asub[numanz] = vIdx;
            mosekSolver.aval[numanz] = 1.0;
            numanz++;
            
            if (!he->edge->isBoundary()) {
                int eIdx = eIntIndices[he->edge->index] +
                           (int)(mesh.faces.size() + mesh.vertices.size() - mesh.boundaries.size());
                mosekSolver.blc[eIdx] -= angle;
                mosekSolver.buc[eIdx] -= angle;
                mosekSolver.asub[numanz] = eIdx;
                mosekSolver.aval[numanz] = 1.0;
                numanz++;
            }
            
            mosekSolver.ptre[a] = numanz;
            
            // ensure positive and non reflex angles
            mosekSolver.bkx[a] = MSK_BK_RA;
            mosekSolver.blx[a] = EPSILON - angle;
            mosekSolver.bux[a] = (M_PI - EPSILON) - angle;
            
            // set objective
            mosekSolver.qsubi[a] = a;
            mosekSolver.qsubj[a] = a;
            mosekSolver.qval[a] = 1.0;
            mosekSolver.c[a] = 0.0;
            
            a++;
        }
    }
}

void CirclePatterns::setThetas()
{
    // set opposite angles
    int a = 0;
    for (HalfEdgeCIter he = mesh.halfEdges.begin(); he != mesh.halfEdges.end(); he++) {
        if (!he->onBoundary) angles[he->index] = he->angle() + mosekSolver.xx[a++];
        else angles[he->index] = 0.0;
    }
    
    // set thetas
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter he = e->he;
        thetas[e->index] = M_PI - angles[he->index] - angles[he->flip->index];
    }
}

bool CirclePatterns::computeAngles()
{
    int variables = 3 * (int)(mesh.faces.size() - mesh.boundaries.size()); 
    int constraints = (int)(mesh.vertices.size() + mesh.edges.size() + mesh.faces.size() -
                            imaginaryHe - mesh.boundaries.size());
    int numanz = 3 * variables - imaginaryHe;
    int numqnz = variables;
    
    // Try MOSEK first (native builds). WASM uses MosekStub and will fail here.
    if (mosekSolver.initialize(variables, constraints, numanz, 0, numqnz)) {
        setupAngleOptProblem();
        bool success = mosekSolver.solve(MosekSolver::QO);
        if (success) {
            setThetas();
            mosekSolver.reset();
            return true;
        }
        mosekSolver.reset();
    }

    // Fallback: solve the same convex QP with a lightweight ADMM-based QP solver.
    // This is particularly important for WASM builds where MOSEK is not available.
    QPsolver::MinNormBoxQPProblem qp;
    qp.n = variables;
    qp.m = constraints;

    const int faceShift = (int)(mesh.faces.size() - mesh.boundaries.size());
    const int edgeShift = faceShift + (int)mesh.vertices.size();

    qp.lowerX = Eigen::VectorXd::Zero(variables);
    qp.upperX = Eigen::VectorXd::Zero(variables);
    qp.lowerY = Eigen::VectorXd::Zero(constraints);
    qp.upperY = Eigen::VectorXd::Zero(constraints);

    // Constraint initialization (triangle sum / vertex sum / local Delaunay range)
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            const int row = f->index;
            if (row < 0 || row >= faceShift) return false;
            qp.lowerY[row] = M_PI;
            qp.upperY[row] = M_PI;
        }
    }
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        const int row = v->index + faceShift;
        if (row < 0 || row >= constraints) return false;
        if (v->isBoundary()) {
            qp.lowerY[row] = M_PI;
            qp.upperY[row] = M_PI;
        } else {
            qp.lowerY[row] = 2.0 * M_PI;
            qp.upperY[row] = 2.0 * M_PI;
        }
    }
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!e->isBoundary()) {
            const int row = eIntIndices[e->index] + edgeShift;
            if (row < 0 || row >= constraints) return false;
            qp.lowerY[row] = EPSILON;
            qp.upperY[row] = M_PI - EPSILON;
        }
    }

    // Build A matrix and adjust bounds by subtracting original angles.
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(std::max(0, numanz));

    int a = 0;
    for (HalfEdgeCIter he = mesh.halfEdges.begin(); he != mesh.halfEdges.end(); he++) {
        if (!he->onBoundary) {
            const double angle = he->angle();

            if (a >= variables) return false;
            qp.lowerX[a] = EPSILON - angle;
            qp.upperX[a] = (M_PI - EPSILON) - angle;

            // face sum contribution: x_a + ... = M_PI - sum(original_angles)
            const int fIdx = he->face->index;
            if (fIdx < 0 || fIdx >= faceShift) return false;
            qp.lowerY[fIdx] -= angle;
            qp.upperY[fIdx] -= angle;
            triplets.emplace_back(fIdx, a, 1.0);

            // vertex sum contribution
            const int vIdx = he->next->next->vertex->index + faceShift;
            if (vIdx < 0 || vIdx >= constraints) return false;
            qp.lowerY[vIdx] -= angle;
            qp.upperY[vIdx] -= angle;
            triplets.emplace_back(vIdx, a, 1.0);

            // local Delaunay constraint contribution (only for interior edges)
            if (!he->edge->isBoundary()) {
                const int eRow = eIntIndices[he->edge->index] + edgeShift;
                if (eRow < 0 || eRow >= constraints) return false;
                qp.lowerY[eRow] -= angle;
                qp.upperY[eRow] -= angle;
                triplets.emplace_back(eRow, a, 1.0);
            }

            a++;
        }
    }
    if (a != variables) return false;

    qp.A.resize(constraints, variables);
    qp.A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd corrections;
    QPsolver::MinNormBoxQPOptions opts;
    opts.admmMaxIters = 250;
    opts.cgMaxIters = 120;
    opts.rhoY = 1.0;
    opts.rhoX = 1.0;

    if (!QPsolver::QPSolver::solveMinNormBoxQPIdentityQ(qp, corrections, opts)) return false;
    if (corrections.size() != variables) return false;

    // Apply corrections back to angles/thetas (same as setThetas(), but without MOSEK memory).
    int b = 0;
    for (HalfEdgeCIter he = mesh.halfEdges.begin(); he != mesh.halfEdges.end(); he++) {
        if (!he->onBoundary) angles[he->index] = he->angle() + corrections[b++];
        else angles[he->index] = 0.0;
    }
    if (b != variables) return false;

    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter he = e->he;
        thetas[e->index] = M_PI - angles[he->index] - angles[he->flip->index];
    }

    return true;
}

double ImLi2Sum(double dp, double theta)
{
    double tStar = M_PI - theta;
    double x = 2*atan(tanh(0.5*dp) * tan(0.5*tStar));
    
    return x*dp + Cl2(x + tStar) + Cl2(-x + tStar) - Cl2(2.0*tStar);
}

double fe(double dp, double theta)
{
    return atan2(sin(theta), exp(dp) - cos(theta));
}

void CirclePatterns::computeEnergy(double& energy, const Eigen::VectorXd& rho)
{
    energy = 0.0;
    
    // sum over edges
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        int fk = e->he->face->index;
        
        if (e->isBoundary()) {
            energy -= 2*(M_PI - thetas[e->index]) * rho[fk];
            
        } else {
            int fl = e->he->flip->face->index;
            energy += ImLi2Sum(rho[fk] - rho[fl], thetas[e->index]) -
                      (M_PI - thetas[e->index])*(rho[fk] + rho[fl]);
        }
    }
    
    // sum over faces
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) energy += 2*M_PI*rho[f->index];
    }
}

void CirclePatterns::computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& rho)
{
    // loop over faces
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fk = f->index;
            gradient[fk] = 2*M_PI;
            
            // sum of adjacent edges
            HalfEdgeCIter he = f->he;
            do {
                EdgeCIter e = he->edge;
                if (e->isBoundary()) {
                    gradient[fk] -= 2*(M_PI - thetas[e->index]);
                    
                } else {
                    HalfEdgeCIter h = e->he;
                    int fl = fk == (int)h->face->index ? h->flip->face->index : h->face->index;
                    gradient[fk] -= 2*fe(rho[fk] - rho[fl], thetas[e->index]);
                }
                
                he = he->next;
            } while (he != f->he);
        }
    }
}

void CirclePatterns::computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& rho)
{
    std::vector<Eigen::Triplet<double>> HTriplets;
    
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!e->isBoundary()) {
            int fk = e->he->face->index;
            int fl = e->he->flip->face->index;
                        
            double hessval = sin(thetas[e->index]) / (cosh(rho(fk) - rho(fl)) - cos(thetas[e->index]));
            HTriplets.push_back(Eigen::Triplet<double>(fk, fk, hessval + 1e-8));
            HTriplets.push_back(Eigen::Triplet<double>(fl, fl, hessval + 1e-8));
            HTriplets.push_back(Eigen::Triplet<double>(fk, fl, -hessval));
            HTriplets.push_back(Eigen::Triplet<double>(fl, fk, -hessval));
        }
    }
    
    hessian.setFromTriplets(HTriplets.begin(), HTriplets.end());
}

void CirclePatterns::setRadii()
{
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) radii[f->index] = exp(solver.x[f->index]);
    }
}

bool CirclePatterns::computeRadii()
{
    MeshHandle handle;
    handle.computeEnergy = std::bind(&CirclePatterns::computeEnergy, this, _1, _2);
    handle.computeGradient = std::bind(&CirclePatterns::computeGradient, this, _1, _2);
    handle.computeHessian = std::bind(&CirclePatterns::computeHessian, this, _1, _2);
    
    solver.handle = &handle;
    if (OptScheme == GRAD_DESCENT) solver.gradientDescent();
    else if (OptScheme == NEWTON) solver.newton();
    else if (OptScheme == TRUST_REGION) solver.trustRegion();
    else solver.lbfgs();
    
    // set radii
    setRadii();
    
    return true;
}

void CirclePatterns::computeAnglesAndEdgeLengths(Eigen::VectorXd& lengths)
{
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter h1 = e->he;
        
        if (e->isBoundary()) {
            angles[h1->index] = M_PI - thetas[e->index];
        
        } else {
            HalfEdgeCIter h2 = h1->flip;
            double dp = log(radii[h1->face->index]) - log(radii[h2->face->index]);
            angles[h1->index] = fe(dp, thetas[e->index]);
            angles[h2->index] = fe(-dp, thetas[e->index]);
        }
        
        lengths[e->index] = 2.0*radii[h1->face->index]*sin(angles[h1->index]);
    }
}

void CirclePatterns::performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir,
                                       Eigen::VectorXd& lengths, std::unordered_map<int, bool>& visited,
                                       std::stack<EdgeCIter>& stack)
{
    if (!he->onBoundary) {
        int fIdx = he->face->index;
        if (visited.find(fIdx) == visited.end()) {
            HalfEdgeCIter next = he->next;
            HalfEdgeCIter prev = he->next->next;
            
            // compute new uv position
            double angle = angles[next->index];
            Eigen::Vector2d newDir = {cos(angle)*dir[0] - sin(angle)*dir[1],
                                      sin(angle)*dir[0] + cos(angle)*dir[1]};
            
            prev->vertex->uv = he->vertex->uv + newDir*lengths[prev->edge->index];
            
            // mark face as visited
            visited[fIdx] = true;
            
            // push edges onto stack
            stack.push(next->edge);
            stack.push(prev->edge);
        }
    }
}

void CirclePatterns::setUVs()
{
    // compute edge lengths
    Eigen::VectorXd lengths(mesh.edges.size());
    computeAnglesAndEdgeLengths(lengths);
    
    // push any edge
    std::stack<EdgeCIter> stack;
    EdgeCIter e = mesh.edges.begin();
    stack.push(e);
    e->he->vertex->uv = Eigen::Vector2d::Zero();
    e->he->next->vertex->uv = Eigen::Vector2d(lengths[e->index], 0);
    
    // perform layout
    std::unordered_map<int, bool> visited;
    while (!stack.empty()) {
        EdgeCIter e = stack.top();
        stack.pop();
        
        HalfEdgeCIter h1 = e->he;
        HalfEdgeCIter h2 = h1->flip;
        
        // compute edge vector
        Eigen::Vector2d dir = h2->vertex->uv - h1->vertex->uv;
        dir.normalize();
        
        performFaceLayout(h1, dir, lengths, visited, stack);
        performFaceLayout(h2, -dir, lengths, visited, stack);
    }
    
    normalize();
}

void CirclePatterns::parameterize()
{
    // set interior edge indices
    int eIdx = 0;
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!e->isBoundary()) eIntIndices[e->index] = eIdx++;
        else {
            eIntIndices[e->index] = -1;
            imaginaryHe++;
        }
    }
    
    // compute angles
    if (!computeAngles()) {
        std::cout << "Unable to compute angles" << std::endl;
        return;
    }
    
    // compute radii
    if (!computeRadii()) {
        std::cout << "Unable to compute radii" << std::endl;
        return;
    }
    
    // set uvs
    setUVs();
}
