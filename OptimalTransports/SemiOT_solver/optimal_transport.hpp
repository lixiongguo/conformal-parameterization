// MongeAmpere++  (geogram port)
// Original Copyright (C) 2014 Quentin Merigot, CNRS (GPLv3+)
//
// Damped Newton solver for semi-discrete optimal transport. CGAL-free: the
// source is a grid density (otgeo::DensityGrid) and the geometry comes from the
// geogram-backed power diagram via kantorovich().

#ifndef MA_OPTIMAL_TRANSPORT_HPP
#define MA_OPTIMAL_TRANSPORT_HPP

#include <Eigen/SparseCholesky>

#include <cassert>
#include <cstddef>
#include <iostream>

#include "kantorovich.hpp"

namespace MA {
struct Statistics {
    size_t niter;
    size_t neval;
};

template <class SparseMatrix, class Vector>
Vector solve_laplacian_matrix(const SparseMatrix& h, const Vector& g,
                              bool verbose = false) {
    size_t N = h.rows();
    assert(N == h.cols());
    auto v = h.diagonal();
    if (v.minCoeff() == 0) {
        size_t i;
        std::cerr << "Error: hessian of Kantorovich's functional "
                  << "is not invertible:\n";
        std::cerr << "diag = " << v.head(10) << " ... in [" << v.minCoeff()
                  << "," << v.minCoeff(&i) << "]\n";
        std::cerr << "minCoeff => " << i << "\n";
    }

    // remove last row and column so that the linear system is invertible
    Vector gs = g.head(N - 1);
    SparseMatrix hs = h.block(0, 0, N - 1, N - 1);  // top-left submatrix
    Eigen::SimplicialLLT<SparseMatrix> solver(hs);
    Vector ds = solver.solve(gs);

    double err = (hs * ds - gs).norm();
    if (err > 1e-7) {
        std::cerr << "WARNING: in solve_laplacian_matrix: err=" << err << "\n";
    }

    Vector d(N);
    d.head(N - 1) = ds;
    d(N - 1) = 0;
    return d;
}

// Solve semi-discrete OT between a grid source density and Dirac sites X with
// prescribed `masses`. The result `x` holds the Kantorovich weights.
template <class Matrix, class Vector>
void ot_solve(const otgeo::DensityGrid& source, const Matrix& X,
              const Vector& masses,
              Vector& x,  // result and initial guess
              double eps_g = 1e-7, size_t maxiter = 100, bool verbose = true,
              struct Statistics* stats = 0) {
    typedef Eigen::SparseMatrix<double> SparseMatrix;

    size_t neval = 0, niter = 0;
    size_t N = X.rows();
    assert(X.cols() == 2);
    assert(static_cast<size_t>(masses.rows()) == N);

    auto f = [&](const Vector& xx,
                 Vector& m,  // masses of the power cells
                 Vector& g,  // gradient
                 SparseMatrix& h) {
        ++neval;
        double r = kantorovich(source, X, xx, g, h);
        m = g;
        g = g - masses;
        return r - masses.dot(xx);
    };

    if (static_cast<size_t>(x.size()) != N) {
        x = Vector::Zero(N);
    }
    Vector g, m;
    SparseMatrix h;
    double fx = f(x, m, g, h);

    double eps0 = std::min(m.minCoeff(), masses.minCoeff()) / 2;
    if (eps0 <= 0) {
        std::cerr << "Error: computed minimum mass is non-positive\n";
        size_t i;
        m.minCoeff(&i);
        std::cerr << "This is because the Laguerre cell for the "
                  << "point [" << X(i, 0) << ", " << X(i, 1) << "],"
                  << " i=" << i << " is empty.\n";
        return;
    }

    while (g.norm() >= eps_g && niter++ <= maxiter) {
        Vector d = -solve_laplacian_matrix(h, g);

        double alpha = 1;
        Vector x0 = x;
        double n0 = g.norm();
        size_t nlinesearch = 0;

        while (1) {
            x = x0 + alpha * d;
            fx = f(x, m, g, h);
            if (m.minCoeff() >= eps0 && g.norm() <= (1 - alpha / 2) * n0) break;
            alpha *= .5;
            if (verbose) {
                std::cerr << "subit " << niter << "." << (nlinesearch++)
                          << ": min(masses)=" << m.minCoeff() << "\n";
            }
        }
        if (verbose) {
            std::cerr << "it " << niter << ":"
                      << " f=" << fx << " |df|=" << g.norm()
                      << " min(m)=" << masses.minCoeff() << " tau = " << alpha
                      << " eval = " << neval << "\n";
        }
    }

    if (stats) {
        stats->niter = niter;
        stats->neval = neval;
    }
}
}  // namespace MA

#endif
