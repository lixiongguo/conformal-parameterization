#pragma once

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <vector>

#include "types.hpp"

namespace aht_ot {

/// 5-point Laplacian + SparseLU — discrete inverse Laplacian (sec 7.3, eq 26-27)
struct PoissonSolver5pt {
    int nx = 0;
    int ny = 0;
    double h2 = 1.0;
    Eigen::SparseMatrix<double> L;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;

    PoissonSolver5pt(int nx_, int ny_, double square_edge_len)
        : nx(nx_), ny(ny_) {
        const double h = square_edge_len / static_cast<double>(nx - 1);
        h2 = h * h;

        const int N = nx * ny;
        L.resize(N, N);
        std::vector<Eigen::Triplet<double>> trips;
        trips.reserve(static_cast<size_t>(N) * 5);

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                const int row = i * ny + j;
                const bool interior =
                    (i > 0 && i < nx - 1 && j > 0 && j < ny - 1);
                if (interior) {
                    trips.emplace_back(row, row, 4.0);
                    trips.emplace_back(row, (i - 1) * ny + j, -1.0);
                    trips.emplace_back(row, (i + 1) * ny + j, -1.0);
                    trips.emplace_back(row, i * ny + (j - 1), -1.0);
                    trips.emplace_back(row, i * ny + (j + 1), -1.0);
                } else {
                    trips.emplace_back(row, row, 1.0);  // Dirichlet f = 0
                }
            }
        }
        L.setFromTriplets(trips.begin(), trips.end());
        L.makeCompressed();
        lu.analyzePattern(L);
        lu.factorize(L);
    }

    /// Solve Δf = rhs with f|∂Ω = 0
    Eigen::VectorXd solve(const Eigen::VectorXd& rhs) const {
        return lu.solve(rhs * h2);
    }
};

}  // namespace aht_ot
