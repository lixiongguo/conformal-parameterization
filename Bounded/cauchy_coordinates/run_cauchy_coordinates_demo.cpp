#include "CauchyCoordinates.hpp"

#include <Eigen/Core>

#include <iostream>

int main() {
    using namespace cauchy;

    VecC cage(4);
    cage << Complex(0.0, 0.0),
            Complex(1.0, 0.0),
            Complex(1.0, 1.0),
            Complex(0.0, 1.0);

    VecC queries(4);
    queries << Complex(0.5, 0.5),
               Complex(0.25, 0.25),
               Complex(0.75, 0.35),
               Complex(0.4, 0.8);

    const MatC C = computeCauchyCoordinates(cage, queries);
    const MatC D = computeCauchyCoordinateDerivatives(cage, queries);
    const PrecisionError precision = checkPrecision(cage, queries);

    VecC target(4);
    target << Complex(0.0, 0.0),
              Complex(1.2, 0.1),
              Complex(0.9, 1.1),
              Complex(-0.1, 0.8);
    const VecC mapped = evaluateCauchyMap(cage, target, queries);

    std::cout << "Cauchy coordinates matrix: " << C.rows() << " x " << C.cols() << '\n';
    std::cout << "Derivative matrix        : " << D.rows() << " x " << D.cols() << '\n';
    std::cout << "constant precision error : " << precision.constant_precision << '\n';
    std::cout << "linear precision error   : " << precision.linear_precision << '\n';

    for (int i = 0; i < queries.size(); ++i) {
        std::cout << "query " << i
                  << " z=(" << queries[i].real() << ", " << queries[i].imag() << ")"
                  << " mapped=(" << mapped[i].real() << ", " << mapped[i].imag() << ")"
                  << " sumC=" << C.row(i).sum()
                  << '\n';
    }

    const bool ok = precision.constant_precision < 1e-10 && precision.linear_precision < 1e-10;
    std::cout << (ok ? "PASS\n" : "FAIL\n");
    return ok ? 0 : 1;
}
