#include "CauchyCoordinates.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cauchy {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
const Complex kCauchyFactor(0.0, -1.0 / (2.0 * kPi));  // 1 / (2 pi i)

int nextIndex(int i, int n) {
    return (i + 1) % n;
}

void validatePolygonAndQueries(const VecC& polygon, const VecC& query_points) {
    if (polygon.size() < 3) {
        throw CauchyError("polygon must contain at least three vertices");
    }
    if (query_points.size() == 0) {
        throw CauchyError("query_points must not be empty");
    }
    for (int i = 0; i < polygon.size(); ++i) {
        const Complex edge = polygon[nextIndex(i, static_cast<int>(polygon.size()))] - polygon[i];
        if (std::abs(edge) <= std::numeric_limits<double>::epsilon()) {
            throw CauchyError("polygon contains a zero-length edge");
        }
    }
}

double pointSegmentDistance(const Complex& p, const Complex& a, const Complex& b) {
    const Complex ab = b - a;
    const double len2 = std::norm(ab);
    if (len2 <= std::numeric_limits<double>::epsilon()) {
        return std::abs(p - a);
    }

    const Complex ap = p - a;
    const double t = std::max(0.0, std::min(1.0, (ap.real() * ab.real() + ap.imag() * ab.imag()) / len2));
    return std::abs(p - (a + t * ab));
}

void validateAwayFromBoundary(const VecC& polygon, const Complex& z, const Options& options) {
    for (int j = 0; j < polygon.size(); ++j) {
        const Complex a = polygon[j];
        const Complex b = polygon[nextIndex(j, static_cast<int>(polygon.size()))];
        if (pointSegmentDistance(z, a, b) <= options.boundary_epsilon) {
            throw CauchyError("query point lies on or too close to polygon boundary");
        }
    }
}

}  // namespace

MatC computeCauchyCoordinates(const VecC& polygon, const VecC& query_points, const Options& options) {
    validatePolygonAndQueries(polygon, query_points);

    const int n = static_cast<int>(polygon.size());
    const int m = static_cast<int>(query_points.size());
    MatC C = MatC::Zero(m, n);

    for (int p = 0; p < m; ++p) {
        const Complex z = query_points[p];
        validateAwayFromBoundary(polygon, z, options);

        for (int j = 0; j < n; ++j) {
            const int jp = nextIndex(j, n);
            const Complex v0 = polygon[j];
            const Complex v1 = polygon[jp];
            const Complex a = v1 - v0;
            const Complex b0 = v0 - z;
            const Complex b1 = v1 - z;
            const Complex R = std::log(b1 / b0);

            C(p, j) += kCauchyFactor * (b1 * R - a) / a;
            C(p, jp) += kCauchyFactor * (a - b0 * R) / a;
        }
    }

    return C;
}

MatC cauchyCoordinates(const VecC& polygon, const VecC& query_points, const Options& options) {
    return computeCauchyCoordinates(polygon, query_points, options);
}

MatC computeCauchyCoordinateDerivatives(const VecC& polygon, const VecC& query_points, const Options& options) {
    validatePolygonAndQueries(polygon, query_points);

    const int n = static_cast<int>(polygon.size());
    const int m = static_cast<int>(query_points.size());
    MatC D = MatC::Zero(m, n);

    for (int p = 0; p < m; ++p) {
        const Complex z = query_points[p];
        validateAwayFromBoundary(polygon, z, options);

        for (int j = 0; j < n; ++j) {
            const int jp = nextIndex(j, n);
            const Complex v0 = polygon[j];
            const Complex v1 = polygon[jp];
            const Complex a = v1 - v0;
            const Complex b0 = v0 - z;
            const Complex b1 = v1 - z;
            const Complex R = std::log(b1 / b0);

            D(p, j) += kCauchyFactor * (-R + a / b0) / a;
            D(p, jp) += kCauchyFactor * (R - a / b1) / a;
        }
    }

    return D;
}

MatC derivativesOfCauchyCoord(const VecC& polygon, const VecC& query_points, const Options& options) {
    return computeCauchyCoordinateDerivatives(polygon, query_points, options);
}

MatC computeCauchyCoordinateSecondDerivatives(const VecC& polygon, const VecC& query_points, const Options& options) {
    validatePolygonAndQueries(polygon, query_points);

    const int n = static_cast<int>(polygon.size());
    const int m = static_cast<int>(query_points.size());
    MatC E = MatC::Zero(m, n);

    for (int p = 0; p < m; ++p) {
        const Complex z = query_points[p];
        validateAwayFromBoundary(polygon, z, options);

        for (int j = 0; j < n; ++j) {
            const int jp = nextIndex(j, n);
            const Complex v0 = polygon[j];
            const Complex v1 = polygon[jp];
            const Complex b0 = v0 - z;
            const Complex b1 = v1 - z;

            E(p, j) += kCauchyFactor * (Complex(1.0, 0.0) / (b0 * b0) - Complex(1.0, 0.0) / (b0 * b1));
            E(p, jp) += kCauchyFactor * (Complex(1.0, 0.0) / (b0 * b1) - Complex(1.0, 0.0) / (b1 * b1));
        }
    }

    return E;
}

MatC secondDerivativesOfCauchyCoord(const VecC& polygon, const VecC& query_points, const Options& options) {
    return computeCauchyCoordinateSecondDerivatives(polygon, query_points, options);
}

VecC evaluateCauchyMap(
    const VecC& polygon,
    const VecC& target_vertices,
    const VecC& query_points,
    const Options& options) {
    if (target_vertices.size() != polygon.size()) {
        throw CauchyError("target_vertices must have one value per polygon vertex");
    }
    return computeCauchyCoordinates(polygon, query_points, options) * target_vertices;
}

PrecisionError checkPrecision(const VecC& polygon, const VecC& query_points, const Options& options) {
    const MatC C = computeCauchyCoordinates(polygon, query_points, options);
    PrecisionError error;
    for (int p = 0; p < query_points.size(); ++p) {
        const Complex constant = C.row(p).sum();
        Complex linear(0.0, 0.0);
        for (int j = 0; j < polygon.size(); ++j) {
            linear += C(p, j) * polygon[j];
        }
        error.constant_precision = std::max(error.constant_precision, std::abs(constant - Complex(1.0, 0.0)));
        error.linear_precision = std::max(error.linear_precision, std::abs(linear - query_points[p]));
    }
    return error;
}

}  // namespace cauchy
