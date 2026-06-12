#pragma once

#include <Eigen/Core>
#include <cmath>
#include <complex>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace vf {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEps = 1e-20;

inline uint64_t edgeKey(int a, int b) {
    const uint32_t lo = static_cast<uint32_t>(std::min(a, b));
    const uint32_t hi = static_cast<uint32_t>(std::max(a, b));
    return (static_cast<uint64_t>(lo) << 32) | hi;
}

struct FaceBasis {
    Eigen::Vector3d x = Eigen::Vector3d::UnitX();
    Eigen::Vector3d y = Eigen::Vector3d::UnitY();
    Eigen::Vector3d n = Eigen::Vector3d::UnitZ();
};

inline FaceBasis buildFaceBasis(
    const Eigen::Vector3d& p0,
    const Eigen::Vector3d& p1,
    const Eigen::Vector3d& p2) {
    FaceBasis b;
    b.x = p1 - p0;
    if (b.x.norm() < kEps) b.x = Eigen::Vector3d::UnitX();
    b.x.normalize();

    b.n = (p1 - p0).cross(p2 - p0);
    if (b.n.norm() < kEps) b.n = Eigen::Vector3d::UnitZ();
    b.n.normalize();

    b.y = b.n.cross(b.x);
    if (b.y.norm() < kEps) b.y = Eigen::Vector3d::UnitY();
    b.y.normalize();
    return b;
}

/// LC transport factor: conj(edge)^N with edge = (pos_vb - pos_va) in face basis.
inline std::complex<double> edgeConjPowN(
    const Eigen::Vector3d& posA,
    const Eigen::Vector3d& posB,
    const FaceBasis& basis,
    int n) {
    Eigen::Vector3d e = posB - posA;
    const double len = e.norm();
    if (len < kEps) e = basis.x;
    else e /= len;

    const std::complex<double> c(e.dot(basis.x), e.dot(basis.y));
    const std::complex<double> conjC = std::conj(c);
    const double ang = std::atan2(conjC.imag(), conjC.real()) * static_cast<double>(n);
    return {std::cos(ang), std::sin(ang)};
}

inline std::complex<double> directionConjPowN(
    const Eigen::Vector3d& dir,
    const FaceBasis& basis,
    int n) {
    Eigen::Vector3d d = dir;
    const double len = d.norm();
    if (len < kEps) d = basis.x;
    else d /= len;

    const std::complex<double> c(d.dot(basis.x), d.dot(basis.y));
    const std::complex<double> conjC = std::conj(c);
    const double ang = std::atan2(conjC.imag(), conjC.real()) * static_cast<double>(n);
    return {std::cos(ang), std::sin(ang)};
}

/// Levi-Civita parallel transport angle from face a to b across edge va->vb.
inline double lcTransportAngle(
    const FaceBasis& a,
    const FaceBasis& b,
    const Eigen::Vector3d& posA,
    const Eigen::Vector3d& posB) {
    Eigen::Vector3d t = posB - posA;
    const double tl = t.norm();
    if (tl < kEps) return 0.0;
    t /= tl;

    Eigen::Vector3d ua = a.n.cross(t);
    Eigen::Vector3d ub = b.n.cross(t);
    const double la = ua.norm();
    const double lb = ub.norm();
    if (la < kEps || lb < kEps) return 0.0;
    ua /= la;
    ub /= lb;

    const double sinb = ua.cross(ub).dot(t);
    const double cosb = ua.dot(ub);
    return std::atan2(sinb, cosb);
}

}  // namespace vf
