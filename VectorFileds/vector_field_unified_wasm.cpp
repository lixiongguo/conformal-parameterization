#include <emscripten.h>

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "vector_field_topology.hpp"

namespace {

constexpr int kRosyN = 4;

struct EdgeUse {
  int face = -1;
  int va = -1;
  int vb = -1;
};

struct InteriorEdgeEq {
  int fi = -1;
  int fj = -1;
  int va_f = -1;
  int vb_f = -1;
  int va_g = -1;
  int vb_g = -1;
};

struct BoundaryConstraint {
  int fi = -1;
  int va = -1;
  int vb = -1;
};

struct Neighbor {
  int other = -1;
  int va = -1;
  int vb = -1;
};

std::vector<float> g_theta;
double g_lastMs = 0.0;
int g_lastAlgorithm = 0;

static std::vector<vf::FaceBasis> buildFaceBases(const float* V_ptr, const int* F_ptr, int F_rows) {
  std::vector<vf::FaceBasis> bases(static_cast<size_t>(F_rows));
  for (int fi = 0; fi < F_rows; ++fi) {
    const int i0 = F_ptr[fi * 3];
    const int i1 = F_ptr[fi * 3 + 1];
    const int i2 = F_ptr[fi * 3 + 2];
    const Eigen::Vector3d p0(V_ptr[i0 * 3], V_ptr[i0 * 3 + 1], V_ptr[i0 * 3 + 2]);
    const Eigen::Vector3d p1(V_ptr[i1 * 3], V_ptr[i1 * 3 + 1], V_ptr[i1 * 3 + 2]);
    const Eigen::Vector3d p2(V_ptr[i2 * 3], V_ptr[i2 * 3 + 1], V_ptr[i2 * 3 + 2]);
    bases[static_cast<size_t>(fi)] = vf::buildFaceBasis(p0, p1, p2);
  }
  return bases;
}

static void buildTopology(
    const int* F_ptr,
    int F_rows,
    std::vector<InteriorEdgeEq>& interiorEdges,
    std::vector<BoundaryConstraint>& boundaryEdges,
    std::vector<std::vector<Neighbor>>& faceAdj) {
  std::unordered_map<uint64_t, std::vector<EdgeUse>> edgeUses;
  edgeUses.reserve(static_cast<size_t>(F_rows) * 3);

  for (int fi = 0; fi < F_rows; ++fi) {
    const int tri[3] = {F_ptr[fi * 3], F_ptr[fi * 3 + 1], F_ptr[fi * 3 + 2]};
    for (int e = 0; e < 3; ++e) {
      const int va = tri[e];
      const int vb = tri[(e + 1) % 3];
      edgeUses[vf::edgeKey(va, vb)].push_back({fi, va, vb});
    }
  }

  interiorEdges.clear();
  boundaryEdges.clear();
  faceAdj.assign(static_cast<size_t>(F_rows), {});
  interiorEdges.reserve(static_cast<size_t>(F_rows) * 2);
  boundaryEdges.reserve(static_cast<size_t>(F_rows));

  for (const auto& kv : edgeUses) {
    const std::vector<EdgeUse>& uses = kv.second;
    if (uses.size() == 2) {
      const EdgeUse& u0 = uses[0];
      const EdgeUse& u1 = uses[1];
      interiorEdges.push_back({u0.face, u1.face, u0.va, u0.vb, u1.va, u1.vb});
      faceAdj[static_cast<size_t>(u0.face)].push_back({u1.face, u0.va, u0.vb});
      faceAdj[static_cast<size_t>(u1.face)].push_back({u0.face, u1.va, u1.vb});
    } else if (uses.size() == 1) {
      const EdgeUse& u = uses[0];
      boundaryEdges.push_back({u.face, u.va, u.vb});
    }
  }
}

static std::complex<double> edgeComplexConjPowN(
    const float* V_ptr,
    const vf::FaceBasis& b,
    int va,
    int vb,
    int n) {
  const Eigen::Vector3d posA(V_ptr[va * 3], V_ptr[va * 3 + 1], V_ptr[va * 3 + 2]);
  const Eigen::Vector3d posB(V_ptr[vb * 3], V_ptr[vb * 3 + 1], V_ptr[vb * 3 + 2]);
  return vf::edgeConjPowN(posA, posB, b, n);
}

}  // namespace

extern "C" {

EMSCRIPTEN_KEEPALIVE
int compute_lc_propagated_field(float* V_ptr, int V_rows, int* F_ptr, int F_rows) {
  (void)V_rows;
  auto t0 = std::chrono::high_resolution_clock::now();
  g_theta.clear();
  g_lastMs = 0.0;
  g_lastAlgorithm = 0;
  if (!V_ptr || !F_ptr || F_rows <= 0) return -1;

  try {
    const std::vector<vf::FaceBasis> bases = buildFaceBases(V_ptr, F_ptr, F_rows);
    std::vector<InteriorEdgeEq> interiorEdges;
    std::vector<BoundaryConstraint> boundaryEdges;
    std::vector<std::vector<Neighbor>> faceAdj;
    buildTopology(F_ptr, F_rows, interiorEdges, boundaryEdges, faceAdj);

    g_theta.assign(static_cast<size_t>(F_rows), 0.0f);
    std::vector<uint8_t> visited(static_cast<size_t>(F_rows), 0);
    std::vector<int> queue;
    queue.reserve(static_cast<size_t>(F_rows));

    for (int seed = 0; seed < F_rows; ++seed) {
      if (visited[static_cast<size_t>(seed)]) continue;
      visited[static_cast<size_t>(seed)] = 1;
      g_theta[static_cast<size_t>(seed)] = 0.0f;
      queue.clear();
      queue.push_back(seed);

      size_t qh = 0;
      while (qh < queue.size()) {
        const int fi = queue[qh++];
        const vf::FaceBasis& fb = bases[static_cast<size_t>(fi)];
        for (const Neighbor& nb : faceAdj[static_cast<size_t>(fi)]) {
          const int gj = nb.other;
          if (visited[static_cast<size_t>(gj)]) continue;
          const vf::FaceBasis& gb = bases[static_cast<size_t>(gj)];
          const Eigen::Vector3d posA(
              V_ptr[nb.va * 3], V_ptr[nb.va * 3 + 1], V_ptr[nb.va * 3 + 2]);
          const Eigen::Vector3d posB(
              V_ptr[nb.vb * 3], V_ptr[nb.vb * 3 + 1], V_ptr[nb.vb * 3 + 2]);
          const double tau = vf::lcTransportAngle(fb, gb, posA, posB);
          g_theta[static_cast<size_t>(gj)] =
              static_cast<float>(g_theta[static_cast<size_t>(fi)] + tau);
          visited[static_cast<size_t>(gj)] = 1;
          queue.push_back(gj);
        }
      }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    g_lastMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    g_lastAlgorithm = 1;
    return 0;
  } catch (...) {
    return -2;
  }
}

EMSCRIPTEN_KEEPALIVE
int compute_trivial_nrosy_field(float* V_ptr, int V_rows, int* F_ptr, int F_rows) {
  return compute_lc_propagated_field(V_ptr, V_rows, F_ptr, F_rows);
}

EMSCRIPTEN_KEEPALIVE
int compute_4rosy_field(float* V_ptr, int V_rows, int* F_ptr, int F_rows) {
  (void)V_rows;
  auto t0 = std::chrono::high_resolution_clock::now();
  g_theta.clear();
  g_lastMs = 0.0;
  g_lastAlgorithm = 0;
  if (!V_ptr || !F_ptr || F_rows <= 0) return -1;

  try {
    const std::vector<vf::FaceBasis> bases = buildFaceBases(V_ptr, F_ptr, F_rows);

    std::vector<InteriorEdgeEq> interiorEdges;
    std::vector<BoundaryConstraint> boundaryEdges;
    std::vector<std::vector<Neighbor>> faceAdj;
    buildTopology(F_ptr, F_rows, interiorEdges, boundaryEdges, faceAdj);

    using SpMat = Eigen::SparseMatrix<std::complex<double>>;
    using Trip = Eigen::Triplet<std::complex<double>>;

    const bool needGauge = boundaryEdges.empty() && F_rows > 0;
    const int rowCount = static_cast<int>(interiorEdges.size()) +
                         static_cast<int>(boundaryEdges.size()) +
                         (needGauge ? 1 : 0);
    std::vector<Trip> trips;
    trips.reserve(static_cast<size_t>(rowCount) * 2);
    Eigen::VectorXcd b = Eigen::VectorXcd::Zero(rowCount);

    int r = 0;
    for (const InteriorEdgeEq& e : interiorEdges) {
      const std::complex<double> ef =
          edgeComplexConjPowN(V_ptr, bases[static_cast<size_t>(e.fi)], e.va_f, e.vb_f, kRosyN);
      const std::complex<double> eg =
          edgeComplexConjPowN(V_ptr, bases[static_cast<size_t>(e.fj)], e.va_g, e.vb_g, kRosyN);
      trips.emplace_back(r, e.fi, ef);
      trips.emplace_back(r, e.fj, -eg);
      ++r;
    }

    for (const BoundaryConstraint& bc : boundaryEdges) {
      const vf::FaceBasis& fb = bases[static_cast<size_t>(bc.fi)];
      Eigen::Vector3d p0(V_ptr[bc.va * 3], V_ptr[bc.va * 3 + 1], V_ptr[bc.va * 3 + 2]);
      Eigen::Vector3d p1(V_ptr[bc.vb * 3], V_ptr[bc.vb * 3 + 1], V_ptr[bc.vb * 3 + 2]);
      Eigen::Vector3d edgeTan = p1 - p0;
      const double dl = edgeTan.norm();
      if (dl < 1e-12) edgeTan = fb.x;
      else edgeTan /= dl;
      Eigen::Vector3d inPlane = edgeTan - edgeTan.dot(fb.n) * fb.n;
      if (inPlane.norm() < 1e-12) inPlane = fb.y;
      else inPlane.normalize();
      const std::complex<double> zn = vf::directionConjPowN(inPlane, fb, kRosyN);
      trips.emplace_back(r, bc.fi, std::complex<double>(1.0, 0.0));
      b[r] = zn;
      ++r;
    }

    if (needGauge) {
      trips.emplace_back(r, 0, std::complex<double>(1.0, 0.0));
      b[r] = std::complex<double>(1.0, 0.0);
      ++r;
    }

    SpMat A(r, F_rows);
    A.setFromTriplets(trips.begin(), trips.end());
    SpMat AH = A.adjoint();
    SpMat ATA = AH * A;
    Eigen::VectorXcd rhs = AH * b.head(r);

    for (int fi = 0; fi < F_rows; ++fi) {
      ATA.coeffRef(fi, fi) += std::complex<double>(1e-8, 0.0);
    }

    Eigen::SimplicialLDLT<SpMat> solver;
    solver.compute(ATA);
    if (solver.info() != Eigen::Success) return -2;

    Eigen::VectorXcd x = solver.solve(rhs);
    if (solver.info() != Eigen::Success) return -3;

    g_theta.resize(static_cast<size_t>(F_rows), 0.0f);
    for (int fi = 0; fi < F_rows; ++fi) {
      const std::complex<double> z = x[fi];
      if (std::isfinite(z.real()) && std::isfinite(z.imag()) && std::abs(z) > vf::kEps) {
        g_theta[static_cast<size_t>(fi)] =
            static_cast<float>(std::atan2(z.imag(), z.real()) / static_cast<double>(kRosyN));
      }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    g_lastMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    g_lastAlgorithm = 2;
    return 0;
  } catch (...) {
    return -4;
  }
}

EMSCRIPTEN_KEEPALIVE
float* get_nrosy_theta() {
  return g_theta.empty() ? nullptr : g_theta.data();
}

EMSCRIPTEN_KEEPALIVE
int get_nrosy_theta_size() {
  return static_cast<int>(g_theta.size());
}

EMSCRIPTEN_KEEPALIVE
double get_nrosy_last_time_ms() {
  return g_lastMs;
}

EMSCRIPTEN_KEEPALIVE
int get_nrosy_last_algorithm() {
  return g_lastAlgorithm;
}

}  // extern "C"
