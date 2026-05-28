#include <emscripten.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
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

namespace {

constexpr double kEps = 1e-20;
constexpr int kRosyN = 4;

struct FaceBasis {
  Eigen::Vector3d x;
  Eigen::Vector3d y;
  Eigen::Vector3d n;
};

struct EdgeUse {
  int face = -1;
  int va = -1;
  int vb = -1;
};

struct InteriorEdgeEq {
  int fi = -1;
  int fj = -1;
  int va = -1;
  int vb = -1;
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
int g_lastAlgorithm = 0;  // 1: trivial, 2: complex-poly

static inline uint64_t edgeKey(int a, int b) {
  const uint32_t lo = static_cast<uint32_t>(std::min(a, b));
  const uint32_t hi = static_cast<uint32_t>(std::max(a, b));
  return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
}

static std::vector<FaceBasis> buildFaceBases(const float* V_ptr, const int* F_ptr, int F_rows) {
  std::vector<FaceBasis> bases(static_cast<size_t>(F_rows));
  for (int fi = 0; fi < F_rows; ++fi) {
    const int i0 = F_ptr[fi * 3];
    const int i1 = F_ptr[fi * 3 + 1];
    const int i2 = F_ptr[fi * 3 + 2];

    const Eigen::Vector3d p0(V_ptr[i0 * 3], V_ptr[i0 * 3 + 1], V_ptr[i0 * 3 + 2]);
    const Eigen::Vector3d p1(V_ptr[i1 * 3], V_ptr[i1 * 3 + 1], V_ptr[i1 * 3 + 2]);
    const Eigen::Vector3d p2(V_ptr[i2 * 3], V_ptr[i2 * 3 + 1], V_ptr[i2 * 3 + 2]);

    FaceBasis b;
    b.x = p1 - p0;
    if (b.x.norm() < 1e-12) b.x = Eigen::Vector3d::UnitX();
    b.x.normalize();

    b.n = (p1 - p0).cross(p2 - p0);
    if (b.n.norm() < 1e-12) b.n = Eigen::Vector3d::UnitZ();
    b.n.normalize();

    b.y = b.n.cross(b.x);
    if (b.y.norm() < 1e-12) b.y = Eigen::Vector3d::UnitY();
    b.y.normalize();

    bases[static_cast<size_t>(fi)] = b;
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
      edgeUses[edgeKey(va, vb)].push_back({fi, va, vb});
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
      const int va = std::min(u0.va, u0.vb);
      const int vb = std::max(u0.va, u0.vb);
      interiorEdges.push_back({u0.face, u1.face, va, vb});
      faceAdj[static_cast<size_t>(u0.face)].push_back({u1.face, va, vb});
      faceAdj[static_cast<size_t>(u1.face)].push_back({u0.face, va, vb});
    } else if (uses.size() == 1) {
      const EdgeUse& u = uses[0];
      boundaryEdges.push_back({u.face, u.va, u.vb});
    }
  }
}

static std::complex<double> unitComplexPowNFromDir(
    const Eigen::Vector3d& dir,
    const FaceBasis& b,
    int n) {
  Eigen::Vector3d d = dir;
  const double len = d.norm();
  if (len < 1e-12) d = b.x;
  else d /= len;

  const std::complex<double> c(d.dot(b.x), d.dot(b.y));
  const double ang = std::atan2(c.imag(), c.real()) * static_cast<double>(n);
  return {std::cos(ang), std::sin(ang)};
}

static std::complex<double> edgeComplexConjPowN(
    const float* V_ptr,
    const FaceBasis& b,
    int va,
    int vb,
    int n) {
  Eigen::Vector3d e(
      V_ptr[va * 3] - V_ptr[vb * 3],
      V_ptr[va * 3 + 1] - V_ptr[vb * 3 + 1],
      V_ptr[va * 3 + 2] - V_ptr[vb * 3 + 2]);
  const double l = e.norm();
  if (l < 1e-12) e = b.x;
  else e /= l;

  std::complex<double> c(e.dot(b.x), e.dot(b.y));
  c = std::conj(c);
  const double ang = std::atan2(c.imag(), c.real()) * static_cast<double>(n);
  return {std::cos(ang), std::sin(ang)};
}

static double parallelTransportTau(
    const FaceBasis& a,
    const FaceBasis& b,
    const float* V_ptr,
    int va,
    int vb) {
  Eigen::Vector3d t(
      V_ptr[vb * 3] - V_ptr[va * 3],
      V_ptr[vb * 3 + 1] - V_ptr[va * 3 + 1],
      V_ptr[vb * 3 + 2] - V_ptr[va * 3 + 2]);
  const double tl = t.norm();
  if (tl < 1e-12) return 0.0;
  t /= tl;

  Eigen::Vector3d ua = a.n.cross(t);
  Eigen::Vector3d ub = b.n.cross(t);
  const double la = ua.norm();
  const double lb = ub.norm();
  if (la < 1e-12 || lb < 1e-12) return 0.0;
  ua /= la;
  ub /= lb;

  const double sinb = ua.cross(ub).dot(t);
  const double cosb = ua.dot(ub);
  return std::atan2(sinb, cosb);
}

}  // namespace

extern "C" {

EMSCRIPTEN_KEEPALIVE
int compute_trivial_nrosy_field(float* V_ptr, int V_rows, int* F_ptr, int F_rows) {
  (void)V_rows;
  auto t0 = std::chrono::high_resolution_clock::now();
  g_theta.clear();
  g_lastMs = 0.0;
  g_lastAlgorithm = 0;
  if (!V_ptr || !F_ptr || F_rows <= 0) return -1;

  try {
    const std::vector<FaceBasis> bases = buildFaceBases(V_ptr, F_ptr, F_rows);
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
        const FaceBasis& fb = bases[static_cast<size_t>(fi)];
        for (const Neighbor& nb : faceAdj[static_cast<size_t>(fi)]) {
          const int gj = nb.other;
          if (visited[static_cast<size_t>(gj)]) continue;
          const FaceBasis& gb = bases[static_cast<size_t>(gj)];
          const double tau = parallelTransportTau(fb, gb, V_ptr, nb.va, nb.vb);
          g_theta[static_cast<size_t>(gj)] = static_cast<float>(g_theta[static_cast<size_t>(fi)] + tau);
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
int compute_4rosy_field(float* V_ptr, int V_rows, int* F_ptr, int F_rows) {
  (void)V_rows;
  auto t0 = std::chrono::high_resolution_clock::now();
  g_theta.clear();
  g_lastMs = 0.0;
  g_lastAlgorithm = 0;
  if (!V_ptr || !F_ptr || F_rows <= 0) return -1;

  try {
    const std::vector<FaceBasis> bases = buildFaceBases(V_ptr, F_ptr, F_rows);

    std::vector<InteriorEdgeEq> interiorEdges;
    std::vector<BoundaryConstraint> boundaryEdges;
    std::vector<std::vector<Neighbor>> faceAdj;
    buildTopology(F_ptr, F_rows, interiorEdges, boundaryEdges, faceAdj);

    using SpMat = Eigen::SparseMatrix<std::complex<double>>;
    using Trip = Eigen::Triplet<std::complex<double>>;

    const int rowCount = static_cast<int>(interiorEdges.size()) +
                         static_cast<int>(boundaryEdges.size()) +
                         ((boundaryEdges.empty() && F_rows > 0) ? 1 : 0);
    std::vector<Trip> trips;
    trips.reserve(static_cast<size_t>(rowCount) * 2);
    Eigen::VectorXcd b = Eigen::VectorXcd::Zero(rowCount);

    int r = 0;
    for (const InteriorEdgeEq& e : interiorEdges) {
      const std::complex<double> ef = edgeComplexConjPowN(V_ptr, bases[static_cast<size_t>(e.fi)], e.va, e.vb, kRosyN);
      const std::complex<double> eg = edgeComplexConjPowN(V_ptr, bases[static_cast<size_t>(e.fj)], e.va, e.vb, kRosyN);
      trips.emplace_back(r, e.fi, ef);
      trips.emplace_back(r, e.fj, -eg);
      ++r;
    }

    for (const BoundaryConstraint& bc : boundaryEdges) {
      const FaceBasis& fb = bases[static_cast<size_t>(bc.fi)];
      Eigen::Vector3d p0(V_ptr[bc.va * 3], V_ptr[bc.va * 3 + 1], V_ptr[bc.va * 3 + 2]);
      Eigen::Vector3d p1(V_ptr[bc.vb * 3], V_ptr[bc.vb * 3 + 1], V_ptr[bc.vb * 3 + 2]);
      Eigen::Vector3d dir = p0 - p1;
      const double dl = dir.norm();
      if (dl < 1e-12) dir = fb.x;
      else dir /= dl;
      const Eigen::Vector3d pp = dir.cross(fb.n);
      const std::complex<double> zn = unitComplexPowNFromDir(pp, fb, kRosyN);
      trips.emplace_back(r, bc.fi, std::complex<double>(1.0, 0.0));
      b[r] = zn;
      ++r;
    }

    if (boundaryEdges.empty() && F_rows > 0) {
      trips.emplace_back(r, 0, std::complex<double>(1.0, 0.0));
      b[r] = std::complex<double>(1.0, 0.0);
      ++r;
    }

    SpMat A(r, F_rows);
    A.setFromTriplets(trips.begin(), trips.end());
    SpMat AH = A.adjoint();
    SpMat ATA = AH * A;
    Eigen::VectorXcd rhs = AH * b.head(r);

    Eigen::SimplicialLDLT<SpMat> solver;
    solver.compute(ATA);
    if (solver.info() != Eigen::Success) return -2;

    Eigen::VectorXcd x = solver.solve(rhs);
    if (solver.info() != Eigen::Success) return -3;

    g_theta.resize(static_cast<size_t>(F_rows), 0.0f);
    for (int fi = 0; fi < F_rows; ++fi) {
      const std::complex<double> z = x[fi];
      if (std::isfinite(z.real()) && std::isfinite(z.imag()) && std::abs(z) > kEps) {
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
