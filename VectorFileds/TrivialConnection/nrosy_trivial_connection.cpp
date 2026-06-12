#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "Mesh.h"

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Args {
  std::string inputObj;
  std::string outputDir;
  std::string singularityCsv;
  int N = 4;
  bool autoSingularity = false;
};

bool parseArgs(int argc, char** argv, Args& args) {
  if (argc < 3) return false;
  args.inputObj = argv[1];
  args.outputDir = argv[2];
  for (int i = 3; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--singularity-csv" && i + 1 < argc) {
      args.singularityCsv = argv[++i];
    } else if (a == "--n" && i + 1 < argc) {
      args.N = std::max(1, std::atoi(argv[++i]));
    } else if (a == "--auto-singularity") {
      args.autoSingularity = true;
    } else {
      std::cerr << "Unknown option: " << a << "\n";
      return false;
    }
  }
  return true;
}

void usage() {
  std::cout << "Usage:\n"
            << "  nrosy_trivial_connection.exe <input.obj> <output_dir> [options]\n\n"
            << "Options:\n"
            << "  --n <int>                 N in N-RoSy (default 4)\n"
            << "  --singularity-csv <path>  CSV: vertex,index  (index can be fractional)\n"
            << "  --auto-singularity        Auto generate singularities with sum index = chi\n";
}

bool isInteriorFace(FaceCIter f) { return !f->isBoundary(); }

int countInteriorFaces(const Mesh& mesh) {
  int cnt = 0;
  for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
    if (isInteriorFace(f)) cnt++;
  }
  return cnt;
}

double vertexCurvature(VertexCIter v) {
  double sumAngles = 0.0;
  HalfEdgeCIter h = v->he;
  do {
    if (!h->onBoundary) sumAngles += h->next->next->angle();
    h = h->flip->next;
  } while (h != v->he);
  return 2.0 * kPi - sumAngles;
}

void buildBasis(FaceCIter f, Eigen::Vector3d& e1, Eigen::Vector3d& e2, Eigen::Vector3d& n) {
  HalfEdgeCIter h = f->he;
  const Eigen::Vector3d p0 = h->vertex->position;
  const Eigen::Vector3d p1 = h->next->vertex->position;
  const Eigen::Vector3d p2 = h->next->next->vertex->position;
  e1 = p1 - p0;
  if (e1.norm() < 1e-12) e1 = Eigen::Vector3d::UnitX();
  e1.normalize();
  n = (p1 - p0).cross(p2 - p0);
  if (n.norm() < 1e-12) n = Eigen::Vector3d::UnitZ();
  n.normalize();
  e2 = n.cross(e1);
  if (e2.norm() < 1e-12) e2 = Eigen::Vector3d::UnitY();
  e2.normalize();
}

double lcTransport(HalfEdgeCIter h) {
  Eigen::Vector3d e1i, e2i, ni;
  Eigen::Vector3d e1j, e2j, nj;
  buildBasis(h->face, e1i, e2i, ni);
  buildBasis(h->flip->face, e1j, e2j, nj);
  const Eigen::Vector3d u = h->flip->vertex->position - h->vertex->position;
  const double thetaIJ = std::atan2(u.dot(e2i), u.dot(e1i));
  const double thetaJI = std::atan2(u.dot(e2j), u.dot(e1j));
  return -thetaIJ + thetaJI;
}

bool loadSingularity(const std::string& path, int N, std::vector<int>& q) {
  std::ifstream in(path);
  if (!in.is_open()) return false;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line.find("vertex") != std::string::npos) continue;
    size_t c = line.find(',');
    if (c == std::string::npos) continue;
    int v = std::stoi(line.substr(0, c));
    double idx = std::stod(line.substr(c + 1));  // index in units of full turn
    int qi = (int)std::llround(idx * N);
    if (v >= 0 && v < (int)q.size()) q[v] = qi;
  }
  return true;
}

void autoSingularity(int N, int chi, const Eigen::VectorXd& K, std::vector<int>& q) {
  int sumQTarget = N * chi;
  if (sumQTarget == 0) return;

  std::vector<std::pair<double, int>> ranked;
  ranked.reserve(q.size());
  for (int i = 0; i < (int)q.size(); ++i) ranked.push_back({std::abs(K[i]), i});
  std::sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

  int sign = sumQTarget > 0 ? 1 : -1;
  int remaining = std::abs(sumQTarget);
  for (size_t i = 0; i < ranked.size() && remaining > 0; ++i) {
    q[ranked[i].second] += sign;
    remaining--;
  }
}

}  // namespace

int main(int argc, char** argv) {
  Args args;
  if (!parseArgs(argc, argv, args)) {
    usage();
    return 1;
  }

  Mesh mesh;
  if (!mesh.read(args.inputObj)) {
    std::cerr << "Failed to read mesh: " << args.inputObj << "\n";
    return 2;
  }
  if (!mesh.boundaries.empty()) {
    std::cerr << "Only closed meshes are supported in this version.\n";
    return 3;
  }

  const int nV = (int)mesh.vertices.size();
  const int nE = (int)mesh.edges.size();
  const int nF = countInteriorFaces(mesh);
  const int chi = nV - nE + nF;
  std::cout << "Mesh V/E/F=" << nV << "/" << nE << "/" << nF
            << " chi=" << chi << " N=" << args.N << "\n";

  // Vertex curvatures (LC holonomy around dual cell)
  Eigen::VectorXd K = Eigen::VectorXd::Zero(nV);
  for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
    K[v->index] = vertexCurvature(v);
  }

  // q = N * index (integer)
  std::vector<int> q(nV, 0);
  if (!args.singularityCsv.empty()) {
    if (!loadSingularity(args.singularityCsv, args.N, q)) {
      std::cerr << "Failed to load singularity CSV: " << args.singularityCsv << "\n";
      return 4;
    }
  } else if (args.autoSingularity) {
    autoSingularity(args.N, chi, K, q);
  } else {
    std::cerr << "Provide --singularity-csv or --auto-singularity.\n";
    return 5;
  }

  int sumQ = 0;
  for (int x : q) sumQ += x;
  if (sumQ != args.N * chi) {
    std::cerr << "Gauss-Bonnet mismatch for N-RoSy: sum(q)=" << sumQ
              << " expected=" << args.N * chi << "\n";
    return 6;
  }

  // Solve minimum-norm connection correction phi on primal edges:
  //   B phi = 2π q/N - K
  std::vector<Eigen::Triplet<double>> Btrips;
  Btrips.reserve(nE * 2);
  for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
    const int i = e->he->vertex->index;
    const int j = e->he->flip->vertex->index;
    Btrips.emplace_back(i, e->index, 1.0);
    Btrips.emplace_back(j, e->index, -1.0);
  }
  Eigen::SparseMatrix<double> B(nV, nE);
  B.setFromTriplets(Btrips.begin(), Btrips.end());

  Eigen::VectorXd rhs(nV);
  for (int i = 0; i < nV; ++i) rhs[i] = 2.0 * kPi * ((double)q[i] / (double)args.N) - K[i];

  const int Nsys = nE + nV;
  std::vector<Eigen::Triplet<double>> KKTtrips;
  KKTtrips.reserve(nE + 2 * B.nonZeros() + 1);
  for (int i = 0; i < nE; ++i) KKTtrips.emplace_back(i, i, 1.0);
  for (int k = 0; k < B.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
      KKTtrips.emplace_back(it.col(), nE + it.row(), it.value());
      KKTtrips.emplace_back(nE + it.row(), it.col(), it.value());
    }
  }
  KKTtrips.emplace_back(nE, nE, 1e-12);
  Eigen::SparseMatrix<double> KKT(Nsys, Nsys);
  KKT.setFromTriplets(KKTtrips.begin(), KKTtrips.end());

  Eigen::VectorXd b = Eigen::VectorXd::Zero(Nsys);
  b.tail(nV) = rhs;

  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(KKT);
  solver.factorize(KKT);
  if (solver.info() != Eigen::Success) {
    std::cerr << "KKT factorization failed.\n";
    return 7;
  }
  Eigen::VectorXd sol = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    std::cerr << "KKT solve failed.\n";
    return 8;
  }
  Eigen::VectorXd phi = sol.head(nE);

  // Period jumps p on edges (integer), and residual continuous connection.
  const double period = 2.0 * kPi / (double)args.N;
  Eigen::VectorXd p = Eigen::VectorXd::Zero(nE);
  Eigen::VectorXd residual = Eigen::VectorXd::Zero(nE);
  for (int ei = 0; ei < nE; ++ei) {
    const double pe = std::round(phi[ei] / period);
    p[ei] = pe;
    residual[ei] = phi[ei] - pe * period;
  }

  // Recover per-face angles α on dual graph:
  std::vector<double> alpha(mesh.faces.size(), 0.0);
  std::vector<char> vis(mesh.faces.size(), 0);
  std::queue<int> qf;
  int root = -1;
  for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
    if (isInteriorFace(f)) {
      root = f->index;
      break;
    }
  }
  if (root < 0) return 9;
  vis[root] = 1;
  qf.push(root);

  while (!qf.empty()) {
    int fi = qf.front();
    qf.pop();
    FaceCIter f = mesh.faces.begin() + fi;
    HalfEdgeCIter h0 = f->he;
    HalfEdgeCIter h = h0;
    do {
      if (!h->edge->isBoundary()) {
        FaceCIter g = h->flip->face;
        if (isInteriorFace(g) && !vis[g->index]) {
          const int ei = h->edge->index;
          const double sgn = (h == h->edge->he) ? 1.0 : -1.0;
          const double tau = lcTransport(h);
          const double conn = sgn * residual[ei];
          alpha[g->index] = alpha[fi] + tau - conn;
          vis[g->index] = 1;
          qf.push(g->index);
        }
      }
      h = h->next;
    } while (h != h0);
  }

  std::filesystem::create_directories(args.outputDir);
  const std::string edgeCsv = (std::filesystem::path(args.outputDir) / "edge_jumps.csv").string();
  const std::string faceCsv = (std::filesystem::path(args.outputDir) / "face_nrosy_field.csv").string();
  const std::string singCsv = (std::filesystem::path(args.outputDir) / "singularity_used.csv").string();

  {
    std::ofstream out(edgeCsv);
    out << "edge,v0,v1,phi,p,residual\n";
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
      const int ei = e->index;
      out << ei << ","
          << e->he->vertex->index << ","
          << e->he->flip->vertex->index << ","
          << phi[ei] << "," << p[ei] << "," << residual[ei] << "\n";
    }
  }

  {
    std::ofstream out(faceCsv);
    out << "face,alpha_mod,dirx,diry,dirz,cx,cy,cz\n";
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
      if (!isInteriorFace(f)) continue;
      Eigen::Vector3d e1, e2, n;
      buildBasis(f, e1, e2, n);
      double a = std::fmod(alpha[f->index], period);
      if (a < 0) a += period;
      Eigen::Vector3d d = std::cos(a) * e1 + std::sin(a) * e2;
      HalfEdgeCIter h = f->he;
      Eigen::Vector3d c = (h->vertex->position + h->next->vertex->position + h->next->next->vertex->position) / 3.0;
      out << f->index << "," << a << ","
          << d.x() << "," << d.y() << "," << d.z() << ","
          << c.x() << "," << c.y() << "," << c.z() << "\n";
    }
  }

  {
    std::ofstream out(singCsv);
    out << "vertex,q,index,K,target_holonomy\n";
    for (int i = 0; i < nV; ++i) {
      const double idx = (double)q[i] / (double)args.N;
      out << i << "," << q[i] << "," << idx << "," << K[i] << "," << 2.0 * kPi * idx << "\n";
    }
  }

  std::cout << "Done.\n"
            << "  edge jumps: " << edgeCsv << "\n"
            << "  face field: " << faceCsv << "\n"
            << "  singularity: " << singCsv << "\n";
  return 0;
}

