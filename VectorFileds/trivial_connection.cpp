#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include "Mesh.h"

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Args {
  std::string inputObj;
  std::string outputDir;
  std::string singularityCsv;
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
    } else if (a == "--auto-singularity") {
      args.autoSingularity = true;
    } else {
      std::cerr << "Unknown option: " << a << "\n";
      return false;
    }
  }
  return true;
}

void printUsage() {
  std::cout
      << "Usage:\n"
      << "  trivial_connection.exe <input.obj> <output_dir> [options]\n\n"
      << "Options:\n"
      << "  --singularity-csv <path>   CSV with lines: vertex,index\n"
      << "  --auto-singularity         Auto place singularities from curvature (debug)\n";
}

bool isInteriorFace(FaceCIter f) { return !f->isBoundary(); }

int countInteriorFaces(const Mesh& mesh) {
  int cnt = 0;
  for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
    if (isInteriorFace(f)) cnt++;
  }
  return cnt;
}

double computeVertexCurvature(VertexCIter v) {
  double sumAngles = 0.0;
  HalfEdgeCIter h = v->he;
  do {
    if (!h->onBoundary) {
      sumAngles += h->next->next->angle();
    }
    h = h->flip->next;
  } while (h != v->he);
  return 2.0 * kPi - sumAngles;
}

void buildFaceBasis(FaceCIter f, Eigen::Vector3d& e1, Eigen::Vector3d& e2, Eigen::Vector3d& n) {
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

double transportNoRotation(HalfEdgeCIter h) {
  Eigen::Vector3d e1i, e2i, ni;
  Eigen::Vector3d e1j, e2j, nj;
  buildFaceBasis(h->face, e1i, e2i, ni);
  buildFaceBasis(h->flip->face, e1j, e2j, nj);

  const Eigen::Vector3d u = h->flip->vertex->position - h->vertex->position;
  const double thetaIJ = std::atan2(u.dot(e2i), u.dot(e1i));
  const double thetaJI = std::atan2(u.dot(e2j), u.dot(e1j));
  return -thetaIJ + thetaJI;
}

bool loadSingularityCsv(const std::string& path, std::vector<int>& s) {
  std::ifstream in(path);
  if (!in.is_open()) return false;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line.find("vertex") != std::string::npos) continue;
    const size_t c = line.find(',');
    if (c == std::string::npos) continue;
    int v = std::stoi(line.substr(0, c));
    int idx = std::stoi(line.substr(c + 1));
    if (v >= 0 && v < static_cast<int>(s.size())) s[v] = idx;
  }
  return true;
}

void autoSingularityFromCurvature(const Mesh& mesh, const Eigen::VectorXd& K, std::vector<int>& s) {
  const int chi = static_cast<int>(mesh.vertices.size()) - static_cast<int>(mesh.edges.size()) + countInteriorFaces(mesh);
  if (chi == 0) return;

  std::vector<std::pair<double, int>> ranked;
  ranked.reserve(mesh.vertices.size());
  for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
    ranked.push_back({std::abs(K[v->index]), v->index});
  }
  std::sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

  const int sign = chi > 0 ? 1 : -1;
  int need = std::abs(chi);
  for (size_t i = 0; i < ranked.size() && need > 0; ++i) {
    s[ranked[i].second] = sign;
    need--;
  }
}

}  // namespace

int main(int argc, char** argv) {
  Args args;
  if (!parseArgs(argc, argv, args)) {
    printUsage();
    return 1;
  }

  Mesh mesh;
  if (!mesh.read(args.inputObj)) {
    std::cerr << "Failed to read mesh: " << args.inputObj << "\n";
    return 2;
  }
  if (!mesh.boundaries.empty()) {
    std::cerr << "Only closed meshes are supported for now.\n";
    return 3;
  }

  const int nV = static_cast<int>(mesh.vertices.size());
  const int nE = static_cast<int>(mesh.edges.size());
  const int nF = countInteriorFaces(mesh);
  const int chi = nV - nE + nF;
  std::cout << "Mesh V/E/F=" << nV << "/" << nE << "/" << nF << " chi=" << chi << "\n";

  Eigen::VectorXd K = Eigen::VectorXd::Zero(nV);
  for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
    K[v->index] = computeVertexCurvature(v);
  }

  std::vector<int> singularity(nV, 0);
  if (!args.singularityCsv.empty()) {
    if (!loadSingularityCsv(args.singularityCsv, singularity)) {
      std::cerr << "Failed to load singularity CSV: " << args.singularityCsv << "\n";
      return 4;
    }
  } else if (args.autoSingularity) {
    autoSingularityFromCurvature(mesh, K, singularity);
  } else {
    std::cerr << "Provide --singularity-csv or --auto-singularity.\n";
    return 5;
  }

  int sumS = 0;
  for (int x : singularity) sumS += x;
  if (sumS != chi) {
    std::cerr << "Gauss-Bonnet violation: sum singularity=" << sumS << " != chi=" << chi << "\n";
    return 6;
  }

  // Build vertex-edge incidence around dual loops: B * x = rhs.
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

  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(nV);
  for (int i = 0; i < nV; ++i) rhs[i] = 2.0 * kPi * singularity[i] - K[i];

  // Solve min ||x||^2 s.t. Bx = rhs using KKT.
  const int N = nE + nV;
  std::vector<Eigen::Triplet<double>> KKTtrips;
  KKTtrips.reserve(nE + 2 * B.nonZeros() + nV);
  for (int i = 0; i < nE; ++i) KKTtrips.emplace_back(i, i, 1.0);
  for (int k = 0; k < B.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
      const int r = it.row();
      const int c = it.col();
      const double v = it.value();
      KKTtrips.emplace_back(c, nE + r, v);
      KKTtrips.emplace_back(nE + r, c, v);
    }
  }
  KKTtrips.emplace_back(nE, nE, 1e-12);  // tiny gauge regularizer

  Eigen::SparseMatrix<double> KKT(N, N);
  KKT.setFromTriplets(KKTtrips.begin(), KKTtrips.end());
  Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
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
  Eigen::VectorXd x = sol.head(nE);  // connection correction per edge

  // Propagate face angles.
  std::vector<double> alpha(mesh.faces.size(), 0.0);
  std::vector<char> visited(mesh.faces.size(), 0);
  std::queue<int> q;
  int root = -1;
  for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
    if (isInteriorFace(f)) {
      root = f->index;
      break;
    }
  }
  if (root < 0) {
    std::cerr << "No interior faces.\n";
    return 9;
  }
  visited[root] = 1;
  q.push(root);

  while (!q.empty()) {
    const int fi = q.front();
    q.pop();
    FaceCIter f = mesh.faces.begin() + fi;
    HalfEdgeCIter h0 = f->he;
    HalfEdgeCIter h = h0;
    do {
      if (!h->edge->isBoundary()) {
        FaceCIter g = h->flip->face;
        if (isInteriorFace(g) && !visited[g->index]) {
          const double tau = transportNoRotation(h);
          const double phi = (h == h->edge->he ? 1.0 : -1.0) * x[h->edge->index];
          alpha[g->index] = alpha[fi] + tau - phi;
          visited[g->index] = 1;
          q.push(g->index);
        }
      }
      h = h->next;
    } while (h != h0);
  }

  std::filesystem::create_directories(args.outputDir);
  const std::string connCsv = (std::filesystem::path(args.outputDir) / "connection.csv").string();
  const std::string faceCsv = (std::filesystem::path(args.outputDir) / "face_field.csv").string();
  const std::string singCsv = (std::filesystem::path(args.outputDir) / "singularity_used.csv").string();

  {
    std::ofstream out(connCsv);
    out << "edge,v0,v1,phi\n";
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
      out << e->index << ","
          << e->he->vertex->index << ","
          << e->he->flip->vertex->index << ","
          << x[e->index] << "\n";
    }
  }

  {
    std::ofstream out(faceCsv);
    out << "face,alpha,vx,vy,vz,cx,cy,cz\n";
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
      if (!isInteriorFace(f)) continue;
      Eigen::Vector3d e1, e2, n;
      buildFaceBasis(f, e1, e2, n);
      const double a = alpha[f->index];
      Eigen::Vector3d v = std::cos(a) * e1 + std::sin(a) * e2;

      HalfEdgeCIter h = f->he;
      Eigen::Vector3d c = (h->vertex->position + h->next->vertex->position + h->next->next->vertex->position) / 3.0;
      out << f->index << "," << a << ","
          << v.x() << "," << v.y() << "," << v.z() << ","
          << c.x() << "," << c.y() << "," << c.z() << "\n";
    }
  }

  {
    std::ofstream out(singCsv);
    out << "vertex,index,K,target_curvature\n";
    for (int i = 0; i < nV; ++i) {
      out << i << "," << singularity[i] << "," << K[i] << "," << (2.0 * kPi * singularity[i]) << "\n";
    }
  }

  std::cout << "Done.\n"
            << "  connection: " << connCsv << "\n"
            << "  face field: " << faceCsv << "\n"
            << "  singularity: " << singCsv << "\n";
  return 0;
}

