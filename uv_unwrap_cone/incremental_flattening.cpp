#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "IncrementalFlatteningWasm.h"
#include "Mesh.h"

namespace {

constexpr double kPi = 3.14159265358979323846;

struct FaceMetricData {
    int v0;
    int v1;
    int v2;
    int e01;
    int e12;
    int e20;
};

struct Params {
    double eps0 = 0.02;
    double epsStep = 0.01;
    double epsMax = 1.2;
    int maxIters = 80;
    int maxCones = 12;
};

double clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

bool computeTriangleAngles(double l01, double l12, double l20,
                           double& a0, double& a1, double& a2) {
    const double eps = 1e-12;
    if (l01 <= eps || l12 <= eps || l20 <= eps) return false;
    if (l01 + l12 <= l20 + eps || l12 + l20 <= l01 + eps || l20 + l01 <= l12 + eps) return false;

    const double c0 = clamp((l01 * l01 + l20 * l20 - l12 * l12) / (2.0 * l01 * l20), -1.0, 1.0);
    const double c1 = clamp((l01 * l01 + l12 * l12 - l20 * l20) / (2.0 * l01 * l12), -1.0, 1.0);
    const double c2 = clamp((l12 * l12 + l20 * l20 - l01 * l01) / (2.0 * l12 * l20), -1.0, 1.0);
    a0 = std::acos(c0);
    a1 = std::acos(c1);
    a2 = std::acos(c2);
    return std::isfinite(a0) && std::isfinite(a1) && std::isfinite(a2);
}

double roundQuarterPi(double x) {
    return std::round(2.0 * x / kPi) * (kPi * 0.5);
}

void buildFaceMetricData(const Mesh& mesh, std::vector<FaceMetricData>& outFaces) {
    outFaces.clear();
    outFaces.reserve(mesh.faces.size());
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        HalfEdgeCIter h = f->he;
        FaceMetricData fm;
        fm.v0 = h->vertex->index;
        fm.v1 = h->next->vertex->index;
        fm.v2 = h->next->next->vertex->index;
        fm.e01 = h->edge->index;
        fm.e12 = h->next->edge->index;
        fm.e20 = h->next->next->edge->index;
        outFaces.push_back(fm);
    }
}

void buildVertexAreas(const Mesh& mesh, Eigen::VectorXd& areas) {
    const int n = static_cast<int>(mesh.vertices.size());
    areas = Eigen::VectorXd::Constant(n, 1e-6);
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
        if (f->isBoundary()) continue;
        const double a = f->area();
        if (a <= 0.0) continue;
        HalfEdgeCIter h = f->he;
        areas[h->vertex->index] += a / 3.0;
        areas[h->next->vertex->index] += a / 3.0;
        areas[h->next->next->vertex->index] += a / 3.0;
    }
}

void buildBoundaryFlags(const Mesh& mesh, std::vector<bool>& isBoundary) {
    const int n = static_cast<int>(mesh.vertices.size());
    isBoundary.assign(n, false);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
        isBoundary[v->index] = v->isBoundary();
    }
}

Eigen::SparseMatrix<double> buildCotanLaplacian(const Mesh& mesh) {
    const int n = static_cast<int>(mesh.vertices.size());
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.edges.size() * 4 + n);

    std::vector<double> diag(n, 0.0);
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
        const int i = e->he->vertex->index;
        const int j = e->he->flip->vertex->index;
        double w = e->cotanWeigth();
        if (!std::isfinite(w)) w = 0.0;
        diag[i] += w;
        diag[j] += w;
        triplets.emplace_back(i, j, -w);
        triplets.emplace_back(j, i, -w);
    }
    for (int i = 0; i < n; ++i) {
        triplets.emplace_back(i, i, diag[i] + 1e-10);
    }

    Eigen::SparseMatrix<double> L(n, n);
    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}

void edgeLengthsFromScale(const Mesh& mesh, const Eigen::VectorXd& u,
                          const std::vector<double>& baseLen, std::vector<double>& edgeLen) {
    edgeLen.resize(mesh.edges.size());
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
        const int i = e->he->vertex->index;
        const int j = e->he->flip->vertex->index;
        edgeLen[e->index] = baseLen[e->index] * std::exp(0.5 * (u[i] + u[j]));
    }
}

void computeMetricCurvature(const std::vector<FaceMetricData>& faces,
                            const std::vector<double>& edgeLen,
                            const std::vector<bool>& isBoundary,
                            Eigen::VectorXd& K) {
    const int n = static_cast<int>(isBoundary.size());
    Eigen::VectorXd angleSums = Eigen::VectorXd::Zero(n);

    for (const FaceMetricData& f : faces) {
        const double l01 = edgeLen[f.e01];
        const double l12 = edgeLen[f.e12];
        const double l20 = edgeLen[f.e20];
        double a0 = 0.0, a1 = 0.0, a2 = 0.0;
        if (!computeTriangleAngles(l01, l12, l20, a0, a1, a2)) {
            continue;
        }
        angleSums[f.v0] += a0;
        angleSums[f.v1] += a1;
        angleSums[f.v2] += a2;
    }

    K = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < n; ++i) {
        const double tgt = isBoundary[i] ? kPi : (2.0 * kPi);
        K[i] = tgt - angleSums[i];
    }
}

bool solveConstrainedMinNorm(
    const Eigen::VectorXd& vertexAreas,
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& Lrow,
    const Eigen::VectorXd& Korig,
    const std::vector<std::pair<int, double>>& constraints,
    Eigen::VectorXd& uOut) {
    const int n = static_cast<int>(vertexAreas.size());
    const int m = static_cast<int>(constraints.size());
    if (m == 0) {
        uOut = Eigen::VectorXd::Zero(n);
        return true;
    }

    const int N = n + m;
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(n + 2 * (Lrow.nonZeros() + m));

    for (int i = 0; i < n; ++i) {
        trips.emplace_back(i, i, std::max(vertexAreas[i], 1e-8));
    }

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(N);
    for (int r = 0; r < m; ++r) {
        const int rowIdx = constraints[r].first;
        const double targetK = constraints[r].second;
        rhs[n + r] = targetK - Korig[rowIdx];

        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Lrow, rowIdx); it; ++it) {
            const int c = it.col();
            const double v = it.value();
            if (std::abs(v) < 1e-14) continue;
            trips.emplace_back(c, n + r, v);
            trips.emplace_back(n + r, c, v);
        }
    }

    Eigen::SparseMatrix<double> KKT(N, N);
    KKT.setFromTriplets(trips.begin(), trips.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(KKT);
    solver.factorize(KKT);
    if (solver.info() != Eigen::Success) return false;
    Eigen::VectorXd x = solver.solve(rhs);
    if (solver.info() != Eigen::Success) return false;
    uOut = x.head(n);
    return true;
}

void writeConeCsv(const std::string& path,
                  const std::vector<bool>& isBoundary,
                  const std::vector<bool>& flattened,
                  const std::set<int>& candidates,
                  const std::map<int, double>& roundedTargets,
                  const Eigen::VectorXd& Korig,
                  const Eigen::VectorXd& Kfinal) {
    std::ofstream out(path);
    out << "vertex,is_boundary,flattened,candidate,K_orig,K_final,K_rounded\n";
    for (int i = 0; i < static_cast<int>(isBoundary.size()); ++i) {
        out << i << ","
            << (isBoundary[i] ? 1 : 0) << ","
            << (flattened[i] ? 1 : 0) << ","
            << (candidates.count(i) ? 1 : 0) << ","
            << std::setprecision(17) << Korig[i] << ","
            << std::setprecision(17) << Kfinal[i] << ",";
        auto it = roundedTargets.find(i);
        if (it != roundedTargets.end()) out << std::setprecision(17) << it->second;
        out << "\n";
    }
}

void writeEdgeCsv(const std::string& path, const Mesh& mesh,
                  const std::vector<double>& baseLen,
                  const std::vector<double>& finalLen) {
    std::ofstream out(path);
    out << "edge,v0,v1,l0,lfinal,scale\n";
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
        const int v0 = e->he->vertex->index;
        const int v1 = e->he->flip->vertex->index;
        const double l0 = baseLen[e->index];
        const double lf = finalLen[e->index];
        out << e->index << "," << v0 << "," << v1 << ","
            << std::setprecision(17) << l0 << ","
            << std::setprecision(17) << lf << ","
            << std::setprecision(17) << (l0 > 0.0 ? lf / l0 : 1.0) << "\n";
    }
}

void writeSummary(const std::string& path, const Params& p, int nVertices, int nEdges, int nFaces,
                  int flattenIters, int nCandidates, int nRounded) {
    std::ofstream out(path);
    out << "Incremental Flattening Reproduction (core stages)\n";
    out << "vertices=" << nVertices << "\n";
    out << "edges=" << nEdges << "\n";
    out << "faces=" << nFaces << "\n";
    out << "flatten_iterations=" << flattenIters << "\n";
    out << "cone_candidates=" << nCandidates << "\n";
    out << "rounded_cones=" << nRounded << "\n";
    out << "eps0=" << p.eps0 << "\n";
    out << "eps_step=" << p.epsStep << "\n";
    out << "eps_max=" << p.epsMax << "\n";
    out << "max_iters=" << p.maxIters << "\n";
    out << "max_cones=" << p.maxCones << "\n";
    out << "\n";
    out << "Notes:\n";
    out << "- Implements paper core stages: flattening + cone rounding.\n";
    out << "- Outputs intrinsic metric and rounded cone curvatures.\n";
    out << "- Holonomy rounding on homology loops and full ARAP seamless solve\n";
    out << "  are intentionally left as next-step extensions.\n";
}

void printUsage() {
    std::cout << "Usage:\n"
              << "  incremental_flattening.exe <input.obj> <output_dir> [options]\n\n"
              << "Options:\n"
              << "  --eps0 <value>        Initial flatten threshold (default 0.02)\n"
              << "  --eps-step <value>    Flatten threshold increment (default 0.01)\n"
              << "  --eps-max <value>     Max flatten threshold (default 1.2)\n"
              << "  --max-iters <int>     Max flatten iterations (default 80)\n"
              << "  --max-cones <int>     Stop flatten when remaining interior vertices <= this (default 12)\n";
}

bool parseArgs(int argc, char** argv, std::string& inputObj, std::string& outDir, Params& p) {
    if (argc < 3) return false;
    inputObj = argv[1];
    outDir = argv[2];

    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        auto needValue = [&](const std::string& name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                return nullptr;
            }
            return argv[++i];
        };
        if (arg == "--eps0") {
            const char* v = needValue(arg);
            if (!v) return false;
            p.eps0 = std::atof(v);
        } else if (arg == "--eps-step") {
            const char* v = needValue(arg);
            if (!v) return false;
            p.epsStep = std::atof(v);
        } else if (arg == "--eps-max") {
            const char* v = needValue(arg);
            if (!v) return false;
            p.epsMax = std::atof(v);
        } else if (arg == "--max-iters") {
            const char* v = needValue(arg);
            if (!v) return false;
            p.maxIters = std::atoi(v);
        } else if (arg == "--max-cones") {
            const char* v = needValue(arg);
            if (!v) return false;
            p.maxCones = std::atoi(v);
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::string inputObj;
    std::string outDir;
    Params params;
    if (!parseArgs(argc, argv, inputObj, outDir, params)) {
        printUsage();
        return 1;
    }

    Mesh mesh;
    if (!mesh.read(inputObj)) {
        std::cerr << "Failed to read mesh: " << inputObj << "\n";
        return 2;
    }

    IncrementalOptions options;
    options.eps0 = params.eps0;
    options.epsStep = params.epsStep;
    options.epsMax = params.epsMax;
    options.maxIters = params.maxIters;
    options.maxCones = params.maxCones;

    IncrementalState state;
    if (!runIncrementalFlattening(mesh, options, state)) {
        std::cerr << "Incremental flattening failed: " << state.error << "\n";
        return 3;
    }

    const int nV = state.nVertices;
    const int nE = state.nEdges;
    const int nF = state.nFaces;
    std::cout << "Loaded mesh: V=" << nV << " E=" << nE << " F=" << nF << "\n";

    std::filesystem::create_directories(outDir);
    const std::string coneCsv = (std::filesystem::path(outDir) / "cones.csv").string();
    const std::string edgeCsv = (std::filesystem::path(outDir) / "edge_metric.csv").string();
    const std::string summaryTxt = (std::filesystem::path(outDir) / "summary.txt").string();
    const std::string scaleTxt = (std::filesystem::path(outDir) / "vertex_log_scale.csv").string();

    writeConeCsv(coneCsv, state.isBoundary, state.flattened, state.candidates, state.roundedTargets, state.Korig, state.Kfinal);
    writeEdgeCsv(edgeCsv, mesh, state.baseLen, state.edgeLen);
    writeSummary(summaryTxt, params, nV, nE, nF, state.flattenIterations,
                 static_cast<int>(state.candidates.size()), static_cast<int>(state.roundedTargets.size()));

    {
        std::ofstream out(scaleTxt);
        out << "vertex,u\n";
        for (int i = 0; i < nV; ++i) {
            out << i << "," << std::setprecision(17) << state.u[i] << "\n";
        }
    }

    std::cout << "\nDone.\n"
              << "  cones: " << coneCsv << "\n"
              << "  edge metric: " << edgeCsv << "\n"
              << "  vertex scale: " << scaleTxt << "\n"
              << "  summary: " << summaryTxt << "\n";
    return 0;
}

