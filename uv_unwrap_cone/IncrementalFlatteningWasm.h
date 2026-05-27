#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

#include <Eigen/Core>

class Mesh;

struct IncrementalOptions {
    double eps0 = 0.02;
    double epsStep = 0.01;
    double epsMax = 1.2;
    int maxIters = 80;
    int maxCones = 12;
};

struct IncrementalState {
    int nVertices = 0;
    int nEdges = 0;
    int nFaces = 0;
    int flattenIterations = 0;

    std::vector<bool> isBoundary;
    std::vector<bool> flattened;
    std::set<int> candidates;
    std::map<int, double> roundedTargets;

    std::vector<double> baseLen;
    std::vector<double> edgeLen;
    Eigen::VectorXd Korig;
    Eigen::VectorXd Kfinal;
    Eigen::VectorXd u;
    std::vector<double> uvFlat;

    std::string error;
};

bool runIncrementalFlattening(Mesh& mesh, const IncrementalOptions& options, IncrementalState& outState);
