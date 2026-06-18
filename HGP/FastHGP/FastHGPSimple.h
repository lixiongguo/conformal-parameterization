#ifndef FAST_HGP_SIMPLE_H
#define FAST_HGP_SIMPLE_H

#include "../Parameterization/Parameterization.h"
#include <Eigen/Dense>
#include <vector>

// WASM subset: BaseMesh + Eigen only (no CGAL / GMM / MATLAB).
// LSCM initial UV, then symmetric Dirichlet minimization with injectivity line search.
// Full desktop FastHGP (KKT / ATP / Newton) is in FastHGP.h and requires a CGAL build.
class FastHGPSimple : public Parameterization {
public:
    explicit FastHGPSimple(Mesh& mesh0, int maxIterations = 80);

    void parameterize() override;

private:
    struct TriData {
        int v[3];
        Eigen::Matrix2d invLocal;
        double area;
    };

    int maxIterations;
    std::vector<TriData> triangles;
    std::vector<int> variableIndex;
    std::vector<int> freeVertices;

    void initializeWithLscm();
    void buildTriangleData();
    void buildFreeVariables();

    Eigen::VectorXd packFreeUvs() const;
    void unpackFreeUvs(const Eigen::VectorXd& x);
    std::vector<Eigen::Vector2d> currentUvs() const;
    std::vector<Eigen::Vector2d> uvsWithStep(const Eigen::VectorXd& direction, double t) const;

    double symmetricDirichletEnergy(const std::vector<Eigen::Vector2d>& uv) const;
    double computeEnergyAndGradient(Eigen::VectorXd& gradient) const;
    double localInjectivityStep(const std::vector<Eigen::Vector2d>& uv,
                                const Eigen::VectorXd& direction) const;
    double decreasingEnergyStep(double energy,
                                const Eigen::VectorXd& gradient,
                                const Eigen::VectorXd& direction,
                                double tMax) const;
};

#endif
