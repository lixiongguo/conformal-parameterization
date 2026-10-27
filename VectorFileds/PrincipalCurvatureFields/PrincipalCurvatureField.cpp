#include "PrincipalCurvatureField.h"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <vector>

using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;

namespace {

void buildLocalFrame(const Vector3d& n, Vector3d& t1, Vector3d& t2) {
    Vector3d ref(1, 0, 0);
    if (std::fabs(n.dot(ref)) > 0.9) ref = Vector3d(0, 1, 0);
    t1 = (ref - ref.dot(n) * n).normalized();
    t2 = n.cross(t1).normalized();
}

// Shape operator in an orthonormal tangent basis: S(dx) = -dn.
// Three edge constraints, symmetric 2x2, so the direction does not depend on
// which triangle corner is listed first.
bool shapeOperator(const Vector3d point[3], const Vector3d normal[3], const Vector3d& t1, const Vector3d& t2,
                   Matrix2d& shape, double& curvatureMax, double& curvatureMin, Vector2d& principal) {
    Matrix3d system = Matrix3d::Zero();
    Vector3d rhs = Vector3d::Zero();
    for (int edge = 0; edge < 3; ++edge) {
        const Vector3d dp = point[(edge + 1) % 3] - point[edge];
        const Vector3d dn = normal[(edge + 1) % 3] - normal[edge];
        const double du = dp.dot(t1);
        const double dv = dp.dot(t2);
        const double ru = -dn.dot(t1);
        const double rv = -dn.dot(t2);
        system(0, 0) += du * du;
        system(0, 1) += du * dv;
        rhs(0) += du * ru;
        system(1, 0) += dv * du;
        system(1, 1) += dv * dv;
        rhs(1) += dv * ru;
        system(1, 1) += du * du;
        system(1, 2) += du * dv;
        rhs(1) += du * rv;
        system(2, 1) += dv * du;
        system(2, 2) += dv * dv;
        rhs(2) += dv * rv;
    }
    const Eigen::LDLT<Matrix3d> solver(system);
    if (solver.info() != Eigen::Success) return false;
    const Vector3d coefficient = solver.solve(rhs);
    if (!coefficient.allFinite()) return false;
    shape << coefficient(0), coefficient(1), coefficient(1), coefficient(2);
    Eigen::SelfAdjointEigenSolver<Matrix2d> eigen(shape);
    if (eigen.info() != Eigen::Success) return false;
    const int dominant = std::fabs(eigen.eigenvalues()(0)) >= std::fabs(eigen.eigenvalues()(1)) ? 0 : 1;
    const int other = 1 - dominant;
    curvatureMax = eigen.eigenvalues()(dominant);
    curvatureMin = eigen.eigenvalues()(other);
    principal = eigen.eigenvectors().col(dominant);
    return principal.allFinite();
}

} // namespace

namespace PrincipalCurvatureField {

void estimateFromWeingartenWithCurvatures(
    const Eigen::MatrixXd& vertPos,
    const Eigen::VectorXi& faces,
    const Eigen::MatrixXd& vertexNormals,
    const Eigen::MatrixXd& faceNormals,
    Eigen::MatrixXd& faceDirs,
    Eigen::MatrixXd& k1k2) {

    const int nF = static_cast<int>(faceNormals.rows());
    const int nV = static_cast<int>(vertPos.rows());
    faceDirs.resize(nF, 3);
    k1k2.resize(nF, 2);

    struct FaceFrame {
        Matrix2d shape = Matrix2d::Zero();
        Vector3d t1 = Vector3d::UnitX();
        Vector3d t2 = Vector3d::UnitY();
        double area = 0;
        bool ok = false;
    };
    std::vector<FaceFrame> frame(static_cast<size_t>(nF));
    std::vector<Matrix2d> vertexShape(static_cast<size_t>(nV), Matrix2d::Zero());
    std::vector<double> vertexWeight(static_cast<size_t>(nV), 0);
    std::vector<Vector3d> vertexT1(static_cast<size_t>(nV), Vector3d::UnitX());
    std::vector<Vector3d> vertexT2(static_cast<size_t>(nV), Vector3d::UnitY());
    for (int vertex = 0; vertex < nV; ++vertex) {
        Vector3d normal = vertexNormals.row(vertex);
        if (normal.norm() < 1e-12) continue;
        buildLocalFrame(normal.normalized(), vertexT1[static_cast<size_t>(vertex)], vertexT2[static_cast<size_t>(vertex)]);
    }

    for (int face = 0; face < nF; ++face) {
        const int corner[3] = {faces[face * 3], faces[face * 3 + 1], faces[face * 3 + 2]};
        const Vector3d point[3] = {vertPos.row(corner[0]), vertPos.row(corner[1]), vertPos.row(corner[2])};
        const Vector3d normal[3] = {vertexNormals.row(corner[0]), vertexNormals.row(corner[1]), vertexNormals.row(corner[2])};
        Vector3d faceNormal = faceNormals.row(face);
        const double area = 0.5 * (point[1] - point[0]).cross(point[2] - point[0]).norm();
        if (faceNormal.norm() < 1e-12 || area < 1e-14) continue;
        faceNormal.normalize();
        FaceFrame& local = frame[static_cast<size_t>(face)];
        buildLocalFrame(faceNormal, local.t1, local.t2);
        local.area = area;
        double curvatureMax = 0;
        double curvatureMin = 0;
        Vector2d principal;
        local.ok = shapeOperator(point, normal, local.t1, local.t2, local.shape, curvatureMax, curvatureMin, principal);
        if (!local.ok) continue;
        for (int cornerId : corner) {
            const Vector3d& s1 = vertexT1[static_cast<size_t>(cornerId)];
            const Vector3d& s2 = vertexT2[static_cast<size_t>(cornerId)];
            Matrix2d rotation;
            rotation << local.t1.dot(s1), local.t2.dot(s1), local.t1.dot(s2), local.t2.dot(s2);
            vertexShape[static_cast<size_t>(cornerId)] += area * (rotation * local.shape * rotation.transpose());
            vertexWeight[static_cast<size_t>(cornerId)] += area;
        }
    }

    for (int face = 0; face < nF; ++face) {
        FaceFrame& local = frame[static_cast<size_t>(face)];
        Vector3d direction = local.t1;
        double curvatureMax = 0;
        double curvatureMin = 0;
        Matrix2d accumulated = Matrix2d::Zero();
        double weight = 0;
        if (local.ok) {
            for (int k = 0; k < 3; ++k) {
                const int vertex = faces[face * 3 + k];
                if (vertexWeight[static_cast<size_t>(vertex)] < 1e-14) continue;
                const Matrix2d vertexOperator = vertexShape[static_cast<size_t>(vertex)] / vertexWeight[static_cast<size_t>(vertex)];
                const Vector3d& s1 = vertexT1[static_cast<size_t>(vertex)];
                const Vector3d& s2 = vertexT2[static_cast<size_t>(vertex)];
                Matrix2d rotation;
                rotation << s1.dot(local.t1), s2.dot(local.t1), s1.dot(local.t2), s2.dot(local.t2);
                accumulated += rotation * vertexOperator * rotation.transpose();
                weight += 1;
            }
        }
        if (weight > 0) {
            const Matrix2d shape = 0.5 * (accumulated + accumulated.transpose()) / weight;
            Eigen::SelfAdjointEigenSolver<Matrix2d> eigen(shape);
            if (eigen.info() == Eigen::Success) {
                const int dominant = std::fabs(eigen.eigenvalues()(0)) >= std::fabs(eigen.eigenvalues()(1)) ? 0 : 1;
                curvatureMax = eigen.eigenvalues()(dominant);
                curvatureMin = eigen.eigenvalues()(1 - dominant);
                const Vector2d principal = eigen.eigenvectors().col(dominant);
                direction = principal(0) * local.t1 + principal(1) * local.t2;
                if (direction.norm() > 1e-12) direction.normalize();
                else direction = local.t1;
            }
        }
        faceDirs.row(face) = direction;
        k1k2(face, 0) = curvatureMax;
        k1k2(face, 1) = curvatureMin;
    }
}

void estimateFromWeingarten(
    const Eigen::MatrixXd& vertPos,
    const Eigen::VectorXi& faces,
    const Eigen::MatrixXd& vertexNormals,
    const Eigen::MatrixXd& faceNormals,
    Eigen::MatrixXd& faceDirs) {

    Eigen::MatrixXd curvature;
    estimateFromWeingartenWithCurvatures(vertPos, faces, vertexNormals, faceNormals, faceDirs, curvature);
}

} // namespace PrincipalCurvatureField
