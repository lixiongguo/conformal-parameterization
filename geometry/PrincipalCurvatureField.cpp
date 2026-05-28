#include "PrincipalCurvatureField.h"
#include <Eigen/Eigenvalues>
#include <cmath>

using Eigen::Matrix2d;
using Eigen::Vector2d;
using Eigen::Vector3d;

namespace {

void buildLocalFrame(const Vector3d& n, Vector3d& t1, Vector3d& t2) {
    Vector3d ref(1, 0, 0);
    if (std::fabs(n.dot(ref)) > 0.9) ref = Vector3d(0, 1, 0);
    t1 = (ref - ref.dot(n) * n).normalized();
    t2 = n.cross(t1).normalized();
}

} // namespace

namespace PrincipalCurvatureField {

void estimateFromWeingarten(
    const Eigen::MatrixXd& vertPos,
    const Eigen::VectorXi& faces,
    const Eigen::MatrixXd& vertexNormals,
    const Eigen::MatrixXd& faceNormals,
    Eigen::MatrixXd& faceDirs) {

    const int nF = (int)faceNormals.rows();
    faceDirs.resize(nF, 3);

    for (int fi = 0; fi < nF; ++fi) {
        const int v0 = faces[fi * 3];
        const int v1 = faces[fi * 3 + 1];
        const int v2 = faces[fi * 3 + 2];
        const Vector3d p0 = vertPos.row(v0);
        const Vector3d p1 = vertPos.row(v1);
        const Vector3d p2 = vertPos.row(v2);

        const Vector3d fn = faceNormals.row(fi).normalized();

        Vector3d t1, t2;
        buildLocalFrame(fn, t1, t2);

        const Vector3d n0 = vertexNormals.row(v0);
        const Vector3d n1 = vertexNormals.row(v1);
        const Vector3d n2 = vertexNormals.row(v2);

        const Vector3d e01 = p1 - p0;
        const Vector3d e02 = p2 - p0;
        const double u1 = e01.dot(t1), v1c = e01.dot(t2);
        const double u2 = e02.dot(t1), v2c = e02.dot(t2);

        Matrix2d I_mat;
        I_mat << u1 * u1 + v1c * v1c, u1 * u2 + v1c * v2c,
                 u1 * u2 + v1c * v2c, u2 * u2 + v2c * v2c;

        const Vector3d dn1 = n1 - n0;
        const Vector3d dn2 = n2 - n0;
        const Vector2d rhs1(-dn1.dot(t1), -dn2.dot(t1));
        const Vector2d rhs2(-dn1.dot(t2), -dn2.dot(t2));

        Matrix2d M;
        M << u1, v1c, u2, v2c;
        const Matrix2d MtM = M.transpose() * M;

        const Vector2d ef = MtM.ldlt().solve(M.transpose() * rhs1);
        const Vector2d fg = MtM.ldlt().solve(M.transpose() * rhs2);
        const double f_sym = 0.5 * (ef(1) + fg(0));

        Matrix2d II_mat;
        II_mat << ef(0), f_sym, f_sym, fg(1);

        const Matrix2d W = I_mat.inverse() * II_mat;
        const Matrix2d W_sym = 0.5 * (W + W.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix2d> es(W_sym);
        Vector2d dir2d = (std::fabs(es.eigenvalues()(0)) >= std::fabs(es.eigenvalues()(1)))
                             ? es.eigenvectors().col(0)
                             : es.eigenvectors().col(1);

        Vector3d dir3d = dir2d(0) * t1 + dir2d(1) * t2;
        const double dlen = dir3d.norm();
        if (dlen > 1e-12) dir3d /= dlen;
        else dir3d = t1;

        faceDirs.row(fi) = dir3d;
    }
}

void estimateFromWeingartenWithCurvatures(
    const Eigen::MatrixXd& vertPos,
    const Eigen::VectorXi& faces,
    const Eigen::MatrixXd& vertexNormals,
    const Eigen::MatrixXd& faceNormals,
    Eigen::MatrixXd& faceDirs,
    Eigen::MatrixXd& k1k2) {

    const int nF = (int)faceNormals.rows();
    faceDirs.resize(nF, 3);
    k1k2.resize(nF, 2);

    for (int fi = 0; fi < nF; ++fi) {
        const int v0 = faces[fi * 3];
        const int v1 = faces[fi * 3 + 1];
        const int v2 = faces[fi * 3 + 2];
        const Vector3d p0 = vertPos.row(v0);
        const Vector3d p1 = vertPos.row(v1);
        const Vector3d p2 = vertPos.row(v2);

        const Vector3d fn = faceNormals.row(fi).normalized();

        Vector3d t1, t2;
        buildLocalFrame(fn, t1, t2);

        const Vector3d n0 = vertexNormals.row(v0);
        const Vector3d n1 = vertexNormals.row(v1);
        const Vector3d n2 = vertexNormals.row(v2);

        const Vector3d e01 = p1 - p0;
        const Vector3d e02 = p2 - p0;
        const double u1 = e01.dot(t1), v1c = e01.dot(t2);
        const double u2 = e02.dot(t1), v2c = e02.dot(t2);

        Matrix2d I_mat;
        I_mat << u1 * u1 + v1c * v1c, u1 * u2 + v1c * v2c,
                 u1 * u2 + v1c * v2c, u2 * u2 + v2c * v2c;

        const Vector3d dn1 = n1 - n0;
        const Vector3d dn2 = n2 - n0;
        const Vector2d rhs1(-dn1.dot(t1), -dn2.dot(t1));
        const Vector2d rhs2(-dn1.dot(t2), -dn2.dot(t2));

        Matrix2d M;
        M << u1, v1c, u2, v2c;
        const Matrix2d MtM = M.transpose() * M;

        const Vector2d ef = MtM.ldlt().solve(M.transpose() * rhs1);
        const Vector2d fg = MtM.ldlt().solve(M.transpose() * rhs2);
        const double f_sym = 0.5 * (ef(1) + fg(0));

        Matrix2d II_mat;
        II_mat << ef(0), f_sym, f_sym, fg(1);

        const Matrix2d W = I_mat.inverse() * II_mat;
        const Matrix2d W_sym = 0.5 * (W + W.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix2d> es(W_sym);
        const double k0 = es.eigenvalues()(0);
        const double k1 = es.eigenvalues()(1);
        k1k2(fi, 0) = k0;
        k1k2(fi, 1) = k1;

        Vector2d dir2d = (std::fabs(k0) >= std::fabs(k1))
                             ? es.eigenvectors().col(0)
                             : es.eigenvectors().col(1);

        Vector3d dir3d = dir2d(0) * t1 + dir2d(1) * t2;
        const double dlen = dir3d.norm();
        if (dlen > 1e-12) dir3d /= dlen;
        else dir3d = t1;

        faceDirs.row(fi) = dir3d;
    }
}

} // namespace PrincipalCurvatureField

