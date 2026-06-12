#ifndef PRINCIPAL_CURVATURE_FIELD_H
#define PRINCIPAL_CURVATURE_FIELD_H

#include <Eigen/Core>

namespace PrincipalCurvatureField {

/**
 * Estimate one principal curvature direction per face on a triangle mesh.
 *
 * Inputs:
 *  - vertPos: nV x 3 vertex positions
 *  - faces:   (nF*3) triangle indices (flattened)
 *  - vertexNormals: nV x 3 precomputed vertex normals (area-weighted suggested)
 *  - faceNormals:   nF x 3 precomputed face normals (unit)
 *
 * Output:
 *  - faceDirs: nF x 3 unit tangent direction per face (principal direction)
 */
void estimateFromWeingarten(
    const Eigen::MatrixXd& vertPos,
    const Eigen::VectorXi& faces,
    const Eigen::MatrixXd& vertexNormals,
    const Eigen::MatrixXd& faceNormals,
    Eigen::MatrixXd& faceDirs);

/**
 * Same as estimateFromWeingarten(), but also returns per-face principal curvatures (k1,k2).
 *
 * Output:
 *  - faceDirs: nF x 3 unit tangent direction per face (dominant principal direction by |k|)
 *  - k1k2:     nF x 2 principal curvature values (eigenvalues of symmetric Weingarten map)
 */
void estimateFromWeingartenWithCurvatures(
    const Eigen::MatrixXd& vertPos,
    const Eigen::VectorXi& faces,
    const Eigen::MatrixXd& vertexNormals,
    const Eigen::MatrixXd& faceNormals,
    Eigen::MatrixXd& faceDirs,
    Eigen::MatrixXd& k1k2);

} // namespace PrincipalCurvatureField

#endif

