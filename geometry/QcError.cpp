#include "geometry/QcError.h"
#include <Eigen/SVD>

double QuasiConformalError::compute(std::vector<Eigen::Vector3d> p, std::vector<Eigen::Vector3d> q)
{
    // p: 3D vertex positions of triangle (p0, p1, p2)
    // q: 2D UV positions of triangle (q0, q1, q2)
    //
    // Conformal distortion K = sigma1 / sigma2,
    // where sigma1 >= sigma2 are the singular values of the Jacobian J.
    //
    // For a triangle, the mapping is piecewise linear.
    // Jacobian J (2x2) maps from the triangle's local 2D coords to UV space.
    //
    // Let e1 = p1-p0, e2 = p2-p0  (3D edge vectors)
    // Let f1 = q1-q0, f2 = q2-q0        (2D UV edge vectors)
    //
    // Build local 2D coordinate system on the triangle:
    //   local_e1 = (|e1|, 0)
    //   local_e2 = (e2·e1/|e1|, |e2 × e1|/|e1|)
    // Then E = [local_e1; local_e2] is a 2x2 matrix.
    //   F = [f1; f2] is a 2x2 matrix.
    // Jacobian J = F * E^{-1}.

    Eigen::Vector3d e1 = p[1] - p[0];
    Eigen::Vector3d e2 = p[2] - p[0];
    Eigen::Vector2d f1(q[1].x() - q[0].x(), q[1].y() - q[0].y());
    Eigen::Vector2d f2(q[2].x() - q[0].x(), q[2].y() - q[0].y());

    double e1_norm = e1.norm();
    if (e1_norm < 1e-12) return 1.0;

    // E: 2x2 matrix mapping local 2D coords to 3D edge coords
    // column 0 = e1 in local coords = (|e1|, 0)
    // column 1 = e2 in local coords = (e2·e1/|e1|, |e2×e1|/|e1|)
    double Ex = e2.dot(e1) / e1_norm;
    double Ey = (e2.cross(e1)).norm() / e1_norm;  // = |e2|*sin(angle)
    if (Ey < 1e-12) return 1.0;  // degenerate triangle

    // E = [e1_norm,  Ex;  0, Ey]
    // E^{-1} = [1/e1_norm, -Ex/(e1_norm*Ey);  0, 1/Ey]
    double inv_e1 = 1.0 / e1_norm;
    double inv_Ey = 1.0 / Ey;

    // J = F * E^{-1}
    // F = [f1.x(), f2.x();  f1.y(), f2.y()]  (2x2, column-major)
    // J(i,j) = sum_k F(i,k) * E^{-1}(k,j)
    double J00 = f1.x()*inv_e1 + f2.x()*(-Ex*inv_e1*inv_Ey);
    double J01 = f2.x() * inv_Ey;
    double J10 = f1.y()*inv_e1 + f2.y()*(-Ex*inv_e1*inv_Ey);
    double J11 = f2.y() * inv_Ey;

    // Compute singular values of J via SVD of 2x2 matrix
    // For 2x2 J, singular values are sqrt(eigenvalues of J^T J)
    double a = J00*J00 + J10*J10;  // (J^T J)_00
    double b = J00*J01 + J10*J11;  // (J^T J)_01 = (J^T J)_10
    double c = J01*J01 + J11*J11;  // (J^T J)_22

    double disc = sqrt((a - c)*(a - c) + 4.0*b*b);
    double sigma1_sq = 0.5 * ((a + c) + disc);
    double sigma2_sq = 0.5 * ((a + c) - disc);
    if (sigma2_sq < 1e-12) sigma2_sq = 1e-12;
    if (sigma1_sq < sigma2_sq) std::swap(sigma1_sq, sigma2_sq);

    double sigma1 = sqrt(sigma1_sq);
    double sigma2 = sqrt(sigma2_sq);

    // Conformal distortion K = sigma1 / sigma2
    // K = 1.0 means perfectly conformal
    double K = sigma1 / (sigma2 + 1e-12);
    return K;
}

Eigen::Vector3d QuasiConformalError::hsv(double h, double s, double v)
{
    double r = 0, g = 0, b = 0;
    
    if (s == 0) {
        r = v;
        g = v;
        b = v;
    
    } else {
        h = (h == 1 ? 0 : h) * 6;
        
        int i = (int)floor(h);
        
        double f = h - i;
        double p = v * (1 - s);
        double q = v * (1 - (s * f));
        double t = v * (1 - s * (1 - f));
        
        switch (i) {
            case 0:
                r = v;
                g = t;
                b = p;
                break;
                
            case 1:
                r = q;
                g = v;
                b = p;
                break;
                
            case 2:
                r = p;
                g = v;
                b = t;
                break;
                
            case 3:
                r = p;
                g = q;
                b = v;
                break;
                
            case 4:
                r = t;
                g = p;
                b = v;
                break;
                
            case 5:
                r = v;
                g = p;
                b = q;
                break;
                
            default:
                break;
        }
    }
    
    return Eigen::Vector3d(r, g, b);
}

Eigen::Vector3d QuasiConformalError::color(double qc)
{
    // clamp to range [1, 1.5]
    qc = std::max(1.0, std::min(1.5, qc));
    
    // compute color
    return hsv((2.0 - 4.0*(qc-1.0))/3.0, 0.7, 0.65);
}
