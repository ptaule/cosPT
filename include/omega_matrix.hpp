#ifndef OMEGA_MATRIX_HPP
#define OMEGA_MATRIX_HPP

#include "interpolation.hpp"
#include "utilities.hpp"

/* How many components to evolve? */
#define COMPONENTS 4

/*Definition of the linear (LHS) part of the ODE system for the kernels, using
 * k, eta = ln(D), kappa, zeta(eta) and xi(eta, k). A factor exp(n eta) is
 * factorized out of the kernels, hence
 * d(y)/d(eta) + omega.y + n.y = non-linear part
 * Can in principle provide multiple kappa/zeta/xi values (hence Vec1D
 * structure) */
inline void update_omega_matrix(
    double eta,
    double k,
    const Vec1D<double>& kappa_vec,
    const Vec1D<Interpolation1D>& zeta_vec,
    const Vec1D<Interpolation2D>& xi_vec,
    Vec2D<double>& omega
    )
{
    double kappa = kappa_vec.at(0);
    double zeta = zeta_vec.at(0)(eta);
    double xi = xi_vec.at(0)(eta, k);

    omega[0][1] = -1;
    omega[1][0] = -1.5*zeta * (1 - kappa);
    omega[1][1] = 1.5*zeta - 1;
    omega[1][2] = -1.5*zeta * kappa;
    omega[2][3] = -1;
    omega[3][0] = -1.5*zeta * (1 - kappa);
    omega[3][2] = -1.5*zeta * (kappa - k * k * xi);
    omega[3][3] = 1.5*zeta - 1;

    /*Alternative definition with only two components and zeta(eta) dependence (EdS
     * relaxation test)*/
    /*UNUSED(eta);*/
    /*UNUSED(k);*/
    /*UNUSED(kappa_vec);*/
    /*UNUSED(xi_vec);*/
    /**/
    /*omega[0][1] = -1; */
    /*omega[1][0] = -1.5 * zeta;*/
    /*omega[1][1] = 1.5 * zeta - 1; */
}

#endif /* #ifndef OMEGA_MATRIX_HPP */
