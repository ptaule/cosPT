#ifndef OMEGA_MATRIX_HPP
#define OMEGA_MATRIX_HPP

/*Definition of the linear (LHS) part of the ODE system for the kernels, using
 * kappa, zeta(eta) and xi(eta, k), where eta=ln(D). A factor exp(n eta) is
 * factorized out of the kernels, hence y' = omega.y - n.y + rhs */
#define UPDATE_OMEGA_MATRIX                                                    \
  omega[0][1] = 1;                                                             \
  omega[1][0] = zeta * (1 - kappa);                                            \
  omega[1][1] = -zeta + 1;                                                     \
  omega[1][2] = zeta * kappa;                                                  \
  omega[2][3] = 1;                                                             \
  omega[3][0] = zeta * (1 - kappa);                                            \
  omega[3][2] = zeta * (kappa - k * k * xi);                                   \
  omega[3][3] = -zeta + 1;

/*Alternative definition with only two components and zeta(eta) dependence (EdS
 * relaxation test)*/
/*#define UPDATE_OMEGA_MATRIX \*/
/*omega[0][1] = 1; \*/
/*omega[1][0] = zeta *: \*/
/*omega[1][1] = -zeta + 1; \*/

#endif /* #ifndef OMEGA_MATRIX_HPP */
