/*
   gls_diff_eq.c

   Created by Petter Taule on 02.04.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/


#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>

#define h 0.72
#define c 299792
#define OMEGA_0_M 0.26
#define OMEGA_0_B 0.044
#define OMEGA_0_G 2.47e-5 / (h*h)
#define M_NU 0.3
#define F_NU 1 / (OMEGA_0_M*h*h) * M_NU / 93.14

typedef struct {
    double k;
    gsl_matrix* omega;
} input_parameters;


void set_omega_matrix(double eta, double k, gsl_matrix** omega) {
    double hubble_factor = OMEGA_0_M * exp(-3*eta) + OMEGA_0_G * exp(-4*eta) + (1 - OMEGA_0_M - OMEGA_0_G);
    double omega_M = OMEGA_0_M * exp(-3*eta)/ hubble_factor;
    double k_FS = 0.908 * sqrt(OMEGA_0_M) * M_NU/3 * exp(eta/2);

    // First row
    gsl_matrix_set(*omega,0,0,  0);
    gsl_matrix_set(*omega,0,1, -1);
    gsl_matrix_set(*omega,0,2,  0);
    gsl_matrix_set(*omega,0,3,  0);
    // Second row
    gsl_matrix_set(*omega,1,0, -3/2.0 * omega_M * (1 - F_NU) );
    gsl_matrix_set(*omega,1,1,  2 - 3/2.0 * omega_M          );
    gsl_matrix_set(*omega,1,2, -3/2.0 * omega_M * F_NU       );
    gsl_matrix_set(*omega,1,3,  0                          );
    // Third row
    gsl_matrix_set(*omega,2,0,  0);
    gsl_matrix_set(*omega,2,1,  0);
    gsl_matrix_set(*omega,2,2,  0);
    gsl_matrix_set(*omega,2,3, -1);
    // Fourth row
    gsl_matrix_set(*omega,3,0, -3/2.0 * omega_M * (1 - F_NU)             );
    gsl_matrix_set(*omega,3,1,  0                                      );
    gsl_matrix_set(*omega,3,2, -3/2.0 * omega_M * (F_NU - pow(k/k_FS,2)) );
    gsl_matrix_set(*omega,3,3,  2 - 3/2.0 * omega_M                      );
}


int func (double eta, const double y[], double f[], void *params) {
    input_parameters input = *(input_parameters*)params;
    double k = input.k;
    gsl_matrix* omega = input.omega;

    set_omega_matrix(eta, k, &omega);

    // Diff eq reads y' = - Omega.y, therefore multiply Omega by -1
    gsl_matrix_scale(omega, -1);

    gsl_vector_const_view y_vec = gsl_vector_const_view_array(y,4);
    gsl_vector_view f_vec = gsl_vector_view_array(f,4);

    gsl_blas_dgemv(CblasNoTrans, 1, omega, &y_vec.vector, 0, &f_vec.vector);

    return GSL_SUCCESS;
}



int main () {
    input_parameters input;
    input.k = 0.1;
    input.omega = gsl_matrix_alloc(4,4);

    double y[4] = {-555.202, -545.305, -4.41303, -8.6141};
    double eta_ini = -log(26);
    double eta0 = 0;

    gsl_odeiv2_system sys = {func, NULL, 4, &input};

    gsl_odeiv2_driver * driver =
        gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
                1e-6, 1e-6, 1e-6);

    int N = 100;
    double eta = eta_ini;

    for (int i = 1; i <= N; i++)
    {
        double eta_i = eta_ini + i * abs(eta0 - eta_ini)/(double)N;
        int status = gsl_odeiv2_driver_apply (driver, &eta, eta_i, y);

        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }

        printf ("%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", eta, y[0], y[1], y[2], y[3]);
    }

    gsl_matrix_free(input.omega);
    gsl_odeiv2_driver_free (driver);
    return 0;
}

/*
int func (double t, const double y[], double f[], void *params)
{
  double mu = *(double *)params;
  f[0] = y[1];
  f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
  return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
  gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

int main () {
  double mu = 10;
  gsl_odeiv2_system sys = {func, jac, 2, &mu};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = 100.0;
  double y[2] = { 1.0, 0.0 };

  for (i = 1; i <= 100; i++)
    {
      double ti = i * t1 / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }

      printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

  gsl_odeiv2_driver_free (d);
  return 0;
}
*/
