/*
   io.h

   Created by Petter Taule on 25.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef IO_H
#define IO_H

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#include "constants.h"

typedef struct {
    const char* input_ps_file;
    const char* input_zeta_file;
    const char* input_redshift_file;
    const char* input_omega_eigvals_file;
    const char** ic_F1_files;
    const char* description;
    double cuba_epsrel;
    double cuba_epsabs;
    double cuba_maxevals;
    double k;
    double lin_ps[INTEGRAND_COMPONENTS];
    double non_lin_ps[INTEGRAND_COMPONENTS];
    double error[INTEGRAND_COMPONENTS];
    double F1_eta_i[COMPONENTS];
} output_t;

void read_and_interpolate(const char* filename, gsl_interp_accel** acc, gsl_spline** spline);
void read_and_interpolate_2d(
        const char* x_grid_file,
        const char* y_grid_file,
        const char* data_file,
        gsl_interp_accel** x_acc,
        gsl_interp_accel** y_acc,
        gsl_spline2d** spline
        );

double get_wavenumber(const char* filename, int a);
void write_PS(const char* filename, const output_t* output);

#endif /* ifndef IO_H */
