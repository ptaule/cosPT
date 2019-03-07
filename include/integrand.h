/*
   integrand.h

   Created by Petter Taule on 24.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef INTEGRAND_H
#define INTEGRAND_H

#include <gsl/gsl_spline.h>
#include "constants.h"

typedef struct {
    short int l;
    short int r;
    short int m;
} diagram_t;

typedef struct {
    vfloat k;
    short int component_a;
    short int component_b;
    gsl_interp_accel* acc;
    gsl_spline* spline;
} integration_input_t;

int diagram_factor(const diagram_t* diagram);
int integrand_symmetrization_factor(const diagram_t* diagram);

void possible_diagrams(diagram_t diagrams[]);

void find_kernel_arguments(
        const diagram_t* diagram,
        const short int rearrangement[],
        const short int signs[],
        short int arguments_l[],
        short int arguments_r[]
        );

vfloat compute_k1(short int m, const vfloat bare_scalar_products[][N_COEFFS]);
int heaviside_theta(short int m, vfloat k1, const vfloat Q_magnitudes[]);

/* vfloat integrand_term( */
/*         const diagram_t* diagram, */
/*         const integration_input_t* data, */
/*         const vfloat Q_magnitudes[], */
/*         const vfloat bare_scalar_products[][N_COEFFS], */
/*         const matrix_vfloat* alpha, */
/*         const matrix_vfloat* beta, */
/*         kernel_value_t* kernels */
/*         ); */

void symmetrized_integrand(const diagram_t* diagram);

vfloat integrand(const integration_input_t* data, const integration_variables_t* vars);

#endif /* ifndef INTEGRAND_H */
