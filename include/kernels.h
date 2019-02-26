/*
   kernels.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef KERNELS_H
#define KERNELS_H

#include <gsl/gsl_matrix.h>
#include "constants.h"

typedef struct {
    vfloat value;
    bool computed;
} kernel_value;

typedef struct {
    vfloat magnitudes[LOOPS]; /* Loop momenta magnitudes                                     */
    vfloat cos_theta[LOOPS];  /* Cosine of polar angles of the loop momenta                  */
    vfloat phi[LOOPS - 1];    /* We assume that the first loop momenta has azimuthal angle 0 */
} integration_variables;

void compute_bare_scalar_products(
        vfloat k,
        const integration_variables* vars,
        vfloat bare_scalar_products[][N_COEFFS]
        );

void compute_scalar_products(
        const vfloat bare_scalar_products[][N_COEFFS],
        matrix_vfloat* scalar_products
        );

void compute_alpha_beta_tables(
        const vfloat bare_scalar_products[][N_COEFFS],
        matrix_vfloat* alpha,
        matrix_vfloat* beta
        );

short int kernel_index_from_fundamental(short int argument);
void kernel_index_from_arguments(const short int arguments[], short int* index, short int* n);
short int combined_kernel_index(short int argument_index,short int component);
void find_kernel_arguments(short int m, short int s, short int arguments[]);

#endif /* ifndef KERNELS_H */
