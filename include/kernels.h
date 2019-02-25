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

void compute_bare_scalar_products(
        vfloat k,
        const vfloat magnitudes[],
        const vfloat cos_theta[],
        const vfloat phi[],
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
void find_kernel_arguments(short int l, short int r, short int m, short int arguments[]);

#endif /* ifndef KERNELS_H */
