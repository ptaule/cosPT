/*
   kernels.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef KERNELS_H
#define KERNELS_H

#include "constants.h"

typedef struct {
    vfloat value;
    bool computed;
} kernel_value_t;

typedef struct {
    vfloat magnitudes[LOOPS]; /* Loop momenta magnitudes                                     */
    vfloat cos_theta[LOOPS];  /* Cosine of polar angles of the loop momenta                  */
    vfloat phi[LOOPS - 1];    /* We assume that the first loop momenta has azimuthal angle 0 */
} integration_variables_t;

void compute_bare_scalar_products(
        vfloat k,
        const integration_variables_t* vars,
        vfloat bare_scalar_products[][N_COEFFS]
        );

void compute_scalar_products(
        const vfloat bare_scalar_products[][N_COEFFS],
        matrix_t* scalar_products
        );

void compute_alpha_beta_tables(
        const vfloat bare_scalar_products[][N_COEFFS],
        matrix_t* alpha,
        matrix_t* beta
        );

short int kernel_index_from_arguments(const short int arguments[]);


inline short int combined_kernel_index(
        short int argument_index,short int component)
{
    return argument_index * COMPONENTS + component;
}

#endif /* ifndef KERNELS_H */
