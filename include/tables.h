/*
   tables.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef KERNELS_H
#define KERNELS_H

#include <stdbool.h>
#include "constants.h"

typedef struct {
    vfloat value;
    bool computed;
} kernel_t;

// Struct storing integration variables
typedef struct {
    vfloat magnitudes[LOOPS]; /* Loop momenta magnitudes                                     */
    vfloat cos_theta[LOOPS];  /* Cosine of polar angles of the loop momenta                  */

    /* We assume that the first loop momenta has azimuthal angle 0 */
#if LOOPS >= 2
    vfloat phi[LOOPS - 1];
#endif
} integration_variables_t;

// Struct storing bare_scalar_products and pointers to various other data
// tables.
typedef struct {
    const vfloat* Q_magnitudes;
    vfloat bare_scalar_products[N_COEFFS][N_COEFFS];
    const short int (*sum_table)[N_CONFIGS];
    matrix_t* alpha;
    matrix_t* beta;
    kernel_t* kernels;
} table_pointers_t;

short int sum_vectors(
        const short int labels[],
        size_t n_vecs,
        const short int (*sum_table)[N_CONFIGS]
        );

void compute_sum_table(short int sum_table[][N_CONFIGS]);

void compute_bare_scalar_products(
        vfloat k,
        const integration_variables_t* vars,
        vfloat bare_scalar_products[][N_COEFFS]
        );

void allocate_tables(table_pointers_t* tables);
void zero_initialize_tables(table_pointers_t* tables);
void gc_tables(table_pointers_t* tables);

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


inline short int combined_kernel_index(short int argument_index, short int component)
{
    return argument_index * COMPONENTS + component;
}

#endif /* ifndef KERNELS_H */
