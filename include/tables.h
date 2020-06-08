/*
   tables.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef TABLES_H
#define TABLES_H

#include <stdbool.h>
#include "constants.h"

typedef struct {
    double** values;    /* TIME_STEPS x COMPONENTS table of kernels      */
    vfloat* spt_values; /* SPT_COMPONENTS array of initial kernels (SPT) */
    bool evolved;       /* Evolved by GSL ODE routine?                   */
    bool spt_computed;  /* (SPT) kernel computed?                        */
} kernel_t;

// Struct storing integration variables
typedef struct {
    vfloat magnitudes[LOOPS]; /* Loop momenta magnitudes                    */
    vfloat cos_theta[LOOPS];  /* Cosine of polar angles of the loop momenta */

    /* We assume that the first loop momenta has azimuthal angle 0 */
#if LOOPS >= 2
    vfloat phi[LOOPS - 1];
#endif
} integration_variables_t;

// Struct storing bare_scalar_products and pointers to various other data
// tables.
typedef struct {
    const vfloat* Q_magnitudes;
    const short int (*sum_table)[N_CONFIGS];
    const double* eta;
    vfloat bare_scalar_products[N_COEFFS][N_COEFFS];
    vfloat scalar_products[N_CONFIGS][N_CONFIGS];
    vfloat alpha[N_CONFIGS][N_CONFIGS];
    vfloat beta[N_CONFIGS][N_CONFIGS];
    kernel_t kernels[N_KERNELS];
} tables_t;

short int sum_vectors(
        const short int labels[],
        size_t n_vecs,
        const short int sum_table[][N_CONFIGS]
        );

void compute_sum_table(short int sum_table[][N_CONFIGS]);

void initialize_timesteps(double eta[], double eta_i, double eta_f, double eta_asymp);

void tables_allocate(tables_t* tables);
void tables_zero_initialize(tables_t* tables);
void tables_gc(tables_t* tables);

void compute_bare_scalar_products(
        vfloat k,
        const integration_variables_t* vars,
        vfloat bare_scalar_products[][N_COEFFS]
        );

void compute_scalar_products(
        const vfloat bare_scalar_products[][N_COEFFS],
        vfloat scalar_products[][N_CONFIGS]
        );

void compute_alpha_beta_tables(
        const vfloat scalar_products[][N_CONFIGS],
        vfloat alpha[N_CONFIGS][N_CONFIGS],
        vfloat beta[N_CONFIGS][N_CONFIGS]
        );

short int kernel_index_from_arguments(const short int arguments[]);

#endif /* ifndef TABLES_H */
