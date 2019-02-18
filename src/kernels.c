/*
   kernels.c

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/kernels.h"


void compute_scalar_products(
        const double k,                /* in, overall k-vector (oriented in z-direction) */
        const double Q,                /* in, loop-momenta, absolute value               */
        const double mu,               /* in, cos of angle between Q and z-axis (k)      */
        matrix_vfloat* scalar_products /* out, computed scalar products                  */
        )
{
    double bare_scalar_products[N_COEFFS][N_COEFFS] = {};

    // Note that the k-vector has the last coefficient in a configuration
    bare_scalar_products[0][0] = k*k;
    bare_scalar_products[0][1] = Q*k*mu;
    bare_scalar_products[1][0] = Q*k*mu;
    bare_scalar_products[1][1] = Q*Q;

    short int a_coeffs[N_COEFFS];
    short int b_coeffs[N_COEFFS];

    // The scalar_products matrix is symmetric, this should be tested for
    // Diagonally symmetric as well as mirrored around middle row/column

    for (int a = 0; a < N_CONFIGS; ++a) {
        for (int b = 0; b < N_CONFIGS; ++b) {
            vfloat product_value = 0;

            label2config(a,a_coeffs,N_COEFFS);
            label2config(b,b_coeffs,N_COEFFS);
            for (int i = 0; i < N_COEFFS; ++i) {
                for (int j = 0; j < N_COEFFS; ++j) {
                    product_value += a_coeffs[i] * b_coeffs[j]
                        * bare_scalar_products[i][j];
                }
            }
            gsl_matrix_set(scalar_products,a,b,product_value);
        }
    }
}



// Potential tests for alpha/beta:
// alpha/beta diagonal should always be 2
void compute_alpha_beta_tables(
        const double k,    // in, overall k-vector (oriented in z-direction)
        const double Q,    // in, loop-momenta, absolute value
        const double mu,   // in, cos of angle between Q and z-axis (k)
        gsl_matrix* alpha, // out, matrix of alpha-func. with possible arguments
        gsl_matrix* beta   // out, matrix of beta-func. with possible arguments
        )
{
    matrix_vfloat* scalar_products = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_scalar_products(k,Q,mu,scalar_products);

    for (int a = 0; a < N_CONFIGS; ++a) {
        for (int b = 0; b < N_CONFIGS; ++b) {
            vfloat alpha_val = 0.0;
            vfloat beta_val  = 0.0;

            // If the first argument is the zero-vector, alpha and beta remains 0
            // If the second argument is the zero-vector, beta remains 0
            if (a != ZERO_LABEL) {
                alpha_val = 1 +
                    (vfloat)(gsl_matrix_get(scalar_products,a,b)) / gsl_matrix_get(scalar_products,a,a);
                if (b != ZERO_LABEL) {
                    beta_val = gsl_matrix_get(scalar_products,a,b) / 2.0
                        * ( 1.0 / gsl_matrix_get(scalar_products,a,a)
                          + 1.0 / gsl_matrix_get(scalar_products,b,b)
                          + 2.0 * gsl_matrix_get(scalar_products,a,b) /
                          (gsl_matrix_get(scalar_products,a,a) * gsl_matrix_get(scalar_products,b,b))
                          );
                }
            }
            gsl_matrix_set(alpha,a,b,alpha_val);
            gsl_matrix_set(beta,a,b,beta_val);
        }
    }
    gsl_matrix_free(scalar_products);
}



// This function computes (addition to) kernel index for fundamental vector
// arguments. (It assumes that the arguments are fundamentals.)
short int kernel_index_from_fundamental(short int argument) {
    short int index = 0;
    short int coeffs[N_COEFFS] = {};
    label2config(argument,coeffs,N_COEFFS);

    // The last coefficient is for k, hence we can skip this (j < N_COEFFS - 1)
    for (int i = 0; i < N_COEFFS - 1; ++i) {
        /* This formula converts fundamental configs to corresponding index
         * if   Q_i is not present, add 0 * 2^(2i + 1/2) = 0        to index
         * if - Q_i is present, add     1 * 2^(2i + 0/2) = 2^(2i)   to index
         * if + Q_i is present, add     1 * 2^(2i + 2/2) = 2^(2i+1) to index
         */
        index += abs(coeffs[i]) * pow(2, 2 * i + (coeffs[i] + 1)/2);
    }
    return index;
}



void kernel_index_from_arguments(
        const short int arguments[], /* in, kernel arguments                            */
        short int* index,            /* out, kernel index                               */
        short int* n                 /* out, number of non-zero arguments/kernel number */
        )
{
    // In DEBUG-mode, check that non-zero arguments (zero_label) are unique
#if DEBUG
    if (!unique_elements(arguments,N_KERNEL_ARGS,ZERO_LABEL))
        warning("Duplicate vector arguments passed to kernel.");
#endif
    //-------------------------------------------//

    // Define temp index used in calculation
    short int temp_index = 0;
    // Define counters for k-type vectors and fundamentals
    short int n_k_vectors = 0;
    short int n_fundamentals = 0;
    // Define block size. A block consists of all fundamental vector argument
    // combinations.
    short int block_size = pow(4,LOOPS);
    // Store zero-label for comparison

    for (int i = 0; i < N_KERNEL_ARGS; ++i) {
        // First, check if argument is a zero vector
        if (arguments[i] == ZERO_LABEL) continue;

        // Argument is a k-type vector (i.e. on the form k + c_i Q_i) if k is
        // present. In our vector-label convention, k is the last coefficient,
        // hence k is present if label >= N_CONFIGS/2
        if (arguments[i] >= N_CONFIGS/2) {
            n_k_vectors++;
            temp_index += (arguments[i] - N_CONFIGS/2 + 1) * block_size;
        }
        else {
            // In DEBUG-mode, check that this is in fact a fundamental vector
#if DEBUG
            if(!is_fundamental(arguments[i]))
                warning("Kernel argument is neither 0, k-type, nor fundamental.");
#endif

            n_fundamentals++;
            temp_index += kernel_index_from_fundamental(arguments[i]);
        }
    }

    // Set out-parameters
    *index = temp_index;
    *n = n_k_vectors + n_fundamentals;

#if DEBUG
    if (n_k_vectors > 1)
        warning("More than one kernel argument is k-type.");
#endif
}



short int combined_kernel_index(short int argument_index,short int component) {
    return argument_index * COMPONENTS + component;
}
