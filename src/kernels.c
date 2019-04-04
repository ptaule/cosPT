/*
   kernels.c

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_matrix.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/kernels.h"


void compute_bare_scalar_products(
        vfloat k,                               /* in, magnitude of k vector         */
        const integration_variables_t* vars,    /* in, loop momenta components       */
        vfloat bare_scalar_products[][N_COEFFS] /* out, scalar product combinations  */
        )
{
    // Diagonal products correspond to Q1*Q1, etc.
    for (int i = 0; i < N_COEFFS - 1; ++i) {
        bare_scalar_products[i][i] = vars->magnitudes[i] * vars->magnitudes[i];
    }
    bare_scalar_products[N_COEFFS - 1][N_COEFFS - 1] = k*k;

    // Products involving k and Q_i has the form k*Q_i*cos(theta_i)
    for (int i = 0; i < N_COEFFS - 1; ++i) {
        vfloat value = k * vars->magnitudes[i] * vars->cos_theta[i];
        bare_scalar_products[N_COEFFS - 1][i] = value;
        bare_scalar_products[i][N_COEFFS - 1] = value;
    }

#if LOOPS >= 2
    // Compute Q_1 * Q_i
    // (This is a special case since phi_1 is chosen to be zero.)
    vfloat cos_theta_1 = vars->cos_theta[0];
    vfloat sin_theta_1 = sqrt(1 - pow(cos_theta_1,2));
    vfloat Q_1 = vars->magnitudes[0];
    for (int i = 1; i < N_COEFFS - 1; ++i) {
        vfloat sin_theta_i = sqrt(1 - pow(vars->cos_theta[i],2));
        vfloat value =
            sin_theta_1 * cos(vars->phi[i-1]) * sin_theta_i
            + cos_theta_1 * vars->cos_theta[i];
        value *= Q_1 * vars->magnitudes[i];

        bare_scalar_products[i][0] = value;
        bare_scalar_products[0][i] = value;
    }

    // Compute Q_i * Q_j for {i,j} != 1
    // (Q_i symbol is 1-indexed while arrays are 0-indexed)
    for (int i = 1; i < N_COEFFS - 1; ++i) {
        for (int j = 1; j < i; ++j) {
            vfloat sin_theta_i = sqrt(1 - pow(vars->cos_theta[i],2));
            vfloat sin_theta_j = sqrt(1 - pow(vars->cos_theta[j],2));

            vfloat value =
                cos(vars->phi[i-1]) * sin_theta_i * cos(vars->phi[j-1]) * sin_theta_j
              + sin(vars->phi[i-1]) * sin_theta_i * sin(vars->phi[j-1]) * sin_theta_j
              + vars->cos_theta[i] * vars->cos_theta[j];

            value *= vars->magnitudes[i] * vars->magnitudes[j];
            bare_scalar_products[i][j] = value;
            bare_scalar_products[j][i] = value;
        }
    }
#endif
}


void compute_scalar_products(
        const vfloat bare_scalar_products[][N_COEFFS], /* in, bare scalar products         */
        matrix_t* scalar_products                      /* out, scalar product combinations */
        )
{
    short int a_coeffs[N_COEFFS];
    short int b_coeffs[N_COEFFS];

    // Scalar product matrix is symmetric, hence compute indices [a,b] and
    // [b,a] simultaneously
    for (int a = 0; a < N_CONFIGS; ++a) {
        for (int b = 0; b <= a; ++b) {
            vfloat product_value = 0;

            label2config(a,a_coeffs,N_COEFFS);
            label2config(b,b_coeffs,N_COEFFS);
            for (int i = 0; i < N_COEFFS; ++i) {
                for (int j = 0; j < N_COEFFS; ++j) {
                    product_value += a_coeffs[i] * b_coeffs[j]
                        * bare_scalar_products[i][j];
                }
            }
            if (a == b) {
                matrix_set(scalar_products,a,a,product_value);
            }
            else {
                matrix_set(scalar_products,a,b,product_value);
                matrix_set(scalar_products,b,a,product_value);
            }
        }
    }
}



void compute_alpha_beta_tables(
        const vfloat bare_scalar_products[][N_COEFFS], /* in, bare scalar products   */
        matrix_t* alpha,                               /* out, matrix of alpha-func. */
        matrix_t* beta                                 /* out, matrix of beta-func.  */
        )
{
    matrix_t* scalar_products = matrix_alloc(N_CONFIGS,N_CONFIGS);
    compute_scalar_products(bare_scalar_products,scalar_products);

    for (int a = 0; a < N_CONFIGS; ++a) {
        for (int b = 0; b < N_CONFIGS; ++b) {
            // Special case when a == b
            if (a == b) {
                matrix_set(alpha,a,b,2.0);
                matrix_set(beta,a,b,2.0);
                continue;
            }

            vfloat alpha_val = 0.0;
            vfloat beta_val  = 0.0;

            // If the first argument is the zero-vector, alpha and beta remains 0
            // If the second argument is the zero-vector, beta remains 0
            if (a != ZERO_LABEL) {
                vfloat product_ab = matrix_get(scalar_products,a,b);
                vfloat product_aa = matrix_get(scalar_products,a,a);

                alpha_val = 1 + product_ab/product_aa;

                if (b != ZERO_LABEL) {
                    vfloat product_bb = matrix_get(scalar_products,b,b);

                    beta_val = product_ab / 2.0
                        * ( 1.0 / product_aa + 1.0 / product_bb
                          + 2.0 * product_ab / (product_aa * product_bb)
                          );
                }
            }
            matrix_set(alpha,a,b,alpha_val);
            matrix_set(beta,a,b,beta_val);
        }
    }
    matrix_free(scalar_products);
}



// This function computes (addition to) kernel index for fundamental vector
// arguments. (It assumes that the arguments are fundamentals.)
inline static short int kernel_index_from_fundamental(short int argument) {
    short int coeffs[N_COEFFS] = {0};
    label2config(argument,coeffs,N_COEFFS);

    // The last coefficient is for k, hence we can skip this (j < N_COEFFS - 1)
    for (int i = 0; i < N_COEFFS - 1; ++i) {
         /* if - Q_i is present, return 2^(2i + 0/2) = 2^(2i)   */
         /* if + Q_i is present, return 2^(2i + 2/2) = 2^(2i+1) */
        if (coeffs[i] != 0) {
            return pow(2, 2 * i + (coeffs[i] + 1)/2);
        }
    }
    return 0;
}



short int kernel_index_from_arguments(const short int arguments[]) {
    // In DEBUG-mode, check that non-zero arguments (zero_label) are unique
#if DEBUG >= 1
    if (!unique_elements(arguments,N_KERNEL_ARGS,ZERO_LABEL))
        warning("Duplicate vector arguments passed to kernel.");
    short int n_k_vectors = 0;
#endif

    short int index = 0;
    // Define block size. A block consists of all fundamental vector argument
    // combinations.
    short int block_size = pow(4,LOOPS);

    for (int i = 0; i < N_KERNEL_ARGS; ++i) {
        // First, check if argument is a zero vector
        if (arguments[i] == ZERO_LABEL) continue;

        // Argument is a k-type vector (i.e. on the form k + c_i Q_i) if k is
        // present. In our vector-label convention, k is the last coefficient,
        // hence k is present if label >= N_CONFIGS/2
        if (arguments[i] >= N_CONFIGS/2) {
            index += (arguments[i] - N_CONFIGS/2 + 1) * block_size;
#if DEBUG >= 1
            n_k_vectors++;
#endif
        }
        else {
            // In DEBUG-mode, check that this is in fact a fundamental vector
#if DEBUG >= 1
            if(!is_fundamental(arguments[i]))
                warning("Kernel argument is neither 0, k-type, nor "
                        "fundamental.");
#endif

            index += kernel_index_from_fundamental(arguments[i]);
        }
    }
#if DEBUG >= 1
    if (n_k_vectors > 1)
        warning("More than one kernel argument is k-type.");
#endif

    return index;
}



extern short int combined_kernel_index(short int argument_index, short int component);
