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


void compute_bare_scalar_products(
        const vfloat magnitudes[],              /* in, vector magnitudes: Q1,Q2,(...),k  */
        const vfloat cos_theta[],               /* in, cosine of polar angles of Q's     */
        const vfloat phi[],                     /* in, azimuthal angles of Q2,Q3,...     */
        vfloat bare_scalar_products[][N_COEFFS] /* out, scalar product combinations      */
        )
{
    // Diagonal products correspond to k*k, Q1*Q1, etc.
    for (int i = 0; i < N_COEFFS; ++i) {
        bare_scalar_products[i][i] = magnitudes[i] * magnitudes[i];
    }

    // Products involving k and Q_i has the form k*Q_i*cos(theta_i)
    vfloat k = magnitudes[N_COEFFS - 1];
    for (int i = 0; i < N_COEFFS - 1; ++i) {
        vfloat value = k * magnitudes[i] * cos_theta[i];
        bare_scalar_products[N_COEFFS - 1][i] = value;
        bare_scalar_products[i][N_COEFFS - 1] = value;
    }

    // Compute Q_1 * Q_i
    // (This is a special case since phi_1 is chosen to be zero.)
    vfloat cos_theta_1 = cos_theta[0];
    vfloat sin_theta_1 = sqrt(1 - pow(cos_theta_1,2));
    vfloat Q_1 = magnitudes[0];

    for (int i = 1; i < N_COEFFS - 1; ++i) {
        vfloat sin_theta_i = sqrt(1 - pow(cos_theta[i],2));
        vfloat value =
            sin_theta_1 * cos(phi[i-1]) * sin_theta_i
            + cos_theta_1 * cos_theta[i];
        value *= Q_1 * magnitudes[i];

        bare_scalar_products[i][0] = value;
        bare_scalar_products[0][i] = value;
    }

    // Compute Q_i * Q_j for {i,j} != 1
    // (Q_i symbol is 1-indexed while arrays are 0-indexed)
    for (int i = 1; i < N_COEFFS - 1; ++i) {
        for (int j = 1; j < i; ++j) {
            vfloat sin_theta_i = sqrt(1 - pow(cos_theta[i],2));
            vfloat sin_theta_j = sqrt(1 - pow(cos_theta[j],2));

            vfloat value =
                cos(phi[i-1]) * sin_theta_i * cos(phi[j-1]) * sin_theta_j
              + sin(phi[i-1]) * sin_theta_i * sin(phi[j-1]) * sin_theta_j
              + cos_theta[i] * cos_theta[j];

            value *= magnitudes[i] * magnitudes[j];
            bare_scalar_products[i][j] = value;
            bare_scalar_products[j][i] = value;
        }
    }
}


void compute_scalar_products(
        const vfloat magnitudes[],     /* in, vector magnitudes: Q1,Q2,(...),k  */
        const vfloat cos_theta[],      /* in, cosine of polar angles of Q's     */
        const vfloat phi[],            /* in, azimuthal angles of Q2,Q3,...     */
        matrix_vfloat* scalar_products /* out, scalar products */
        )
{
    vfloat bare_scalar_products[N_COEFFS][N_COEFFS] = {};

    compute_bare_scalar_products(magnitudes,cos_theta,phi,bare_scalar_products);

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
        const vfloat magnitudes[], /* in, vector magnitudes: Q1,Q2,(...),k               */
        const vfloat cos_theta[],  /* in, cosine of polar angles of Q's                  */
        const vfloat phi[],        /* in, azimuthal angles of Q2,Q3,...                  */
        matrix_vfloat* alpha,      /* out, matrix of alpha-func. with possible arguments */
        matrix_vfloat* beta        /* out, matrix of beta-func. with possible arguments  */
        )
{
    matrix_vfloat* scalar_products = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);
    compute_scalar_products(magnitudes,cos_theta,phi,scalar_products);

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



// For debugging purposes
void kernel_index_to_arguments(short int index, short int arguments[]) {

    // Initialize arguments to zero
    for (int i = 0; i < N_KERNEL_ARGS; ++i) {
        arguments[i] = 0;
    }

    // Define block size. A block consists of all fundamental vector argument
    // combinations.
    short int block_size = pow(4,LOOPS);
    // Argument number we are working with
    short int j = 0;

    // If kernel index is >= block_size, a k-vector is present. We choose the
    // k-vector be the first argument without loss of generality, since the
    // arguments-array is symmetric w.r.t exchanges of arguments.
    if (index >= block_size) {
        arguments[0] += N_CONFIGS/2;
        index -= block_size;
        arguments[0] += index/block_size;

        j++;
        index %= block_size;
    }


    // Check which fundamentals are present
    for (int i = 0; i < 2 * LOOPS; ++i) {
        if (index % 2 == 1) {
            short int fundamental[N_COEFFS] = {};
            fundamental[i/2] = pow(-1,i);

            arguments[j] = config2label(fundamental,N_COEFFS);
            j++;
        }
        index /= 2;
    }
}


short int combined_kernel_index(short int argument_index,short int component) {
    return argument_index * COMPONENTS + component;
}
