/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>

#include "utilities.h"

static void compute_scalar_products(
        const double k,             // in, overall k-vector (oriented in z-direction)
        const double Q,             // in, loop-momenta, absolute value
        const double mu,            // in, cos of angle between Q and z-axis (k)
        gsl_matrix* scalar_products // out, computed scalar products
        )
{
    double bare_scalar_products[N_COEFFS][N_COEFFS] = {};

    bare_scalar_products[0][0] = k*k;
    bare_scalar_products[0][1] = k*Q*mu;
    bare_scalar_products[1][0] = k*Q*mu;
    bare_scalar_products[1][1] = Q*Q;

    short int a_coeffs[N_COEFFS];
    short int b_coeffs[N_COEFFS];

    // Can potentially also use that only +k appears (not -k)

    // The scalar_products matrix is symmetric, this should be tested for
    // Diagonally symmetric as well as mirrored around middle row/column

    for (int a = 0; a < N_CONFIGS; ++a) {
        for (int b = 0; b < N_CONFIGS; ++b) {
            double product_value = 0;

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

/* // Potential tests for alpha/beta: */
/* // alpha/beta diagonal should always be 2 */

static void compute_alpha_beta_tables(
        const double k,    // in, overall k-vector (oriented in z-direction)
        const double Q,    // in, loop-momenta, absolute value
        const double mu,   // in, cos of angle between Q and z-axis (k)
        gsl_matrix* alpha, // out, matrix of alpha-func. with possible arguments
        gsl_matrix* beta   // out, matrix of beta-func. with possible arguments
        )
{
    gsl_matrix* scalar_products = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_scalar_products(k,Q,mu,scalar_products);

    for (int a = 0; a < N_CONFIGS; ++a) {
        for (int b = 0; b < N_CONFIGS; ++b) {
            double alpha_val = 1 + (double)(gsl_matrix_get(scalar_products,a,b)) / gsl_matrix_get(scalar_products,a,a);
            double beta_val = gsl_matrix_get(scalar_products,a,b) / 2.0
                * ( 1.0 / gsl_matrix_get(scalar_products,a,a)
                  + 1.0 / gsl_matrix_get(scalar_products,b,b)
                  + 2.0 * gsl_matrix_get(scalar_products,a,b) /
                    ( gsl_matrix_get(scalar_products,a,a) * gsl_matrix_get(scalar_products,b,b) )
                  );
            gsl_matrix_set(alpha,a,b,alpha_val);
            gsl_matrix_set(beta,a,b,beta_val);
        }
    }
    gsl_matrix_free(scalar_products);
}



// This function computes (additive) kernel index for fundamental vector arguments
// (It assumes that the arguments are fundamentals)
short int kernel_index_fundamental_vectors(short int arguments[], short int length) {
    short int index = 0;
    short int coeffs[N_COEFFS] = {};
    for (int i = 0; i < length; ++i) {
        label2config(arguments[i],coeffs,N_COEFFS);
        // The last coefficient is for k, hence we can skip this (j < N_COEFFS - 1)
        for (int j = 0; j < N_COEFFS - 1; ++j) {
            /* This formula converts fundamental configs to corresponding index
             * if   Q_i is not present, add 0 * 2^(2j + 1/2) = 0        to index
             * if - Q_i is present, add     1 * 2^(2j + 0/2) = 2^(2i)   to index
             * if + Q_i is present, add     1 * 2^(2j + 2/2) = 2^(2i+1) to index
             */
            index += abs(coeffs[j]) * pow(2, 2 * j + (coeffs[j] + 1)/2);
        }
    }
    return index;
}



short int kernel_index(short int arguments[], short int length) {
    // First, do some checks if DEBUG-mode is on:
    // - The most possible kernel arguments is (2 LOOPS + 1)
    // - No argument should appear twice
    // - Arguments 2-end should be fundamental vectors (no k present, only one loop vector)
    // - No argument should be the 0 (partly covered by the above point, since
    //      a 0 argument is no fundamental)
#if DEBUG
    if (length > 2 * LOOPS + 1) {
        fprintf(stderr, "WARNING: Number of arguments exceeds what is possible "
                "for %d-LOOP.\n",LOOPS);
    }
    if (!unique_elements(arguments,length))
        fprintf(stderr,"WARNING: Duplicate vector arguments passed to kernel_index().\n");

    for (int i = 1; i < length; ++i) {
        if (!is_fundamental(arguments[i]))
            fprintf(stderr,"WARNING: One of the 2nd to last kernel arguments is "
                    "not a fundamental.\n");
    }
    short int zero_config[N_COEFFS] = {};
    short int zero_label = config2label(zero_config,N_COEFFS);
    if (arguments[0] == zero_label)
        fprintf(stderr,"WARNING: kernel_index() given a 0 argument.\n");
#endif
    //-------------------------------------------//

    short int index = 0;

    // A block consists of all fundamental vector argument combinations
    short int block_size = pow(4,LOOPS);

    // Two major cases: is a k-vector-combination present? It is assumed that
    // this can only occur in the first argument. In our vector-label
    // convention, k is the last coefficient, hence k is not present if
    // label < N_CONFIGS/2
    if (arguments[0] < N_CONFIGS/2) {
        // k-vector config is _not_ present
        index = N_KERNELS - block_size;

        // In debug-mode, check that if k=0, then the first argument should be fundamental
#if DEBUG
        if (!is_fundamental(arguments[0]))
            fprintf(stderr,"WARNING: kernel_index(): k not present in first "
                    "argument and the argument is not a fundamental.\n");
#endif
        index += kernel_index_fundamental_vectors(arguments,length);
    }
    else {
        // k-vector config _is_ present
        index = block_size * (arguments[0] - N_CONFIGS/2);
        // Compute addition to index from fundamental vectors (not including
        // first argument, which is a k-vector-combination)
        index += kernel_index_fundamental_vectors(arguments+1,length-1);
    }
    return index;
}



void kernel_index2arguments(short int kernel_index) { };


int main () {
    debug_print("LOOPS      = %d\n", LOOPS);
    debug_print("N_CONFIGS  = %d\n", N_CONFIGS);
    debug_print("N_KERNELS  = %d\n", N_KERNELS);
    debug_print("COMPONENTS = %d\n", COMPONENTS);

    short int arguments[5] = {4,7,1,5,3};
    printf("kernel_index = %d\n",kernel_index(arguments,5));
}




void print_configs() {
    short int configs[N_COEFFS];

    for (short int label = 0; label < N_CONFIGS; ++label) {
        label2config(label,configs,N_COEFFS);
        for (int i = 0; i < N_COEFFS; ++i) {
            printf("%d,",configs[i]);
        }
        printf("\n");
    }
}

void printAlphaBetaTest() {
    double k = 1;
    double Q = 1;
    double mu = 0.5;

    gsl_matrix* alpha = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);
    gsl_matrix* beta = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_alpha_beta_tables(k,Q,mu,alpha,beta);

    printf("alpha=\n");
    print_gsl_matrix(alpha,N_CONFIGS,N_CONFIGS);
    printf("beta=\n");
    print_gsl_matrix(beta,N_CONFIGS,N_CONFIGS);

    gsl_matrix_free(alpha);
    gsl_matrix_free(beta);
}


void test_kernel_index_fundamental_vectors() {
    short int label1 = 7;
    short int label2 = 1;
    short int label3 = 5;
    short int label4 = 3;

    short int arguments[4] = {label1,label2,label3,label4};
    short int config1[N_COEFFS];
    short int config2[N_COEFFS];
    short int config3[N_COEFFS];
    short int config4[N_COEFFS];

    label2config(label1,config1,N_COEFFS);
    label2config(label2,config2,N_COEFFS);
    label2config(label3,config3,N_COEFFS);
    label2config(label4,config4,N_COEFFS);

    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config1[i]);
    }
    printf("\n");
    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config2[i]);
    }
    printf("\n");
    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config3[i]);
    }
    printf("\n");
    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config4[i]);
    }
    printf("\n");

    printf("kernel_index_fundamental_vectors = %d\n",kernel_index_fundamental_vectors(arguments,4));
}

