/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
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

short int kernelIndex(short int arguments[], short int n) { }

int main (void) {

    debug_print("LOOPS = %d\n", LOOPS);
    debug_print("N_COEFFS = %d\n", N_COEFFS);
    debug_print("N_CONFIGS = %d\n", N_CONFIGS);
    debug_print("COMPONENTS = %d\n", COMPONENTS);

    /* double k = 1; */
    /* double Q = 1; */
    /* double mu = 0.5; */

    /* gsl_matrix* alpha = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS); */
    /* gsl_matrix* beta = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS); */


    /* compute_alpha_beta_tables(k,Q,mu,alpha,beta); */

    /* printf("alpha=\n"); */
    /* print_gsl_matrix(alpha,N_CONFIGS,N_CONFIGS); */
    /* printf("beta=\n"); */
    /* print_gsl_matrix(beta,N_CONFIGS,N_CONFIGS); */


    /* gsl_matrix_free(alpha); */
    /* gsl_matrix_free(beta); */

    return 0;
}
