/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_combination.h>

#include "utilities.h"

static void compute_scalar_products(
        const double k,                /* in, overall k-vector (oriented in z-direction) */
        const double Q,                /* in, loop-momenta, absolute value               */
        const double mu,               /* in, cos of angle between Q and z-axis (k)      */
        matrix_vfloat* scalar_products /* out, computed scalar products                  */
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

// Potential tests for alpha/beta:
// alpha/beta diagonal should always be 2

static void compute_alpha_beta_tables(
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
    if (!unique_elements(arguments,N_KERNEL_ARGS,zero_vector_label()))
        fprintf(stderr,"%s:%d:\tWARNING: Duplicate vector arguments passed to "
                "kernel.\n",__FILE__,__LINE__);
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
    short int zero_label = zero_vector_label();

    for (int i = 0; i < N_KERNEL_ARGS; ++i) {
        // First, check if argument is a zero vector
        if (arguments[i] == zero_label) continue;

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
                fprintf(stderr,"%s:%d:\tWARNING: Kernel argument is neither 0, "
                        "k-type, nor fundamental.\n",__FILE__,__LINE__);
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
        fprintf(stderr,"%s:%d:\tWARNING: More than one kernel argument is "
                "k-type.\n",__FILE__,__LINE__);
#endif
}



short int combined_kernel_index(short int argument_index,short int component) {
    return argument_index * COMPONENTS + component;
}



vfloat compute_SPT_kernel(
        const short int arguments[], /* kernel arguments                                  */
        short int component,         /* component to compute, NB: assumed to be 0-indexed */
        const gsl_matrix* alpha,     /* table of alpha function values for various input  */
        const gsl_matrix* beta,      /* table of beta function values for various input   */
        kernel_value* kernels        /* kernel table                                      */
        )
{
    // Compute kernel index, this depends on arguments (argument_index) and
    // which component is to be computed
    short int argument_index = 0;
    short int n = 0;
    kernel_index_from_arguments(arguments,&argument_index,&n);
    short int index = combined_kernel_index(argument_index,component);

    // First check if the kernel is already computed
    if (kernels[index].computed) return kernels[index].value;

    // For SPT kernels, F_1 = G_1 = ... = 1
    if (n == 1) {
        kernels[index].computed = true;
        return kernels[index].value = 1.0;
    }

    // Define some factors dependent on component to compute
    /* short int a,b; */
    /* if (component == 1) { */
    /*     a = 2 * n + 1; */
    /*     b = 2; */
    /* } */
    /* else { */
    /*     a = 3; */
    /*     b = 2 * n; */
    /* } */

    /* vfloat value = 0.0; */

    for (int m = 1; m < n; ++m) {
        // - comb_l starts at {0,1,...,m} and in the while-loop goes over all
        //   combinations of m elements from {0,...,n} (n choose m possibilities)
        // - comb_r starts at {m+1,...,n} and in the while-loop goes
        //   ("backwards") over all combinations of (n-m) elements from {0,...,n}
        //   (n choose (n-m) possibilities)

        gsl_combination* comb_l = gsl_combination_alloc(n,m);
        gsl_combination* comb_r = gsl_combination_alloc(n,n-m);

        gsl_combination_init_first(comb_l);
        gsl_combination_init_first(comb_r);
    }
    return 0.0;
}



void testKernelComputer() {
    vfloat k = 1;
    vfloat Q = 1;
    vfloat mu = 0.5;

    matrix_vfloat* alpha = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);
    matrix_vfloat* beta = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_alpha_beta_tables(k,Q,mu,alpha,beta);

    // Allocate space for kernels
    kernel_value* kernels = (kernel_value*)malloc(COMPONENTS * N_KERNELS * sizeof(kernel_value));

    for (int i = 0; i < COMPONENTS * N_KERNELS; ++i) {
        kernels[i].value = 0.0;
        kernels[i].computed = false;
    }


    // Free allocated memory
    free(kernels);
    gsl_matrix_free(alpha);
    gsl_matrix_free(beta);
}



int main () {
    debug_print("LOOPS      = %d\n", LOOPS);
    debug_print("N_CONFIGS  = %d\n", N_CONFIGS);
    debug_print("N_KERNELS  = %d\n", N_KERNELS);
    debug_print("COMPONENTS = %d\n", COMPONENTS);

    /* test_kernel_index_from_arguments(); */
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

    matrix_vfloat* alpha = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);
    matrix_vfloat* beta = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_alpha_beta_tables(k,Q,mu,alpha,beta);

    printf("alpha=\n");
    print_gsl_matrix(alpha,N_CONFIGS,N_CONFIGS);
    printf("beta=\n");
    print_gsl_matrix(beta,N_CONFIGS,N_CONFIGS);

    gsl_matrix_free(alpha);
    gsl_matrix_free(beta);
}


void test_kernel_index_from_arguments() {
    short int label1 = 17;
    short int label2 = 5;
    short int label3 = 1;
    short int label4 = 3;
    short int label5 = 7;

    short int arguments[5] = {label1,label2,label3,label4,label5};
    short int config1[N_COEFFS];
    short int config2[N_COEFFS];
    short int config3[N_COEFFS];
    short int config4[N_COEFFS];
    short int config5[N_COEFFS];

    label2config(label1,config1,N_COEFFS);
    label2config(label2,config2,N_COEFFS);
    label2config(label3,config3,N_COEFFS);
    label2config(label4,config4,N_COEFFS);
    label2config(label5,config5,N_COEFFS);

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
    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config5[i]);
    }
    printf("\n");

    short int index = 0;
    short int n = 0;
    kernel_index_from_arguments(arguments,&index,&n);

    printf("kernel_index_from_arguments = %d\n",index);
}

