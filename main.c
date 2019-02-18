/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

#include "include/constants.h"
#include "include/utilities.h"
#include "include/kernels.h"
#include "include/spt_kernels.h"

void testVectorSum();
void testAlphaBeta();
void testKernelComputer();
void print_configs();
void test_kernel_index_from_arguments();



int main () {
    debug_print("LOOPS      = %d\n", LOOPS);
    debug_print("N_CONFIGS  = %d\n", N_CONFIGS);
    debug_print("N_KERNELS  = %d\n", N_KERNELS);
    debug_print("COMPONENTS = %d\n", COMPONENTS);
    debug_print("ZERO_LABEL = %d\n", ZERO_LABEL);

    testKernelComputer();
}



void testKernelComputer() {
    vfloat k = 2;
    vfloat Q = 1;
    vfloat mu = 0.5;

    matrix_vfloat* alpha = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);
    matrix_vfloat* beta = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_alpha_beta_tables(k,Q,mu,alpha,beta);

    // Allocate space for kernels (calloc also initializes values to 0)
    kernel_value* kernels = (kernel_value*)calloc(COMPONENTS * N_KERNELS, sizeof(kernel_value));

    short int args[N_KERNEL_ARGS] = {2,0,4};
    vfloat value = compute_SPT_kernel(args,0,alpha,beta,kernels);
    printf("value = %f\n",value);

    // Free allocated memory
    free(kernels);
    gsl_matrix_free(alpha);
    gsl_matrix_free(beta);
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


void testAlphaBeta() {
    double k = 2;
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



void testVectorSum() {
    short int args[3] = {9,7,5};

    for (int j = 0; j < 3; ++j) {
        short int config[N_COEFFS];
        label2config(args[j],config,N_COEFFS);
        for (int i = 0; i < N_COEFFS; ++i) {
            printf("%d,",config[i]);
        }
        printf(" + ");
    }
    printf("=\n");

    short int sum = sum_vectors(args,3);

    short int sumConfig[N_COEFFS];
    label2config(sum,sumConfig,N_COEFFS);

    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",sumConfig[i]);
    }
    printf("\n");
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

