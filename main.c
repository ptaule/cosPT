/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_sf.h>

#include "include/constants.h"
#include "include/utilities.h"
#include "include/kernels.h"
#include "include/spt_kernels.h"
#include "include/integrand.h"
#include "include/power_spectrum_io.h"

// Min/max wavenumber
#define K_MIN 1e-4
#define K_MAX 1e2

// Number of evaluation points between k=K_MIN and k=K_MAX
#define N_POINTS 139

// Input power spectrum file
#define INPUT_FILE "input/PS_linear_z000_pk.dat"
// Output power spectrum to file
#define OUTPUT_FILE "output_" TOSTRING(LOOPS) "loop.dat"

void testKernelComputer(vfloat k, const integration_variables_t* vars);


int main () {
    gsl_interp_accel* acc;
    gsl_spline* spline;

    read_PS(INPUT_FILE,&acc,&spline);

    integration_variables_t vars = {
        .magnitudes = {0.3,0.2},
        .cos_theta  = {1,1},
        .phi        = {0}
    };

    integration_input_t data = {
        .k = 100,
        .component_a = 0,
        .component_b = 0,
        .acc = acc,
        .spline = spline
    };

    printf("k     = %Lf\n",data.k);
    printf("Q1    = %Lf\n",vars.magnitudes[0]);
    printf("Q2    = %Lf\n",vars.magnitudes[1]);
    printf("cosk1 = %Lf\n",vars.cos_theta[0]);
    printf("cosk2 = %Lf\n",vars.cos_theta[1]);
    printf("phi   = %Lf\n",vars.phi[0]);

    printf("====================================\n");

    double result = integrand(&data,&vars);
    printf("result  = %e\n", result );

    testKernelComputer(data.k,&vars);

    return 0;
}



void testKernelComputer(vfloat k, const integration_variables_t* vars) {
    printf("testKernelComputer\n");

    short int component = 0;

    short int rearrangement[LOOPS] = {0,1};
    short int signs[LOOPS] = {1,1};

    diagram_t diagram = {.m = 1, .l = 1, .r = 1};

    short int args_l[N_KERNEL_ARGS] = {};
    short int args_r[N_KERNEL_ARGS] = {};

    find_kernel_arguments(&diagram,rearrangement,signs,args_l,args_r);

    table_pointers_t data_tables;
    data_tables.alpha = matrix_alloc(N_CONFIGS,N_CONFIGS);
    data_tables.beta  = matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_sum_table(data_tables.sum_table);
    compute_bare_scalar_products(k,vars,data_tables.bare_scalar_products);
    compute_alpha_beta_tables(data_tables.bare_scalar_products,data_tables.alpha,data_tables.beta);

    // Allocate space for kernels (calloc also initializes values to 0)
    data_tables.kernels = (kernel_value_t*)calloc(COMPONENTS * N_KERNELS, sizeof(kernel_value_t));

    vfloat k1 = compute_k1(diagram.m,rearrangement,signs,data_tables.bare_scalar_products);
    printf("k1 = %Le\n",k1);

    vfloat value_l = compute_SPT_kernel(args_l,diagram.m + 2*diagram.l,component,&data_tables);
    vfloat value_r = compute_SPT_kernel(args_r,diagram.m + 2*diagram.r,component,&data_tables);


    printf("F(k");
    short int config[N_COEFFS];
    label2config(args_l[0],config,N_COEFFS);
    for (int i = 0; i < LOOPS; ++i) {
        if (config[i] == 0) continue;
        else if (config[i] == -1) printf("-Q%d",i+1);
        else if (config[i] == 1)  printf("+Q%d",i+1);
    }
    printf(", ");

    for (int i = 1; i < N_KERNEL_ARGS; ++i) {
        if (args_l[i] == ZERO_LABEL) break;
        label2config(args_l[i],config,N_COEFFS);
        for (int j = 0; j < LOOPS; ++j) {
            if (config[j] == 0) continue;
            else if (config[j] == -1) printf("-Q%d, ",j+1);
            else if (config[j] == 1)  printf("+Q%d, ",j+1);
        }
    }
    printf(") = %Le\n",value_l);
    printf("F(k");
    label2config(args_r[0],config,N_COEFFS);
    for (int i = 0; i < LOOPS; ++i) {
        if (config[i] == 0) continue;
        else if (config[i] == -1) printf("-Q%d",i+1);
        else if (config[i] == 1)  printf("+Q%d",i+1);
    }
    printf(", ");
    for (int i = 1; i < N_KERNEL_ARGS; ++i) {
        if (args_r[i] == ZERO_LABEL) break;
        label2config(args_r[i],config,N_COEFFS);
        for (int j = 0; j < LOOPS; ++j) {
            if (config[j] == 0) continue;
            else if (config[j] == -1) printf("-Q%d, ",j+1);
            else if (config[j] == 1)  printf("+Q%d, ",j+1);
        }
    }
    printf(") = %Lf\n",value_r);

    printf("bare_scalar_products:\n");
    for (int i = 0; i < N_COEFFS; ++i) {
        for (int j = 0; j < N_COEFFS; ++j) {
            printf("%Le  ",data_tables.bare_scalar_products[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    // Free allocated memory
    free(data_tables.kernels);
    matrix_free(data_tables.alpha);
    matrix_free(data_tables.beta);
}
