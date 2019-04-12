/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <math.h>

#include "include/constants.h"
#include "include/tables.h"
#include "include/integrand.h"
#include "include/power_spectrum_io.h"
#include "include/spt_kernels.h"
#include "include/evolve_kernels.h"

// Min/max wavenumber
#define K_MIN 1e-4
#define K_MAX 1e2

// Number of evaluation points between k=K_MIN and k=K_MAX
#define N_POINTS 139

// Input power spectrum file
#define INPUT_FILE "input/PS_linear_z000_pk.dat"
// Output power spectrum to file
#define OUTPUT_FILE "output_" TOSTRING(LOOPS) "loop.dat"

void test_kernel_evolution(
        vfloat k,
        gsl_matrix* omega,
        const parameters_t* params,
        const integration_variables_t* vars
        );


int main () {
    gsl_interp_accel* acc;
    gsl_spline* spline;

    read_PS(INPUT_FILE,&acc,&spline);

    const parameters_t params = {
        .omega_m0 = 0.26,
        .f2 = 1,
        .f_nu = 0,
        .eta_i = - log(25 + 1),
        .eta_f = 0,
    };

    integration_variables_t vars = {
        .magnitudes = {10,1},
        .cos_theta  = {0.5,0.5},
        .phi        = {0}
    };

    integration_input_t input = {
        .k = 0.001,
        .component_a = 0,
        .component_b = 0,
        .acc = acc,
        .spline = spline,
        .params = &params,
        .omega = gsl_matrix_alloc(COMPONENTS, COMPONENTS)
    };

    test_kernel_evolution(input.k,input.omega,input.params,&vars);

    gsl_matrix_free(input.omega);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return 0;
}



void test_kernel_evolution(
        vfloat k,
        gsl_matrix* omega,
        const parameters_t* params,
        const integration_variables_t* vars
        )
{
    printf("testKernelComputer\n");

    short int rearrangement[LOOPS] = {0,1};
    short int signs[LOOPS] = {1,1};

    diagram_t diagram = {.m = 1, .l = 1, .r = 1};

    short int args_l[N_KERNEL_ARGS] = {0};
    short int args_r[N_KERNEL_ARGS] = {0};

    find_kernel_arguments(&diagram,rearrangement,signs,args_l,args_r);

    table_pointers_t data_tables;
    data_tables.alpha = matrix_alloc(N_CONFIGS,N_CONFIGS);
    data_tables.beta  = matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_sum_table(data_tables.sum_table);
    compute_bare_scalar_products(k,vars,data_tables.bare_scalar_products);
    compute_alpha_beta_tables((const vfloat (*)[])data_tables.bare_scalar_products,data_tables.alpha,data_tables.beta);

    // Allocate space for kernels (calloc also initializes values to 0)
    data_tables.kernels = (kernel_t*)calloc(N_KERNELS, sizeof(kernel_t));

    // Allocate time/component dimensions of kernels
    for (int i = 0; i < N_KERNELS; ++i) {
        data_tables.kernels[i].values = (double**)calloc(TIME_STEPS, sizeof(double*));
        data_tables.kernels[i].spt_values = (vfloat*)calloc(COMPONENTS, sizeof(vfloat));
        for (int j = 0; j < TIME_STEPS; ++j) {
            data_tables.kernels[i].values[j] = (double*)calloc(COMPONENTS, sizeof(double));
        }
    }

    short int index_l = compute_SPT_kernels(args_l,diagram.m + 2*diagram.l,0,&data_tables);
    short int index_r = compute_SPT_kernels(args_r,diagram.m + 2*diagram.r,0,&data_tables);

    // Then, evolve kernels
    kernel_evolution(args_l, index_l, diagram.m + 2*diagram.l, omega, params, &data_tables);
    kernel_evolution(args_r, index_r, diagram.m + 2*diagram.r, omega, params, &data_tables);

    print_evolved_kernel(args_r,index_r,diagram.m + 2*diagram.r, &data_tables);
    print_evolved_kernel(args_l,index_l,diagram.m + 2*diagram.l, &data_tables);

    // Free allocated memory
    for (int i = 0; i < N_KERNELS; ++i) {
        for (int j = 0; j < TIME_STEPS; ++j) {
            free(data_tables.kernels[i].values[j]);
        }
        free(data_tables.kernels[i].values);
    }
    free(data_tables.kernels);
    matrix_free(data_tables.alpha);
    matrix_free(data_tables.beta);
}
