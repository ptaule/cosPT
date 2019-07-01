/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_sf.h>
#include <cuba.h>

#include "include/constants.h"
#include "include/utilities.h"
#include "include/tables.h"
#include "include/io.h"
#include "include/integrand.h"
#include "include/spt_kernels.h"
#include "include/evolve_kernels.h"


void draw_momenta(double* k, integration_variables_t* vars) {
    double ratio = Q_MAX/(double)Q_MIN;
    double rand_f = (double)(rand()) / RAND_MAX;

    *k = Q_MIN * pow(ratio, rand_f);
    printf("k          = %e \n", *k);

    for (int i = 0; i < LOOPS; ++i) {
        rand_f = (double)(rand()) / RAND_MAX;
        vars->magnitudes[i] = Q_MIN * pow(ratio,rand_f);
        vars->cos_theta[i] = 2 * (double)(rand())/RAND_MAX - 1;

        printf("Q%d         = %Le \n", i + 1, vars->magnitudes[0]);
        printf("mu_%d       = %Lf \n", i + 1, vars->cos_theta[0]);
    }

#if LOOPS >= 2
    for (int i = 0; i < LOOPS - 1; ++i) {
        vars->phi[i] = 2 * PI * (rand() / (double) RAND_MAX);

        printf("phi_%d      = %Lf \n", i + 1, vars->phi[0]);
    }
#endif
}



int main () {
    srand(time(0));

    char* input_zeta_file = "input/zeta.dat";

    printf("LOOPS      = %d\n", LOOPS);
    printf("COMPONENTS = %d\n", COMPONENTS);
    printf("TIME STEPS = %d\n", TIME_STEPS);

    gsl_interp_accel *zeta_acc;
    gsl_spline *zeta_spline;

    read_and_interpolate(input_zeta_file,&zeta_acc,&zeta_spline);

    const evolution_params_t params = {
        .zeta_acc = zeta_acc,
        .zeta_spline = zeta_spline,
        .omega = gsl_matrix_alloc(COMPONENTS, COMPONENTS)
    };

    integration_variables_t vars;
    double k;

    draw_momenta(&k, &vars);

    // Store pointers to the computed tables in struct for convenience
    tables_t tables;
    tables.Q_magnitudes = vars.magnitudes;

    // Initialize time steps in eta
    double eta[TIME_STEPS];
    initialize_timesteps(eta, ETA_I, ETA_F);
    tables.eta = eta;

    short int sum_table[N_CONFIGS][N_CONFIGS];
    compute_sum_table(sum_table);
    tables.sum_table = (const short int (*)[])sum_table;

    tables_allocate(&tables);

    // Initialize bare_scalar_products-, alpha- and beta-tables
    compute_bare_scalar_products(k, &vars, tables.bare_scalar_products);
    // Cast bare_scalar_products to const vfloat 2D-array
    compute_alpha_beta_tables((const vfloat (*)[])tables.bare_scalar_products,
            tables.alpha, tables.beta);

    tables_zero_initialize(&tables);

    tables.Q_magnitudes = vars.magnitudes;

    // Initialize sum-, bare_scalar_products-, alpha- and beta-tables
    compute_bare_scalar_products(k, &vars,
            tables.bare_scalar_products);
    // Cast bare_scalar_products to const vfloat 2D-array
    compute_alpha_beta_tables(
            (const vfloat (*)[])tables.bare_scalar_products,
            tables.alpha, tables.beta);

    diagram_t diagrams[N_DIAGRAMS];
    initialize_diagrams(diagrams);


    short int component_a = 0;
    short int component_b = 0;

    for (int i = 0; i < N_DIAGRAMS; ++i) {
        // Pointer alias for convenience
        const diagram_t* const dg = &(diagrams[i]);

        short int m = dg->m;
        short int l = dg->l;
        short int r = dg->r;

        for (int a = 0; a < dg->n_rearrangements; ++a) {
            for (int b = 0; b < dg->n_sign_configs; ++b) {
                const short int* const arguments_l =
                    dg->argument_configs_l[a][b];
                const short int* const arguments_r =
                    dg->argument_configs_r[a][b];
                short int kernel_index_l = dg->kernel_indices_l[a][b];
                short int kernel_index_r = dg->kernel_indices_r[a][b];

                compute_SPT_kernels(arguments_l, kernel_index_l, 2*l + m, &tables);
                compute_SPT_kernels(arguments_r, kernel_index_r, 2*r + m, &tables);
                // Then, evolve kernels
                kernel_evolution(arguments_l, kernel_index_l, 2*l + m, &params,
                        &tables);
                kernel_evolution(arguments_r, kernel_index_r, 2*r + m, &params,
                        &tables);

                printf(ANSI_COLOR_MAGENTA "(m,l,r) = (%d,%d,%d)\n" ANSI_COLOR_RESET,m,l,r);
                printf("%c%d", component_a + 'F', m + 2*l);
                print_labels(arguments_l, N_KERNEL_ARGS);
                printf(" = %Le (SPT) \t|\t %e (LCDM)\n",
                        tables.kernels[kernel_index_l].spt_values[component_a],
                        tables.kernels[kernel_index_l].values[TIME_STEPS - 1][component_a]);

                printf("%c%d", component_b + 'F', m + 2*r);
                print_labels(arguments_r, N_KERNEL_ARGS);
                printf(" = %Le (SPT) \t|\t %e (LCDM)\n",
                        tables.kernels[kernel_index_r].spt_values[component_b],
                        tables.kernels[kernel_index_r].values[TIME_STEPS - 1][component_b]);
            }
        }
    }

    diagrams_gc(diagrams);

    tables_gc(&tables);

    gsl_matrix_free(params.omega);

    gsl_spline_free(zeta_spline);
    gsl_interp_accel_free(zeta_acc);

    return 0;
}
