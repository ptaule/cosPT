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

#include "include/constants.h"
#include "include/utilities.h"
#include "include/tables.h"
#include "include/io.h"
#include "include/diagrams.h"
#include "include/integrand.h"
#include "include/evolve_kernels.h"


#define INPUT "/home/t30/ben/ge52sir/non_linear_PS/input/"


int main (int argc, char* argv[]) {
    const char* input_ps_file       = INPUT "massive_nu_0.1eV_z2_pk_cb.dat";
    /* const char* input_zeta_file     = INPUT "zeta.dat"; */
    const char* input_redshift_file = INPUT "redshift.dat";
    const char* output_ps_file      = "cb_cb_L" TOSTRING(LOOPS) ".dat";

    char* ic_perturbations_files[2];
    ic_perturbations_files[0] = INPUT "delta_nu_over_delta_cb_lin_z25.dat";
    ic_perturbations_files[1] = INPUT "theta_nu_over_aHf_delta_cb_lin_z25.dat";

    if (argc == 2) {
        input_ps_file = argv[1];
    }
    else if (argc == 3) {
        input_ps_file = argv[1];
        output_ps_file = argv[2];
    }

    printf("LOOPS                 = %d\n", LOOPS);
    printf("COMPONENTS            = %d\n", COMPONENTS);
    printf("TIME STEPS            = %d\n", TIME_STEPS);
    printf("Reading input power spectrum from %s.\n",input_ps_file);
    /* printf("Reading input zeta function from %s.\n",input_zeta_file); */
    printf("Results will be written to %s.\n",output_ps_file);

#if N_CORES >= 0
    cubacores(N_CORES, 10000);
    printf("Using %d cores.\n",N_CORES);
#endif

    // Array of table_ptrs, one for each worker (thread)
    tables_t table;

    // Initialize time steps in eta
    double eta[TIME_STEPS];
    initialize_timesteps(eta, ETA_I, ETA_F);

    // Sum table can be computed right away
    short int sum_table[N_CONFIGS][N_CONFIGS];
    compute_sum_table(sum_table);
        table.sum_table = (const short int (*)[])sum_table;
        table.eta = eta;

    // Initialize diagrams to compute at this order in PT
    diagram_t diagrams[N_DIAGRAMS];
    initialize_diagrams(diagrams);

    evolution_params_t params = {
        /* .zeta_acc = NULL, */
        /* .zeta_spline = NULL, */
        .redshift_acc = NULL,
        .redshift_spline = NULL,
        .omega = gsl_matrix_alloc(COMPONENTS, COMPONENTS)
    };

    integration_input_t input = {
        .k = 0.0,
        .component_a = 0,
        .component_b = 0,
        .ps_acc = NULL,
        .ps_spline = NULL,
        .diagrams = diagrams,
        .params = &params,
        .worker_mem = &table
    };

    integration_variables_t vars;

    // Read input files and interpolate
    read_and_interpolate(input_ps_file,&input.ps_acc,&input.ps_spline);
    read_and_interpolate(input_redshift_file,&params.redshift_acc,&params.redshift_spline);
    /* read_and_interpolate(input_zeta_file,&params.zeta_acc,&params.zeta_spline); */

    for (int i = 0; i < 2; ++i) {
        read_and_interpolate(ic_perturbations_files[i],
                &params.ic_perturb_accs[i], &params.ic_perturb_splines[i]);
    }

    double delta_factor = pow((double)K_MAX/K_MIN, 1.0/N_POINTS);
    double k = K_MIN;

    for (int i = 0; i < N_POINTS; ++i) {
        k *= delta_factor;

        table.Q_magnitudes = vars.magnitudes;

        // Initialize sum-, bare_scalar_products-, alpha- and beta-table
        compute_bare_scalar_products(input.k, &vars,
                table.bare_scalar_products);
        compute_scalar_products((const vfloat (*)[])table.bare_scalar_products,
                table.scalar_products);
        compute_alpha_beta_tables((const vfloat (*)[])table.scalar_products,
                table.alpha, table.beta);

        integrand(&input, &table);
    }

    diagrams_gc(diagrams);

    gsl_matrix_free(params.omega);

    gsl_spline_free(input.ps_spline);
    gsl_interp_accel_free(input.ps_acc);
    gsl_spline_free(params.redshift_spline);
    gsl_interp_accel_free(params.redshift_acc);
    /* gsl_spline_free(params.zeta_spline); */
    /* gsl_interp_accel_free(params.zeta_acc); */

    for (int i = 0; i < 2; ++i) {
        gsl_interp_accel_free(params.ic_perturb_accs[i]);
        gsl_spline_free(params.ic_perturb_splines[i]);
    }

    return 0;
}
