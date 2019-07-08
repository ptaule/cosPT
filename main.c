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
#include "include/integrand.h"
#include "include/io.h"



int main (int argc, char* argv[]) {
    char* input_ps_file   = "input/PS_linear_z000_pk.dat";
    char* input_zeta_file = "input/zeta.dat";
    /* char* output_ps_file  = "output/PS_" TOSTRING(LOOPS) "loop.dat"; */

    /* if (argc == 2) { */
    /*     input_ps_file = argv[1]; */
    /* } */
    /* else if (argc == 3) { */
    /*     input_ps_file = argv[1]; */
    /*     output_ps_file = argv[2]; */
    /* } */

    printf("LOOPS                 = %d\n", LOOPS);
    printf("COMPONENTS            = %d\n", COMPONENTS);
    printf("TIME STEPS            = %d\n", TIME_STEPS);
    printf("MONTE CARLO MAX EVALS = %.2e\n", CUBA_MAXEVAL);
    printf("Reading input power spectrum from %s.\n",input_ps_file);
    printf("Reading input zeta function from %s.\n",input_zeta_file);
    /* printf("Results will be written to %s.\n",output_ps_file); */

    gsl_interp_accel *ps_acc, *zeta_acc;
    gsl_spline *ps_spline, *zeta_spline;

    read_and_interpolate(input_ps_file,&ps_acc,&ps_spline);
    read_and_interpolate(input_zeta_file,&zeta_acc,&zeta_spline);

    const evolution_params_t params = {
        .zeta_acc = zeta_acc,
        .zeta_spline = zeta_spline,
        .omega = gsl_matrix_alloc(COMPONENTS, COMPONENTS)
    };

    // Initialize time steps in eta
    double eta[TIME_STEPS];
    initialize_timesteps(eta, ETA_I, ETA_F);

    table_pointers_t data_tables = {.eta = eta};
    allocate_tables(&data_tables);

    double test = data_tables.partial_rhs_sum[1][0][1];
    printf("test = %f\n", test);

    integration_input_t input = {
        .k = 0.0,
        .component_a = 0,
        .component_b = 0,
        .ps_acc = ps_acc,
        .ps_spline = ps_spline,
        .params = &params,
        .data_tables = &data_tables
    };

    double* const wavenumbers    = (double*)calloc(N_POINTS, sizeof(double));

    // Overall factors:
    // - Only integrating over cos_theta_i between 0 and 1, multiply by 2 to
    //   obtain [-1,1] (for each loop momenta)
    // - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
    // - Phi integration of first loop momenta gives a factor 2pi
    vfloat overall_factor = pow(2,LOOPS) * gsl_sf_fact(LOOPS) * TWOPI;

    /* int nregions, neval, fail; */
    cubareal result[1];

    double delta_factor = pow((double)K_MAX/K_MIN, 1.0/N_POINTS);
    double k = K_MIN;

    integration_variables_t vars = {
        .magnitudes = {0.3},
        .cos_theta  = {1},
        /* .phi        = {0} */
    };


    for (int i = 0; i < N_POINTS; ++i) {
        k *= delta_factor;
        wavenumbers[i] = k;
        input.k = k;

        result[0]= integrand(&input,&vars);
        result[0] *= overall_factor;

        printf("k = %f, result = %e \n", k, (double)*result);
    }


    gc_tables(&data_tables);

    free(wavenumbers);
    gsl_matrix_free(params.omega);

    gsl_spline_free(ps_spline);
    gsl_interp_accel_free(ps_acc);

    gsl_spline_free(zeta_spline);
    gsl_interp_accel_free(zeta_acc);

    return 0;
}
