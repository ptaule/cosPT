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

// Min/max wavenumber
#define K_MIN 1e-4
#define K_MAX 1e2

// Number of evaluation points between k=K_MIN and k=K_MAX
#define N_POINTS 139

// Input power spectrum file
#define INPUT_FILE "input/PS_linear_z000_pk.dat"
// Output power spectrum to file
#define OUTPUT_FILE "output_" TOSTRING(LOOPS) "loop.dat"



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
        .magnitudes = {0.7,0.2},
        .cos_theta  = {0.5,1},
        .phi        = {0}
    };

    integration_input_t input = {
        .k = 100,
        .component_a = 0,
        .component_b = 0,
        .acc = acc,
        .spline = spline,
        .params = &params,
        .omega = gsl_matrix_alloc(COMPONENTS, COMPONENTS)
    };


    printf("====================================\n");

    double result = integrand(&input,&vars);
    printf("result  = %e\n", result );

    gsl_matrix_free(input.omega);

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return 0;
}
