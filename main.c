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
#define OUTPUT_FILE \
    "output_" TOSTRING(LOOPS) "loop_fixed_momenta.dat"



int main () {
    printf("LOOPS         = %d\n", LOOPS);
    printf("N_CONFIGS     = %d\n", N_CONFIGS);
    printf("N_KERNELS     = %d\n", N_KERNELS);
    printf("N_KERNEL_ARGS = %d\n", N_KERNEL_ARGS);
    printf("ZERO_LABEL    = %d\n", ZERO_LABEL);
    printf("COMPONENTS    = %d\n", COMPONENTS);

    gsl_interp_accel* acc;
    gsl_spline* spline;

    read_PS(INPUT_FILE,&acc,&spline);

    integration_input_t data = {
        .k = 0.0,
        .component_a = 0,
        .component_b = 0,
        .acc = acc,
        .spline = spline
    };

    integration_variables_t vars = {
        .magnitudes = {0.3,0.2},
        .cos_theta  = {1,1},
        .phi        = {0}
    };

    double* wavenumbers    = (double*)calloc(N_POINTS, sizeof(double));
    double* power_spectrum = (double*)calloc(N_POINTS, sizeof(double));

    double delta_logk = log(K_MAX/K_MIN) / N_POINTS;
    double k = K_MIN;

    for (int i = 0; i < N_POINTS; ++i) {
        k *= exp(delta_logk);
        wavenumbers[i] = k;
        data.k = k;

        power_spectrum[i] = integrand(&data,&vars) / gsl_spline_eval(data.spline,k,data.acc);

        printf("k  = %f, result = %e,\n", k, power_spectrum[i]);
    }

    write_PS(OUTPUT_FILE,N_POINTS,wavenumbers,power_spectrum);

    free(wavenumbers);
    free(power_spectrum);

    return 0;
}
