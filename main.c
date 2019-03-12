/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
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

    integration_input_t input = {
        .k = 0.0,
        .component_a = 0,
        .component_b = 0,
        .acc = acc,
        .spline = spline
    };
    integration_variables_t vars;

    vfloat ratio = K_MAX/K_MIN;

    double xx[5] = {0};
    vfloat jacobian = 0.0;

    // Overall factors:
    // - Only integrating over cos_theta_i between 0 and 1, multiply by 2 to
    //   obtain [-1,1]
    // - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
    // - Phi integration of first loop momenta gives a factor 2pi
    // - Conventionally divide by (2pi)^3
    vfloat overall_factor = pow(2,LOOPS) * gsl_sf_fact(LOOPS) * TWOPI;

    double delta_logk = log(K_MAX/K_MIN) / N_POINTS;
    double k = K_MIN;

    for (int i = 0; i < N_POINTS; ++i) {
        double result = 0;
        k *= exp(delta_logk);
        input.k = k;

        for (int i = 0; i < (int)1e3; ++i) {
            for (int i = 0; i < 5; ++i) {
                xx[i] = (double)rand() / (double)RAND_MAX;
            }
            switch (LOOPS) {
                case 1:
                        vars.magnitudes[0] = K_MIN * pow(ratio,xx[0]);
                        vars.cos_theta[0] = xx[1];
                        jacobian = K_MIN * log(ratio);
                        break;
                case 2:
                    vars.magnitudes[0] = K_MIN * pow(ratio,xx[0]);
                    vars.magnitudes[1] = K_MIN * pow(ratio, xx[0] * xx[1]);
                    vars.cos_theta[0] = xx[2];
                    vars.cos_theta[1] = xx[3];
                    vars.phi[0] = xx[4] * TWOPI;
                    jacobian = TWOPI * pow(K_MIN * log(ratio),2)
                        * xx[0] * pow(K_MIN/K_MAX,xx[0]*(1 + xx[1]));
                    break;
                default:
                    warning_verbose("No jacobian for LOOPS = %d",LOOPS);
            }

            result += jacobian * integrand(&input,&vars);

        }

        result *= overall_factor/1e3;

        printf("k  = %f, result = %f\n", k, result);
    }

    return 0;
}
