/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_sf.h>

#include <cuba.h>

#include "include/constants.h"
#include "include/utilities.h"
#include "include/kernels.h"
#include "include/spt_kernels.h"
#include "include/integrand.h"
#include "include/power_spectrum_io.h"

// Input power spectrum file
#define INPUT_FILE "input/PS_linear_z000_pk.dat"
// Output power spectrum to file
#define OUTPUT_FILE "output_" TOSTRING(LOOPS) "loop.dat"


int cuba_integrand(
        __attribute__((unused)) const int *ndim,
        const cubareal xx[],
        __attribute__((unused)) const int *ncomp,
        cubareal ff[],
        void *userdata
        )
{
    integration_input_t* input = (integration_input_t*)userdata;
    integration_variables_t vars;

    vfloat ratio = Q_MAX/Q_MIN;

    vfloat jacobian = 0.0;
#if LOOPS == 1
    vars.magnitudes[0] = Q_MIN * pow(ratio,xx[0]);
    vars.cos_theta[0] = xx[1];
    jacobian = log(ratio) * pow(vars.magnitudes[0],3);
#elif LOOPS == 2
    vars.magnitudes[0] = Q_MIN * pow(ratio,xx[0]);
    vars.magnitudes[1] = Q_MIN * pow(ratio,xx[0] * xx[1]);
    vars.cos_theta[0] = xx[2];
    vars.cos_theta[1] = xx[3];
    vars.phi[0] = xx[4] * TWOPI;
    jacobian = TWOPI * xx[0]
        * pow(log(ratio),2)
        * pow(vars.magnitudes[0],3)
        * pow(vars.magnitudes[1],3);
#else
    warning_verbose("Monte-carlo integration not implemented for LOOPS = %d.",LOOPS);
#endif

    ff[0] = jacobian * integrand(input,&vars);
    return 0;
}



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

    double* const wavenumbers    = (double*)calloc(N_POINTS, sizeof(double));
    double* const power_spectrum = (double*)calloc(N_POINTS, sizeof(double));
    double* const errors         = (double*)calloc(N_POINTS, sizeof(double));

    // Overall factors:
    // - Only integrating over cos_theta_i between 0 and 1, multiply by 2 to
    //   obtain [-1,1]
    // - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
    // - Phi integration of first loop momenta gives a factor 2pi
    // - Conventionally divide by (2pi)^3
    vfloat overall_factor = pow(2,LOOPS) * gsl_sf_fact(LOOPS) * TWOPI;

    int nregions, neval, fail;
    cubareal result[1], error[1], prob[1];

    double delta_logk = log(K_MAX/K_MIN) / N_POINTS;
    double k = K_MIN;

    for (int i = 0; i < N_POINTS; ++i) {
        k *= exp(delta_logk);
        wavenumbers[i] = k;
        input.k = k;

        Suave(N_DIMS, 1, cuba_integrand, &input, CUBA_NVEC, CUBA_EPSREL,
                CUBA_EPSABS, CUBA_VERBOSE | CUBA_LAST, CUBA_SEED, CUBA_MINEVAL,
                CUBA_MAXEVAL, CUBA_NNEW, CUBA_NMIN, CUBA_FLATNESS,
                CUBA_STATEFILE, CUBA_SPIN, &nregions, &neval, &fail, result,
                error, prob);

        result[0] *= overall_factor;
        error[0] *= overall_factor;

        power_spectrum[i] = (double)result[0];
        errors[i]         = (double)error[0];

        printf("k  = %f, result = %e, error = %f, prob = %f\n",
                k, (double)*result, (double)error[0], (double)prob[0]);
    }

    write_PS(OUTPUT_FILE,N_POINTS,wavenumbers,power_spectrum, errors);

    free(wavenumbers);
    free(power_spectrum);

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return 0;
}
