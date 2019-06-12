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
#include "include/tables.h"
#include "include/integrand.h"
#include "include/power_spectrum_io.h"



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



int main (int argc, char* argv[]) {
    char* input_ps_file   = "input/PS_linear_z000_pk.dat";
    char* input_zeta_file = "input/zeta.dat";
    char* output_ps_file  = "output/PS_" TOSTRING(LOOPS) "loop.dat";

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
    printf("MONTE CARLO MAX EVALS = %.2e\n", CUBA_MAXEVAL);
    printf("Reading input power spectrum from %s.\n",input_ps_file);
    printf("Reading input zeta function from %s.\n",input_zeta_file);
    printf("Results will be written to %s.\n",output_ps_file);

    gsl_interp_accel *ps_acc, *zeta_acc;
    gsl_spline *ps_spline, *zeta_spline;

    read_and_interpolate(input_ps_file,&ps_acc,&ps_spline);
    read_and_interpolate(input_zeta_file,&zeta_acc,&zeta_spline);

    const evolution_params_t params = {
        .eta_i = -3.18,
        .eta_f = 0,
        .zeta_acc = zeta_acc,
        .zeta_spline = zeta_spline,
        .omega = gsl_matrix_alloc(COMPONENTS, COMPONENTS)
    };

    integration_input_t input = {
        .k = 0.0,
        .component_a = 0,
        .component_b = 0,
        .ps_acc = ps_acc,
        .ps_spline = ps_spline,
        .params = &params,
    };

    double* const wavenumbers    = (double*)calloc(N_POINTS, sizeof(double));
    double* const power_spectrum = (double*)calloc(N_POINTS, sizeof(double));
    double* const errors         = (double*)calloc(N_POINTS, sizeof(double));

    // Overall factors:
    // - Only integrating over cos_theta_i between 0 and 1, multiply by 2 to
    //   obtain [-1,1] (for each loop momenta)
    // - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
    // - Phi integration of first loop momenta gives a factor 2pi
    vfloat overall_factor = pow(2,LOOPS) * gsl_sf_fact(LOOPS) * TWOPI;

    int nregions, neval, fail;
    cubareal result[1], error[1], prob[1];

    double delta_factor = pow(K_MAX/K_MIN, 1.0/N_POINTS);
    double k = K_MIN;

    // Timing
    time_t beginning, end;

    for (int i = 0; i < N_POINTS; ++i) {
        k *= delta_factor;
        wavenumbers[i] = k;
        input.k = k;

        time(&beginning);

        Suave(N_DIMS, 1, cuba_integrand, &input, CUBA_NVEC, CUBA_EPSREL,
                CUBA_EPSABS, CUBA_VERBOSE | CUBA_LAST, CUBA_SEED, CUBA_MINEVAL,
                CUBA_MAXEVAL, CUBA_NNEW, CUBA_NMIN, CUBA_FLATNESS,
                CUBA_STATEFILE, CUBA_SPIN,
                &nregions, &neval, &fail, result, error, prob);

        time(&end);

        result[0] *= overall_factor;
        error[0] *= overall_factor;

        power_spectrum[i] = (double)result[0];
        errors[i]         = (double)error[0];

        printf("k = %f, result = %e, error = %e, prob = %f, "
                "elapsed time = %.0fs\n",
                k, (double)*result, (double)error[0], (double)prob[0],
                difftime(end,beginning));
    }

    write_PS(output_ps_file, N_POINTS, wavenumbers, power_spectrum, errors);

    free(wavenumbers);
    free(power_spectrum);
    gsl_matrix_free(params.omega);

    gsl_spline_free(ps_spline);
    gsl_interp_accel_free(ps_acc);

    gsl_spline_free(zeta_spline);
    gsl_interp_accel_free(zeta_acc);

    return 0;
}
