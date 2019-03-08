/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_matrix.h>

#include <cuba.h>

#include "include/constants.h"
#include "include/utilities.h"
#include "include/kernels.h"
#include "include/spt_kernels.h"
#include "include/input_PS.h"
#include "include/integrand.h"



int cuba_integrand(
        const int *ndim,
        const cubareal xx[],
        const int *ncomp,
        cubareal ff[],
        void *userdata
        )
{
    integration_input_t* data = (integration_input_t*)userdata;
    integration_variables_t vars;

    vfloat log_interval = log(K_MAX/K_MIN);

    vars.magnitudes[0] = K_MIN * exp(log_interval * xx[0]);
    vars.cos_theta[0] = xx[1];

    vfloat jacobian = K_MIN * log_interval;

    // Only integrating over cos_theta between 0 and 1, multiply by 2 to obtain
    // [-1,1] (integrand symmetric w.r.t spatial inversion of loop momentum)
    int sym_factor = 2;

    ff[0] = sym_factor * jacobian * integrand(data,&vars);
    return 0;
}

/* CUBA settings */
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 10e6

#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.



int main () {
    debug_print("LOOPS         = %d\n", LOOPS);
    debug_print("N_CONFIGS     = %d\n", N_CONFIGS);
    debug_print("N_KERNELS     = %d\n", N_KERNELS);
    debug_print("N_KERNEL_ARGS = %d\n", N_KERNEL_ARGS);
    debug_print("ZERO_LABEL    = %d\n", ZERO_LABEL);
    debug_print("COMPONENTS    = %d\n", COMPONENTS);

    gsl_interp_accel* acc;
    gsl_spline* spline;
    const char* filename = "/home/pettertaule/Dropbox/Mathematica/simple00_pk.dat";

    read_input_PS(filename,&acc,&spline);

    integration_input_t data = {
        .k = 0.1,
        .component_a = 0,
        .component_b = 0,
        .acc = acc,
        .spline = spline
    };

    int nregions, neval, fail;
    cubareal result[1], error[1], prob[1];

    Suave(2, 1, cuba_integrand, &data,
            NVEC, EPSREL, EPSABS, VERBOSE | LAST, SEED,
            MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
            STATEFILE, SPIN,
            &nregions, &neval, &fail, result, error, prob);

    printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
            nregions, neval, fail);

    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
            (double)*result, (double)*error, (double)*prob);

    return 0;
}
