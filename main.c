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



int main () {
    char* input_ps_file = "/home/t30/all/ge52sir/non_linear_PS/input/PS_linear_z000_pk.dat";
    char* output_ps_file = "/space/ge52sir/non_linear_PS/output/spt_" TOSTRING(LOOPS) "loop.dat";

    printf("LOOPS         = %d\n", LOOPS);
    printf("COMPONENTS    = %d\n", COMPONENTS);
    printf("Reading input power spectrum from %s.\n",input_ps_file);
    printf("Results will be written to %s.\n",output_ps_file);

    gsl_interp_accel* acc;
    gsl_spline* spline;

    read_PS(input_ps_file,&acc,&spline);

    // Compute table of sums of two and two vector labels
    short int sum_table[N_CONFIGS][N_CONFIGS];
    compute_sum_table(sum_table);

    // Possible diagrams (m,l,r) at this loop order
    diagram_t diagrams[N_DIAGRAMS];
    possible_diagrams(diagrams);

    integration_input_t input = {
        .k = 0.0,
        .component_a = 0,
        .component_b = 0,
        .acc = acc,
        .spline = spline,
        .sum_table = sum_table,
        .diagrams = diagrams
    };

    // Overall factors:
    // - Only integrating over cos_theta_i between 0 and 1, multiply by 2 to
    //   obtain [-1,1]
    // - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
    // - Phi integration of first loop momenta gives a factor 2pi
    // - Conventionally divide by (2pi)^3
    vfloat overall_factor = pow(2,LOOPS) * gsl_sf_fact(LOOPS) * TWOPI;

    /* int nregions, neval, fail; */
    cubareal result[1];

    double delta_logk = log(K_MAX/K_MIN) / N_POINTS;
    double k = K_MIN;

    integration_variables_t vars = {
        .magnitudes = {0.3},
        .cos_theta  = {1},
        /* .phi        = {0} */
    };


    for (int i = 0; i < N_POINTS; ++i) {
        k *= exp(delta_logk);
        input.k = k;

        result[0]= integrand(&input,&vars);
        result[0] *= overall_factor;

        printf("k = %f, result = %e \n", k, (double)*result);
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return 0;
}
