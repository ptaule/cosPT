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
#include "include/diagrams.h"
#include "include/integrand.h"

int cuba_integrand(const int *ndim, const cubareal xx[], const int *ncomp,
        cubareal ff[], void *userdata, const int *nvec, const int *core);
void init_worker(table_ptrs_t* worker_mem, const int* core);
void exit_worker(table_ptrs_t* worker_mem, const int* core);



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

    table_ptrs_t worker_mem[CUBA_MAXCORES];

    short int sum_table[N_CONFIGS][N_CONFIGS];
    compute_sum_table(sum_table);

    for (int i = 0; i < CUBA_MAXCORES; ++i) {
        worker_mem[i].sum_table = (const short int (*)[])sum_table;
    }

    cubainit(init_worker, worker_mem);
    cubaexit(exit_worker, worker_mem);

    // Initialize diagrams to compute at this order in PT
    diagram_t diagrams[N_DIAGRAMS];
    initialize_diagrams(diagrams);

    integration_input_t input = {
        .k = 0.0,
        .component_a = 0,
        .component_b = 0,
        .diagrams = diagrams,
        .acc = acc,
        .spline = spline,
        .worker_mem = (table_ptrs_t*)worker_mem
    };

    double* const wavenumbers    = (double*)calloc(N_POINTS, sizeof(double));
    double* const power_spectrum = (double*)calloc(N_POINTS, sizeof(double));
    double* const errors         = (double*)calloc(N_POINTS, sizeof(double));

    // Overall factors:
    // - Only integrating over cos_theta_i between 0 and 1, multiply by 2 to
    //   obtain [-1,1]
    // - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
    // - Phi integration of first loop momenta gives a factor 2pi
    // - Conventionally divide by ((2pi)^3)^(LOOPS)
    vfloat overall_factor =
        pow(2,LOOPS) * gsl_sf_fact(LOOPS) * pow(TWOPI, 1 - 3*LOOPS);

    int nregions, neval, fail;
    cubareal result[1], error[1], prob[1];

    double delta_factor = pow((vfloat)K_MAX/K_MIN, 1.0/N_POINTS);
    double k = K_MIN;

    for (int i = 0; i < N_POINTS; ++i) {
        k *= delta_factor;
        wavenumbers[i] = k;
        input.k = k;

        Suave(N_DIMS, 1, (integrand_t)cuba_integrand, &input, CUBA_NVEC,
                CUBA_EPSREL, CUBA_EPSABS, CUBA_VERBOSE | CUBA_LAST, CUBA_SEED,
                CUBA_MINEVAL, CUBA_MAXEVAL, CUBA_NNEW, CUBA_NMIN,
                CUBA_FLATNESS, CUBA_STATEFILE, CUBA_SPIN, &nregions, &neval,
                &fail, result, error, prob);

        result[0] *= overall_factor;
        error[0] *= overall_factor;

        power_spectrum[i] = (double)result[0];
        errors[i]         = (double)error[0];

        printf("k  = %f, result = %e, error = %f, prob = %f\n",
                k, (double)*result, (double)error[0], (double)prob[0]);
    }

    write_PS(output_ps_file, N_POINTS, wavenumbers, power_spectrum, errors);

    diagrams_gc(diagrams);

    free(wavenumbers);
    free(power_spectrum);
    free(errors);

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return 0;
}



void init_worker(table_ptrs_t* worker_mem, const int* core) {
    // Master core has number 2^15 = 32768
    if (*core == 32768) {
        allocate_tables(&worker_mem[0]);
        return;
    }
    if (*core + 1 >= CUBA_MAXCORES) {
        error_verbose("Tried to start worker %d (in addition to the master "
                "fork), which exceeds MAXCORES = %d.", *core + 1, CUBA_MAXCORES);
    }
    allocate_tables(&worker_mem[*core + 1]);
}



void exit_worker(table_ptrs_t* worker_mem, const int* core) {
    // Master core has number 2^15 = 32768
    if (*core == 32768) {
        gc_tables(&worker_mem[0]);
    }
    else {
        gc_tables(&worker_mem[*core + 1]);
    }
}



int cuba_integrand(
        __attribute__((unused)) const int *ndim,
        const cubareal xx[],
        __attribute__((unused)) const int *ncomp,
        cubareal ff[],
        void *userdata,
        __attribute__((unused)) const int *nvec,
        const int *core
        )
{
    integration_input_t* input = (integration_input_t*)userdata;
    integration_variables_t vars;

    vfloat ratio = (vfloat)Q_MAX/Q_MIN;

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

    // tables points to memory allocated for worker number <*core + 1>
    // (index 0 is reserved for master(
    table_ptrs_t* tables = &input->worker_mem[*core + 1];
    // Set tables to zero
    zero_initialize_tables(tables);

    tables->Q_magnitudes = vars.magnitudes,

    // Initialize sum-, bare_scalar_products-, alpha- and beta-tables
    compute_bare_scalar_products(input->k, &vars,
            tables->bare_scalar_products);
    // Cast bare_scalar_products to const vfloat 2D-array
    compute_alpha_beta_tables((const vfloat (*)[])tables->bare_scalar_products,
            tables->alpha, tables->beta);

    ff[0] = jacobian * integrand(input,tables);
    return 0;
}
