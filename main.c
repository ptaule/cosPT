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
#include "include/io.h"
#include "include/diagrams.h"
#include "include/integrand.h"
#include "include/evolve_kernels.h"

void init_worker(tables_t* worker_mem, const int* core);
void exit_worker(tables_t* worker_mem, const int* core);

int cuba_integrand(const int *ndim, const cubareal xx[], const int *ncomp,
        cubareal ff[], void *userdata, const int *nvec, const int *core);



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

#if N_CORES >= 0
    cubacores(N_CORES, 10000);
    printf("Using %d cores.\n",N_CORES);
#endif

    // Array of table_ptrs, one for each worker (thread)
    tables_t* worker_mem = (tables_t*)malloc(MAXCORES * sizeof(tables_t));

    // Set routines to be run before process forking (new worker)
    cubainit(init_worker, worker_mem);
    cubaexit(exit_worker, worker_mem);

    // Initialize time steps in eta
    double eta[TIME_STEPS];
    initialize_timesteps(eta, ETA_I, ETA_F);

    // Sum table can be computed right away
    short int sum_table[N_CONFIGS][N_CONFIGS];
    compute_sum_table(sum_table);
    for (int i = 0; i < MAXCORES; ++i) {
        worker_mem[i].sum_table = (const short int (*)[])sum_table;
        worker_mem[i].eta = eta;
    }

    // Initialize diagrams to compute at this order in PT
    diagram_t diagrams[N_DIAGRAMS];
    initialize_diagrams(diagrams);

    evolution_params_t params = {
        .zeta_acc = NULL,
        .zeta_spline = NULL,
    };

    integration_input_t input = {
        .k = 0.0,
        .component_a = 0,
        .component_b = 0,
        .ps_acc = NULL,
        .ps_spline = NULL,
        .diagrams = diagrams,
        .params = &params,
        .worker_mem = worker_mem
    };

    read_and_interpolate(input_ps_file,&input.ps_acc,&input.ps_spline);
    read_and_interpolate(input_zeta_file,&params.zeta_acc,&params.zeta_spline);

    output_t output = {
        .input_ps_file = input_ps_file,
        .wavenumbers = (double*)calloc(N_POINTS, sizeof(double)),
        .lin_ps      = (double*)calloc(N_POINTS, sizeof(double)),
        .non_lin_ps  = (double*)calloc(N_POINTS, sizeof(double)),
        .errors      = (double*)calloc(N_POINTS, sizeof(double))
    };

    // Overall factors:
    // - Only integrating over cos_theta_i between 0 and 1, multiply by 2 to
    //   obtain [-1,1] (for each loop momenta)
    // - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
    // - Phi integration of first loop momenta gives a factor 2pi
    // - Conventionally divide by ((2pi)^3)^(LOOPS)
    vfloat overall_factor =
        pow(2,LOOPS) * gsl_sf_fact(LOOPS) * pow(TWOPI, 1 - 3*LOOPS);

    int nregions, neval, fail;
    cubareal result[1], error[1], prob[1];

    double delta_factor = pow((double)K_MAX/K_MIN, 1.0/N_POINTS);
    double k = K_MIN;

    // Timing
    time_t beginning, end;

    for (int i = 0; i < N_POINTS; ++i) {
        k *= delta_factor;
        output.wavenumbers[i] = k;
        input.k = k;


        time(&beginning);

        Suave(N_DIMS, 1, (integrand_t)cuba_integrand, &input, CUBA_NVEC,
                CUBA_EPSREL, CUBA_EPSABS, CUBA_VERBOSE | CUBA_LAST, CUBA_SEED,
                CUBA_MINEVAL, CUBA_MAXEVAL, CUBA_NNEW, CUBA_NMIN,
                CUBA_FLATNESS, CUBA_STATEFILE, CUBA_SPIN, &nregions, &neval,
                &fail, result, error, prob);

        time(&end);

        output.non_lin_ps[i] = (double)result[0] * overall_factor;
        output.errors[i]     = (double)error[0]  * overall_factor;

        /* (F1(z_0) */
        double F1[COMPONENTS];
        compute_F1(k, input.params, eta, F1);

        output.lin_ps[i] = gsl_spline_eval(input.ps_spline, k, input.ps_acc)
            * F1[input.component_a] * F1[input.component_b];

        printf("k = %f, lin_ps = %e, %d-loop = %e, error = %e, prob = %f, "
                "elapsed time = %.0fs\n", k, output.lin_ps[i], LOOPS,
                output.non_lin_ps[i], output.errors[i], (double)prob[0],
                difftime(end,beginning));
    }

    write_PS(output_ps_file, N_POINTS, &output);

    diagrams_gc(diagrams);

    free(worker_mem);

    free(output.wavenumbers);
    free(output.lin_ps);
    free(output.non_lin_ps);
    free(output.errors);

    gsl_spline_free(input.ps_spline);
    gsl_interp_accel_free(input.ps_acc);
    gsl_spline_free(params.zeta_spline);
    gsl_interp_accel_free(params.zeta_acc);

    return 0;
}



void init_worker(tables_t* worker_mem, const int* core) {
    // Master core has number 2^15 = 32768
    if (*core == 32768) {
        tables_allocate(&worker_mem[0]);
        return;
    }
    if (*core + 1 >= MAXCORES) {
        error_verbose("Tried to start worker %d (in addition to the master "
                "fork), which exceeds MAXCORES = %d.", *core + 1, MAXCORES);
    }
    tables_allocate(&worker_mem[*core + 1]);
}



void exit_worker(tables_t* worker_mem, const int* core) {
    // Master core has number 2^15 = 32768
    if (*core == 32768) {
        tables_gc(&worker_mem[0]);
    }
    else {
        tables_gc(&worker_mem[*core + 1]);
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
    tables_t* tables = &input->worker_mem[*core + 1];
    // Set tables to zero
    tables_zero_initialize(tables);

    tables->Q_magnitudes = vars.magnitudes;

    // Initialize sum-, bare_scalar_products-, alpha- and beta-tables
    compute_bare_scalar_products(input->k, &vars,
            tables->bare_scalar_products);
    compute_scalar_products((const vfloat (*)[])tables->bare_scalar_products,
            tables->scalar_products);
    compute_alpha_beta_tables((const vfloat (*)[])tables->scalar_products,
            tables->alpha, tables->beta);

    ff[0] = jacobian * integrand(input,tables);
    return 0;
}
