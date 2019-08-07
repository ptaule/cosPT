/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "include/worker_mem.h"

void init_worker(worker_mem_t* worker_mem, const int* core);
void exit_worker(worker_mem_t* worker_mem, const int* core);

int cuba_integrand(const int *ndim, const cubareal xx[], const int *ncomp,
        cubareal ff[], void *userdata, const int *nvec, const int *core);
void set_output_filepaths(char output_ps_file[], char cuba_statefile[], double k);


// System specific paths
#define INPUT_PATH  "/home/t30/all/ge52sir/non_linear_PS/input/"
#define OUTPUT_PATH "/space/ge52sir/non_linear_PS/output/"
#define CUBA_STATEPATH "/space/ge52sir/non_linear_PS/output/CUBA_statefiles/"


int main (int argc, char* argv[]) {
    char* input_ps_file     = INPUT_PATH "PS_linear_z000_pk.dat";
    char* input_zeta_file   = INPUT_PATH "zeta.dat";
    char* input_wavenumbers = INPUT_PATH "wavenumbers.dat";

    // Wavenumber index
    int a = 0;

    if (argc == 2) {
        a = atoi(argv[1]);
    }
    else {
        error("Please provide wavenumber index as command argument.");
    }

    double k = get_wavenumber(input_wavenumbers, a);

    char output_ps_file[100];
    char cuba_statefile[100];
    set_output_filepaths(output_ps_file, cuba_statefile, k);

    printf("k                     = %e\n", k);
    printf("LOOPS                 = %d\n", LOOPS);
    printf("COMPONENTS            = %d\n", COMPONENTS);
    printf("TIME STEPS            = %d\n", TIME_STEPS);
    printf("MONTE CARLO MAX EVALS = %.2e\n", CUBA_MAXEVAL);

    printf("Reading input power spectrum from %s.\n", input_ps_file);
    printf("Reading input zeta function from %s.\n", input_zeta_file);
    printf("Results will be written to %s.\n", output_ps_file);
    printf("Cuba statefile: %s.\n", cuba_statefile);

#if N_CORES >= 0
    cubacores(N_CORES, 10000);
    printf("Using %d cores.\n",N_CORES);
#endif

    gsl_interp_accel *ps_acc, *zeta_acc;
    gsl_spline *ps_spline, *zeta_spline;

    read_and_interpolate(input_ps_file, &ps_acc, &ps_spline);
    read_and_interpolate(input_zeta_file, &zeta_acc, &zeta_spline);

    const evolution_params_t params = {
        .zeta_acc = zeta_acc,
        .zeta_spline = zeta_spline,
        .omega = gsl_matrix_alloc(COMPONENTS, COMPONENTS)
    };

    // Initialize time steps in eta
    double eta[TIME_STEPS];
    initialize_timesteps(eta, ETA_I, ETA_F);

    // Sum table can be computed right away
    short int sum_table[N_CONFIGS][N_CONFIGS];
    compute_sum_table(sum_table);

    // Array of table_ptrs, one for each worker (thread) Allocate initially for
    // 10 workers, if more processes are forked, this array is reallocated
    worker_mem_t worker_mem;
    init_worker_mem(&worker_mem, 1, (const short int (*)[])sum_table, eta);

    // Set routines to be run before process forking (new worker)
    cubainit(init_worker, &worker_mem);
    cubaexit(exit_worker, &worker_mem);

    // Initialize diagrams to compute at this order in PT
    diagram_t diagrams[N_DIAGRAMS];
    initialize_diagrams(diagrams);

    integration_input_t input = {
        .k = k,
        .component_a = 0,
        .component_b = 0,
        .ps_acc = ps_acc,
        .ps_spline = ps_spline,
        .diagrams = diagrams,
        .params = &params,
        .worker_mem = &worker_mem
    };

    double* const wavenumbers    = (double*)calloc(1, sizeof(double));
    double* const power_spectrum = (double*)calloc(1, sizeof(double));
    double* const errors         = (double*)calloc(1, sizeof(double));

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

    // Timing
    time_t beginning, end;
    time(&beginning);

    Suave(N_DIMS, 1, (integrand_t)cuba_integrand, &input, CUBA_NVEC,
            CUBA_EPSREL, CUBA_EPSABS, CUBA_VERBOSE | CUBA_LAST, CUBA_SEED,
            CUBA_MINEVAL, CUBA_MAXEVAL, CUBA_NNEW, CUBA_NMIN, CUBA_FLATNESS,
            cuba_statefile, CUBA_SPIN, &nregions, &neval, &fail, result,
            error, prob);

    time(&end);

    result[0] *= overall_factor;
    error[0] *= overall_factor;

    wavenumbers[0]    = k;
    power_spectrum[0] = (double)result[0];
    errors[0]         = (double)error[0];

    printf("k = %f, result = %e, error = %e, prob = %f, "
            "elapsed time = %.0fs\n",
            k, (double)*result, (double)error[0], (double)prob[0],
            difftime(end,beginning));

    /* write_PS(output_ps_file, 1, wavenumbers, power_spectrum, errors); */

    diagrams_gc(diagrams);

    worker_mem_gc(&worker_mem);

    free(wavenumbers);
    free(power_spectrum);
    free(errors);

    gsl_matrix_free(params.omega);

    gsl_spline_free(ps_spline);
    gsl_interp_accel_free(ps_acc);
    gsl_spline_free(zeta_spline);
    gsl_interp_accel_free(zeta_acc);

    return 0;
}



void init_worker(worker_mem_t* worker_mem, const int* core) {
    size_t core_num;
    // Master core has number 2^15 = 32768
    if (*core == 32768) {
        core_num = 0;
    }
    else {
        core_num = (size_t)*core + 1;
    }

    printf("core_num  = %li\n", core_num );

    if (core_num >= worker_mem->size) {
        resize_worker_mem(worker_mem, (core_num * 3)/2);
    }

    allocate_tables(&worker_mem->data[core_num]);
}



void exit_worker(worker_mem_t* worker_mem, const int* core) {
    size_t core_num;
    // Master core has number 2^15 = 32768
    if (*core == 32768) {
        core_num = 0;
    }
    else {
        core_num = (size_t)*core + 1;
    }

    gc_tables(&worker_mem->data[core_num]);
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
    // (index 0 is reserved for master)
    table_ptrs_t* tables = &input->worker_mem->data[*core + 1];
    // Set tables to zero
    zero_initialize_tables(tables);

    // Set sum_table and eta pointers from worker_mem
    tables->sum_table = input->worker_mem->sum_table;
    tables->eta = input->worker_mem->eta;

    tables->Q_magnitudes = vars.magnitudes,

    // Initialize sum-, bare_scalar_products-, alpha- and beta-tables
    compute_bare_scalar_products(input->k, &vars,
            tables->bare_scalar_products);
    // Cast bare_scalar_products to const vfloat 2D-array
    compute_alpha_beta_tables(
            (const vfloat (*)[])tables->bare_scalar_products,
            tables->alpha, tables->beta);

    ff[0] = jacobian * integrand(input,tables);
    return 0;
}



void set_output_filepaths(char output_ps_file[], char cuba_statefile[], double k) {
    char num[20];
    sprintf(num,"%e",k);

    // The name of the output directory is defined through parameter values
    strcpy(output_ps_file, OUTPUT_PATH "lcdm_L" TOSTRING(LOOPS) "_nT"
            TOSTRING(TIME_STEPS) "_N" TOSTRING(CUBA_MAXEVAL));

    // Check that directory exists
    if (!does_directory_exist(output_ps_file)) {
        error_verbose("Directory %s does not exist.", output_ps_file);
    }

    // Name of output file is defined by k-value
    strcat(output_ps_file, "/k_");
    strcat(output_ps_file, num);
    strcat(output_ps_file, ".dat");

    // CUBA statefile
    strcpy(cuba_statefile, CUBA_STATEPATH "lcdm_L" TOSTRING(LOOPS) "_nT"
            TOSTRING(TIME_STEPS) "_N" TOSTRING(CUBA_MAXEVAL) "_k_");
    strcat(cuba_statefile, num);
}
