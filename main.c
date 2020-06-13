/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <string.h>

#include <gsl/gsl_sf.h>

#include <cuba.h>

#include "include/constants.h"
#include "include/utilities.h"
#include "include/tables.h"
#include "include/io.h"
#include "include/diagrams.h"
#include "include/integrand.h"
#include "include/evolve_kernels.h"

int max_n_threads;

void init_worker(tables_t* worker_mem, const int* core);
void exit_worker(tables_t* worker_mem, const int* core);

int cuba_integrand(const int *ndim, const cubareal xx[], const int *ncomp,
        cubareal ff[], void *userdata, const int *nvec, const int *core);
void set_output_filepaths(
        char output_ps_file[],
        char cuba_statefile[],
        char description[],
        char output_path[],
        char cuba_statefile_path[],
        int a);

void print_help();
void print_compilation_settings();



int main (int argc, char* argv[]) {
    // Input files
    const char* input_ps_file            = CLASS_PATH "z10_pk_cb.dat";
    const char* input_zeta_file          = CLASS_PATH "zeta_of_etaD.dat";
    const char* input_redshift_file      = CLASS_PATH "redshift_of_etaD.dat";
    const char* input_omega_eigvals_file =
        "input/m_nu_" M_NU_STRING "eV/cg2_growing_mode_eigenvalues_etaD_ini.dat";
    const char* input_wavenumbers   = "input/wavenumbers_bao_zoom.dat";

    const char* ic_F1_files[COMPONENTS];
    ic_F1_files[0] = "input/m_nu_" M_NU_STRING "eV/cg2_F1_growing_mode_etaD_-10.dat";
    ic_F1_files[1] = "input/m_nu_" M_NU_STRING "eV/cg2_F2_growing_mode_etaD_-10.dat";
    ic_F1_files[2] = "input/m_nu_" M_NU_STRING "eV/cg2_F3_growing_mode_etaD_-10.dat";
    ic_F1_files[3] = "input/m_nu_" M_NU_STRING "eV/cg2_F4_growing_mode_etaD_-10.dat";

    // Constants fixed by command line options (and default values)
    max_n_threads        = 100;
    double cuba_epsabs   = 1e-12;
    double cuba_epsrel   = 1e-3;
    double cuba_maxevals = 1e6;
    int cuba_verbose     = 1;

    char description[100] = "";
    char output_path[200] = "/space/ge52sir/non_linear_PS/output/";
    char cuba_statefile_path[200] =
        "/space/ge52sir/non_linear_PS/output/CUBA_statefiles/";

    int c = 0;
    while ((c = getopt(argc, argv, "a:c:C:d:ho:N:r:s:v:")) != -1) {
        switch (c) {
        case 'a':
            cuba_epsabs = atof(optarg);
            break;
        case 'c': {
                      int n_threads = atoi(optarg);
                      printf("Using %d cores.\n", n_threads);
                      cubacores(n_threads, 10000);
                      max_n_threads = n_threads + 1;
                      break;
                  }
        case 'C':
            max_n_threads = atoi(optarg);
            break;
        case 'd':
            strcpy(description, optarg);
            break;
        case 'h':
            print_help();
            return 0;
        case 'o':
            strcpy(output_path, optarg);
            break;
        case 'N':
            cuba_maxevals = atof(optarg);
            break;
        case 'r':
            cuba_epsrel = atof(optarg);
            break;
        case 's':
            strcpy(cuba_statefile_path, optarg);
            break;
        case 'v':
            cuba_verbose = atoi(optarg);
            break;
        case '?':
            if (optopt == 'a' || optopt == 'c' ||
                    optopt == 'C' || optopt == 'd' ||
                    optopt == 'o' || optopt == 'N' ||
                    optopt == 'r' || optopt == 's' ||
                    optopt == 'v')
            {
                error_verbose("Option %c requires a keyword as argument, see manual "
                        "(-h) for more information.", optopt);
            }
            else {
                error("Unknown option.");
            }
        default:
            exit(EXIT_FAILURE);
        }
    }

    int argOffset = optind;
    if (argc != argOffset + 1) {
        error("The number of options and arguments given is not correct. Please "
                "see manual (-h) for more information.");
    }
    // Get wavenumber from wavenumber index
    int wavenumber_index = atoi(argv[argOffset]);
    double k = get_wavenumber(input_wavenumbers, wavenumber_index);

    char output_ps_file[200];
    char cuba_statefile[200];
    set_output_filepaths(output_ps_file, cuba_statefile, description,
            output_path, cuba_statefile_path, wavenumber_index);

    printf("Compilation settings:\n");
    print_compilation_settings();
    printf("\nRuntime settings:\n");
    printf("Monte Carlo max evals = %.2e\n", cuba_maxevals);
    printf("Input power spectrum  = %s.\n", input_ps_file);
    printf("Output file           = %s.\n", output_ps_file);
    printf("Cuba_statefile        = %s.\n", cuba_statefile);

    // Array of table_ptrs, one for each worker (thread)
    tables_t* worker_mem = (tables_t*)malloc(max_n_threads * sizeof(tables_t));

    // Set routines to be run before process forking (new worker)
    cubainit(init_worker, worker_mem);
    cubaexit(exit_worker, worker_mem);

    // Initialize time steps in eta
    double eta[TIME_STEPS];
    initialize_timesteps(eta, ETA_I, ETA_F, ETA_ASYMP);

    // Sum table can be computed right away
    short int sum_table[N_CONFIGS][N_CONFIGS];
    compute_sum_table(sum_table);
    for (int i = 0; i < max_n_threads; ++i) {
        worker_mem[i].sum_table = (const short int (*)[])sum_table;
        worker_mem[i].eta = eta;
    }

    // Initialize diagrams to compute at this order in PT
    diagram_t diagrams[N_DIAGRAMS];
    initialize_diagrams(diagrams);

    evolution_params_t params = {
        .zeta_acc        = NULL,
        .zeta_spline     = NULL,
        .redshift_acc    = NULL,
        .redshift_spline = NULL,
    };

    integration_input_t input = {
        .k           = k,
        .ps_acc      = NULL,
        .ps_spline   = NULL,
        .diagrams    = diagrams,
        .params      = &params,
        .worker_mem  = worker_mem
    };

    // Read input files and interpolate
    read_and_interpolate(input_ps_file, &input.ps_acc, &input.ps_spline);
    read_and_interpolate(input_redshift_file, &params.redshift_acc,
            &params.redshift_spline);
    read_and_interpolate(input_zeta_file, &params.zeta_acc, &params.zeta_spline);
    read_and_interpolate(input_omega_eigvals_file, &params.omega_eigvals_acc,
            &params.omega_eigvals_spline);

    for (int i = 0; i < COMPONENTS; ++i) {
        read_and_interpolate(ic_F1_files[i], &params.ic_F1_accs[i],
                &params.ic_F1_splines[i]);
    }

    output_t output = {
        .input_ps_file            = input_ps_file,
        .input_zeta_file          = input_zeta_file,
        .input_redshift_file      = input_redshift_file,
        .input_omega_eigvals_file = input_omega_eigvals_file,
        .ic_F1_files              = ic_F1_files,
        .description              = description,
        .cuba_epsrel              = cuba_epsrel,
        .cuba_epsabs              = cuba_epsabs,
        .cuba_maxevals            = cuba_maxevals,
        .k                        = k,
        .lin_ps                   = {0.0},
        .non_lin_ps               = {0.0},
        .error                    = {0.0},
        .F1_eta_i                 = {0.0}
    };

    /* Linear evolution */
    double F1_eta_f[COMPONENTS];
    compute_F1(k, input.params, eta, output.F1_eta_i, F1_eta_f);

    output.lin_ps[0] = gsl_spline_eval(input.ps_spline, k, input.ps_acc)
        * F1_eta_f[COMPONENT_A] * F1_eta_f[COMPONENT_A];
    output.lin_ps[1] = gsl_spline_eval(input.ps_spline, k, input.ps_acc)
        * F1_eta_f[COMPONENT_A] * F1_eta_f[COMPONENT_B];
    output.lin_ps[2] = gsl_spline_eval(input.ps_spline, k, input.ps_acc)
        * F1_eta_f[COMPONENT_B] * F1_eta_f[COMPONENT_B];

    /* Non-linear evolution */
    // Overall factors:
    // - Only integrating over cos_theta_i between 0 and 1, multiply by 2 to
    //   obtain [-1,1] (for each loop momenta)
    // - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
    // - Phi integration of first loop momenta gives a factor 2pi
    // - Conventionally divide by ((2pi)^3)^(LOOPS)
    vfloat overall_factor =
        pow(2,LOOPS) * gsl_sf_fact(LOOPS) * pow(TWOPI, 1 - 3*LOOPS);

    int nregions, neval, fail;
    cubareal result[INTEGRAND_COMPONENTS];
    cubareal error[INTEGRAND_COMPONENTS];
    cubareal prob[INTEGRAND_COMPONENTS];

    // Timing
    time_t beginning, end;
    time(&beginning);

    // CUBA settings
#define CUBA_NVEC 1
#define CUBA_LAST 4
#define CUBA_RETAIN_STATEFILE 16
#define CUBA_SEED 0
#define CUBA_MINEVAL 0
#define CUBA_SPIN NULL
#define CUBA_NNEW 1000
#define CUBA_NMIN 2
#define CUBA_FLATNESS 25.
    Suave(N_DIMS, INTEGRAND_COMPONENTS, (integrand_t)cuba_integrand, &input,
            CUBA_NVEC, cuba_epsrel, cuba_epsabs,
            (cuba_verbose | CUBA_LAST | CUBA_RETAIN_STATEFILE), CUBA_SEED,
            CUBA_MINEVAL, cuba_maxevals, CUBA_NNEW, CUBA_NMIN, CUBA_FLATNESS,
            cuba_statefile, CUBA_SPIN, &nregions, &neval, &fail, result, error,
            prob);

    time(&end);

    for (int i = 0; i < INTEGRAND_COMPONENTS; ++i) {
        output.non_lin_ps[i] = (double)result[i] * overall_factor;
        output.error[i]      = (double)error[i]  * overall_factor;
    }

    printf("\nElapsed time: %.0fs\n", difftime(end, beginning));
    printf("k = %e\n", k);

    char* corr_strings[INTEGRAND_COMPONENTS] = {"<AA>","<AB>","<BB>"};
    for (int i = 0; i < INTEGRAND_COMPONENTS; ++i) {
        printf("%s:\tP_lin = % .6e, P_%d-loop = % .6e, err_%d-loop = % .6e, "
                "prob = %f\n", corr_strings[i], output.lin_ps[i], LOOPS,
                output.non_lin_ps[i], LOOPS, output.error[i], (double)prob[i]);
    }

    write_PS(output_ps_file, &output);

    /* Free allocated memory */
    diagrams_gc(diagrams);

    free(worker_mem);

    gsl_spline_free(input.ps_spline);
    gsl_interp_accel_free(input.ps_acc);
    gsl_spline_free(params.redshift_spline);
    gsl_interp_accel_free(params.redshift_acc);
    gsl_spline_free(params.zeta_spline);
    gsl_interp_accel_free(params.zeta_acc);
    gsl_spline_free(params.omega_eigvals_spline);
    gsl_interp_accel_free(params.omega_eigvals_acc);

    for (int i = 0; i < COMPONENTS; ++i) {
        gsl_interp_accel_free(params.ic_F1_accs[i]);
        gsl_spline_free(params.ic_F1_splines[i]);
    }

    return 0;
}



void init_worker(tables_t* worker_mem, const int* core) {
    // Master core has number 2^15 = 32768
    if (*core == 32768) {
        tables_allocate(&worker_mem[0]);
        return;
    }
    if (*core + 1 >= max_n_threads) {
        error_verbose("Tried to start worker %d (in addition to the master "
                "fork), which exceeds MAXCORES = %d.", *core + 1,
                max_n_threads);
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

    double* results = ff;
    for (int i = 0; i < INTEGRAND_COMPONENTS; ++i) results[i] = 0;

    integrand(input, tables, results);

    for (int i = 0; i < INTEGRAND_COMPONENTS; ++i) results[i] *= jacobian;

    return 0;
}



void set_output_filepaths(
        char output_ps_file[], /* Out */
        char cuba_statefile[], /* Out */
        char description[],
        char output_path[],
        char cuba_statefile_path[],
        int a
        )
{
    char a_string[4];
    sprintf(a_string,"%02d",a);

    // The name of the output directory is defined through parameter values
    strcpy(output_ps_file, output_path);
    strcat(output_ps_file, description);
    strcat(output_ps_file, "_L" TOSTRING(LOOPS) "_nT" TOSTRING(TIME_STEPS) "/");

    // Check that directory exists
    if (!does_directory_exist(output_ps_file)) {
        error_verbose("Directory %s does not exist.", output_ps_file);
    }

    // Basename of output file is defined by a-int
    strcat(output_ps_file, a_string);
    strcat(output_ps_file, ".dat");

    // CUBA statefile
    strcpy(cuba_statefile, cuba_statefile_path);
    strcat(cuba_statefile, description);
    strcat(cuba_statefile,  "_L" TOSTRING(LOOPS) "_nT" TOSTRING(TIME_STEPS)
            "_k_index_");
    strcat(cuba_statefile, a_string);
}



void print_help() {
    printf( \
ANSI_COLOR_MAGENTA"cosPT " ANSI_COLOR_RESET "computes 1- and 2-loop \
corrections to various correlation functions in cosmological perturbation \
theory.\nThe program takes one argument: the wavenumber index corresponding to \
the wavenumber at which to compute the power spectrum correction. See "
ANSI_COLOR_RED "input/wavenumbers.dat" ANSI_COLOR_RESET ".\nIn addition, the \
program takes these options (which must be given before the arguments):\n\n"
ANSI_COLOR_BLUE "-a " ANSI_COLOR_RESET
"Specify absolute tolerance for Monte Carlo integrator. Default: 1e-12.\n"
ANSI_COLOR_BLUE "-c " ANSI_COLOR_RESET
"Specify number of threads to spawn. Default: CUBA decides.\n"
ANSI_COLOR_BLUE "-C " ANSI_COLOR_RESET
"Specify maximum number of threads CUBA can spawn. Default: 100.\n"
ANSI_COLOR_BLUE "-d " ANSI_COLOR_RESET
"Give a description. This partly defines the path of output file.\n"
ANSI_COLOR_BLUE "-h " ANSI_COLOR_RESET
"Print help and exit.\n"
ANSI_COLOR_BLUE "-o " ANSI_COLOR_RESET
"Specify output_path (to which description is appended).  Default: \
/space/ge52sir/non_linear_PS/output/.\n"
ANSI_COLOR_BLUE "-N " ANSI_COLOR_RESET
"Specify maximum number of CUBA Monte Carlo evaluations. Default: 1e6.\n"
ANSI_COLOR_BLUE "-r " ANSI_COLOR_RESET
"Specify relative tolerance for Monte Carlo integrator. Default: 1e-3.\n"
ANSI_COLOR_BLUE "-s " ANSI_COLOR_RESET
"Specify path for CUBA statefile. Default: \
/space/ge52sir/non_linear_PS/output/CUBA_statefiles/.\n"
ANSI_COLOR_BLUE "-v " ANSI_COLOR_RESET
"CUBA verbosity level (1, 2 or 3). Default: 1.\n\n"
"Example usage:\n\n\
$ cosPT -C 100 -d 2fluid_m_nu_0.07_IC1 -N 1e5 4\n\n\
Current compilation settings:\n\n");
    print_compilation_settings();
}



void print_compilation_settings() {
    printf("Loops                     = %d\n", LOOPS);
    printf("Components                = %d\n", COMPONENTS);
    printf("Time steps                = %d\n", TIME_STEPS);
    printf("Integration limits        = [%e,%e]\n", Q_MIN, Q_MAX);
    printf("Initial/final times       = [%e,%e]\n", ETA_I, ETA_F);
    printf("Neutrino mass             = %f\n", M_NU);
}
