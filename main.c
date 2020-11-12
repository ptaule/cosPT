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
#include "include/spt_kernels.h"
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
    const char* input_ps_file        = CLASS_PATH "z10_pk_cb.dat";
    const char* zeta_file            = CLASS_PATH "zeta_of_etaD.dat";
    const char* redshift_file        = CLASS_PATH "redshift_of_etaD.dat";

#if SOUND_SPEED == CG2
    const char* effcs2_etaD_grid_file = "";
    const char* effcs2_k_grid_file    = "";
    const char* effcs2_file           = "";
    const char* omega_eigvals_file =
        "input/m_nu_" M_NU_STRING "/cg2/growing_mode_eigenvalues_etaD_ini.dat";
    const char* ic_F1_files[COMPONENTS];
    ic_F1_files[0] = "input/m_nu_" M_NU_STRING "/cg2/F1_growing_mode_etaD_-10.dat";
    ic_F1_files[1] = "input/m_nu_" M_NU_STRING "/cg2/F2_growing_mode_etaD_-10.dat";
    ic_F1_files[2] = "input/m_nu_" M_NU_STRING "/cg2/F3_growing_mode_etaD_-10.dat";
    ic_F1_files[3] = "input/m_nu_" M_NU_STRING "/cg2/F4_growing_mode_etaD_-10.dat";
#elif SOUND_SPEED == EFFCS2
    const char* effcs2_etaD_grid_file = EFFCS2_PATH "m_nu_" M_NU_STRING
        "/etaD_grid.dat";
    const char* effcs2_k_grid_file    = EFFCS2_PATH "k_grid.dat";
    const char* effcs2_file           = EFFCS2_PATH "m_nu_" M_NU_STRING
        "/effcs2_exact.dat";

    const char* omega_eigvals_file =
        "input/m_nu_" M_NU_STRING "/effcs2_exact/growing_mode_eigenvalues_etaD_ini.dat";
    const char* ic_F1_files[COMPONENTS];
    ic_F1_files[0] = "input/m_nu_" M_NU_STRING "/effcs2_exact/F1_growing_mode_etaD_-10.dat";
    ic_F1_files[1] = "input/m_nu_" M_NU_STRING "/effcs2_exact/F2_growing_mode_etaD_-10.dat";
    ic_F1_files[2] = "input/m_nu_" M_NU_STRING "/effcs2_exact/F3_growing_mode_etaD_-10.dat";
    ic_F1_files[3] = "input/m_nu_" M_NU_STRING "/effcs2_exact/F4_growing_mode_etaD_-10.dat";
#endif


    // Constants fixed by command line options (and default values)
    max_n_threads        = 100;
    double cuba_epsabs   = 1e-12;
    double cuba_epsrel   = 1e-4;
    double cuba_maxevals = 1e5;
    int cuba_verbose     = 0;

    char output_ps_file[] = "old_m_nu_0.07ev.dat";
    char* cuba_statefile = NULL;
    /* set_output_filepaths(output_ps_file, cuba_statefile, description, */
    /*         output_path, cuba_statefile_path, wavenumber_index); */

    printf("Compilation settings:\n");
    print_compilation_settings();
    printf("\nRuntime settings:\n");
    printf("Monte Carlo max evals = %.2e\n", cuba_maxevals);
    printf("Input power spectrum  = %s.\n", input_ps_file);
    /* printf("Output file           = %s.\n", output_ps_file); */
    /* printf("Cuba_statefile        = %s.\n", cuba_statefile); */

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

    evolution_params_t params = {
        .zeta_acc             = NULL,
        .zeta_spline          = NULL,
        .redshift_acc         = NULL,
        .redshift_spline      = NULL,
        .omega_eigvals_acc    = NULL,
        .omega_eigvals_spline = NULL,
        .ic_F1_accs           = {NULL},
        .ic_F1_splines        = {NULL},
        .effcs2_x_acc         = NULL,
        .effcs2_y_acc         = NULL,
        .effcs2_spline        = NULL
    };

    integration_input_t input = {
        .k           = 0.0,
        .q           = 0.0,
        .ps_acc      = NULL,
        .ps_spline   = NULL,
        .diagrams    = NULL,
        .params      = &params,
        .worker_mem  = worker_mem
    };

    // Read input files and interpolate
    read_and_interpolate(input_ps_file, &input.ps_acc, &input.ps_spline);
    read_and_interpolate(redshift_file, &params.redshift_acc,
            &params.redshift_spline);
    read_and_interpolate(zeta_file, &params.zeta_acc, &params.zeta_spline);
    read_and_interpolate(omega_eigvals_file, &params.omega_eigvals_acc,
            &params.omega_eigvals_spline);

#if SOUND_SPEED == EFFCS2
    read_and_interpolate_2d(effcs2_etaD_grid_file, effcs2_k_grid_file,
            effcs2_file, &params.effcs2_x_acc, &params.effcs2_y_acc,
            &params.effcs2_spline);
#endif

    for (int i = 0; i < COMPONENTS; ++i) {
        read_and_interpolate(ic_F1_files[i], &params.ic_F1_accs[i],
                &params.ic_F1_splines[i]);
    }

    vfloat k_min = 1e-3;
    vfloat k_max = 1e1;
#define n_points 25

    vfloat k = 0;
    vfloat q = 0;

    vfloat F2s[n_points][n_points] = {0};
    vfloat G2s[n_points][n_points] = {0};

    cubacores(4,10000);

    for (int i = 0; i < n_points; ++i) {
        printf("i = %i\n",i);
        k = k_min * pow(k_max/k_min, i/(vfloat)(n_points - 1));
        input.k = k;
        for (int j = 0; j < n_points; ++j) {
            q = k_min * pow(k_max/k_min, j/(vfloat)(n_points - 1));
            input.q = q;

            int nregions, neval, fail;
            cubareal result[INTEGRAND_COMPONENTS] = {0};
            cubareal error[INTEGRAND_COMPONENTS] = {0};
            cubareal prob[INTEGRAND_COMPONENTS] = {0};
    // CUBA settings
#define CUBA_NVEC 1
#define CUBA_LAST 4
#define CUBA_RETAIN_STATEFILE 0
#define CUBA_SEED 0
#define CUBA_MINEVAL 0
#define CUBA_SPIN NULL
#define CUBA_NNEW 1000
#define CUBA_NMIN 2
#define CUBA_FLATNESS 25.

#define N_DIMS 1
            Suave(N_DIMS, INTEGRAND_COMPONENTS, (integrand_t)cuba_integrand, &input,
                    CUBA_NVEC, cuba_epsrel, cuba_epsabs,
                    (cuba_verbose | CUBA_LAST | CUBA_RETAIN_STATEFILE), CUBA_SEED,
                    CUBA_MINEVAL, cuba_maxevals, CUBA_NNEW, CUBA_NMIN, CUBA_FLATNESS,
                    cuba_statefile, CUBA_SPIN, &nregions, &neval, &fail, result, error,
                    prob);

            vfloat overall_factor = 1;

            F2s[i][j] = overall_factor * result[0];
            G2s[i][j] = overall_factor * result[1];
        }
    }
    printf("F2:\n");
    printf("===================================================\n");
    for (int i = 0; i < n_points; ++i) {
        for (int j = 0; j < n_points; ++j) {
            printf("%Le\t", F2s[i][j]);
        }
        printf("\n");
    }
    printf("===================================================\n");
    printf("G2:\n");
    printf("===================================================\n");
    for (int i = 0; i < n_points; ++i) {
        for (int j = 0; j < n_points; ++j) {
            printf("%Le\t", G2s[i][j]);
        }
        printf("\n");
    }

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

    vars.magnitudes[0] = input->q;
    vars.cos_theta[0] = 2.0 * xx[0] - 1.0;

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

    short int args[] = {4, 2, ZERO_LABEL};
    short int kernel_index = kernel_evolution(args, -1, 2, input->params, tables);

    ff[0] = vars.cos_theta[0] * tables->kernels[kernel_index].values[TIME_STEPS - 1][0];
    ff[1] = vars.cos_theta[0] * tables->kernels[kernel_index].values[TIME_STEPS - 1][1];

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
