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

#include "include/wavenumbers.h"

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

void draw_momenta(double* k, integration_variables_t* vars);



int main (int argc, char* argv[]) {
    srand(time(0));

    // Input files
    const char* input_zeta_file     = INPUT_PATH "zeta.dat";
    const char* input_redshift_file = INPUT_PATH "redshift.dat";

    const char* ic_perturbations_files[3];
    ic_perturbations_files[0] = INPUT_PATH "z10_theta_cb_over_aHf_delta_cb.dat";
    ic_perturbations_files[1] = INPUT_PATH "z10_delta_nu_over_delta_cb.dat";
    ic_perturbations_files[2] = INPUT_PATH "z10_theta_nu_over_aHf_delta_cb.dat";

    tables_t tables;

    // Initialize time steps in eta
    double eta[TIME_STEPS];
    initialize_timesteps(eta, ETA_I, ETA_F);
    tables.eta = eta;

    // Sum table can be computed right away
    short int sum_table[N_CONFIGS][N_CONFIGS];
    compute_sum_table(sum_table);
    tables.sum_table = (const short int (*)[])sum_table;

    tables_allocate(&tables);

    integration_variables_t vars = {
        .magnitudes = {0.1,0.2},
        .cos_theta = {-0.5,0.6},
        .phi = {2}
    };

#if LOOPS == 2
    double cos12 = cos(vars.phi[0]) *
        sin(acos(vars.cos_theta[0])) *
        sin(acos(vars.cos_theta[1]))
        + vars.cos_theta[0] * vars.cos_theta[1];

    printf("cos12 = %f\n", cos12);
#endif

    tables.Q_magnitudes = vars.magnitudes;

    diagram_t diagrams[N_DIAGRAMS];
    initialize_diagrams(diagrams);

    short int component_a = 0;
    short int component_b = 0;

    evolution_params_t params = {
        .zeta_acc        = NULL,
        .zeta_spline     = NULL,
        .redshift_acc    = NULL,
        .redshift_spline = NULL,
    };

    // Read input files and interpolate
    read_and_interpolate(input_redshift_file,&params.redshift_acc,&params.redshift_spline);
    read_and_interpolate(input_zeta_file,&params.zeta_acc,&params.zeta_spline);

    for (int i = 0; i < 3; ++i) {
        read_and_interpolate(ic_perturbations_files[i],
                &params.ic_perturb_accs[i], &params.ic_perturb_splines[i]);
    }

    double k = 0;

    for (int i = 0; i < NUM_WAVENUMBERS; ++i) {
        k = wavenumbers[i];

        tables_zero_initialize(&tables);

        // Initialize sum-, bare_scalar_products-, alpha- and beta-tables
        compute_bare_scalar_products(k, &vars, tables.bare_scalar_products);
        compute_scalar_products((const vfloat (*)[])tables.bare_scalar_products,
                tables.scalar_products);
        compute_alpha_beta_tables((const vfloat (*)[])tables.scalar_products,
                tables.alpha, tables.beta);

        short int m = 1;
        short int l = 2;
        short int r = 0;
        short int arguments_l[] = {13,5,3,7,1};
        short int arguments_r[] = {13,ZERO_LABEL,ZERO_LABEL,ZERO_LABEL,ZERO_LABEL};

        short int kernel_index_l = compute_SPT_kernels(arguments_l, -1, 2*l + m, &tables);
        short int kernel_index_r = compute_SPT_kernels(arguments_r, -1, 2*r + m, &tables);
        // Then, evolve kernels
        kernel_evolution(arguments_l, kernel_index_l, 2*l + m, &params, &tables);
        kernel_evolution(arguments_r, kernel_index_r, 2*r + m, &params, &tables);

        /* printf(ANSI_COLOR_MAGENTA "(m,l,r) = (%d,%d,%d)\n" ANSI_COLOR_RESET,m,l,r); */
        /* printf("%c%d", component_a + 'F', m + 2*l); */
        /* print_labels(arguments_l, N_KERNEL_ARGS); */
        printf("%f\t%Le\t%e",
                k,
                tables.kernels[kernel_index_l].spt_values[component_a],
                tables.kernels[kernel_index_l].values[TIME_STEPS - 1][component_a]);

        /* printf("%c%d", component_b + 'F', m + 2*r); */
        /* print_labels(arguments_r, N_KERNEL_ARGS); */
        printf("\t%Le\t%e\n",
                tables.kernels[kernel_index_r].spt_values[component_b],
                tables.kernels[kernel_index_r].values[TIME_STEPS - 1][component_b]);
    }

    diagrams_gc(diagrams);
    tables_gc(&tables);

    gsl_spline_free(params.redshift_spline);
    gsl_interp_accel_free(params.redshift_acc);
    gsl_spline_free(params.zeta_spline);
    gsl_interp_accel_free(params.zeta_acc);

    for (int i = 0; i < 3; ++i) {
        gsl_interp_accel_free(params.ic_perturb_accs[i]);
        gsl_spline_free(params.ic_perturb_splines[i]);
    }

    return 0;
}



void draw_momenta(double* k, integration_variables_t* vars) {
    double ratio = Q_MAX/(double)Q_MIN;
    double rand_f = (double)(rand()) / RAND_MAX;

    *k = Q_MIN * pow(ratio, rand_f);
    printf("values = {k -> %f", *k);

    for (int i = 0; i < LOOPS; ++i) {
        rand_f = (double)(rand()) / RAND_MAX;
        vars->magnitudes[i] = Q_MIN * pow(ratio,rand_f);
        vars->cos_theta[i] = 2 * (double)(rand())/RAND_MAX - 1;

        printf(", Q%d -> %Lf", i + 1, vars->magnitudes[i]);
        printf(", cosk%d -> %Lf", i + 1, vars->cos_theta[i]);
    }

#if LOOPS == 2
    vars->phi[0] = 2 * PI * (rand() / (double) RAND_MAX);

    double cos12 = cos(vars->phi[0]) *
        sin(acos(vars->cos_theta[0])) *
        sin(acos(vars->cos_theta[1]))
        + vars->cos_theta[0] * vars->cos_theta[1];

    printf(", cos12 -> %f", cos12);
#endif
    printf(" }\n");
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
    printf("Neutrino kernels (n>1) IC = %d\n", NEUTRINO_KERNEL_IC);
}
