/*
   evolve_kernels.c

   Created by Petter Taule on 04.04.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <string.h>
#include <math.h>

#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/tables.h"
#include "../include/evolve_kernels.h"


static void vertex(
        short int m_l,
        short int m_r,
        const short int args_l[],
        const short int args_r[],
        short int sum_l,
        short int sum_r,
        double (*partial_rhs_sum)[TIME_STEPS],
        const evolution_params_t* params,
        tables_t* tables
        )
{
    vfloat alpha_lr = tables->alpha[sum_l][sum_r];
    /* vfloat alpha_rl = tables->alpha[sum_r][sum_l]; */
    vfloat beta = tables->beta[sum_l][sum_r];

    short int index_l = kernel_evolution(args_l, -1, m_l, params, tables);
    short int index_r = kernel_evolution(args_r, -1, m_r, params, tables);

    short int a,b,c;

    // Note: all components are zero-indexed
    for (int i = 0; i < TIME_STEPS; ++i) {

        switch (COMPONENTS) {
            case 4:
                // Component a = 2
                a = 2, b = 3, c = 2;
                // gamma_223 = alpha_lr
                partial_rhs_sum[a][i] += alpha_lr
                    * tables->kernels[index_l].values[i][b]
                    * tables->kernels[index_r].values[i][c];

                /* The term below is redundant; due to momentum symmetrization
                 * of the evolution eq. RHS, it will be covered by the term
                 * above. */
                /* b = 2, c = 3; */
                /* // gamma_232 = alpha_rl */
                /* partial_rhs_sum[a][i] += 0.5 * alpha_rl */
                /*     * tables->kernels[index_l].values[i][b] */
                /*     * tables->kernels[index_r].values[i][c]; */

                // Component a = 3
                a = 3, b = 3, c = 3;
                // gamma_444 = beta
                partial_rhs_sum[a][i] += beta
                    * tables->kernels[index_l].values[i][b]
                    * tables->kernels[index_r].values[i][c];

            // Use switch fallthrough; for case 4 also code in case 2 should be
            // executed
            __attribute__((fallthrough));
            case 2:
                // Component a = 0
                a = 0, b = 1, c = 0;
                // gamma_001 = alpha_lr
                partial_rhs_sum[a][i] += alpha_lr
                    * tables->kernels[index_l].values[i][b]
                    * tables->kernels[index_r].values[i][c];

                /* The term below is redundant; due to momentum symmetrization
                 * of the evolution eq. RHS, it will be covered by the term
                 * above. */
                /* b = 0, c = 1; */
                /* // gamma_010 = alpha_rl */
                /* partial_rhs_sum[a][i] += 0.5 * alpha_rl */
                /*     * tables->kernels[index_l].values[i][b] */
                /*     * tables->kernels[index_r].values[i][c]; */

                // Component a = 1
                a = 1, b = 1, c = 1;
                // gamma_111 = beta
                partial_rhs_sum[a][i] += beta
                    * tables->kernels[index_l].values[i][b]
                    * tables->kernels[index_r].values[i][c];

                break;
            default:
                warning_verbose("No vertex implemented for COMPONENTS = %d.",COMPONENTS);
        }
    }
}



void compute_RHS_sum(
        const short int arguments[],
        short int n,
        const evolution_params_t* params,
        tables_t* tables,
        gsl_spline* rhs_splines[],   /* out, interpolated RHS sum for each component      */
        gsl_interp_accel* rhs_accs[] /* out, gsl_interpolation accelerated lookup objects */
        )
{
    double rhs_sum[COMPONENTS][TIME_STEPS] = {0};
    double partial_rhs_sum[COMPONENTS][TIME_STEPS] = {0};

    short int args_l[N_KERNEL_ARGS] = {0};
    short int args_r[N_KERNEL_ARGS] = {0};

    for (int m = 1; m <= n/2; ++m) {
        // Initialize partial_rhs_sum to 0
        for (int i = 0; i < COMPONENTS; ++i) {
            for (int j = 0; j < TIME_STEPS; ++j) {
                partial_rhs_sum[i][j] = 0;
            }
        }

        // Initialize args_l and args_r
        for (int i = 0; i < N_KERNEL_ARGS; ++i) {
            args_l[i] = ZERO_LABEL;
            args_r[i] = ZERO_LABEL;
        }

        // - comb_l starts at {0,1,...,m} and in the while-loop goes over all
        //   combinations of m elements from {0,...,n} (n choose m possibilities)
        // - comb_r starts at {m+1,...,n} and in the while-loop goes
        //   ("backwards") over all combinations of (n-m) elements from {0,...,n}
        //   (n choose (n-m) possibilities)

        gsl_combination* comb_l = gsl_combination_alloc(n,m);
        gsl_combination* comb_r = gsl_combination_alloc(n,n-m);

        gsl_combination_init_first(comb_l);
        gsl_combination_init_last(comb_r);

        do {
            // Use comb_l and comb_r to find argument combination
            for (int i = 0; i < m; ++i) {
                args_l[i] = arguments[gsl_combination_get(comb_l,i)];
            }
            for (int i = 0; i < n - m; ++i) {
                args_r[i] = arguments[gsl_combination_get(comb_r,i)];
            }

            short int sum_l = sum_vectors(args_l,N_KERNEL_ARGS,tables->sum_table);
            short int sum_r = sum_vectors(args_r,N_KERNEL_ARGS,tables->sum_table);

            vertex(m, n-m, args_l, args_r, sum_l, sum_r, partial_rhs_sum,
                    params, tables);

            // When m != (n - m), we may additionally compute the (n-m)-term by
            // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
            // compute_RHS_sum() only needs to sum up to (including) floor(n/2).
            if (m != n - m) {
                vertex(n-m, m, args_r, args_l, sum_r, sum_l, partial_rhs_sum,
                        params, tables);
            }
        } while (gsl_combination_next(comb_l) == GSL_SUCCESS &&
                gsl_combination_prev(comb_r) == GSL_SUCCESS
                );

        // Devide through by symmetrization factor (n choose m)
        int n_choose_m = gsl_sf_choose(n,m);
        for (int i = 0; i < COMPONENTS; ++i) {
            for (int j = 0; j < TIME_STEPS; ++j) {
                rhs_sum[i][j] += partial_rhs_sum[i][j] / n_choose_m;
            }
        }

        gsl_combination_free(comb_l);
        gsl_combination_free(comb_r);
    }

    // Interpolate rhs_sum's
    for (int i = 0; i < COMPONENTS; ++i) {
        rhs_accs[i] = gsl_interp_accel_alloc();
        rhs_splines[i] = gsl_spline_alloc(INTERPOL_TYPE, TIME_STEPS);
        gsl_spline_init(rhs_splines[i], tables->eta, rhs_sum[i],
                TIME_STEPS);
    }
}



typedef struct {
    short int n;
    const double k;
    gsl_spline** rhs_splines;
    gsl_interp_accel** rhs_accs;
    const evolution_params_t* parameters;
} ode_input_t;



int kernel_gradient(double eta, const double y[], double f[], void *ode_input) {
    ode_input_t input = *(ode_input_t*)ode_input;
    short int n = input.n;
    double k = input.k;
    const evolution_params_t* params = input.parameters;

    // Between eta_asymp and eta_i, set eta = eta_i
    if (eta < ETA_I) {
        eta = ETA_I;
    }

    double zeta = gsl_spline_eval(params->zeta_spline, eta, params->zeta_acc);

#define FS_factor 0.907778 * M_NU * SQRT_OMEGA_M

    double k_FS2 = FS_factor * FS_factor
        * pow(1 + gsl_spline_eval(params->redshift_spline, eta, params->redshift_acc), -1);
    /*
    // With eta parametrization:
    double k_FS2 = FS_factor * FS_factor * exp(eta);
    */

    double rhs[COMPONENTS] = {0};
    // If n == 1, rhs = 0
    if (n > 1) {
        for (int i = 0; i < COMPONENTS; ++i) {
            rhs[i] = gsl_spline_eval(input.rhs_splines[i], eta, input.rhs_accs[i]);
        }
    }

    // etaD parametrization with zeta(etaD) from CLASS
    f[0] = rhs[0] - n * y[0] + y[1];
    f[1] = rhs[1] + 1.5 * zeta * (1 - F_NU) * y[0] + (-1.5 * zeta + 1 - n) * y[1] + 1.5 * zeta * F_NU * y[2];
    f[2] = rhs[2] - n * y[2] + y[3];
    f[3] = rhs[3] + 1.5 * zeta * (1 - F_NU) * y[0] + 1.5 * zeta * (F_NU - k*k/k_FS2) * y[2] +
        (-1.5 * zeta + 1 - n) * y[3];
    /*
    // eta parametrization with omegaM(eta) from CLASS
    f[0] = rhs[0] - n * y[0] + y[1];
    f[1] = rhs[1] + 1.5 * omegaM * (1 - F_NU) * y[0] + (-2 + 1.5 * omegaM - n) * y[1] +
        1.5 * omegaM *  F_NU * y[2];
    f[2] = rhs[2] - n * y[2] + y[3];
    f[3] = rhs[3] + 1.5 * omegaM * (1 - F_NU) * y[0] + 1.5 * omegaM * (F_NU - k*k/k_FS2) * y[2] +
        (-2 + 1.5 * omegaM - n) * y[3];
    */

#undef FS_factor

    return GSL_SUCCESS;
}



void evolve_kernels(
        ode_input_t* input,
        const double* eta,
        double** kernels  /* TIME_STEPS*COMPONENTS table of kernels */
        )
{
    gsl_odeiv2_system sys = {kernel_gradient, NULL, COMPONENTS, input};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys,
            gsl_odeiv2_step_rkf45, ODE_HSTART, ODE_RTOL, ODE_ATOL);

    double eta_current = eta[0];

    // Evolve system
    for (int i = 1; i < TIME_STEPS; i++) {
        // y is a pointer to the kernel table at time i
        double* y = kernels[i];
        // First copy previous values (i-1) to y
        memcpy(y, kernels[i-1], COMPONENTS * sizeof(double));
        // Then evolve y to time index i
        int status = gsl_odeiv2_driver_apply(driver, &eta_current, eta[i], y);

        if (status != GSL_SUCCESS) {
            warning_verbose("GLS ODE driver gave error value = %d.", status);
            break;
        }
    }

    // Free GSL ODE driver
    gsl_odeiv2_driver_free(driver);
}



static void kernel_initial_conditions(
        short int kernel_index,
        short int n,
        double k,
        const evolution_params_t* params,
        tables_t* tables
        )
{
    // Use growing mode IC for F1 at eta_asymp
    if (n == 1) {
        for (int i = 0; i < COMPONENTS; ++i) {
            tables->kernels[kernel_index].values[0][i] =
                gsl_spline_eval(params->ic_F1_splines[i], k,
                        params->ic_F1_accs[i]);
        }
    }
    else {
        // Set F(n>1)-kernels to zero at eta_asymp
        for (int i = 0; i < COMPONENTS; ++i) {
            tables->kernels[kernel_index].values[0][i] = 0;
        }
    }
}



short int kernel_evolution(
        const short int arguments[],
        short int kernel_index,             /* -1 indicates not known */
        short int n,
        const evolution_params_t* params,
        tables_t* tables
        )
{
#if DEBUG >= 1
    // DEBUG: check that the number of non-zero arguments is in fact n
    int n_args = 0;
    for (int i = 0; i < N_KERNEL_ARGS; ++i) {
        if (arguments[i] != ZERO_LABEL) n_args++;
    }
    if (n_args != n)
        warning_verbose("Number of arguments is %d, while n is %d.", n_args,n);

    // DEBUG: check that if an kernel_index is provided (kernel_index != -1),
    // it is in fact equivalent to arguments
    if (kernel_index != -1) {
        short int argument_index = kernel_index_from_arguments(arguments);
        if (argument_index != kernel_index) {
            warning_verbose("Index computed from kernel arguments (%d) does not "
                    "equal kernel_index argument (%d).", argument_index,
                    kernel_index);
        }
    }
#endif

    if (kernel_index == -1) {
        kernel_index = kernel_index_from_arguments(arguments);
    }

    // If the kernel is already computed, return kernel_index
    if (tables->kernels[kernel_index].evolved) return kernel_index;

    // Compute k (sum of kernel arguments)
    short int sum = sum_vectors(arguments, N_KERNEL_ARGS, tables->sum_table);
    double k = sqrt(tables->scalar_products[sum][sum]);

    // Set initial conditions
    kernel_initial_conditions(kernel_index, n, k, params, tables);

    // GSL interpolation variables for interpolated RHS
    gsl_spline*       rhs_splines[COMPONENTS] = {NULL};
    gsl_interp_accel* rhs_accs   [COMPONENTS] = {NULL};

    // Compute RHS sum in evolution equation if n > 1. If n == 1, the RHS
    // equals 0, which is implemented in evolve_kernels().
    if (n > 1) {
        compute_RHS_sum(arguments, n, params, tables, rhs_splines, rhs_accs);
    }

    // Set up ODE input and system
    ode_input_t input = {
        .n = n,
        .k = k,
        .parameters = params,
        .rhs_splines = rhs_splines,
        .rhs_accs = rhs_accs,
    };

    evolve_kernels(&input, tables->eta, tables->kernels[kernel_index].values);

    // Free GSL interpolation objects
    for (int i = 0; i < COMPONENTS; ++i) {
        gsl_spline_free(rhs_splines[i]);
        gsl_interp_accel_free(rhs_accs[i]);
    }

    tables->kernels[kernel_index].evolved = true;
    return kernel_index;
}



void compute_F1(
        double k,
        const evolution_params_t* params,
        const double* eta,
        double* F1_eta_i, /* out */
        double* F1_eta_f  /* out */
        )
{
    double** values = (double**)calloc(TIME_STEPS, sizeof(double*));
    for (int i = 0; i < TIME_STEPS; ++i) {
        values[i] = (double*)calloc(COMPONENTS, sizeof(double));
    }

    for (int i = 0; i < COMPONENTS; ++i) {
        values[0][i] = gsl_spline_eval(params->ic_F1_splines[i], k,
                params->ic_F1_accs[i]);
    }

    // Set up ODE input and system
    ode_input_t input = {
        .n = 1,
        .k = k,
        .parameters = params,
        .rhs_splines = NULL,
        .rhs_accs = NULL,
    };

    evolve_kernels(&input, eta, values);

    for (int i = 0; i < COMPONENTS; ++i) {
        F1_eta_i[i] = values[PRE_TIME_STEPS][i];
        F1_eta_f[i] = values[TIME_STEPS - 1][i];
    }

    for (int i = 0; i < TIME_STEPS; ++i) {
        free(values[i]);
    }
    free(values);
}
