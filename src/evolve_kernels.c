/*
   evolve_kernels.c

   Created by Petter Taule on 04.04.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <gsl/gsl_combination.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/kernels.h"
#include "../include/evolve_kernels.h"


inline static void initialize_timesteps(vfloat eta[], vfloat eta_i, vfloat eta_f) {
    // Linear time step:
    vfloat d_eta = abs(eta_f - eta_i)/TIME_STEPS;
    for (int i = 0; i < TIME_STEPS; ++i) {
        eta[i] = eta_i + i*d_eta;
    }
}


inline static void vertex(
        short int m1,
        short int m2,
        short int component,
        const short int args_l[],
        const short int args_r[],
        short int sum_l,
        short int sum_r,
        vfloat rhs_sum[],
        const parameters_t* params,
        const table_pointers_t* data_tables
        )
{
    // Compute RHS expression depending on component
    short int index_l, index_r, b, c;
    vfloat gamma_value = 0;
    switch (component) {
        case 0: {
                    // Note that component-variables (b,c) are zero indexed
                    b = 0;
                    c = 1;
                    // gamma_112 = alpha(k1,k2)
                    gamma_value = matrix_get(data_tables->alpha,sum_l,sum_r);
                    index_l = compute_time_dependent_kernel(args_l, m1, b,
                            params, data_tables);
                    index_r = compute_time_dependent_kernel( args_r, m2, c,
                            params, data_tables);

                    for (int i = 0; i < TIME_STEPS; ++i) {
                        rhs_sum[i] += gamma_value
                            * data_tables->kernels[index_l].values[i]
                            * data_tables->kernels[index_r].values[i];
                    }

                    b = 1;
                    c = 0;
                    // gamma_121 = alpha(k2,k1)
                    gamma_value = matrix_get(data_tables->alpha,sum_r,sum_l);
                    index_l = compute_time_dependent_kernel(args_l, m1, b,
                            params, data_tables);
                    index_r = compute_time_dependent_kernel(args_r, m2, c,
                            params, data_tables);

                    for (int i = 0; i < TIME_STEPS; ++i) {
                        rhs_sum[i] += gamma_value
                            * data_tables->kernels[index_l].values[i]
                            * data_tables->kernels[index_r].values[i];
                    }
                    break;
                }
        case 1: {
                    b = 1;
                    c = 1;
                    // gamma_222 = beta(k1,k2)
                    gamma_value = matrix_get(data_tables->beta,sum_l,sum_r);
                    index_l = compute_time_dependent_kernel(args_l, m1, b,
                            params,data_tables);
                    index_r = compute_time_dependent_kernel(args_r, m2, c,
                            params,data_tables);

                    for (int i = 0; i < TIME_STEPS; ++i) {
                        rhs_sum[i] += gamma_value
                            * data_tables->kernels[index_l].values[i]
                            * data_tables->kernels[index_r].values[i];
                    }
                    break;
                }
        case 2: {
                    // Note that component-variables (b,c) are zero indexed
                    b = 2;
                    c = 3;
                    // gamma_334 = alpha(k1,k2)
                    gamma_value = matrix_get(data_tables->alpha,sum_l,sum_r);
                    index_l = compute_time_dependent_kernel(args_l, m1, b,
                            params, data_tables);
                    index_r = compute_time_dependent_kernel(args_r, m2, c,
                            params, data_tables);

                    for (int i = 0; i < TIME_STEPS; ++i) {
                        rhs_sum[i] += gamma_value
                            * data_tables->kernels[index_l].values[i]
                            * data_tables->kernels[index_r].values[i];
                    }

                    b = 3;
                    c = 2;
                    // gamma_343 = alpha(k2,k1)
                    gamma_value = matrix_get(data_tables->alpha,sum_r,sum_l);
                    index_l = compute_time_dependent_kernel(args_l, m1, b,
                            params, data_tables);
                    index_r = compute_time_dependent_kernel(args_r, m2, c,
                            params, data_tables);

                    for (int i = 0; i < TIME_STEPS; ++i) {
                        rhs_sum[i] += gamma_value
                            * data_tables->kernels[index_l].values[i]
                            * data_tables->kernels[index_r].values[i];
                    }
                    break;
                }
        case 3: {
                    b = 3;
                    c = 3;
                    // gamma_444 = beta(k1,k2)
                    gamma_value = matrix_get(data_tables->beta,sum_l,sum_r);
                    index_l = compute_time_dependent_kernel(args_l, m1, b,
                            params,data_tables);
                    index_r = compute_time_dependent_kernel(args_r, m2, c,
                            params,data_tables);

                    for (int i = 0; i < TIME_STEPS; ++i) {
                        rhs_sum[i] += gamma_value
                            * data_tables->kernels[index_l].values[i]
                            * data_tables->kernels[index_r].values[i];
                    }
                    break;
                }
        default:
                warning_verbose("Cannot compute vertex for component a = %d.",component);
    }
}



void compute_RHS_sum(
        const short int arguments[],
        short int n,
        short int component,
        const parameters_t* params,
        const table_pointers_t* data_tables,
        const vfloat eta[],
        gsl_spline** spline, /* out, interpolated RHS sum */
        gsl_interp_accel** acc /* out, gsl_interpolation accelerated lookup object */
        )
{
    vfloat rhs_sum[TIME_STEPS];
    for (int i = 0; i < TIME_STEPS; ++i) rhs_sum[i] = 0.0;

    short int args_l[N_KERNEL_ARGS] = {0};
    short int args_r[N_KERNEL_ARGS] = {0};

    for (int m = 1; m <= n/2; ++m) {
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

            short int sum_l = sum_vectors(args_l,N_KERNEL_ARGS,data_tables->sum_table);
            short int sum_r = sum_vectors(args_r,N_KERNEL_ARGS,data_tables->sum_table);

            vertex(m, n-m, component, args_l, args_r, sum_l, sum_r, rhs_sum,
                    params, data_tables);

            // When m != (n - m), we may additionally compute the (n-m)-term by
            // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
            // compute_SPT_kernel() only needs to sum up to (including) floor(n/2).
            if (m != n - m) {
                vertex(n-m, m, component, args_l, args_r, sum_l, sum_r, rhs_sum,
                        params, data_tables);
            }
        } while (gsl_combination_next(comb_l) == GSL_SUCCESS &&
                gsl_combination_prev(comb_r) == GSL_SUCCESS
                );

        // Devide through by symmetrization factor (n choose m)
        gsl_combination_free(comb_l);
        gsl_combination_free(comb_r);
    }

    *acc = gsl_interp_accel_alloc();
    *spline = gsl_spline_alloc(INTERPOL_TYPE, TIME_STEPS);
    gsl_spline_init(*spline, eta, rhs_sum, TIME_STEPS);
}



typedef struct {
    short int n;
    short int component;
    const parameters_t* parameters;
} ode_input_t;



int evolve_kernel(double eta, const double y[], double f[], void *ode_input) {
    ode_input_t input = *(ode_input_t*)ode_input;
    short int n = input.n;
    short int component = input.component;
    const parameters_t* parameters = input.parameters;

    double hubble_factor = parameters->omega_m0 * exp(-3*eta) + (1 - parameters->omega_m0);
    double omega_M = parameters->omega_m0 * exp(-3*eta)/ hubble_factor;


}



short int compute_time_dependent_kernel(
        const short int arguments[],
        short int n,
        short int component,
        const parameters_t* parameters,
        const table_pointers_t* data_tables
        )
{
    // DEBUG: check that the number of non-zero arguments is in fact n
#if DEBUG >= 1
    int n_args = 0;
    for (int i = 0; i < N_KERNEL_ARGS; ++i) {
        if (arguments[i] != ZERO_LABEL) n_args++;
    }
    if (n_args != n)
        warning_verbose("Number of arguments is %d, while n is %d.", n_args,n);
#endif

    // Compute kernel index, this depends on arguments (argument_index), the
    // component to be computed and the time_step
    short int argument_index = kernel_index_from_arguments(arguments);
    short int index          = combined_kernel_index(argument_index,component);

    // const pointer alias to data_tables->kernels
    kernel_values_t* const kernels = data_tables->kernels;

    // Check if the kernel is already computed
    if (kernels[index].computed) return index;

    // GSL interpolation variables for interpolated RHS
    gsl_spline* spline;
    gsl_interp_accel* acc;
    // Initialize time steps in eta
    vfloat eta[TIME_STEPS];
    initialize_timesteps(eta, parameters->eta_i, parameters->eta_f);
    // Compute RHS sum in evolution equation
    compute_RHS_sum(arguments, n, component, parameters, data_tables, eta, &spline, &acc);

    ode_input_t input = {
        .n = n,
        .component = component,
        .parameters = parameters,
    };

    gsl_odeiv2_system sys = {evolve_kernel, NULL, COMPONENTS, &input};

    /* gsl_odeiv2_driver * driver = */
    /*     gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, */
    /*             1e-6, 1e-6, 1e-6); */

/*     for (int i = 1; i <= TIME_STEPS; i++) */
/*     { */
/*         double eta_i = eta_ini + i * abs(eta0 - eta_ini)/(double)N; */
/*         int status = gsl_odeiv2_driver_apply (driver, &eta, eta_i, y); */

/*         if (status != GSL_SUCCESS) */
/*         { */
/*             printf ("error, return value=%d\n", status); */
/*             break; */
/*         } */

/*         printf ("%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", eta, y[0], y[1], y[2], y[3]); */
/*     } */

    gsl_matrix_free (input.omega);
    gsl_odeiv2_driver_free (driver);

    // Update kernel table
    /* kernels[index].value = value; */

    kernels[index].computed = true;
    return index;
}
