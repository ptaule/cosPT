/*
   evolve_kernels.c

   Created by Petter Taule on 04.04.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <gsl/gsl_combination.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/kernels.h"
#include "../include/evolve_kernels.h"


inline static void initialize_timesteps(vfloat eta[], vfloat eta_i, vfloat eta_f) {
    // Linear time step:
    vfloat d_eta = abs(eta_f - eta_i)/(TIME_STEPS - 1);
    for (int i = 0; i < TIME_STEPS; ++i) {
        eta[i] = eta_i + i*d_eta;
    }
}


static void vertex(
        short int m1,
        short int m2,
        const short int args_l[],
        const short int args_r[],
        short int sum_l,
        short int sum_r,
        vfloat** rhs_sum,
        const parameters_t* params,
        const table_pointers_t* data_tables
        )
{
    // Compute RHS expression depending on component

    vfloat alpha_lr = matrix_get(data_tables->alpha,sum_l,sum_r);
    vfloat alpha_rl = matrix_get(data_tables->alpha,sum_r,sum_l);
    vfloat beta     = matrix_get(data_tables->beta ,sum_r,sum_l);

    short int index_l = compute_time_dependent_kernels(args_l, m1, params,
            data_tables);
    short int index_r = compute_time_dependent_kernels(args_r, m2, params,
            data_tables);

    short int a,b,c;

    // Note: all components are zero-indexed
    for (int i = 0; i < TIME_STEPS; ++i) {
        // Component a = 0, two contributing terms:
        {
            a = 0, b = 0, c = 1;
            // gamma_001 = alpha_lr
            rhs_sum[a][i] += 0.5 * alpha_lr
                * data_tables->kernels[index_l].values[i][b]
                * data_tables->kernels[index_r].values[i][c];

            b = 1, c = 0;
            // gamma_010 = alpha_rl
            rhs_sum[0][i] += 0.5 * alpha_rl
                * data_tables->kernels[index_l].values[i][b]
                * data_tables->kernels[index_r].values[i][c];
        }
        // Component a = 1, one contributing term
        {
            a = 1, b = 1, c = 1;
            // gamma_111 = beta
            rhs_sum[a][i] += beta
                * data_tables->kernels[index_l].values[i][b]
                * data_tables->kernels[index_r].values[i][c];
        }
        // Component a = 2, two contributing terms
        {
            a = 2, b = 2, c = 3;
            // gamma_223 = alpha_lr
            rhs_sum[a][i] += 0.5 * alpha_lr
                * data_tables->kernels[index_l].values[i][b]
                * data_tables->kernels[index_r].values[i][c];

            b = 3, c = 2;
            // gamma_232 = alpha_rl
            rhs_sum[a][i] += 0.5 * alpha_rl
                * data_tables->kernels[index_l].values[i][b]
                * data_tables->kernels[index_r].values[i][c];
        }
        // Component a = 3, one contributing term
        {
            a = 2, b = 2, c = 3;
            // gamma_444 = beta
            rhs_sum[a][i] += beta
                * data_tables->kernels[index_l].values[i][b]
                * data_tables->kernels[index_r].values[i][c];
        }
    }
}



void compute_RHS_sum(
        const short int arguments[],
        short int n,
        const parameters_t* params,
        const table_pointers_t* data_tables,
        const vfloat eta[],
        gsl_spline* spline[],   /* out, interpolated RHS sum for each component     */
        gsl_interp_accel* acc[] /* out, gsl_interpolation accelerated lookup object */
        )
{
    // Allocate memory for rhs_sum
    vfloat** rhs_sum = (vfloat**)calloc(COMPONENTS, sizeof(vfloat*));
    for (int i = 0; i < COMPONENTS; ++i) {
        rhs_sum[i] = (vfloat*)calloc(TIME_STEPS, sizeof(vfloat));
    }

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

            vertex(m, n-m, args_l, args_r, sum_l, sum_r, rhs_sum, params,
                    data_tables);

            // When m != (n - m), we may additionally compute the (n-m)-term by
            // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
            // compute_SPT_kernel() only needs to sum up to (including) floor(n/2).
            if (m != n - m) {
                vertex(n-m, m, args_l, args_r, sum_l, sum_r, rhs_sum, params,
                        data_tables);
            }
        } while (gsl_combination_next(comb_l) == GSL_SUCCESS &&
                gsl_combination_prev(comb_r) == GSL_SUCCESS
                );

        // Devide through by symmetrization factor (n choose m)
        gsl_combination_free(comb_l);
        gsl_combination_free(comb_r);
    }

    for (int component = 0; component < COMPONENTS; ++component) {
        acc[component] = gsl_interp_accel_alloc();
        spline[component] = gsl_spline_alloc(INTERPOL_TYPE, TIME_STEPS);
        gsl_spline_init(spline[component], eta, rhs_sum[component], TIME_STEPS);

        // Free allocated memory
        free(rhs_sum[component]);
    }
    free(rhs_sum);
}



typedef struct {
    short int n;
    gsl_spline** splines;
    gsl_interp_accel** accs;
    const parameters_t* parameters;
} ode_input_t;



inline static void set_omega_matrix(vfloat eta, const parameters_t* params) {
    matrix_t* omega = params->omega;
    double hubble_factor = params->omega_m0 * exp(-3*eta) + (1 - params->omega_m0);
    double omega_M = params->omega_m0 * exp(-3*eta)/ hubble_factor;

    /* double k_FS = 0.908 * sqrt(params->omega_m0) * params->m_nu/3 * exp(eta/2); */

    if (COMPONENTS != 2) {
        warning_verbose("No implementation for COMPONENTS = %d (yet).",COMPONENTS);
    }

    // First row
    matrix_set(omega,0,0,  0);
    matrix_set(omega,0,1, -1);
    // Second row
    matrix_set(omega,1,0, -3/2.0 * omega_M/params->f2    );
    matrix_set(omega,1,1,  3/2.0 * omega_M/params->f2 - 1);

/*
    // First row
    matrix_set(omega,0,0,  0);
    matrix_set(omega,0,1, -1);
    matrix_set(omega,0,2,  0);
    matrix_set(omega,0,3,  0);
    // Second row
    matrix_set(omega,1,0, -3/2.0 * omega_M * (1 - params->f_nu));
    matrix_set(omega,1,1,  2 - 3/2.0 * omega_M                 );
    matrix_set(omega,1,2, -3/2.0 * omega_M * params->f_nu      );
    matrix_set(omega,1,3,  0                                   );
    // Third row
    matrix_set(omega,2,0,  0);
    matrix_set(omega,2,1,  0);
    matrix_set(omega,2,2,  0);
    matrix_set(omega,2,3, -1);
    // Fourth row
    matrix_set(omega,3,0, -3/2.0 * omega_M * (1 - params->f_nu)     );
    matrix_set(omega,3,1,  0                                        );
    matrix_set(omega,3,2, -3/2.0 * omega_M * (F_NU - pow(k/k_FS,2)) );
    matrix_set(omega,3,3,  2 - 3/2.0 * omega_M                      );
*/

}


int evolve_kernels(double eta, const double y[], double f[], void *ode_input) {
    ode_input_t input = *(ode_input_t*)ode_input;
    short int n = input.n;
    const parameters_t* params = input.parameters;

    set_omega_matrix(eta, params);

    matrix_t* omega = params->omega;

    gsl_vector_const_view y_vec = gsl_vector_const_view_array(y,COMPONENTS);
    gsl_vector_view f_vec = gsl_vector_view_array(f,COMPONENTS);

    for (int i = 0; i < COMPONENTS; ++i) {
        vfloat value = matrix_get(omega,i,i);
        matrix_set(omega,i,i,value + n);

        vfloat rhs = gsl_spline_eval(input.splines[i], eta, input.accs[i]);
        gsl_vector_set(&f_vec.vector, i, rhs);
    }

    gsl_blas_dgemv(CblasNoTrans, 1, omega, &y_vec.vector, 1, &f_vec.vector);
    return GSL_SUCCESS;
}



short int compute_time_dependent_kernels(
        const short int arguments[],
        short int n,
        const parameters_t* params,
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
    short int index = kernel_index_from_arguments(arguments);

    // Check if the kernel is already computed
    if (data_tables->kernels[index].evolved) return index;

    // GSL interpolation variables for interpolated RHS
    gsl_spline*       splines[COMPONENTS];
    gsl_interp_accel* accs   [COMPONENTS];

    // Initialize time steps in eta
    vfloat eta[TIME_STEPS];
    initialize_timesteps(eta, params->eta_i, params->eta_f);

    // Compute RHS sum in evolution equation
    compute_RHS_sum(arguments, n, params, data_tables, eta, splines, accs);

    // Set up ODE system
    ode_input_t input = {
        .n = n,
        .parameters = params,
        .splines = splines,
        .accs = accs
    };

    gsl_odeiv2_system sys = {evolve_kernels, NULL, COMPONENTS, &input};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(
            &sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 1e-6);

    vfloat eta_temp = params->eta_i;
    for (int i = 1; i < TIME_STEPS; i++)
    {
        vfloat* y = data_tables->kernels[index].values[i];
        int status = gsl_odeiv2_driver_apply(driver, &eta_temp, eta[i], y);

        if (status != GSL_SUCCESS) {
            warning_verbose("GLS ODE driver gave error value = %d.", status);
            break;
        }
    }

    gsl_odeiv2_driver_free(driver);

    data_tables->kernels[index].evolved = true;
    return index;
}
