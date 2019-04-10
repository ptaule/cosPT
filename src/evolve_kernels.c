/*
   evolve_kernels.c

   Created by Petter Taule on 04.04.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <string.h>

#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>
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
    vfloat d_eta = fabs(eta_f - eta_i)/(TIME_STEPS - 1);
    for (int i = 0; i < TIME_STEPS; ++i) {
        eta[i] = eta_i + i*d_eta;
    }
}


static void vertex(
        short int m_l,
        short int m_r,
        const short int args_l[],
        const short int args_r[],
        short int sum_l,
        short int sum_r,
        vfloat** rhs_sum,
        matrix_t* omega,
        const parameters_t* params,
        const table_pointers_t* data_tables
        )
{
    vfloat alpha_lr = matrix_get(data_tables->alpha,sum_l,sum_r);
    vfloat alpha_rl = matrix_get(data_tables->alpha,sum_r,sum_l);
    vfloat beta     = matrix_get(data_tables->beta ,sum_l,sum_r);

    short int index_l = kernel_evolution(args_l, -1, m_l, omega, params, data_tables);
    short int index_r = kernel_evolution(args_r, -1, m_r, omega, params, data_tables);

    short int a,b,c;

    // Note: all components are zero-indexed
    for (int i = 0; i < TIME_STEPS; ++i) {

        switch (COMPONENTS) {
            case 4:
                // Component a = 2, two contributing terms
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

                // Component a = 3, one contributing term
                a = 3, b = 3, c = 3;
                // gamma_444 = beta
                rhs_sum[a][i] += beta
                    * data_tables->kernels[index_l].values[i][b]
                    * data_tables->kernels[index_r].values[i][c];

            // Use switch fallthrough; for case 4 also code in case 2 should be
            // executed
            __attribute__((fallthrough));
            case 2:
                // Component a = 0, two contributing terms:
                a = 0, b = 0, c = 1;
                // gamma_001 = alpha_lr
                rhs_sum[a][i] += 0.5 * alpha_lr
                    * data_tables->kernels[index_l].values[i][b]
                    * data_tables->kernels[index_r].values[i][c];

                b = 1, c = 0;
                // gamma_010 = alpha_rl
                rhs_sum[a][i] += 0.5 * alpha_rl
                    * data_tables->kernels[index_l].values[i][b]
                    * data_tables->kernels[index_r].values[i][c];

                // Component a = 1, one contributing term
                a = 1, b = 1, c = 1;
                // gamma_111 = beta
                rhs_sum[a][i] += beta
                    * data_tables->kernels[index_l].values[i][b]
                    * data_tables->kernels[index_r].values[i][c];

                break;
            default:
                warning_verbose("No vertex implemented for COMPONENTS = %d.",COMPONENTS)
        }
    }
}



void compute_RHS_sum(
        const short int arguments[],
        short int n,
        matrix_t* omega,
        const parameters_t* params,
        const table_pointers_t* data_tables,
        const vfloat eta[],
        gsl_spline* splines[],   /* out, interpolated RHS sum for each component      */
        gsl_interp_accel* accs[] /* out, gsl_interpolation accelerated lookup objects */
        )
{
    // Allocate memory for rhs_sum
    vfloat** const rhs_sum = (vfloat**)calloc(COMPONENTS, sizeof(vfloat*));
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

            vertex(m, n-m, args_l, args_r, sum_l, sum_r, rhs_sum, omega, params,
                    data_tables);

            // When m != (n - m), we may additionally compute the (n-m)-term by
            // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
            // compute_RHS_sum() only needs to sum up to (including) floor(n/2).
            if (m != n - m) {
                vertex(n-m, m, args_r, args_l, sum_r, sum_l, rhs_sum, omega, params,
                        data_tables);
            }
        } while (gsl_combination_next(comb_l) == GSL_SUCCESS &&
                gsl_combination_prev(comb_r) == GSL_SUCCESS
                );

        // Devide through by symmetrization factor (n choose m)
        int n_choose_m = gsl_sf_choose(n,m);
        for (int i = 0; i < COMPONENTS; ++i) {
            for (int j = 0; j < TIME_STEPS; ++j) {
                rhs_sum[i][j] /= n_choose_m;
            }
        }

        gsl_combination_free(comb_l);
        gsl_combination_free(comb_r);
    }

    // Interpolate rhs_sum's
    for (int component = 0; component < COMPONENTS; ++component) {
        accs[component] = gsl_interp_accel_alloc();
        splines[component] = gsl_spline_alloc(INTERPOL_TYPE, TIME_STEPS);
        gsl_spline_init(splines[component], eta, rhs_sum[component], TIME_STEPS);

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
    matrix_t* omega;
} ode_input_t;



inline static void set_omega_matrix(matrix_t* omega, vfloat eta, const parameters_t* params) {
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

    matrix_t* omega = input.omega;

    set_omega_matrix(omega, eta, params);

    // Scale omega by -1 (i.e. move it to RHS of ODE)
    gsl_matrix_scale(omega,-1);

    gsl_vector_const_view y_vec = gsl_vector_const_view_array(y,COMPONENTS);
    gsl_vector_view f_vec = gsl_vector_view_array(f,COMPONENTS);

    for (int i = 0; i < COMPONENTS; ++i) {
        // Subtract n*I from omega
        vfloat value = matrix_get(omega,i,i);
        matrix_set(omega,i,i,value - n);

        vfloat rhs = gsl_spline_eval(input.splines[i], eta, input.accs[i]);
        gsl_vector_set(&f_vec.vector, i, rhs);
    }

    // Compute omega*y + f and store in f
    gsl_blas_dgemv(CblasNoTrans, 1, omega, &y_vec.vector, 1, &f_vec.vector);
    return GSL_SUCCESS;
}



short int kernel_evolution(
        const short int arguments[],
        short int index,             /* index, if known, else -1 */
        short int n,
        matrix_t* omega,
        const parameters_t* params,
        const table_pointers_t* data_tables
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

    // DEBUG: check that if an index is provided (index != -1), it is in fact
    // equivalent to arguments
    if (index != -1) {
        short int argument_index = kernel_index_from_arguments(arguments);
        if (argument_index != index) {
            warning_verbose("Index computed from kernel arguments (%d) does not "
                    "equal index argument (%d).", argument_index, index);
        }
    }
#endif

    if (index == -1) {
        index = kernel_index_from_arguments(arguments);
    }

    // Check if the kernel is already computed
    if (data_tables->kernels[index].evolved) return index;

    // GSL interpolation variables for interpolated RHS
    gsl_spline*       splines[COMPONENTS];
    gsl_interp_accel* accs   [COMPONENTS];

    // Initialize time steps in eta
    vfloat eta[TIME_STEPS];
    initialize_timesteps(eta, params->eta_i, params->eta_f);

    // Compute RHS sum in evolution equation
    compute_RHS_sum(arguments, n, omega, params, data_tables, eta, splines, accs);

    // Set up ODE system
    ode_input_t input = {
        .n = n,
        .parameters = params,
        .splines = splines,
        .accs = accs,
        .omega = omega
    };

    gsl_odeiv2_system sys = {evolve_kernels, NULL, COMPONENTS, &input};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(
            &sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 1e-6);

    vfloat eta_temp = params->eta_i;
    for (int i = 1; i < TIME_STEPS; i++)
    {
        // y is a pointer to the kernel table at time i
        vfloat* y = data_tables->kernels[index].values[i];
        // First copy previous values (i-1) to y
        memcpy(y, data_tables->kernels[index].values[i-1], COMPONENTS * sizeof(vfloat));
        // Then evolve y to time i
        int status = gsl_odeiv2_driver_apply(driver, &eta_temp, eta[i], y);

        if (status != GSL_SUCCESS) {
            warning_verbose("GLS ODE driver gave error value = %d.", status);
            break;
        }
    }


    // Free GSL ODE driver
    gsl_odeiv2_driver_free(driver);
    // Free GSL interpolation objects
    for (int i = 0; i < COMPONENTS; ++i) {
        gsl_spline_free(splines[i]);
        gsl_interp_accel_free(accs[i]);
    }


    data_tables->kernels[index].evolved = true;
    return index;
}
