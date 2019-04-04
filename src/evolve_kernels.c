/*
   evolve_kernels.c

   Created by Petter Taule on 04.04.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <gsl/gsl_combination.h>

#include "../include/constants.h"
#include "../include/kernels.h"
#include "../include/evolve_kernels.h"

vfloat partial_kernel(
        const short int n,
        const short int m,
        vfloat eta_i,
        vfloat eta_f,
        const table_pointers_t* data_tables
        )
{
    vfloat value = 0;

    short int args_l[N_KERNEL_ARGS] = {0};
    short int args_r[N_KERNEL_ARGS] = {0};

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

        // F_n <-> component 0; G_n <-> component 1


        // When m != (n - m), we may additionally compute the (n-m)-term by
        // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
        // compute_SPT_kernel() only needs to sum up to (including) floor(n/2).


    } while (gsl_combination_next(comb_l) == GSL_SUCCESS &&
             gsl_combination_prev(comb_r) == GSL_SUCCESS
            );

    // Devide through by symmetrization factor (n choose m)
    value /= gsl_sf_choose(n,m);

    gsl_combination_free(comb_l);
    gsl_combination_free(comb_r);

    return value;

}

vfloat evolve_kernel(
        const short int arguments[],
        short int n,
        short int component,
        vfloat eta_i,
        vfloat eta_f,
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
    short int index          = combined_kernel_index(argument_index,0,component);

    // const pointer alias to data_tables->kernels
    kernel_value_t* const kernels = data_tables->kernels;

    // Check if the kernel is already computed
    if (kernels[index].computed) return kernels[index].value;

    vfloat value = 0.0;

    // Only sum up to (including) floor(n/2), since partial_SPT_sum()
    // simultaneously computes terms m and (n-m)
    for (int m = 1; m <= n/2; ++m) {
        value += partial_kernel(arguments,n,m,data_tables);
    }

    // Evolution comes here?


    // Divide by overall factor in SPT recursion relation
    value /= (2*n + 3) * (n - 1);

    // Update kernel table
    kernels[index].value = value;
    kernels[index].computed = true;

    return value;



}
