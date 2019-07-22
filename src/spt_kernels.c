/*
   spt_kernels.c

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/


#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/tables.h"
#include "../include/spt_kernels.h"


vfloat partial_SPT_sum(
        const short int arguments[], /* kernel arguments                                  */
        short int n,                 /* kernel number                                     */
        short int m,                 /* sum index in kernel recursion relation            */
        short int a,                 /* coefficient: (2n+1) for F, 2 for G                */
        short int b,                 /* coefficient: 3 for F, 2n for G                    */
        const table_ptrs_t* tables
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

        short int sum_l = sum_vectors(args_l,N_KERNEL_ARGS,tables->sum_table);
        short int sum_r = sum_vectors(args_r,N_KERNEL_ARGS,tables->sum_table);

        // Send -1 as (unknown) kernel_index argument
        value += compute_SPT_kernel(args_l, -1, m, 1, tables) *
            (  a * matrix_get(tables->alpha, sum_l, sum_r)
               * compute_SPT_kernel(args_r, -1, n-m, 0, tables)
             + b * matrix_get(tables->beta, sum_l, sum_r)
               * compute_SPT_kernel(args_r, -1, n-m, 1, tables)
            );

        // When m != (n - m), we may additionally compute the (n-m)-term by
        // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
        // compute_SPT_kernel() only needs to sum up to (including) floor(n/2).
        if (m != n - m) {
            value += compute_SPT_kernel(args_r, -1, n-m, 1, tables) *
                (  a * matrix_get(tables->alpha, sum_r, sum_l)
                   * compute_SPT_kernel(args_l, -1, m, 0, tables)
                   + b * matrix_get(tables->beta, sum_r, sum_l)
                   * compute_SPT_kernel(args_l, -1, m, 1, tables)
                );
        }

    } while (gsl_combination_next(comb_l) == GSL_SUCCESS &&
             gsl_combination_prev(comb_r) == GSL_SUCCESS
            );

    // Devide through by symmetrization factor (n choose m)
    value /= gsl_sf_choose(n,m);

    gsl_combination_free(comb_l);
    gsl_combination_free(comb_r);

    return value;
}



vfloat compute_SPT_kernel(
        const short int arguments[], /* kernel arguments                                   */
        short int kernel_index,      /* index for kernel table (to be combined with comp.) */
        short int n,                 /* order in perturbation theory expansion            */
        short int component,         /* component to compute, NB: assumed to be 0-indexed */
        const table_ptrs_t* tables
        )
{
    // DEBUG: check that the number of non-zero arguments is in fact n, and
    // that kernel_index is in fact equivalent to arguments
#if DEBUG >= 1
    short int argument_index = kernel_index_from_arguments(arguments);
    if (kernel_index != -1 && argument_index != kernel_index) {
        warning_verbose("Index computed from kernel arguments (%d) does not "
                "equal kernel_index (%d).", argument_index, kernel_index);
    }

    int n_args = 0;
    for (int i = 0; i < N_KERNEL_ARGS; ++i) {
        if (arguments[i] != ZERO_LABEL) n_args++;
    }
    if (n_args != n)
        warning_verbose("Number of arguments is %d, while n is %d.", n_args,n);
#endif

    // For SPT kernels, F_1 = G_1 = ... = 1
    if (n == 1) {
        return 1.0;
    }

    // If kernel_index is not known, -1 is sent as argument
    if (kernel_index == -1) {
        kernel_index = kernel_index_from_arguments(arguments);
    }
    // Combine kernel_index with component
    short int index = combined_kernel_index(kernel_index,component);

    // Check if the kernel is already computed
    if (tables->kernels[index].computed) {
        return tables->kernels[index].value;
    }

    // Define some factors dependent on component to compute
    short int a,b;
    if (component == 0) {
        a = 2 * n + 1;
        b = 2;
    }
    else {
        a = 3;
        b = 2 * n;
    }

    vfloat value = 0.0;

    // Only sum up to (including) floor(n/2), since partial_SPT_sum()
    // simultaneously computes terms m and (n-m)
    for (int m = 1; m <= n/2; ++m) {
        value += partial_SPT_sum(arguments,n,m,a,b,tables);
    }

    // Divide by overall factor in SPT recursion relation
    value /= (2*n + 3) * (n - 1);

    // Update kernel table
    tables->kernels[index].value = value;
    tables->kernels[index].computed = true;

    return value;
}
