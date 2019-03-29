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
#include "../include/spt_kernels.h"


/* This file contains compute_SPT_kernel(), which computes the SPT kernel at
 * order n. If component == 0, F kernel is computed, while G kernel is computed
 * if component == 1.
 *
 * To reduce the number of calls to kernel_index_from_arguments(), the
 * compute_SPT_kernel() function takes a kernel index as input (the kernel
 * indices of the first calls to the recursive SPT-computation can be
 * precomputed), hence it does not need to call kernel_index_from_arguments().
 * However, kernel_index_from_arguments() must be called in the recursive
 * computation, hence the recursion uses a second function,
 * compute_SPT_kernel_() to accomplish this.
 *
 * Since compute_SPT_kernel() and compute_SPT_kernel_() duplicate much of the
 * same code, these code instructions is moved to a helper function:
 * compute_SPT_kernel_helper().
 *
 * partial_SPT_sum() computes one term of the SPT sum.
 *
 * If DEBUG >= 1, various checks are performed to make sure that the parameters
 * kernel index and arguments are equivalent (they have a one-to-one
 * correspondance), and that n is in fact the number of arguments that is not
 * the zero-label. */



static vfloat compute_SPT_kernel_(
        const short int arguments[],
        short int n,
        short int component,
        const table_pointers_t* data_tables
        );
static inline vfloat compute_SPT_kernel_helper(
        short int combined_index,
        const short int arguments[],
        short int n,
        short int component,
        const table_pointers_t* data_tables
        );



static vfloat partial_SPT_sum(
        const short int arguments[], /* kernel arguments                                  */
        short int component,         /* component to compute, NB: assumed to be 0-indexed */
        short int n,                 /* kernel number                                     */
        short int m,                 /* sum index in kernel recursion relation            */
        short int a,                 /* coefficient: (2n+1) for F, 2 for G                */
        short int b,                 /* coefficient: 3 for F, 2n for G                    */
        const table_pointers_t* data_tables
        )
{
    vfloat value = 0;

    short int args_l[N_KERNEL_ARGS] = {};
    short int args_r[N_KERNEL_ARGS] = {};

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
        value += compute_SPT_kernel_(args_l,m,1,data_tables) *
            (  a * matrix_get(data_tables->alpha,sum_l,sum_r)
               * compute_SPT_kernel_(args_r,n-m,0,data_tables)
             + b * matrix_get(data_tables->beta ,sum_l,sum_r)
               * compute_SPT_kernel_(args_r,n-m,1,data_tables)
            );

    } while (gsl_combination_next(comb_l) == GSL_SUCCESS &&
             gsl_combination_prev(comb_r) == GSL_SUCCESS
            );

    // Devide through by symmetrization factor (n choose m)
    value /= gsl_sf_choose(n,m);

    gsl_combination_free(comb_l);
    gsl_combination_free(comb_r);

    return value;
}



static inline vfloat compute_SPT_kernel_helper(
        short int combined_index,    /* kernel table index                                 */
        const short int arguments[], /* kernel arguments                                   */
        short int n,                 /* order in perturbation theory expansion             */
        short int component,         /* component to compute, NB: assumed to be 0-indexed  */
        const table_pointers_t* data_tables
        )
{
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

    for (int m = 1; m < n; ++m) {
        value += partial_SPT_sum(arguments,component,n,m,a,b,data_tables);
    }

    // Divide by overall factor in SPT recursion relation
    value /= (2*n + 3) * (n - 1);

    // Update kernel table
    data_tables->kernels[combined_index].value = value;
    data_tables->kernels[combined_index].computed = true;

    return value;
}



static vfloat compute_SPT_kernel_(
        const short int arguments[], /* kernel arguments                                   */
        short int n,                 /* order in perturbation theory expansion             */
        short int component,         /* component to compute, NB: assumed to be 0-indexed  */
        const table_pointers_t* data_tables
        )
{
    // DEBUG: check that the number of non-zero arguments is in fact n, and
    // that kernel_index is in fact equivalent to arguments
#if DEBUG >= 1
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

    // Compute kernel index, this depends on arguments (argument_index) and
    // which component is to be computed
    short int argument_index = kernel_index_from_arguments(arguments);
    short int index          = combined_kernel_index(argument_index,component);

    // Check if the kernel is already computed
    if (data_tables->kernels[index].computed) {
        return data_tables->kernels[index].value;
    }

    return compute_SPT_kernel_helper(index,arguments,n,component,data_tables);
}



vfloat compute_SPT_kernel(
        short int kernel_index,      /* index for kernel table (to be combined with comp.) */
        const short int arguments[], /* kernel arguments                                   */
        short int n,                 /* order in perturbation theory expansion             */
        short int component,         /* component to compute, NB: assumed to be 0-indexed  */
        const table_pointers_t* data_tables
        )
{
    // DEBUG: check that the number of non-zero arguments is in fact n, and
    // that kernel_index is in fact equivalent to arguments
#if DEBUG >= 1
    short int argument_index = kernel_index_from_arguments(arguments);
    if (argument_index != kernel_index) {
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

    // Compute index in kernel-table, this depends on arguments (kernel_index)
    // and which component is to be computed
    short int index = combined_kernel_index(kernel_index,component);

    // Check if the kernel is already computed
    if (data_tables->kernels[index].computed) {
        return data_tables->kernels[index].value;
    }

    return compute_SPT_kernel_helper(index,arguments,n,component,data_tables);
}
