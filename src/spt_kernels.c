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


vfloat partial_SPT_sum(
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
        value += compute_SPT_kernel(args_l,m,1,data_tables) *
            (  a * matrix_get(data_tables->alpha,sum_l,sum_r)
               * compute_SPT_kernel(args_r,n-m,0,data_tables)
             + b * matrix_get(data_tables->beta ,sum_l,sum_r)
               * compute_SPT_kernel(args_r,n-m,1,data_tables)
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



// For debuggin purposes
__attribute__((unused))
static void print_kernel_info(
        const short int arguments[],
        short int n,
        short int component,
        vfloat value
        )
{
    printf("F%d(k",n);
    short int config[N_COEFFS];
    label2config(arguments[0],config,N_COEFFS);
    for (int i = 0; i < LOOPS; ++i) {
        if (config[i] == 0) continue;
        else if (config[i] == -1) printf("-Q%d",i+1);
        else if (config[i] == 1)  printf("+Q%d",i+1);
    }
    printf(", ");

    for (int i = 1; i < N_KERNEL_ARGS; ++i) {
        if (arguments[i] == ZERO_LABEL) break;
        label2config(arguments[i],config,N_COEFFS);
        for (int j = 0; j < LOOPS; ++j) {
            if (config[j] == 0) continue;
            else if (config[j] == -1) printf("-Q%d, ",j+1);
            else if (config[j] == 1)  printf("+Q%d, ",j+1);
        }
    }
    printf(") (comp. %d) \t\t\t= " vfloat_fmt "\n",component, value);
}



vfloat compute_SPT_kernel(
        const short int arguments[], /* kernel arguments                                  */
        short int n,                 /* order in perturbation theory expansion            */
        short int component,         /* component to compute, NB: assumed to be 0-indexed */
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

    // For SPT kernels, F_1 = G_1 = ... = 1
    if (n == 1) {
        return 1.0;
    }

    // Compute kernel index, this depends on arguments (argument_index) and
    // which component is to be computed
    short int argument_index = kernel_index_from_arguments(arguments);
    short int index          = combined_kernel_index(argument_index,component);

    // const pointer alias to data_tables->kernels
    kernel_value_t* const kernels = data_tables->kernels;

    // Check if the kernel is already computed
    if (kernels[index].computed) return kernels[index].value;

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
    kernels[index].value = value;
    kernels[index].computed = true;

    return value;
}
