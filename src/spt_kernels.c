/*
   spt_kernels.c

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_combination.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/spt_kernels.h"


vfloat partial_SPT_sum(
        const short int arguments[], /* kernel arguments                                  */
        short int component,         /* component to compute, NB: assumed to be 0-indexed */
        const gsl_matrix* alpha,     /* table of alpha function values for various input  */
        const gsl_matrix* beta,      /* table of beta function values for various input   */
        kernel_value* kernels,       /* kernel table                                      */
        const short int n,
        const short int m,
        const short int a,
        const short int b
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

    // DEBUG-mode: check that groupings are equal in size
#if DEBUG
    if (gsl_combination_n(comb_l) != gsl_combination_n(comb_r))
        warning("Left/right grouping of arguments does not have equal sizes.");
#endif

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
        /* printf("args_l = "); */
        /* for (int i = 0; i < N_KERNEL_ARGS; ++i) { */
        /*     printf("%d, ",args_l[i]); */
        /* } */
        /* printf("args_r = "); */
        /* for (int i = 0; i < N_KERNEL_ARGS; ++i) { */
        /*     printf("%d, ",args_r[i]); */
        /* } */
        /* printf("\t|\t"); */

        short int sum_l = sum_vectors(args_l,N_KERNEL_ARGS);
        short int sum_r = sum_vectors(args_r,N_KERNEL_ARGS);
        /* printf("sum_l  = %i, ", sum_l ); */
        /* printf("sum_r  = %i\n", sum_r ); */

        // F_n <-> component 0; G_n <-> component 1
        value += compute_SPT_kernel(args_l,1,alpha,beta,kernels) *
            (  a * gsl_matrix_get(alpha,sum_l,sum_r) * compute_SPT_kernel(args_r,0,alpha,beta,kernels)
             + b * gsl_matrix_get(beta ,sum_l,sum_r) * compute_SPT_kernel(args_r,1,alpha,beta,kernels)
            );

    } while (gsl_combination_next(comb_l) == GSL_SUCCESS &&
             gsl_combination_prev(comb_r) == GSL_SUCCESS
            );

    value /= gsl_combination_n(comb_l);

    gsl_combination_free(comb_l);
    gsl_combination_free(comb_r);

    return value;
}



vfloat compute_SPT_kernel(
        const short int arguments[], /* kernel arguments                                  */
        short int component,         /* component to compute, NB: assumed to be 0-indexed */
        const gsl_matrix* alpha,     /* table of alpha function values for various input  */
        const gsl_matrix* beta,      /* table of beta function values for various input   */
        kernel_value* kernels        /* kernel table                                      */
        )
{
    // Compute kernel index, this depends on arguments (argument_index) and
    // which component is to be computed
    short int argument_index = 0;
    short int n = 0;
    kernel_index_from_arguments(arguments,&argument_index,&n);
    short int index = combined_kernel_index(argument_index,component);

    // For SPT kernels, F_1 = G_1 = ... = 1
    if (n == 1) {
        return 1.0;
    }

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
        value += partial_SPT_sum(arguments,component,alpha,beta,kernels,n,m,a,b);
        debug_print("(n,m) = (%d,%d), \tvalue += %f\n",
            n,m,
            (partial_SPT_sum(arguments,component,alpha,beta,kernels,n,m,a,b)/((2*n + 3)*(n - 1))));
    }

    // Divide by overall factor in SPT recursion relation
    value /= (2*n + 3) * (n - 1);

    // Update kernel table
    kernels[index].value = value;
    kernels[index].computed = true;

    return value;
}
