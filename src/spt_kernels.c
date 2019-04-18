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



static inline vfloat SPT_term(
        short int m_l,
        short int m_r,
        short int component,
        short int time_step,
        short int args_l[],
        short int args_r[],
        short int sum_l,
        short int sum_r,
        const table_pointers_t* data_tables
        )
{
    short int n = m_l + m_r;

    short int index_l = compute_SPT_kernels(args_l,m_l,time_step,data_tables);
    short int index_r = compute_SPT_kernels(args_r,m_r,time_step,data_tables);

    short int a,b;
    switch (component) {
        case 0:
            a = 2 * n + 1;
            b = 2;
            break;
        case 1:
            a = 3;
            b = 2 * n;
            break;
        default:
            warning("SPT_term() does not accept argument 'component' which does "
                    "not equal 0 or 1.")
    }

    return data_tables->kernels[index_l].spt_values[1] *
        (    a * matrix_get(data_tables->alpha,sum_l,sum_r)
           * data_tables->kernels[index_r].spt_values[0]
           + b * matrix_get(data_tables->beta, sum_l,sum_r)
           * data_tables->kernels[index_r].spt_values[1]
        );
}



static void partial_SPT_sum(
        const short int arguments[], /* kernel arguments                                  */
        short int n,                 /* kernel number                                     */
        short int m,                 /* sum index in kernel recursion relation            */
        short int kernel_index,
        short int time_step,
        const table_pointers_t* data_tables
        )
{
    vfloat kernel_values[COMPONENTS] = {0};

    short int args_l[N_KERNEL_ARGS] = {0};
    short int args_r[N_KERNEL_ARGS] = {0};

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

        short int sum_l = sum_vectors(args_l,N_KERNEL_ARGS,data_tables->sum_table);
        short int sum_r = sum_vectors(args_r,N_KERNEL_ARGS,data_tables->sum_table);

        for (int i = 0; i < COMPONENTS; ++i) {
            kernel_values[i] += SPT_term(m,n-m,i,time_step,args_l,args_r,sum_l,sum_r,data_tables);
        }

        // When m != (n - m), we may additionally compute the (n-m)-term by
        // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
        // compute_SPT_kernels() only needs to sum up to (including) floor(n/2).
        if (m != n - m) {
            for (int i = 0; i < COMPONENTS; ++i) {
                kernel_values[i] +=
                    SPT_term(n-m,m,i,time_step,args_r,args_l,sum_r,sum_l,data_tables);
            }
        }
    } while (gsl_combination_next(comb_l) == GSL_SUCCESS &&
             gsl_combination_prev(comb_r) == GSL_SUCCESS
            );

    // Devide through by symmetrization factor (n choose m)
    for (int i = 0; i < COMPONENTS; ++i) {
        kernel_values[i] /= gsl_sf_choose(n,m);
    }

    // Add calculated term for each component to kernel table
    for (int i = 0; i < COMPONENTS; ++i) {
        data_tables->kernels[kernel_index].spt_values[i]
            += kernel_values[i];
    }

    gsl_combination_free(comb_l);
    gsl_combination_free(comb_r);
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



short int compute_SPT_kernels(
        const short int arguments[], /* kernel arguments                                    */
        short int n,                 /* order in perturbation theory expansion              */
        short int time_step,         /* time step where SPT kernels are stored, should be 0 */
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

    short int kernel_index = kernel_index_from_arguments(arguments);

    // Alias pointer to kernel (TIME_STEPS x COMPONENTS table) we are working
    // with for convenience/readability
    kernel_t* const kernel = &data_tables->kernels[kernel_index];

    // Check if the SPT kernels are already computed
    if (kernel->ic_computed) return kernel_index;

    // For SPT kernels, F_1 = G_1 = ... = 1
    if (n == 1) {
        for (int i = 0; i < COMPONENTS; ++i) {
            kernel->spt_values[i] = 1.0;
        }
        return kernel_index;
    }

    // Only sum up to (including) floor(n/2), since partial_SPT_sum()
    // simultaneously computes terms m and (n-m)
    for (int m = 1; m <= n/2; ++m) {
        partial_SPT_sum(arguments,n,m,kernel_index,time_step,data_tables);
    }

    // Divide by overall factor in SPT recursion relation
    for (int i = 0; i < COMPONENTS; ++i) {
        kernel->spt_values[i] /= (2*n + 3) * (n - 1);
    }

    // Update kernel table
    kernel->ic_computed = true;

    return kernel_index;
}
