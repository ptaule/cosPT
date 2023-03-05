/*
   rsd.cpp

   Created by Petter Taule on 23.09.2022
   Copyright (c) 2022 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <cstddef>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

#include "../include/combinatorics.hpp"
#include "../include/kernel_evolution.hpp"
#include "../include/integrand.hpp"
#include "../include/parameters.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/tables.hpp"
#include "../include/rsd.hpp"

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif


struct RSDIntegrationParameters {
    double k;
    const InputPowerSpectrum& ps;
};



double rsd_tree_level_integrand_l0(double mu, void* parameters)
{
    RSDIntegrationParameters params =
        *static_cast<RSDIntegrationParameters*>(parameters);
    return params.ps(params.k, mu);
}



double rsd_tree_level_integrand_l2(double mu, void* parameters)
{
    RSDIntegrationParameters params =
        *static_cast<RSDIntegrationParameters*>(parameters);
    return 0.5 * (3*mu*mu - 1) * params.ps(params.k, mu);
}



double rsd_tree_level_integrand_l4(double mu, void* parameters)
{
    RSDIntegrationParameters params =
        *static_cast<RSDIntegrationParameters*>(parameters);
    return params.ps(params.k, mu)
        * 0.125 * (35 * POW4(mu) - 30 * mu * mu + 3);
}



Vec1D<double> rsd_tree_level(
    double k,
    const InputPowerSpectrum& ps,
    size_t integrate_sub_regions,
    double integrate_atol,
    double integrate_rtol,
    int integrate_key
)
{
    Vec1D<double> results(3);
    double rsd_growth_f = ps.rsd_growth_f();

    if (!ps.ir_resum()) {
        for (auto& el : results) el = ps(k, 0);
        /* Linear monopole */
        results.at(0) *= (
            1 + 2.0/3.0 * rsd_growth_f +
            1.0/5.0 * SQUARE(rsd_growth_f)
        );
        /* Linear quadrupole */
        results.at(1) *= (
            4.0/3.0 * rsd_growth_f
            + 4.0/7.0 * SQUARE(rsd_growth_f)
        );
        /* Linear hexadecapole */
        results.at(2) *= (8.0/35.0 * SQUARE(rsd_growth_f));
    }
    else {
        gsl_integration_workspace* workspace =
            gsl_integration_workspace_alloc(integrate_sub_regions);
        gsl_function F;
        double error;
        int status;

        RSDIntegrationParameters params = { .k = k, .ps = ps };

        /* l=0 integration */
        F.function = rsd_tree_level_integrand_l0;
        F.params = static_cast<void*>(&params);

        status = gsl_integration_qag(&F, 0, 1, integrate_atol,
                                         integrate_rtol,
                                         integrate_sub_regions, integrate_key,
                                         workspace, &results.at(0), &error);

        if (status != 0) {
            throw std::runtime_error("RSD tree-level integration l=0 failed \
                    with error code" + std::to_string(status));
        }

        /* l=2 integration */
        F.function = rsd_tree_level_integrand_l2;
        status = gsl_integration_qag(&F, 0, 1, integrate_atol,
                                         integrate_rtol,
                                         integrate_sub_regions, integrate_key,
                                         workspace, &results.at(1), &error);

        if (status != 0) {
            throw std::runtime_error("RSD tree-level integration l=2 failed \
                    with error code" + std::to_string(status));
        }

        /* l=4 integration */
        F.function = rsd_tree_level_integrand_l4;
        status = gsl_integration_qag(&F, 0, 1, integrate_atol,
                                         integrate_rtol,
                                         integrate_sub_regions, integrate_key,
                                         workspace, &results.at(2), &error);

        if (status != 0) {
            throw std::runtime_error("RSD tree-level integration l=4 failed \
                    with error code" + std::to_string(status));
        }

        /* Multiply by two for integration mu=[-1,0] (integrands are even) */
        for (auto& el : results) el *= 2;

    }
    return results;
}



/* RSD factor from coordinate transformation */
inline double rsd_coord_transformation(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        )
{
#if DEBUG >= 1
    kernel_computer_validate_n(arguments, n, tables);
    kernel_computer_validate_kernel_index(arguments, kernel_index, tables);
#endif
    /* Compute sum of arguments */
    int sum = tables.sum_table.sum_labels(arguments,
            tables.loop_params.n_kernel_args());
    /* If sum of arguments is zero, the PT kernels are zero */
    if (sum == tables.loop_params.zero_label()) return 0;

    /* k = absolute value of sum of arguments */
    double k = std::sqrt(tables.comp_dot_products()
            .at(static_cast<size_t>(sum))
            .at(static_cast<size_t>(sum)));
    /* Sum of arguments dotted with L.o.S., normalized */
    double mu = tables.comp_los_projection().at(static_cast<size_t>(sum)) / k;

    kernel_index = compute_SPT_kernels(arguments, kernel_index, n, tables);

    double result =  (
            tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[0] +
            tables.rsd_growth_f() * SQUARE(mu) *
            tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[1]
           );
    return result;
}



/* RSD factor from jacobian: one factor from taylor expansion of exponential */
inline double rsd_jac_transformation(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        )
{
#if DEBUG >= 1
    kernel_computer_validate_n(arguments, n, tables);
    kernel_computer_validate_kernel_index(arguments, kernel_index, tables);
#endif
    /* Compute sum of arguments */
    int sum = tables.sum_table.sum_labels(arguments,
            tables.loop_params.n_kernel_args());
    /* If sum of arguments is zero, the PT kernels are zero */
    if (sum == tables.loop_params.zero_label()) return 0;

    /* k = absolute value of sum of arguments */
    double k = std::sqrt(tables.comp_dot_products()
            .at(static_cast<size_t>(sum))
            .at(static_cast<size_t>(sum)));
    /* Sum of arguments dotted with L.o.S. */
    double mu = tables.comp_los_projection().at(static_cast<size_t>(sum)) / k;

    kernel_index = compute_SPT_kernels(arguments, kernel_index, n, tables);

    return mu / k *
        tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[1];
}



int rsd_velocity_power(
        const int arguments[],
        int kernel_index,
        int n,
        int N,                 /* Number of velocity powers */
        IntegrandTables& tables
        )
{
#if DEBUG >= 1
    kernel_computer_validate_n(arguments, n, tables);
    kernel_computer_validate_kernel_index(arguments, kernel_index, tables);
#endif
    // If kernel_index is not known, -1 is sent as argument
    if (kernel_index == -1) {
        kernel_index = tables.loop_params.args_2_kernel_index(arguments);
    }

    /* If N < 1 or there are more factors N than wavenumbers, do nothing */
    if (N < 1 || n < N) return kernel_index;

    // Alias reference to kernel we are working with for convenience/readability
    RSDKernel& kernel = tables.vel_power_kernels.at(static_cast<size_t>(kernel_index)).
        at(static_cast<size_t>(N-1));

    // Check if the kernels are already computed
    if (kernel.computed) return kernel_index;

    size_t n_kernel_args = tables.loop_params.n_kernel_args();

    /* Base case N == 1 */
    if (N == 1) {
        double result = rsd_jac_transformation(arguments, kernel_index, n, tables);

        /* Update vel_power_kernels table, N is 1-indexed, therefore subtract 1 */
        tables.vel_power_kernels.at(static_cast<size_t>(kernel_index)).
            at(static_cast<size_t>(N-1)).value = result;
        tables.vel_power_kernels.at(static_cast<size_t>(kernel_index)).
            at(static_cast<size_t>(N-1)).computed = true;
        return kernel_index;
    }
    /* Special case N==n where each velocity factor only contains one argument */
    else if (N == n) {
        /* Temp args variable meant to contain *one* non-zero-label argument */
        int args[N_KERNEL_ARGS_MAX];
        std::fill(&args[0], &args[n_kernel_args],
                tables.loop_params.zero_label());

        double product = 1;
        for (int i = 0; i < n; ++i) {
            args[0] = arguments[i];
            product *= rsd_jac_transformation(args, -1, 1, tables);
        }
        /* Update vel_power_kernels table, N is 1-indexed, therefore subtract 1 */
        tables.vel_power_kernels.at(static_cast<size_t>(kernel_index)).
            at(static_cast<size_t>(N-1)).value = product;
        tables.vel_power_kernels.at(static_cast<size_t>(kernel_index)).
            at(static_cast<size_t>(N-1)).computed = true;
        return kernel_index;
    }
    /* Otherwise, recursive calculation: */

    double result = 0;
    double result_for_m = 0;

    int args_l[N_KERNEL_ARGS_MAX] = {0};
    int args_r[N_KERNEL_ARGS_MAX] = {0};

    for (int m = 1; m <= n/2; ++m) {
        result_for_m = 0;

        /* Initialize args_l and args_r */
        std::fill(&args_l[0], &args_l[n_kernel_args],
                tables.loop_params.zero_label());
        std::fill(&args_r[0], &args_r[n_kernel_args],
                tables.loop_params.zero_label());

        /* Go through all ways to pick m (unordered) elements from group of n */
        Combinations comb(n, m);
        do {
            /* Set args_l and args_r from current combination and complement,
             * respectively */
            comb.rearrange_from_current_combination(arguments, args_l,
                                                    static_cast<size_t>(m));
            comb.rearrange_from_current_complement(arguments, args_r,
                                                   static_cast<size_t>(n - m));

            /* Recursive step: the N-th power equals one power times the (N-1)-th power */
            int kernel_idx_l = rsd_velocity_power(args_l, -1, m, N-1, tables);

            result_for_m += (
                    rsd_jac_transformation(args_r, -1, n-m, tables)
                    * tables.vel_power_kernels.at(static_cast<size_t>(kernel_idx_l)).
                        at(static_cast<size_t>(N-2)).value
                    );

            /* When m != (n - m), we may additionally compute the (n-m)-term by
             * swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then the
             * for loop only needs to sum up to (including) floor(n/2). */
            if (m != n - m) {
                int kernel_idx_r = rsd_velocity_power(args_r, -1, n-m, N-1, tables);

                result_for_m += (
                        rsd_jac_transformation(args_l, -1, m, tables)
                        * tables.vel_power_kernels.at(static_cast<size_t>(kernel_idx_r)).
                            at(static_cast<size_t>(N-2)).value
                        );
            }
        } while (comb.next());

        /* Devide through by symmetrization factor (n choose m) */
        int n_choose_m = static_cast<int>(
                gsl_sf_choose(static_cast<unsigned int>(n),
                    static_cast<unsigned int>(m))
                );
        result += result_for_m / n_choose_m;
    }

    /* Update vel_power_kernels table */
    tables.vel_power_kernels.at(static_cast<size_t>(kernel_index)).
        at(static_cast<size_t>(N-1)).value = result;
    tables.vel_power_kernels.at(static_cast<size_t>(kernel_index)).
        at(static_cast<size_t>(N-1)).computed = true;
    return kernel_index;
}



double compute_rsd_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        )
{
#if DEBUG >= 1
    kernel_computer_validate_n(arguments, n, tables);
    kernel_computer_validate_kernel_index(arguments, kernel_index, tables);
#endif

    // If kernel_index is not known, -1 is sent as argument
    if (kernel_index == -1) {
        kernel_index = tables.loop_params.args_2_kernel_index(arguments);
    }

    // Alias reference to kernel we are working with for convenience/readability
    RSDKernel& kernel = tables.rsd_kernels.at(static_cast<size_t>(kernel_index));

    // Check if the kernels are already computed
    if (kernel.computed) {
        return kernel.value;
    }

    size_t n_kernel_args = tables.loop_params.n_kernel_args();

    /* Compute sum of arguments, and its absolute value k */
    int sum = tables.sum_table.sum_labels(arguments,
            tables.loop_params.n_kernel_args());
    double k = std::sqrt(tables.comp_dot_products()
            .at(static_cast<size_t>(sum))
            .at(static_cast<size_t>(sum)));
    /* RSD growth factor f and L.o.S angle */
    double mu_los = tables.vars.mu_los;
    double rsd_f = tables.rsd_growth_f();

    int args_l[N_KERNEL_ARGS_MAX] = {0};
    int args_r[N_KERNEL_ARGS_MAX] = {0};

    double result = 0;
    double result_for_m = 0;

    for (int m = 1; m <= n/2; ++m) {
        result_for_m = 0;

        /* Initialize args_l and args_r */
        std::fill(&args_l[0], &args_l[n_kernel_args],
                tables.loop_params.zero_label());
        std::fill(&args_r[0], &args_r[n_kernel_args],
                tables.loop_params.zero_label());

        /* Go through all ways to pick m (unordered) elements from group of n */
        Combinations comb(n, m);
        do {
            /* Set args_l and args_r from current combination and complement,
             * respectively */
            comb.rearrange_from_current_combination(arguments, args_l,
                                                    static_cast<size_t>(m));
            comb.rearrange_from_current_complement(arguments, args_r,
                                                   static_cast<size_t>(n - m));

            double result_N_loop = 0;
            for (int N = 1; N <= m; N++) {
                int kernel_idx_l = rsd_velocity_power(args_l, -1, m, N, tables);
                result_N_loop += (
                        pow(rsd_f * mu_los * k, N) / gsl_sf_fact(static_cast<unsigned int>(N)) *
                        tables.vel_power_kernels.at(static_cast<size_t>(kernel_idx_l)).
                            at(static_cast<size_t>(N-1)).value
                        );
            }
            result_for_m += result_N_loop * rsd_coord_transformation(args_r, -1, n-m, tables);

            /* When m != (n - m), we may additionally compute the (n-m)-term by
             * swapping args_l and m with args_r and (n-m). Then the
             * for loop only needs to sum up to (including) floor(n/2). */
            if (m != n - m) {
                result_N_loop = 0;
                for (int N = 1; N <= n-m; N++) {
                    int kernel_idx_r = rsd_velocity_power(args_r, -1, n-m, N, tables);
                    result_N_loop += (
                            pow(rsd_f * mu_los * k, N) / gsl_sf_fact(static_cast<unsigned int>(N)) *
                            tables.vel_power_kernels.at(static_cast<size_t>(kernel_idx_r)).
                                at(static_cast<size_t>(N-1)).value
                            );
                }
                result_for_m += result_N_loop * rsd_coord_transformation(args_l, -1, m, tables);
            }
        } while (comb.next());

        /* Devide through by symmetrization factor (n choose m) */
        int n_choose_m = static_cast<int>(
                gsl_sf_choose(static_cast<unsigned int>(n),
                    static_cast<unsigned int>(m))
                );
        result += result_for_m / n_choose_m;
    }

    /* Add contribution with no velocity powers */
    result += rsd_coord_transformation(arguments, kernel_index, n, tables);

    // Update kernel table
    kernel.computed = true;
    kernel.value = result;
    return result;
}
