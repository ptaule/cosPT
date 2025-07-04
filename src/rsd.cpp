#include <cmath>
#include <cstddef>

extern "C" {
    #include <gsl/gsl_sf.h>
}

#include "../include/combinatorics.hpp"
#include "../include/parameters.hpp"
#include "../include/rsd.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/tables.hpp"

using std::size_t;

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

inline double k_of_args(size_t label, IntegrandTables& tables) {
    return std::sqrt(tables.composite_dot_products()(label, label));
}



inline double mu_of_args(size_t label, IntegrandTables& tables) {
   return
       tables.composite_los_projection().at(label) / k_of_args(label, tables);
}



inline double mu_of_args(size_t label, IntegrandTables& tables, double k) {
   return
       tables.composite_los_projection().at(label) / k;
}



inline size_t get_label(const int arguments[], IntegrandTables& tables) {
    return static_cast<size_t>(
        tables.sum_table(arguments,tables.loop_structure.n_kernel_args()));
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
    size_t sum = get_label(arguments, tables);
    /* If sum of arguments is zero, the PT kernels are zero */
    if (sum == static_cast<size_t>(tables.loop_structure.zero_label())) return 0;

    double mu = mu_of_args(sum, tables);

    kernel_index = compute_SPT_kernels(arguments, kernel_index, n, tables);
    auto& spt_values =
        tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values;
    return spt_values.at(0) + tables.rsd_growth_f() * mu * mu * spt_values.at(1);
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
    size_t sum = get_label(arguments, tables);
    /* If sum of arguments is zero, the PT kernels are zero */
    if (sum == static_cast<size_t>(tables.loop_structure.zero_label())) return 0;

    double k = k_of_args(static_cast<size_t>(sum), tables);
    double mu = mu_of_args(static_cast<size_t>(sum), tables, k);

    kernel_index = compute_SPT_kernels(arguments, kernel_index, n, tables);
    auto& spt_value =
        tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values.at(1);
    return mu / k * spt_value;
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
        kernel_index = tables.loop_structure.args_to_kernel_index(arguments);
    }

    /* If N < 1 or there are more factors N than wavenumbers, do nothing */
    if (N < 1 || n < N) return kernel_index;

    int zero_label = tables.loop_structure.zero_label();
    size_t n_kernel_args = tables.loop_structure.n_kernel_args();
    /* Convert N and kernel_index to size_t for convenience since they are used
     * as indices */
    size_t k_idx_t = static_cast<size_t>(kernel_index);
    size_t N_m1_t = static_cast<size_t>(N-1);
    size_t N_m2_t = static_cast<size_t>(N-2);
    // Alias reference to kernel we are working with for convenience/readability
    RSDKernel& kernel = tables.vel_power_kernels(k_idx_t, N_m1_t);

    // Check if the kernels are already computed
    if (kernel.computed) return kernel_index;

    /* Base case N == 1 */
    if (N == 1) {
        kernel.value = rsd_jac_transformation(arguments, kernel_index, n, tables);
        kernel.computed = true;
        return kernel_index;
    }
    /* Special case N==n where each velocity factor only contains one argument */
    else if (N == n) {
        double product = 1;
        int args[N_KERNEL_ARGS_MAX];
        std::fill(&args[0], &args[n_kernel_args], zero_label);

        for (int i = 0; i < n; ++i) {
            args[0] = arguments[i];
            product *= rsd_jac_transformation(args, -1, 1, tables);
        }
        kernel.value = product;
        kernel.computed = true;
        return kernel_index;
    }
    /* Otherwise, recursive calculation: */

    double result = 0;
    double result_for_m = 0;

    int args_l[N_KERNEL_ARGS_MAX];
    int args_r[N_KERNEL_ARGS_MAX];

    for (int m = 1; m <= n/2; ++m) {
        result_for_m = 0;

        /* Initialize args_l and args_r */
        std::fill(&args_l[0], &args_l[n_kernel_args], zero_label);
        std::fill(&args_r[0], &args_r[n_kernel_args], zero_label);

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
            int kernel_idx_l = rsd_velocity_power(args_l,
                                                  -1, m,
                                                  N-1, tables);
            size_t k_idx_l_t = static_cast<size_t>(kernel_idx_l);

            result_for_m += (
                    rsd_jac_transformation(args_r, -1, n-m, tables)
                    * tables.vel_power_kernels(k_idx_l_t, N_m2_t).value
                    );

            /* When m != (n - m), we may additionally compute the (n-m)-term by
             * swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then the
             * for loop only needs to sum up to (including) floor(n/2). */
            if (m != n - m) {
                int kernel_idx_r = rsd_velocity_power(args_r,
                                                      -1, n-m,
                                                      N-1, tables);
                size_t k_idx_r_t = static_cast<size_t>(kernel_idx_r);

                result_for_m += (
                        rsd_jac_transformation(args_l,
                                               -1, m,
                                               tables)
                        * tables.vel_power_kernels(k_idx_r_t, N_m2_t).value
                        );
            }
        } while (comb.next());

        /* Devide through by symmetrization factor (n choose m) */
        size_t n_t = static_cast<size_t>(n);
        size_t m_t = static_cast<size_t>(m);
        result += result_for_m / binomial_coeffs[n_t][m_t];
    }

    /* Update vel_power_kernels table */
    kernel.value = result;
    kernel.computed = true;
    return kernel_index;
}



void compute_rsd_kernels(
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
        kernel_index = tables.loop_structure.args_to_kernel_index(arguments);
    }

    // Alias reference to kernel we are working with for convenience/readability
    RSDKernel& kernel = tables.rsd_kernels.at(static_cast<size_t>(kernel_index));

    // Check if the kernels are already computed
    if (kernel.computed) {
        return;
    }

    size_t n_kernel_args = tables.loop_structure.n_kernel_args();
    int zero_label = tables.loop_structure.zero_label();

    /* Compute sum of arguments, and its absolute value k */
    size_t sum = get_label(arguments, tables);
    double k = std::sqrt(tables.composite_dot_products()(sum, sum));
    /* RSD growth factor f and L.o.S angle */
    double mu_los = tables.vars.mu_los;
    double rsd_f = tables.rsd_growth_f();

    int args_l[N_KERNEL_ARGS_MAX] = {0};
    int args_r[N_KERNEL_ARGS_MAX] = {0};

    double result = 0;
    double result_for_m = 0;

    // Precompute powers of rsd_f * mu_los * k / (N-1)! for N up to n
    size_t n_t = static_cast<size_t>(n);
    std::vector<double> powers(n_t > 1 ? n_t : 1);
    double fmuk = rsd_f * mu_los * k;
    powers[0] = fmuk;
    for (size_t N = 1; N < powers.size(); ++N) {
        powers[N] = powers[N-1] * fmuk / static_cast<double>(N+1);
    }

    auto get_velocity_power= [&](const int* args, int m) -> double {
        double result = 0.0;
        for (size_t N = 1; N <= static_cast<size_t>(m); ++N) {
            int kernel_idx = rsd_velocity_power(args,
                                                -1, m,
                                                static_cast<int>(N),
                                                tables);
            size_t k_idx_t = static_cast<size_t>(kernel_idx);
            size_t N_m1_t = N - 1;
            result += powers[N_m1_t] * tables.vel_power_kernels(k_idx_t, N_m1_t).value;
        }
        return result;
    };

    for (int m = 1; m <= n/2; ++m) {
        result_for_m = 0;

        /* Initialize args_l and args_r */
        std::fill(&args_l[0], &args_l[n_kernel_args], zero_label);
        std::fill(&args_r[0], &args_r[n_kernel_args], zero_label);

        /* Go through all ways to pick m (unordered) elements from group of n */
        Combinations comb(n, m);
        do {
            /* Set args_l and args_r from current combination and complement,
             * respectively */
            comb.rearrange_from_current_combination(arguments, args_l,
                                                    static_cast<size_t>(m));
            comb.rearrange_from_current_complement(arguments, args_r,
                                                   static_cast<size_t>(n - m));

            double vel_power_result = get_velocity_power(args_l, m);
            result_for_m += vel_power_result *
                rsd_coord_transformation(args_r, -1,
                                         n-m, tables);

            if (m != n - m) {
                vel_power_result = get_velocity_power(args_r, n-m);
                result_for_m += vel_power_result *
                    rsd_coord_transformation(args_l,
                                             -1, m,
                                             tables);
            }
        } while (comb.next());

        /* Devide through by symmetrization factor (n choose m) */
        size_t m_t = static_cast<size_t>(m);
        result += result_for_m / binomial_coeffs[n_t][m_t];
    }

    /* Add contribution with no velocity powers */
    result += rsd_coord_transformation(arguments, kernel_index, n, tables);

    // Update kernel table
    kernel.computed = true;
    kernel.value = result;
}
