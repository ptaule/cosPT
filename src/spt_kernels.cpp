#include <algorithm>
#include <stdexcept>

extern "C" {
    #include <gsl/gsl_sf.h>
}

#include "../include/utilities.hpp"
#include "../include/combinatorics.hpp"
#include "../include/parameters.hpp"
#include "../include/tables.hpp"
#include "../include/spt_kernels.hpp"

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif


const std::array<double, EDS_SPT_COMPONENTS> vertex(
    int m_l,
    int m_r,
    const int args_l[],
    const int args_r[],
    int sum_l,
    int sum_r,
    IntegrandTables& tables
) {
    int n = m_l + m_r;

    int index_l = compute_SPT_kernels(args_l, -1, m_l, tables);
    int index_r = compute_SPT_kernels(args_r, -1, m_r, tables);

    const auto& alpha_val =
        tables.alpha().at(static_cast<size_t>(sum_l)).at(static_cast<size_t>(sum_r));
    const auto& beta_val  =
        tables.beta().at(static_cast<size_t>(sum_l)).at(static_cast<size_t>(sum_r));

    const auto& spt_l = tables.spt_kernels.at(static_cast<size_t>(index_l)).values;
    const auto& spt_r = tables.spt_kernels.at(static_cast<size_t>(index_r)).values;

    std::array<double, EDS_SPT_COMPONENTS> results;
    // component 0
    {
        int a = 2 * n + 1;
        int b = 2;
        results.at(0) = spt_l.at(1) * (a * alpha_val * spt_r.at(0)
            + b * beta_val * spt_r.at(1));
    }
    // component 1
    {
        int a = 3;
        int b = 2 * n;
        results.at(1) = spt_l.at(1) * (a * alpha_val * spt_r.at(0)
            + b * beta_val * spt_r.at(1));
    }
    return results;
}



void partial_SPT_sum(
        const int arguments[], /* kernel arguments                       */
        int n,                 /* kernel number                          */
        int m,                 /* sum index in kernel recursion relation */
        int kernel_index,
        IntegrandTables& tables
        )
{
    std::array<double, EDS_SPT_COMPONENTS> partial_kernel_values = {0};

    int zero_label = tables.loop_params.zero_label();

    size_t n_kernel_args = tables.loop_params.n_kernel_args();
    int args_l[N_KERNEL_ARGS_MAX] = {0};
    int args_r[N_KERNEL_ARGS_MAX] = {0};

    // Initialize args_l and args_r
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

        int sum_l = tables.sum_table.sum_labels(args_l, n_kernel_args);
        int sum_r = tables.sum_table.sum_labels(args_r, n_kernel_args);

        partial_kernel_values += vertex(m, n-m, args_l, args_r,
                                        sum_l, sum_r, tables);

        // When m != (n - m), we may additionally compute the (n-m)-term by
        // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
        // compute_SPT_kernels() only needs to sum up to (including) floor(n/2).
        if (m != n - m) {
            partial_kernel_values += vertex(n-m, m, args_r, args_l, sum_r, sum_l, tables);
        }
    } while (comb.next());

    // Devide through by symmetrization factor (n choose m)
    int n_choose_m = static_cast<int>(
            gsl_sf_choose(static_cast<unsigned int>(n),
                          static_cast<unsigned int>(m))
            );
    for (auto& val : partial_kernel_values) {
        val /= n_choose_m;
    }

    // Add calculated term for each component to kernel table
    auto& kernel_values =
        tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values;
    for (size_t i = 0; i < kernel_values.size(); ++i) {
        kernel_values.at(i) += partial_kernel_values.at(i);
    }
}



int compute_SPT_kernels(
        const int arguments[], /* kernel arguments             */
        int kernel_index,      /* index for kernel table       */
        int n,                 /* order in perturbation theory */
        IntegrandTables& tables
        )
{
    /* DEBUG: check that the number of non-zero arguments is in fact n, and
     * that kernel_index is in fact equivalent to arguments */
#if DEBUG >= 1
    kernel_computer_validate_n(arguments, n, tables);
    kernel_computer_validate_kernel_index(arguments, kernel_index, tables);
#endif

    // If kernel_index is not known, -1 is sent as argument
    if (kernel_index == -1) {
        kernel_index = tables.loop_params.args_2_kernel_index(arguments);
    }

    // Alias reference to kernel we are working with for convenience/readability
    SPTKernel& kernel = tables.spt_kernels.at(static_cast<size_t>(kernel_index));

    // Check if the SPT kernels are already computed
    if (kernel.computed) return kernel_index;

    // For SPT kernels, F_1 = G_1 = ... = 1
    if (n == 1) {
        for (size_t i = 0; i < EDS_SPT_COMPONENTS; ++i) {
            kernel.values.at(i) = 1.0;
        }
        return kernel_index;
    }

    // Only sum up to (including) floor(n/2), since partial_SPT_sum()
    // simultaneously computes terms m and (n-m)
    for (int m = 1; m <= n/2; ++m) {
        partial_SPT_sum(arguments, n, m, kernel_index, tables);
    }

    // Divide by overall factor in SPT recursion relation
    for (size_t i = 0; i < EDS_SPT_COMPONENTS; ++i) {
        kernel.values.at(i) /= (2*n + 3) * (n - 1);
    }

    // Update kernel table
    kernel.computed = true;

    return kernel_index;
}



void kernel_computer_validate_n(
        const int arguments[],
        int n,
        IntegrandTables& tables
        )
{
    int n_args = 0;
    for (size_t i = 0; i < tables.loop_params.n_kernel_args(); ++i) {
        if (arguments[i] != tables.loop_params.zero_label()){
            n_args++;
        }
    }
    if (n_args != n) {
        throw(std::logic_error(
            "number of non-zero-label arguments does not equal n."));
    }
}



void kernel_computer_validate_kernel_index(
        const int arguments[],
        int kernel_index,
        IntegrandTables& tables
        )
{
    int argument_index =
        tables.loop_params.args_2_kernel_index(arguments);
    if (kernel_index != -1 && argument_index != kernel_index) {
        throw(std::logic_error("Index computed from kernel arguments does not "
                               "equal kernel index."));
    }
}
