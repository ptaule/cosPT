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

double spt_term(
        int m_l,
        int m_r,
        int component,
        int args_l[],
        int args_r[],
        int sum_l,
        int sum_r,
        IntegrandTables& tables
        )
{
    int n = m_l + m_r;

    int index_l = compute_SPT_kernels(args_l, -1, m_l, tables);
    int index_r = compute_SPT_kernels(args_r, -1, m_r, tables);

    int a,b;
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
            throw(std::invalid_argument(
                "SPT_term() does not accept argument 'component' which does "
                "not equal 0 or 1."));
        }

        return tables.spt_kernels.at(static_cast<size_t>(index_l)).values[1] *
               (a * tables.alpha()
                        .at(static_cast<size_t>(sum_l))
                        .at(static_cast<size_t>(sum_r)) *
                    tables.spt_kernels.at(static_cast<size_t>(index_r))
                        .values[0] +
                b * tables.beta()
                        .at(static_cast<size_t>(sum_l))
                        .at(static_cast<size_t>(sum_r)) *
                    tables.spt_kernels.at(static_cast<size_t>(index_r))
                        .values[1]);
}



void partial_SPT_sum(
        const int arguments[], /* kernel arguments                       */
        int n,                 /* kernel number                          */
        int m,                 /* sum index in kernel recursion relation */
        int kernel_index,
        IntegrandTables& tables
        )
{
    double partial_kernel_values[EDS_SPT_COMPONENTS] = {0};

    size_t n_kernel_args = tables.loop_params.n_kernel_args();
    int args_l[N_KERNEL_ARGS_MAX] = {0};
    int args_r[N_KERNEL_ARGS_MAX] = {0};

    // Initialize args_l and args_r
    std::fill(&args_l[0], &args_l[n_kernel_args], tables.loop_params.zero_label());
    std::fill(&args_r[0], &args_r[n_kernel_args], tables.loop_params.zero_label());

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

        for (int i = 0; i < EDS_SPT_COMPONENTS; ++i) {
            partial_kernel_values[i] +=
                spt_term(m, n-m, i, args_l, args_r, sum_l, sum_r, tables);
        }

        // When m != (n - m), we may additionally compute the (n-m)-term by
        // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
        // compute_SPT_kernels() only needs to sum up to (including) floor(n/2).
        if (m != n - m) {
            for (int i = 0; i < EDS_SPT_COMPONENTS; ++i) {
                partial_kernel_values[i] +=
                    spt_term(n-m, m, i, args_r, args_l, sum_r, sum_l, tables);
            }
        }
    } while (comb.next());

    // Devide through by symmetrization factor (n choose m)
    int n_choose_m = static_cast<int>(
            gsl_sf_choose(static_cast<unsigned int>(n),
                          static_cast<unsigned int>(m))
            );
    for (int i = 0; i < EDS_SPT_COMPONENTS; ++i) {
        partial_kernel_values[i] /= n_choose_m;
    }

    // Add calculated term for each component to kernel table
    for (size_t i = 0; i < EDS_SPT_COMPONENTS; ++i) {
        tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[i] +=
            partial_kernel_values[i];
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
        for (int i = 0; i < EDS_SPT_COMPONENTS; ++i) {
            kernel.values[i] = 1.0;
        }
        return kernel_index;
    }

    // Only sum up to (including) floor(n/2), since partial_SPT_sum()
    // simultaneously computes terms m and (n-m)
    for (int m = 1; m <= n/2; ++m) {
        partial_SPT_sum(arguments, n, m, kernel_index, tables);
    }

    // Divide by overall factor in SPT recursion relation
    for (int i = 0; i < EDS_SPT_COMPONENTS; ++i) {
        kernel.values[i] /= (2*n + 3) * (n - 1);
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
