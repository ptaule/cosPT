/* Expressions for biased tracers in redshift-space are taken from 2004.10607 */
/* Notice however that the authors adopt a renormalization of the biases in which */
/* two b2-terms in Z3 do not contribute in the end. This is also the */
/* implementation here. */

#include <array>
#include <cmath>
#include <stdexcept>

#include "../include/biased_tracers.hpp"
#include "../include/combinatorics.hpp"
#include "../include/kernel_evolution.hpp"
#include "../include/parameters.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/tables.hpp"

using std::size_t;

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

inline double n1(
    double b1F1, /* b1 * F1 */
    double fG1, /* f * G1 */
    double mu
    )
{
    return b1F1 + fG1 * SQUARE(mu);
}

inline double n2(
        const int arguments[],
        int kernel_index,
        IntegrandTables& tables,
        double k,
        double mu,
        double f
        )
{
    double b1 = tables.bias_parameters.at(0);
    double b2 = tables.bias_parameters.at(1);
    double bG2 = tables.bias_parameters.at(2);

    double q1 = std::sqrt(tables.composite_dot_products()
            .at(static_cast<size_t>(arguments[0]))
            .at(static_cast<size_t>(arguments[0])));
    double q2 = std::sqrt(tables.composite_dot_products()
            .at(static_cast<size_t>(arguments[1]))
            .at(static_cast<size_t>(arguments[1])));

    /* cosine(q1,q2) */
    double cos12 = tables.composite_dot_products()
            .at(static_cast<size_t>(arguments[0]))
            .at(static_cast<size_t>(arguments[1]))
            / (q1*q2);

    double mu1 = tables.composite_los_projection().at(static_cast<size_t>(arguments[0])) / q1;
    double mu2 = tables.composite_los_projection().at(static_cast<size_t>(arguments[1])) / q2;

    double f_mu_k = f * mu * k;

    std::array<double, 2> F1; /* F1(q1), F1(q2) */
    std::array<double, 2> G1; /* G1(q1), G1(q2) */
    std::array<double, 2> Z1; /* Z1(q1), Z1(q2) */
    double F2 = 0;
    double G2 = 0;

    /* Initialise F1/G1 to SPT values */
    for (auto& el : F1) el = 1;
    for (auto& el : G1) el = 1;

    if (tables.loop_params.dynamics() == EDS_SPT) {
        compute_SPT_kernels(arguments, kernel_index, 2, tables);
        F2 = tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[0];
        G2 = tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[1];
    }
    else {
        KernelEvolver kernel_evolver(tables);
        kernel_evolver.compute(arguments, kernel_index, 2);
        F2 = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(0);
        G2 = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(1);

        int n1_args[N_KERNEL_ARGS_MAX];
        std::fill(&n1_args[1],
                &n1_args[tables.loop_params.n_kernel_args()],
                tables.loop_params.zero_label());

        for (size_t i = 0; i < 2; ++i) {
            n1_args[0] = arguments[i];
            kernel_index = kernel_evolver.compute(n1_args, -1, 1);
            F1[i] = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(0);
            G1[i] = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(1);
        }
    }

    Z1.at(0) = n1(b1 * F1.at(0), f * G1.at(0), mu1);
    Z1.at(1) = n1(b1 * F1.at(1), f * G1.at(1), mu2);

    return (
        (0.5 * b2 + bG2 * (SQUARE(cos12) - 1))
            * F1.at(0) * F1.at(1)
        + b1 * F2 + f * SQUARE(mu) * G2
        + 0.5 * f_mu_k * (
            mu1/q1 * G1.at(0) * Z1.at(1)
            + mu2/q2 * G1.at(1) * Z1.at(0)
        )
    );
}



/* Third order kernel for biased tracers in redshift-space, NOT symmetrized
 * over arguments */
inline double n3(
        const int arguments[],
        int kernel_index,
        IntegrandTables& tables,
        double k,
        double mu,
        double f
        )
{
    double b1 = tables.bias_parameters.at(0);
    double bG2 = tables.bias_parameters.at(2);
    double bGamma3 = tables.bias_parameters.at(3);

    double q1 = std::sqrt(tables.composite_dot_products()
            .at(static_cast<size_t>(arguments[0]))
            .at(static_cast<size_t>(arguments[0])));
    double q2 = std::sqrt(tables.composite_dot_products()
            .at(static_cast<size_t>(arguments[1]))
            .at(static_cast<size_t>(arguments[1])));
    double q3 = std::sqrt(tables.composite_dot_products()
            .at(static_cast<size_t>(arguments[2]))
            .at(static_cast<size_t>(arguments[2])));

    /* cosine(q1,q2) */
    double cos12 = tables.composite_dot_products()
            .at(static_cast<size_t>(arguments[0]))
            .at(static_cast<size_t>(arguments[1]))
            / (q1*q2);

    double mu1 = tables.composite_los_projection().at(static_cast<size_t>(arguments[0])) / q1;
    double mu2 = tables.composite_los_projection().at(static_cast<size_t>(arguments[1])) / q2;
    double mu3 = tables.composite_los_projection().at(static_cast<size_t>(arguments[2])) / q3;

    double f_mu_k = f * mu * k;

    std::array<double, 3> F1; /* F1(q1), F1(q2), F1(q3) */
    std::array<double, 3> G1; /* G1(q1), F1(q2), F1(q3) */
    std::array<double, 3> Z1; /* G1(q1), F1(q2), F1(q3) */

    double F3 = 0; /* F3(q1,q2,q3) */
    double G3 = 0; /* G3(q1,q2,q3) */

    /* Initialise F1/G1 to SPT values */
    for (auto& el : F1) el = 1;
    for (auto& el : G1) el = 1;

    if (tables.loop_params.dynamics() == EDS_SPT) {
        compute_SPT_kernels(arguments, kernel_index, 3, tables);
        F3 = tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[0];
        G3 = tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[1];
    }
    else {
        KernelEvolver kernel_evolver(tables);
        kernel_evolver.compute(arguments, kernel_index, 3);

        F3 = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(0);
        G3 = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(1);

        int n1_args[N_KERNEL_ARGS_MAX];
        std::fill(&n1_args[1],
                &n1_args[tables.loop_params.n_kernel_args()],
                tables.loop_params.zero_label());

        for (size_t i = 0; i < 3; ++i) {
            n1_args[0] = arguments[i];
            kernel_index = kernel_evolver.compute(n1_args, -1, 1);
            F1[i] = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(0);
            G1[i] = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(1);
        }
    }

    Z1.at(0) = n1(b1 * F1.at(0), f * G1.at(0), mu1);
    Z1.at(1) = n1(b1 * F1.at(1), f * G1.at(1), mu2);
    Z1.at(2) = n1(b1 * F1.at(2), f * G1.at(2), mu3);

    double result = (
        + b1 * F3 + f * SQUARE(mu) * G3
        + 0.5 * SQUARE(f_mu_k) * Z1.at(2) * mu1 * mu2 / (q1 * q2) * G1.at(0) * G1.at(1)
        + bG2 * f_mu_k * mu3/q3 * G1.at(2) * (SQUARE(cos12) - 1) * F1.at(0) * F1.at(1)
    );

    /* Get label corresponding to q1 + q2 */
    int n2_args[N_KERNEL_ARGS_MAX];
    n2_args[0] = arguments[0];
    n2_args[1] = arguments[1];
    std::fill(&n2_args[2],
              &n2_args[tables.loop_params.n_kernel_args()],
              tables.loop_params.zero_label());

    int q12_label = tables.sum_table.sum_labels(
        n2_args, tables.loop_params.n_kernel_args());

    if (q12_label != tables.loop_params.zero_label()) {
        double q12 = std::sqrt(tables.composite_dot_products()
                .at(static_cast<size_t>(q12_label))
                .at(static_cast<size_t>(q12_label)));

        /* cosine(q1+q2,q3) */
        double cos3_12 = tables.composite_dot_products()
                .at(static_cast<size_t>(q12_label))
                .at(static_cast<size_t>(arguments[2]))
                / (q12 * q3);
        /*  cosine(q1+q2,q3)^2 - 1 */
        double cos3_12_sqm1 = SQUARE(cos3_12) - 1;

        double mu12 = tables.composite_los_projection().at(static_cast<size_t>(q12_label)) / q12;

        int n2_kernel_index = compute_SPT_kernels(n2_args, -1, 2, tables); /* F2(q1,q2) and G2(q1,q2) */
        double F2 = tables.spt_kernels.at(static_cast<size_t>(n2_kernel_index)).values[0];
        double G2 = tables.spt_kernels.at(static_cast<size_t>(n2_kernel_index)).values[1];

        result += (
            + f_mu_k * mu3/q3 * G1.at(2) * (b1 * F2 + f * SQUARE(mu12) * G2)
            + f_mu_k * Z1.at(2) * mu12/q12 * G2
            + 2 * bG2 * cos3_12_sqm1 * F1.at(2) * F2
            + 2 * bGamma3 * cos3_12_sqm1 * F1.at(2) * (F2 - G2)
        );
    }
    return result;
}



int compute_rsd_biased_kernels(
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
        return kernel_index;
    }

    /* Compute sum of arguments, and its absolute value k */
    int sum = tables.sum_table.sum_labels(arguments,
                                          tables.loop_params.n_kernel_args());
    double k = std::sqrt(tables.composite_dot_products()
            .at(static_cast<size_t>(sum))
            .at(static_cast<size_t>(sum)));
    /* RSD growth factor f and L.o.S angle */
    double mu_los = tables.vars.mu_los;
    double rsd_f = tables.rsd_growth_f();

    double result = 0;

    if (n == 1) {
        double F1 = 1;
        double G1 = 1;
        if (tables.loop_params.dynamics() != EDS_SPT) {
            KernelEvolver kernel_evolver(tables);
            kernel_evolver.compute(arguments, kernel_index, 1);
            F1 = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(0);
            G1 = tables.kernels.at(static_cast<size_t>(kernel_index)).values.back().at(1);
        }
        result = tables.bias_parameters.at(0) * F1 + rsd_f * SQUARE(mu_los) * G1;
    }
    else if (n == 2) {
        result = n2(arguments, kernel_index, tables, k, mu_los, rsd_f);
    }
    else if (n == 3) {
        int args_rearranged[N_KERNEL_ARGS_MAX] = {};
        Combinations comb(3, 2);
        do {
            comb.rearrange_from_current_combination(arguments, args_rearranged,
                                                    static_cast<size_t>(2));
            comb.rearrange_from_current_complement(arguments, args_rearranged + 2,
                                                   static_cast<size_t>(1));
            result += n3(args_rearranged, kernel_index, tables, k, mu_los, rsd_f);
        } while (comb.next());
        result /= 3; /* 3 choose 2 */
    }
    else {
        throw std::runtime_error("n = 4 not implemented for biased tracers");
    }

    // Update kernel table
    kernel.computed = true;
    kernel.value = result;
    return kernel_index;
}
