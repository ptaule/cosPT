/*
   tables.cpp

   Created by Petter Taule on 29.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "../include/utilities.hpp"
#include "../include/parameters.hpp"
#include "../include/tables.hpp"

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

using std::size_t;

int SumTable::sum_two_labels(int a, int b)
{
    int n_coeffs = params.n_coeffs;

    int a_coeffs[N_COEFFS_MAX]   = {0};
    int b_coeffs[N_COEFFS_MAX]   = {0};
    int res_coeffs[N_COEFFS_MAX] = {0};

    label2config(a, a_coeffs, n_coeffs);
    label2config(b, b_coeffs, n_coeffs);

    for (int i = 0; i < n_coeffs; ++i) {
        res_coeffs[i] = a_coeffs[i] + b_coeffs[i];
    }
    return config2label(res_coeffs, n_coeffs);
}



SumTable::SumTable(const Parameters& params) : params(params)
{
    int n_configs = params.n_configs;
    int zero_label = params.zero_label;

    sum_table.resize(n_configs);
    for (int a = 0; a < n_configs; ++a) {
        sum_table[a].resize(n_configs);
        for (int b = 0; b < n_configs; ++b) {
            if (a == zero_label) {
                sum_table[a][b] = b;
            }
            else if (b == zero_label) {
                sum_table[a][b] = a;
            }
            else {
                sum_table[a][b] = sum_two_labels(a,b);
            }
        }
    }
}



int SumTable::sum_labels(const int labels[], size_t size) const
{
    int zero_label = params.zero_label;

    if (size == 1) return labels[0];

    int result = labels[0];

    for (size_t i = 1; i < size; ++i) {
        if (labels[i] == zero_label) continue;
        result = sum_table[result][labels[i]];
    }

    // If DEBUG>=1, check that sum is an appropriate vector configuration,
    // i.e. that Q-coefficients are elements of (-1,0,1) and k-coefficient is
    // an element of (0,1)
#if DEBUG >= 1
    int n_coeffs = params.n_coeffs;
    int res_coeffs[N_COEFFS_MAX];
    label2config(result, res_coeffs, n_coeffs);
    for (int i = 0; i < n_coeffs; ++i) {
        int c = res_coeffs[i];
        if (!(c == -1 || c == 0 || c == 1))
            throw(std::logic_error(
                "SumTable::sum_labels(): Sum of labels does not correspond to "
                "an appropriate configuration."));
    }
#endif

    return result;
}



IntegrandTables::IntegrandTables(
        double k_a,
        double k_b,
        const Parameters& params,
        const SumTable& sum_table,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid
        ) :
    k_a(k_a), k_b(k_b), params(params), sum_table(sum_table),
    ev_params(ev_params), eta_grid(eta_grid),
    vars(IntegrationVariables(params.n_loops))
{
    int n_coeffs = params.n_coeffs;
    int n_configs = params.n_configs;
    int n_kernels = params.n_kernels;

    if (params.spectrum == POWERSPECTRUM) {
        kernel_index_from_arguments = ps::kernel_index_from_arguments;
    }
    else if (params.spectrum == BISPECTRUM) {
        kernel_index_from_arguments = bs::kernel_index_from_arguments;
    }
    else {
        throw(std::invalid_argument("IntegrandTables::IntegrandTables(): Unknown "
                                    "dynamics in Parameters params."));
    }

    bare_scalar_products.resize(n_coeffs);
    for (int i = 0; i < n_coeffs; ++i) {
        bare_scalar_products.at(i).resize(n_coeffs);
    }

    scalar_products.resize(n_configs);
    alpha.resize(n_configs);
    beta.resize(n_configs);
    for (int i = 0; i < n_configs; ++i) {
        scalar_products.at(i).resize(n_configs);
        alpha.at(i).resize(n_configs);
        beta.at(i).resize(n_configs);
    }

    if (params.dynamics == SPT || params.dynamics == EVOLVE_SPT_IC) {
        spt_kernels.resize(params.n_kernels);
    }
    if (params.dynamics == EVOLVE_SPT_IC ||
        params.dynamics == EVOLVE_ASYMP_IC
        ) {
        spt_kernels.resize(n_kernels);
        kernels.resize(n_kernels);

        for (int i = 0; i < n_kernels; ++i) {
            kernels.at(i).values.assign(eta_grid.get_time_steps(),
                                        Vec1D<double>(COMPONENTS));
        }
    }
}



IntegrandTables::IntegrandTables(
        double k_a,
        const Parameters& params,
        const SumTable& sum_table,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid
        ) :
    IntegrandTables(k_a, 0, params, sum_table, ev_params, eta_grid)
{
    if (params.spectrum == BISPECTRUM) {
        throw(std::invalid_argument(
            "IntegrandTables::IntegrandTables(): this constructor can only be "
            "used for spetrum = POWERSPECTRUM."));
    }
}



void IntegrandTables::reset()
{
    // bare_scalar_products, alpha, beta tables etc. are completely rewritten
    // by their respective compute-functions, hence no need to zero initialize

    if (params.dynamics == SPT ||
        params.dynamics == EVOLVE_SPT_IC) {
        reset_spt_kernels();
    }
    if (params.dynamics == EVOLVE_SPT_IC ||
        params.dynamics == EVOLVE_ASYMP_IC) {
        reset_kernels();
    }
}



void IntegrandTables::reset_spt_kernels()
{
    for (int i = 0; i < params.n_kernels; ++i) {
        spt_kernels.at(i).computed = false;
        spt_kernels.at(i).values[0] = 0;
        spt_kernels.at(i).values[1] = 0;
    }
}



void IntegrandTables::reset_kernels()
{
    for (int i = 0; i < params.n_kernels; ++i) {
        kernels.at(i).computed = false;

        for (int j = 0; j < eta_grid.get_time_steps(); ++j) {
          std::fill(kernels.at(i).values.at(j).begin(),
                    kernels.at(i).values.at(j).end(), 0);
        }
    }
}



void IntegrandTables::ps_compute_bare_scalar_products()
{
    int n_loops = params.n_loops;
    int n_coeffs = params.n_coeffs;

    // Diagonal products correspond to Q1*Q1, etc.
    for (int i = 0; i < n_loops; ++i) {
        bare_scalar_products.at(i).at(i) =
            vars.magnitudes.at(i) * vars.magnitudes.at(i);
    }
    bare_scalar_products.at(n_coeffs - 1).at(n_coeffs - 1) = k_a * k_a;

    // Products involving k and Q_i has the form k*Q_i*cos(theta_i)
    for (int i = 0; i < n_coeffs - 1; ++i) {
        double value = k_a * vars.magnitudes.at(i) * vars.cos_theta.at(i);
        bare_scalar_products.at(n_coeffs - 1).at(i) = value;
        bare_scalar_products.at(i).at(n_coeffs - 1) = value;
    }

    if (n_loops > 1) {
        // Compute Q_1 * Q_i
        // (This is a special case since phi_1 is chosen to be zero.)
        double cos_theta_1 = vars.cos_theta.at(0);
        double sin_theta_1 = sqrt(1 - pow(cos_theta_1,2));
        double Q_1 = vars.magnitudes.at(0);
        for (int i = 1; i < n_loops; ++i) {
            double sin_theta_i = sqrt(1 - pow(vars.cos_theta.at(i),2));
            double value =
                sin_theta_1 * cos(vars.phi.at(i-1)) * sin_theta_i
                + cos_theta_1 * vars.cos_theta.at(i);
            value *= Q_1 * vars.magnitudes.at(i);

            bare_scalar_products.at(i).at(0) = value;
            bare_scalar_products.at(0).at(i) = value;
        }

        // Compute Q_i * Q_j for {i,j} != 1
        // (Q_i symbol is 1-indexed while arrays are 0-indexed)
        for (int i = 1; i < n_loops; ++i) {
            for (int j = 1; j < i; ++j) {
                double sin_theta_i = sqrt(1 - pow(vars.cos_theta.at(i),2));
                double sin_theta_j = sqrt(1 - pow(vars.cos_theta.at(j),2));

                double value =
                    cos(vars.phi.at(i-1)) * sin_theta_i * cos(vars.phi.at(j-1)) * sin_theta_j
                    + sin(vars.phi.at(i-1)) * sin_theta_i * sin(vars.phi.at(j-1)) * sin_theta_j
                    + vars.cos_theta.at(i) * vars.cos_theta.at(j);

                value *= vars.magnitudes.at(i) * vars.magnitudes.at(j);
                bare_scalar_products.at(i).at(j) = value;
                bare_scalar_products.at(j).at(i) = value;
            }
        }
    }
}



void IntegrandTables::bs_compute_bare_scalar_products()
{
}



/* Computes table of scalar_products given bare_scalar_products table */
void IntegrandTables::compute_scalar_products()
{
    int n_coeffs = params.n_coeffs;
    int n_configs = params.n_configs;

    int a_coeffs[N_COEFFS_MAX];
    int b_coeffs[N_COEFFS_MAX];

    // Scalar product matrix is symmetric, hence compute indices [a,b] and
    // [b,a] simultaneously
    for (int a = 0; a < n_configs; ++a) {
        for (int b = 0; b <= a; ++b) {
            double product_value = 0;

            label2config(a, a_coeffs, n_coeffs);
            label2config(b, b_coeffs, n_coeffs);
            for (int i = 0; i < n_coeffs; ++i) {
                for (int j = 0; j < n_coeffs; ++j) {
                    product_value += a_coeffs[i] * b_coeffs[j]
                        * bare_scalar_products.at(i).at(j);
                }
            }
            if (a == b) {
                scalar_products.at(a).at(a) = product_value;
            }
            else {
                scalar_products.at(a).at(b) = product_value;
                scalar_products.at(b).at(a) = product_value;
            }
        }
    }
}



void IntegrandTables::compute_alpha_beta()
{
    int n_configs = params.n_configs;
    for (int a = 0; a < n_configs; ++a) {
        for (int b = 0; b < n_configs; ++b) {
            // Special case when a == b
            if (a == b) {
                alpha.at(a).at(b) = 2.0;
                beta.at(a).at(b) = 2.0;
                continue;
            }

            double alpha_val = 0.0;
            double beta_val  = 0.0;

            // If the first argument is the zero-vector, alpha and beta remains 0
            // If the second argument is the zero-vector, beta remains 0
            if (a != params.zero_label) {
                double product_ab = scalar_products.at(a).at(b);
                double product_aa = scalar_products.at(a).at(a);

                alpha_val = 1 + product_ab/product_aa;

                if (b != params.zero_label) {
                    double product_bb = scalar_products.at(b).at(b);

                    beta_val = product_ab / 2.0
                        * ( 1.0 / product_aa + 1.0 / product_bb
                          + 2.0 * product_ab / (product_aa * product_bb)
                          );
                }
            }
            alpha.at(a).at(b) = alpha_val;
            beta.at(a).at(b) = beta_val;
        }
    }
}



void IntegrandTables::compute_tables()
{
    if (params.spectrum == POWERSPECTRUM) {
        ps_compute_bare_scalar_products();
    }
    else if (params.spectrum == BISPECTRUM) {
        bs_compute_bare_scalar_products();
    }
    compute_scalar_products();
    compute_alpha_beta();
}



int ps::kernel_index_from_arguments(
        const int arguments[],
        const Parameters& params
        )
{
   /* Precompute powers of two for speedup */
    int pow2[] = {1,2,4,8,16,32,64,128};

    int zero_label                 = params.zero_label;
    int n_kernel_args              = params.n_kernel_args;
    int single_loop_block_size     = params.single_loop_block_size;
    int single_loop_label_max      = params.single_loop_label_max;

    const Vec1D<int>& single_loops = params.single_loops;

    // In DEBUG-mode, check that non-zero arguments (zero_label) are unique
#if DEBUG >= 1
    int n_coeffs              = params.n_coeffs;
    int single_loop_label_min = params.single_loop_label_min;

    if (!unique_elements(arguments, n_kernel_args, zero_label))
        throw(std::logic_error("ps::kernel_index_from_arguments(): duplicate "
                               "vector arguments passed."));
    int n_k_labels = 0;
#endif

    int index = 0;

    for (int i = 0; i < n_kernel_args; ++i) {
        // First, check if argument is a zero vector
        if (arguments[i] == zero_label) continue;

        // Argument is a k-type vector (i.e. on the form k + c_i Q_i) if k is
        // present. In our vector-label convention, k is the last coefficient,
        // hence +k is present if label > single_loop_label_max
        if (arguments[i] > single_loop_label_max) {
            index += (arguments[i] - single_loop_label_max) * single_loop_block_size;
#if DEBUG >= 1
            /* Count k-type labels */
            ++n_k_labels;
#endif
        }
#if DEBUG >= 1
        /* We should not get -k in power spectrum computation */
        else if (arguments[i] < single_loop_label_min) {
            throw(std::logic_error(
                "ps::kernel_index_from_arguments(): got argument with -k."));
        }
#endif
        else {
            /* Single loop */
            for (size_t j = 0; j < single_loops.size(); ++j) {
                if (arguments[i] == single_loops[j]) {
                    index += pow2[j];
                    break;
                }
            }
#if DEBUG >= 1
            /* Check that this is in fact a fundamental vector */
            if(!single_loop_label(arguments[i], n_coeffs, params.spectrum))
                throw(std::logic_error(
                    "ps::kernel_index_from_arguments(): argument is neither 0, "
                    "k-type, or fundamental."));
#endif
        }
    }
#if DEBUG >= 1
    if (n_k_labels > 1)
        throw(std::logic_error("ps::kernel_index_from_arguments(): more than one "
                               "argument is of k-type."));
#endif

    return index;
}



int bs::kernel_index_from_arguments(
        const int arguments[],
        const Parameters& params
        )
{
   /* Precompute powers of two for speedup */
    int pow2[] = {1,2,4,8,16,32,64,128};

    int n_configs                  = params.n_configs;
    int zero_label                 = params.zero_label;
    int n_kernel_args              = params.n_kernel_args;
    int single_loop_block_size     = params.single_loop_block_size;

    int single_loop_label_min      = params.single_loop_label_min;
    int single_loop_label_max      = params.single_loop_label_max;
    int first_composite_block_size = params.first_composite_block_size;

    const Vec1D<int>& single_loops = params.single_loops;

    // In DEBUG-mode, check that non-zero arguments (zero_label) are unique
#if DEBUG >= 1
    if (!unique_elements(arguments, n_kernel_args, zero_label))
        throw(std::logic_error("bs::kernel_index_from_arguments(): duplicate "
                               "vector arguments passed."));
#endif

    int index = 0;

    /* Counter of composite arguments, i.e. not single loop momenta */
    int n_composite = 0;

    for (int i = 0; i < n_kernel_args; ++i) {
        // First, check if argument is a zero vector
        if (arguments[i] == zero_label) continue;

        // Argument is not single loop if label < single_loop_label_min or
        // label > single_loop_label_max */
        if (arguments[i] < single_loop_label_min) {
            index += (arguments[i] + 1) *
                (n_composite > 0 ? single_loop_block_size :
                 first_composite_block_size);

            ++n_composite;
        }
        else if (arguments[i] > single_loop_label_max) {
            int block = n_composite > 0 ? single_loop_block_size : first_composite_block_size;
            index += (arguments[i] - n_configs/9 + 1) * block;

            ++n_composite;
        }
        else {
            /* Single loop */
            for (size_t j = 0; j < single_loops.size(); ++j) {
                if (arguments[i] == single_loops[j]) {
                    index += pow2[j];
                    goto found_single_loop;
                }
            }
            /* Went through loop, which means that argument is composite, hence
             * connecting line with overall loop */
            index += (arguments[i] - single_loop_label_min + 1) *
                (n_composite > 0 ? single_loop_block_size :
                 first_composite_block_size);
            ++n_composite;

found_single_loop: ;
        }
    }
#if DEBUG >= 1
    if (n_composite > 2)
        throw(std::logic_error("bs::kernel_index_from_arguments(): more than two "
                               "arguments is of composite type."));
#endif

    return index;
}
