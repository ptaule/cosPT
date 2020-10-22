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



SumTable::SumTable(const LoopParameters& loop_params) :
    zero_label(loop_params.get_zero_label()),
    n_coeffs(loop_params.get_n_coeffs())
{
    int n_configs = loop_params.get_n_configs();
    int zero_label = loop_params.get_zero_label();

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
        double cos_ab,
        const LoopParameters& loop_params,
        const SumTable& sum_table,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid
        ) :
    k_a(k_a), k_b(k_b), cos_ab(cos_ab), loop_params(loop_params),
    sum_table(sum_table), ev_params(ev_params), eta_grid(eta_grid),
    vars(IntegrationVariables(loop_params.get_n_loops()))
{
    int n_coeffs = loop_params.get_n_coeffs();
    int n_configs = loop_params.get_n_configs();
    int n_kernels = loop_params.get_n_kernels();

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

    if (loop_params.get_dynamics() == SPT ||
        loop_params.get_dynamics() == EVOLVE_SPT_IC) {
        spt_kernels.resize(loop_params.get_n_kernels());
    }
    if (loop_params.get_dynamics() == EVOLVE_SPT_IC ||
        loop_params.get_dynamics() == EVOLVE_ASYMP_IC
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
        const LoopParameters& loop_params,
        const SumTable& sum_table,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid
        ) :
    IntegrandTables(k_a, 0, 0, loop_params, sum_table, ev_params, eta_grid)
{
    if (loop_params.get_spectrum() == BISPECTRUM) {
        throw(std::invalid_argument(
            "IntegrandTables::IntegrandTables(): this constructor can only be "
            "used for spetrum = POWERSPECTRUM."));
    }
}



void IntegrandTables::reset()
{
    // bare_scalar_products, alpha, beta tables etc. are completely rewritten
    // by their respective compute-functions, hence no need to zero initialize

    if (loop_params.get_dynamics() == SPT ||
        loop_params.get_dynamics() == EVOLVE_SPT_IC) {
        reset_spt_kernels();
    }
    if (loop_params.get_dynamics() == EVOLVE_SPT_IC ||
        loop_params.get_dynamics() == EVOLVE_ASYMP_IC) {
        reset_kernels();
    }
}



void IntegrandTables::reset_spt_kernels()
{
    for (int i = 0; i < loop_params.get_n_kernels(); ++i) {
        spt_kernels.at(i).computed = false;
        spt_kernels.at(i).values[0] = 0;
        spt_kernels.at(i).values[1] = 0;
    }
}



void IntegrandTables::reset_kernels()
{
    for (int i = 0; i < loop_params.get_n_kernels(); ++i) {
        kernels.at(i).computed = false;

        for (int j = 0; j < eta_grid.get_time_steps(); ++j) {
          std::fill(kernels.at(i).values.at(j).begin(),
                    kernels.at(i).values.at(j).end(), 0);
        }
    }
}



void IntegrandTables::ps_compute_bare_scalar_products()
{
    int n_loops = loop_params.get_n_loops();
    int n_coeffs = loop_params.get_n_coeffs();

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
    int n_coeffs = loop_params.get_n_coeffs();
    int n_configs = loop_params.get_n_configs();

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
    int n_configs = loop_params.get_n_configs();
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
            if (a != loop_params.get_zero_label()) {
                double product_ab = scalar_products.at(a).at(b);
                double product_aa = scalar_products.at(a).at(a);

                alpha_val = 1 + product_ab/product_aa;

                if (b != loop_params.get_zero_label()) {
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
    if (loop_params.get_spectrum() == POWERSPECTRUM) {
        ps_compute_bare_scalar_products ();
    }
    else if (loop_params.get_spectrum() == BISPECTRUM) {
        bs_compute_bare_scalar_products();
    }
    compute_scalar_products();
    compute_alpha_beta();
}
