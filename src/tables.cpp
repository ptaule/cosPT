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
#include "../include/tables.hpp"

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

using std::size_t;

Settings::Settings(
        short int n_loops,
        Spectrum spectrum,
        Dynamics dynamics,
        short int time_steps,
        short int pre_time_steps,
        short int components,
        double eta_i,
        double eta_f,
        double eta_asymp
        ) :
    dynamics(dynamics), spectrum(spectrum), n_loops(n_loops),
    components(components), time_steps(time_steps),
    pre_time_steps(pre_time_steps), eta_i(eta_i), eta_f(eta_f),
    eta_asymp(eta_asymp)
{
    if (n_loops < 1 || n_loops > 2) {
        throw(std::invalid_argument("Settings::Settings(): n_loops should be between 1 and 2."));
    }
    if (spectrum == POWERSPECTRUM) {
        n_coeffs = n_loops + 1;
        n_configs = pow(3, n_coeffs);
        n_kernels = (n_configs/3 + 1) * pow(4,n_loops);
        n_kernel_args = 2 * n_loops + 1;
        zero_label = get_zero_label(n_coeffs);

        /* Define block size for kernel indexing. A block consists of all */
        /* fundamental vector argument combinations. */
        kernel_index_block_size = pow(4, n_loops);
    }
    else if (spectrum == BISPECTRUM) {
        n_coeffs = n_loops + 2;
        n_configs = pow(3, n_coeffs);
        /* n_kernels = (n_configs/2 + 1) * pow(4, n_loops); */
        n_kernel_args = 2 * n_loops + 2;
        zero_label = get_zero_label(n_coeffs);
    }
}



Settings::Settings(
        short int n_loops,
        Spectrum spectrum,
        Dynamics dynamics,
        short int time_steps,
        short int components,
        double eta_i,
        double eta_f
        ) : Settings(n_loops, spectrum, dynamics, time_steps, 0, components,
            eta_i, eta_f, 0.0)
{
    if (dynamics != EVOLVE_SPT_IC) {
        throw(std::invalid_argument("Settings::Settings(): This constructor is only used for EVOLVE_SPT_IC dynamics."));
    }
}



Settings::Settings(
        short int n_loops,
        Spectrum spectrum,
        Dynamics dynamics
        ) :
    Settings(n_loops, spectrum, dynamics, 0, 0, 0, 0.0, 0.0, 0.0)
{
    if (dynamics != SPT) {
        throw(std::invalid_argument("This constructor is only used for SPT dynamics."));
    }
}



short int SumTable::sum_two_labels(short int a, short int b)
{
    short int n_coeffs = settings.n_coeffs;

    short int a_coeffs[N_COEFFS_MAX]   = {0};
    short int b_coeffs[N_COEFFS_MAX]   = {0};
    short int res_coeffs[N_COEFFS_MAX] = {0};

    label2config(a, a_coeffs, n_coeffs);
    label2config(b, b_coeffs, n_coeffs);

    for (int i = 0; i < n_coeffs; ++i) {
        res_coeffs[i] = a_coeffs[i] + b_coeffs[i];
    }
    return config2label(res_coeffs, n_coeffs);
}



SumTable::SumTable(const Settings& settings) : settings(settings)
{
    short int n_configs = settings.n_configs;
    short int zero_label = settings.zero_label;

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



short int SumTable::sum_labels(const short int labels[], size_t size) const
{
    short int zero_label = settings.zero_label;

    if (size == 1) return labels[0];

    short int result = labels[0];

    for (size_t i = 1; i < size; ++i) {
        if (labels[i] == zero_label) continue;
        result = sum_table[result][labels[i]];
    }

    // If DEBUG>=1, check that sum is an appropriate vector configuration,
    // i.e. that Q-coefficients are elements of (-1,0,1) and k-coefficient is
    // an element of (0,1)
#if DEBUG >= 1
    short int n_coeffs = settings.n_coeffs;
    short int res_coeffs[N_COEFFS_MAX];
    label2config(result, res_coeffs, n_coeffs);
    for (int i = 0; i < n_coeffs; ++i) {
        short int c = res_coeffs[i];
        if (!(c == -1 || c == 0 || c == 1))
            throw(std::logic_error("SumTable::sum_labels(): Sum of labels does not correspond to an appropriate configuration."));
    }
#endif

    return result;
}



Vec1D<double> initialize_eta_grid(const Settings& settings)
{
    Vec1D<double> eta_grid(settings.time_steps, 0.0);

    double eta_i = settings.eta_i;
    double eta_f = settings.eta_f;
    double time_steps = settings.time_steps;

    Dynamics dynamics = settings.dynamics;

    if (dynamics == EVOLVE_SPT_IC) {
        // Linear time step (including endpoints)
        double d_eta = std::abs(eta_f - eta_i)/(time_steps - 1);
        for (int i = 0; i < time_steps; ++i) {
            eta_grid.at(i) = eta_i + i * d_eta;
        }
    }
    else if (dynamics == EVOLVE_ASYMP_IC) {
        short int pre_time_steps = settings.pre_time_steps;
        double eta_asymp = settings.eta_asymp;
        // Linear time step (including endpoints)
        // Between eta_asymp and eta_i
        double d_eta = std::abs(eta_i - eta_asymp)/(pre_time_steps);
        for (int i = 0; i < pre_time_steps; ++i) {
            eta_grid.at(i) = eta_asymp + i*d_eta;
        }
        // Between eta_i and eta_f
        d_eta = std::abs(eta_f - eta_i)/(time_steps - pre_time_steps - 1);
        for (int i = pre_time_steps; i < time_steps; ++i) {
            eta_grid.at(i) = eta_i + (i - pre_time_steps) * d_eta;
        }
    }
    else {
        throw(std::invalid_argument("initialize_eta_grid(): incorrect dynamics."));
    }
    return eta_grid;
}



IntegrandTables::IntegrandTables(
        double k_a,
        double k_b,
        const Settings& settings,
        const SumTable& sum_table,
        const Vec1D<double>& eta_grid
        ) :
    k_a(k_a), k_b(k_b), settings(settings), sum_table(sum_table),
    vars(IntegrationVariables(settings.n_loops)), eta_grid(eta_grid)
{
    short int n_coeffs = settings.n_coeffs;
    short int n_configs = settings.n_configs;
    short int n_kernels = settings.n_kernels;
    short int time_steps = settings.time_steps;
    short int components = settings.components;

    bare_scalar_products.resize(n_coeffs);
    for (short int i = 0; i < n_coeffs; ++i) {
        bare_scalar_products.at(i).resize(n_coeffs);
    }

    scalar_products.resize(n_configs);
    alpha.resize(n_configs);
    beta.resize(n_configs);
    for (short int i = 0; i < n_configs; ++i) {
        scalar_products.at(i).resize(n_configs);
        alpha.at(i).resize(n_configs);
        beta.at(i).resize(n_configs);
    }

    if (settings.dynamics == SPT || settings.dynamics == EVOLVE_SPT_IC) {
        spt_kernels.resize(settings.n_kernels);
    }
    if (settings.dynamics == EVOLVE_SPT_IC ||
        settings.dynamics == EVOLVE_ASYMP_IC
        ) {
        spt_kernels.resize(n_kernels);
        kernels.resize(n_kernels);

        for (short int i = 0; i < n_kernels; ++i) {
            kernels.at(i).values.resize(time_steps);
            for (short int j = 0; j < time_steps; ++j) {
                kernels.at(i).values.at(j).resize(components);
            }
        }
    }
}



IntegrandTables::IntegrandTables(
        double k_a,
        const Settings& settings,
        const SumTable& sum_table,
        const Vec1D<double>& eta_grid
        ) :
    IntegrandTables(k_a, 0, settings, sum_table, eta_grid)
{
    if (settings.spectrum == BISPECTRUM) {
        throw(std::invalid_argument("IntegrandTables::IntegrandTables(): \
this constructor can only be used for spetrum = POWERSPECTRUM."));
    }
}



void IntegrandTables::reset()
{
    // bare_scalar_products, alpha, beta tables etc. are completely rewritten
    // by their respective compute-functions, hence no need to zero initialize

    if (settings.dynamics == SPT ||
        settings.dynamics == EVOLVE_SPT_IC) {
        reset_spt_kernels();
    }
    if (settings.dynamics == EVOLVE_SPT_IC ||
        settings.dynamics == EVOLVE_ASYMP_IC) {
        reset_kernels();
    }
}



void IntegrandTables::reset_spt_kernels()
{
    for (int i = 0; i < settings.n_kernels; ++i) {
        spt_kernels.at(i).computed = false;
        spt_kernels.at(i).values[0] = 0;
        spt_kernels.at(i).values[1] = 0;
    }
}



void IntegrandTables::reset_kernels()
{
    for (int i = 0; i < settings.n_kernels; ++i) {
        kernels.at(i).computed = false;

        for (short int j = 0; j < settings.time_steps; ++j) {
            std::fill(kernels.at(i).values.at(j).begin(),
                    kernels.at(i).values.at(j).end(), 0);
        }
    }
}



void IntegrandTables::ps_compute_bare_scalar_products()
{
    short int n_loops = settings.n_loops;
    short int n_coeffs = settings.n_coeffs;

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



/* Computes table of scalar_products given bare_scalar_products table */
void IntegrandTables::ps_compute_scalar_products()
{
    short int n_coeffs = settings.n_coeffs;
    short int n_configs = settings.n_configs;

    short int a_coeffs[N_COEFFS_MAX];
    short int b_coeffs[N_COEFFS_MAX];

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



void IntegrandTables::ps_compute_alpha_beta()
{
    short int n_configs = settings.n_configs;
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
            if (a != settings.zero_label) {
                double product_ab = scalar_products.at(a).at(b);
                double product_aa = scalar_products.at(a).at(a);

                alpha_val = 1 + product_ab/product_aa;

                if (b != settings.zero_label) {
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

void IntegrandTables::bs_compute_bare_scalar_products() {}
void IntegrandTables::bs_compute_scalar_products() {}
void IntegrandTables::bs_compute_alpha_beta() {}


void IntegrandTables::compute_tables()
{
    if (settings.spectrum == POWERSPECTRUM) {
        ps_compute_bare_scalar_products();
        ps_compute_scalar_products();
        ps_compute_alpha_beta();
    }
    else if (settings.spectrum == BISPECTRUM) {
        bs_compute_bare_scalar_products();
        bs_compute_scalar_products();
        bs_compute_alpha_beta();
    }
}



short int kernel_index_from_arguments(
        const short int arguments[],
        const Settings& settings
        )
{
    short int pow2[] = {1,2,4,8,16,32,64,128};

    short int n_coeffs = settings.n_coeffs;
    short int n_configs = settings.n_configs;
    short int zero_label = settings.zero_label;
    short int n_kernel_args = settings.n_kernel_args;
    short int kernel_index_block_size = settings.kernel_index_block_size;

    // In DEBUG-mode, check that non-zero arguments (zero_label) are unique
#if DEBUG >= 1
    if (!unique_elements(arguments, n_kernel_args,zero_label))
        throw(std::logic_error("kernel_index_from_arguments(): duplicate vector arguments passed."));
    short int n_k_labels = 0;
#endif

    short int index = 0;

    for (int i = 0; i < n_kernel_args; ++i) {
        // First, check if argument is a zero vector
        if (arguments[i] == zero_label) continue;

        // Argument is a k-type vector (i.e. on the form k + c_i Q_i) if k is
        // present. In our vector-label convention, k is the last coefficient,
        // hence +k is present if label >= 2 * n_configs/3
        short int k_threshold = 2 * n_configs/3;
        if (arguments[i] >= k_threshold) {
            index += (arguments[i] - k_threshold + 1) * kernel_index_block_size;
#if DEBUG >= 1
            /* Count k-type labels */
            ++n_k_labels;
#endif
        }
#if DEBUG >= 1
        /* We should not get -k in power spectrum computation */
        else if (arguments[i] < n_configs/3) {
            throw(std::logic_error("kernel_index_from_arguments(): got argument with -k."));
        }
#endif
        else {
#if DEBUG >= 1
            /* Check that this is in fact a fundamental vector */
            if(!pure_loop_label(arguments[i], n_coeffs, settings.spectrum))
                throw(std::logic_error("kernel_index_from_arguments(): argument is neither 0, k-type, or fundamental."));
#endif

            /* Convert argument_label to coefficient array */
            short int coeffs[N_COEFFS_MAX] = {0};
            label2config(arguments[i], coeffs, n_coeffs);

            /* The last coefficient is for k, hence we can skip this */
            for (int j = 0; j < n_coeffs - 1; ++j) {
                /* if - Q_j is present, add 2^(2j + 0/2) = 2^(2j)   */
                /* if + Q_j is present, add 2^(2j + 2/2) = 2^(2j+1) */
                if (coeffs[j] != 0) {
                    index += pow2[2 * j + (coeffs[j] + 1)/2];
                }
            }
        }
    }
#if DEBUG >= 1
    if (n_k_labels > 1)
        throw(std::logic_error("kernel_index_from_arguments(): more than one argument is of k-type."));
#endif

    return index;
}
