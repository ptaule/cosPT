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

using std::size_t;

int SumTable::sum_two_labels(int a, int b)
{
    Vec1D<int> a_coeffs(n_coeffs, 0);
    Vec1D<int> b_coeffs(n_coeffs, 0);
    Vec1D<int> res_coeffs(n_coeffs, 0);

    label2config(a, a_coeffs);
    label2config(b, b_coeffs);

    for (size_t i = 0; i < n_coeffs; ++i) {
        res_coeffs.at(i) = a_coeffs.at(i) + b_coeffs.at(i);
    }
    return config2label(res_coeffs);
}



SumTable::SumTable(const LoopParameters& loop_params) :
    zero_label(loop_params.zero_label()),
    n_coeffs(loop_params.n_coeffs())
{
    size_t n_configs = loop_params.n_configs();

    sum_table.resize(n_configs);
    for (size_t a = 0; a < n_configs; ++a) {
        sum_table.at(a).resize(n_configs);
        for (size_t b = 0; b < n_configs; ++b) {
            if (a == static_cast<size_t>(zero_label)) {
                sum_table.at(a).at(b) = static_cast<int>(b);
            }
            else if (b == static_cast<size_t>(zero_label)) {
                sum_table.at(a).at(b) = static_cast<int>(a);
            }
            else {
                sum_table.at(a).at(b) =
                    sum_two_labels(static_cast<int>(a),static_cast<int>(b));
            }
        }
    }
}



/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

int SumTable::sum_labels(const int labels[], size_t size) const
{
    if (size == 1) return labels[0];

    int result = labels[0];

    for (size_t i = 1; i < size; ++i) {
        if (labels[i] == zero_label) continue;
        result = sum_table.at(static_cast<size_t>(result))
                     .at(static_cast<size_t>(labels[i]));
    }

    /* If DEBUG>=1, check that sum is an appropriate vector configuration, i.e.
     * that Q-coefficients are elements of (-1,0,1) and k-coefficient is an
     * element of (0,1) */
#if DEBUG >= 1
    Vec1D<int> res_coeffs(n_coeffs, 0);
    label2config(result, res_coeffs);
    for (size_t i = 0; i < n_coeffs; ++i) {
        int c = res_coeffs.at(i);
        if (!(c == -1 || c == 0 || c == 1))
            throw(std::logic_error(
                "SumTable::sum_labels(): Sum of labels does not correspond to "
                "an appropriate configuration."));
    }
#endif

    return result;
}



EtaGrid::EtaGrid(
        size_t pre_time_steps,
        size_t time_steps,
        double eta_ini,
        double eta_fin,
        double eta_asymp
       ) :
    pre_time_steps_(pre_time_steps), time_steps_(time_steps),
    eta_ini_(eta_ini), eta_fin_(eta_fin), eta_asymp_(eta_asymp)
{
    grid_.resize(time_steps_);

    // Linear time step (including endpoints)
    // Between eta_asymp and eta_ini
    double d_eta = std::abs(eta_ini_ - eta_asymp_) /
                   static_cast<double>(pre_time_steps_);
    for (size_t i = 0; i < pre_time_steps_; ++i) {
        grid_.at(i) = eta_asymp_ + static_cast<double>(i)*d_eta;
    }
    // Between eta_ini and eta_f
    d_eta = std::abs(eta_fin_ - eta_ini_) /
            static_cast<double>(time_steps_ - pre_time_steps_ - 1);
    for (size_t i = pre_time_steps_; i < time_steps_; ++i) {
        grid_.at(i) = eta_ini_ + static_cast<double>(i - pre_time_steps_) * d_eta;
    }
}



EtaGrid::EtaGrid(
        size_t time_steps,
        double eta_ini,
        double eta_fin
       ) :
    pre_time_steps_(0), time_steps_(time_steps),
    eta_ini_(eta_ini), eta_fin_(eta_fin), eta_asymp_(0.0)
{
    grid_.resize(time_steps_);

    // Linear time step (including endpoints)
    double d_eta =
        std::abs(eta_fin_ - eta_ini_) / static_cast<double>(time_steps_ - 1);
    for (size_t i = 0; i < time_steps_; ++i) {
        grid_.at(i) = eta_ini_ + static_cast<double>(i) * d_eta;
    }
}



std::ostream& operator<<(std::ostream& out, const EtaGrid& eta_grid) {
    for (size_t i = 0; i < eta_grid.time_steps(); ++i) {
        out << eta_grid[i] << std::endl;
    }
    return out;
}



IntegrandTables::IntegrandTables(
        double k_a,
        double k_b,
        double cos_ab,
        double mu_los,
        double rsd_growth_f,
        const LoopParameters& loop_params,
        const SumTable& sum_table,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid
        ) :
    k_a(k_a), k_b(k_b), cos_ab(cos_ab), mu_los_(mu_los), rsd_f(rsd_growth_f),
    loop_params(loop_params), sum_table(sum_table), ev_params(ev_params),
    eta_grid(eta_grid),
    vars(IntegrationVariables(static_cast<size_t>(loop_params.n_loops())))
{
    int n_loops = loop_params.n_loops();
    size_t n_coeffs = loop_params.n_coeffs();
    size_t n_configs = loop_params.n_configs();
    size_t n_kernels = loop_params.n_kernels();

    bare_dot_prod.resize(n_coeffs);
    for (size_t i = 0; i < n_coeffs; ++i) {
        bare_dot_prod.at(i).resize(n_coeffs);
    }

    comp_dot_prod.resize(n_configs);
    alpha_.resize(n_configs);
    beta_.resize(n_configs);
    for (size_t i = 0; i < n_configs; ++i) {
        comp_dot_prod.at(i).resize(n_configs);
        alpha_.at(i).resize(n_configs);
        beta_.at(i).resize(n_configs);
    }
    a_coeffs.resize(n_coeffs);
    b_coeffs.resize(n_coeffs);

    if (loop_params.dynamics() == EDS_SPT ||
        loop_params.dynamics() == EVOLVE_IC_EDS) {
        spt_kernels.resize(n_kernels);
    }
    if (loop_params.dynamics() == EVOLVE_IC_EDS ||
        loop_params.dynamics() == EVOLVE_IC_ASYMP
        ) {
        kernels.resize(n_kernels);

        for (size_t i = 0; i < n_kernels; ++i) {
            kernels.at(i).values.assign(eta_grid.time_steps(),
                Vec1D<double>(COMPONENTS));
        }
    }
    if (loop_params.rsd()) {
        bare_los_projection_.resize(n_coeffs);
        comp_los_projection_.resize(n_configs);

        rsd_kernels.resize(n_kernels);
        vel_power_kernels.resize(n_kernels);

        for (size_t i = 0; i < n_kernels; ++i) {
            vel_power_kernels.at(i).resize(2*static_cast<size_t>(n_loops));
        }
    }
}



IntegrandTables::IntegrandTables(
        double k_a,
        double mu_los,
        double rsd_growth_f,
        const LoopParameters& loop_params,
        const SumTable& sum_table,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid
        ) :
    IntegrandTables(k_a, 0, 0, mu_los, rsd_growth_f, loop_params, sum_table, ev_params, eta_grid)
{
    if (loop_params.spectrum() == BISPECTRUM) {
        throw(std::invalid_argument(
            "IntegrandTables::IntegrandTables(): this constructor can only be "
            "used for spetrum = POWERSPECTRUM."));
    }
}



void IntegrandTables::reset()
{
    // bare_dot_prod, alpha_, beta_ tables etc. are completely rewritten
    // by their respective compute-functions, hence no need to zero initialize

    if (loop_params.dynamics() == EDS_SPT) {
        reset_spt_kernels();
    }
    else if (loop_params.dynamics() == EVOLVE_IC_ASYMP ) {
        reset_kernels();
    }
    else {
        /* loop_params.dynamics() == EVOLVE_IC_EDS */
        reset_spt_kernels();
        reset_kernels();
    }
    if (loop_params.rsd()) {
        reset_rsd_kernels();
    }
}



void IntegrandTables::reset_spt_kernels()
{
    for (size_t i = 0; i < loop_params.n_kernels(); ++i) {
        spt_kernels.at(i).computed = false;
        spt_kernels.at(i).values[0] = 0;
        spt_kernels.at(i).values[1] = 0;
    }
}



void IntegrandTables::reset_kernels()
{
    for (size_t i = 0; i < loop_params.n_kernels(); ++i) {
        kernels.at(i).computed = false;

        for (size_t j = 0; j < eta_grid.time_steps(); ++j) {
          std::fill(kernels.at(i).values.at(j).begin(),
                    kernels.at(i).values.at(j).end(), 0);
        }
    }
}



void IntegrandTables::reset_rsd_kernels()
{
    for (size_t i = 0; i < loop_params.n_kernels(); ++i) {
        rsd_kernels.at(i).computed = false;
        rsd_kernels.at(i).value = 0;

        for (auto& el : vel_power_kernels.at(i)) {
            el.value    = 0;
            el.computed = false;
        }
    }
}



void IntegrandTables::ps_compute_bare_dot_prod()
{
    int n_loops = loop_params.n_loops();
    size_t n_coeffs = loop_params.n_coeffs();

    size_t k_a_idx = n_coeffs - 1;

    // Diagonal products correspond to Q1*Q1, etc.
    for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
        bare_dot_prod.at(i).at(i) = SQUARE(vars.magnitudes.at(i));
    }
    bare_dot_prod.at(k_a_idx).at(k_a_idx) = SQUARE(k_a);

    // Products involving k_a and Q_i has the form k*Q_i*cos(theta_i)
    for (size_t i = 0; i < n_coeffs - 1; ++i) {
        double value = k_a * vars.magnitudes.at(i) * vars.cos_theta.at(i);
        bare_dot_prod.at(k_a_idx).at(i) = value;
        bare_dot_prod.at(i).at(k_a_idx) = value;
    }

    if (n_loops > 1) {
        // Compute Q_1 * Q_i
        // (This is a special case since phi_1 is chosen to be zero.)
        double cos_theta_1 = vars.cos_theta.at(0);
        double sin_theta_1 = sqrt(1 - SQUARE(cos_theta_1));
        double Q_1 = vars.magnitudes.at(0);
        for (size_t i = 1; i < static_cast<size_t>(n_loops); ++i) {
            double sin_theta_i = sqrt(1 - SQUARE(vars.cos_theta.at(i)));
            double value = sin_theta_1 * cos(vars.phi.at(i)) * sin_theta_i +
                           cos_theta_1 * vars.cos_theta.at(i);
            value *= Q_1 * vars.magnitudes.at(i);

            bare_dot_prod.at(i).at(0) = value;
            bare_dot_prod.at(0).at(i) = value;
        }

        // Compute Q_i * Q_j for {i,j} != 1
        for (size_t i = 1; i < static_cast<size_t>(n_loops); ++i) {
            for (size_t j = 1; j < i; ++j) {
                double sin_theta_i = sqrt(1 - SQUARE(vars.cos_theta.at(i)));
                double sin_theta_j = sqrt(1 - SQUARE(vars.cos_theta.at(j)));

                double value = cos(vars.phi.at(i)) * sin_theta_i *
                                   cos(vars.phi.at(j)) * sin_theta_j +
                               sin(vars.phi.at(i)) * sin_theta_i *
                                   sin(vars.phi.at(j)) * sin_theta_j +
                               vars.cos_theta.at(i) * vars.cos_theta.at(j);

                value *= vars.magnitudes.at(i) * vars.magnitudes.at(j);
                bare_dot_prod.at(i).at(j) = value;
                bare_dot_prod.at(j).at(i) = value;
            }
        }
    }
}



void IntegrandTables::bs_compute_bare_dot_prod()
{
    int n_loops = loop_params.n_loops();
    size_t n_coeffs = loop_params.n_coeffs();

    size_t k_a_idx = n_coeffs - 1;
    size_t k_b_idx = n_coeffs - 2;

    /* Diagonal products correspond to Q1*Q1, etc. */
    for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
        bare_dot_prod.at(i).at(i) = SQUARE(vars.magnitudes.at(i));
    }
    bare_dot_prod.at(k_a_idx).at(k_a_idx) = SQUARE(k_a);
    bare_dot_prod.at(k_b_idx).at(k_b_idx) = SQUARE(k_b);

    double value = k_a * k_b * cos_ab;
    bare_dot_prod.at(k_a_idx).at(k_b_idx) = value;
    bare_dot_prod.at(k_b_idx).at(k_a_idx) = value;

    /* Products involving k_a and Q_i has the form k_a*Q_i*cos(theta_i) */
    for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
        value = k_a * vars.magnitudes.at(i) * vars.cos_theta.at(i);
        bare_dot_prod.at(k_a_idx).at(i) = value;
        bare_dot_prod.at(i).at(k_a_idx) = value;
    }
    /* Products involving k_b and Q_i (special case since k_b has azimuthal
     * angle 0)*/
    double sin_ab = sqrt(1 - SQUARE(cos_ab));
    for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
        double sin_theta_i = sqrt(1 - SQUARE(vars.cos_theta.at(i)));

        value = sin_ab * cos(vars.phi.at(i)) * sin_theta_i +
                cos_ab * vars.cos_theta.at(i);
        value *= k_b * vars.magnitudes.at(i);
        bare_dot_prod.at(k_b_idx).at(i) = value;
        bare_dot_prod.at(i).at(k_b_idx) = value;
    }

    if (n_loops > 1) {
        // Compute Q_i * Q_j
        for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
            for (size_t j = 0; j < i; ++j) {
                double sin_theta_i = sqrt(1 - SQUARE(vars.cos_theta.at(i)));
                double sin_theta_j = sqrt(1 - SQUARE(vars.cos_theta.at(j)));

                value = cos(vars.phi.at(i)) * sin_theta_i *
                            cos(vars.phi.at(j)) * sin_theta_j +
                        sin(vars.phi.at(i)) * sin_theta_i *
                            sin(vars.phi.at(j)) * sin_theta_j +
                        vars.cos_theta.at(i) * vars.cos_theta.at(j);

                value *= vars.magnitudes.at(i) * vars.magnitudes.at(j);
                bare_dot_prod.at(i).at(j) = value;
                bare_dot_prod.at(j).at(i) = value;
            }
        }
    }
}



void IntegrandTables::compute_bare_los_proj() {
    int n_loops = loop_params.n_loops();
    size_t n_coeffs = loop_params.n_coeffs();

    size_t k_a_idx = n_coeffs - 1;

    double sin_los_angle = sqrt(1 - SQUARE(mu_los_));

    /* k_a along L.o.S. is simply k_a * mu_los_ */
    bare_los_projection_.at(k_a_idx) = k_a * mu_los_;

    /* Products involving k_a and Q_i has the form k_a*Q_i*cos(theta_i) */
    for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
        double sin_theta_i = sqrt(1 - SQUARE(vars.cos_theta.at(i)));
        bare_los_projection_.at(i) = vars.magnitudes.at(i) * (
                sin_los_angle * sin_theta_i * cos(vars.phi.at(i)) +
                mu_los_ * vars.cos_theta.at(i)
                );
    }
}



void IntegrandTables::compute_comp_los_proj() {
    size_t n_coeffs = loop_params.n_coeffs();
    size_t n_configs = loop_params.n_configs();

    for (size_t a = 0; a < n_configs; ++a) {
        double product_value = 0;

        label2config(static_cast<int>(a), a_coeffs);
        for (size_t i = 0; i < n_coeffs; ++i) {
            product_value += a_coeffs[i] * bare_los_projection_.at(i);
        }
        comp_los_projection_.at(a) = product_value;
    }
}



/* Computes table of composite dot products given bare_dot_prod table */
void IntegrandTables::compute_comp_dot_prod()
{
    size_t n_coeffs = loop_params.n_coeffs();
    size_t n_configs = loop_params.n_configs();

    // Scalar product matrix is symmetric, hence compute indices [a,b] and
    // [b,a] simultaneously
    for (size_t a = 0; a < n_configs; ++a) {
        for (size_t b = 0; b <= a; ++b) {
            double product_value = 0;

            label2config(static_cast<int>(a), a_coeffs);
            label2config(static_cast<int>(b), b_coeffs);
            for (size_t i = 0; i < n_coeffs; ++i) {
                for (size_t j = 0; j < n_coeffs; ++j) {
                    product_value += a_coeffs[i] * b_coeffs[j]
                        * bare_dot_prod.at(i).at(j);
                }
            }
            if (a == b) {
                comp_dot_prod.at(a).at(a) = product_value;
            }
            else {
                comp_dot_prod.at(a).at(b) = product_value;
                comp_dot_prod.at(b).at(a) = product_value;
            }
        }
    }
}



void IntegrandTables::compute_alpha_beta()
{
    size_t n_configs = loop_params.n_configs();
    for (size_t a = 0; a < n_configs; ++a) {
        for (size_t b = 0; b < n_configs; ++b) {
            // Special case when a == b
            if (a == b) {
                alpha_.at(a).at(b) = 2.0;
                beta_.at(a).at(b)  = 2.0;
                continue;
            }

            double alpha_val = 0.0;
            double beta_val  = 0.0;

            // If the first argument is the zero-vector, alpha and beta remains 0
            // If the second argument is the zero-vector, beta remains 0
            if (a != static_cast<size_t>(loop_params.zero_label())) {
                double product_ab = comp_dot_prod.at(a).at(b);
                double product_aa = comp_dot_prod.at(a).at(a);

                alpha_val = 1 + product_ab/product_aa;

                if (b != static_cast<size_t>(loop_params.zero_label())) {
                    double product_bb = comp_dot_prod.at(b).at(b);

                    beta_val = product_ab / 2.0
                        * ( 1.0 / product_aa + 1.0 / product_bb
                          + 2.0 * product_ab / (product_aa * product_bb)
                          );
                }
            }
            alpha_.at(a).at(b) = alpha_val;
            beta_.at(a).at(b) = beta_val;
        }
    }
}



void IntegrandTables::compute_tables()
{
    if (loop_params.spectrum() == POWERSPECTRUM) {
        ps_compute_bare_dot_prod();
        if (loop_params.rsd()) {
            compute_bare_los_proj();
            compute_comp_los_proj();
        }
    }
    else if (loop_params.spectrum() == BISPECTRUM) {
        bs_compute_bare_dot_prod();
    }
    compute_comp_dot_prod();
    compute_alpha_beta();
}
