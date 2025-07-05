#include <iostream>
#include <algorithm>
#include <cmath>

#include "../include/parameters.hpp"
#include "../include/tables.hpp"

using std::size_t;

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif


/* Helper function for constructing the sum table */
int SumTable::convert_and_sum(int a, int b) {
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



/* Check that sum is an appropriate vector configuration, i.e. that
 * Q-coefficients are elements of (-1,0,1) and k-coefficient is an element
 * of (0,1) */
void SumTable::check_result(int res) const {
    Vec1D<int> res_coeffs(n_coeffs, 0);
    label2config(res, res_coeffs);
    for (size_t i = 0; i < n_coeffs; ++i) {
        int c = res_coeffs.at(i);
        if (!(c == -1 || c == 0 || c == 1)) {
            throw(std::logic_error(
                "SumTable::check_result(): Sum of labels does not correspond to "
                "an appropriate configuration."));
        }
    }
}



SumTable::SumTable(const LoopStructure& loop_structure) :
    zero_label(loop_structure.zero_label()),
    n_coeffs(loop_structure.n_coeffs())
{
    size_t n_configs = loop_structure.n_configs();

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
                    convert_and_sum(static_cast<int>(a),static_cast<int>(b));
            }
        }
    }
}




int SumTable::operator()(const int labels[], size_t size) const
{
    if (size == 1) return labels[0];

    int result = labels[0];
    for (size_t i = 1; i < size; ++i) {
        if (labels[i] == zero_label) continue;
        result = sum_table.at(static_cast<size_t>(result))
                     .at(static_cast<size_t>(labels[i]));
    }
#if DEBUG >= 1
    check_result(result);
#endif
    return result;
}



EtaGrid::EtaGrid(
        double eta_ini,
        double eta_fin,
        size_t time_steps,
        size_t pre_time_steps,
        double eta_asymp
       ) :
    eta_ini_(eta_ini), eta_fin_(eta_fin),
    time_steps_(time_steps), pre_time_steps_(pre_time_steps),
    eta_asymp_(eta_asymp)
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
        bool rsd,
        double rsd_growth_f,
        bool biased_tracers,
        const Vec1D<double>& bias_parameters,
        const Dynamics dynamics,
        const LoopStructure& loop_structure,
        const SumTable& sum_table,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid,
        const OmegaEigenspace& omega_eigenspace
        ) :
    k_a(k_a), k_b(k_b), cos_ab(cos_ab), rsd_f(rsd_growth_f),
    rsd(rsd), biased_tracers(biased_tracers), bias_parameters(bias_parameters),
    dynamics(dynamics), loop_structure(loop_structure),
    sum_table(sum_table), ev_params(ev_params), eta_grid(eta_grid),
    omega_eigenspace(omega_eigenspace),
    vars(IntegrationVariables(static_cast<size_t>(loop_structure.n_loops())))
{
    if (biased_tracers && !rsd) {
        throw std::logic_error(
            "IntegrandTables::IntegrandTables(): Biased tracers only implemented "
            "for rsd = true.");
    }

    int n_loops = loop_structure.n_loops();
    size_t n_coeffs = loop_structure.n_coeffs();
    size_t n_configs = loop_structure.n_configs();
    size_t n_kernels = loop_structure.n_kernels();

    a_coeffs.resize(n_coeffs);
    b_coeffs.resize(n_coeffs);

    bare_dot_prod.resize(n_coeffs, n_coeffs);
    composite_dot_prod.resize(n_configs, n_configs);
    alpha_.resize(n_configs, n_configs);
    beta_.resize(n_configs, n_configs);

    spt_kernels.resize(n_kernels, SPTKernel());
    kernels.resize(n_kernels,
                   Kernel(eta_grid.time_steps()));
    if (rsd) {
        bare_los_projection_.resize(n_coeffs);
        composite_los_projection_.resize(n_configs);

        rsd_kernels.resize(n_kernels, RSDKernel());
        vel_power_kernels.resize(n_kernels,
                                 2*static_cast<size_t>(n_loops) ,
                                 RSDKernel());
    }
}



/* Generic function for resetting kernels */
template <typename KernelType, typename ResetFn>
void reset_kernels(KernelType& kernels, ResetFn reset_fn) {
    for (auto& kernel : kernels) {
        kernel.computed = false;
        reset_fn(kernel);
    }
}



void IntegrandTables::reset() {
    switch (dynamics) {
        case EDS_SPT:
            reset_kernels(spt_kernels, [](auto &k) {
                std::fill(k.values.begin(), k.values.end(), 0.0);
            });
            break;
        case EVOLVE_ASYMPTOTIC_ICS:
            reset_kernels(kernels, [](auto &k) {
                std::fill(k.values.begin(), k.values.end(), 0.0);
            });
            break;
        case EVOLVE_EDS_ICS:
            reset_kernels(spt_kernels, [](auto &k) {
                std::fill(k.values.begin(), k.values.end(), 0.0);
            });
            reset_kernels(kernels, [](auto &k) {
                std::fill(k.values.begin(), k.values.end(), 0.0);
            });
            break;
        default:
            throw std::runtime_error("IntegrandTables::reset(): Invalid dynamics type.");
    }

    if (rsd) {
        reset_kernels(rsd_kernels, [](auto &k) {
            k.value = 0.0;
        });
        reset_kernels(vel_power_kernels, [](auto &k) {
            k.value = 0.0;
        });
    }
}



void IntegrandTables::compute_bare_dot_prod()
{
    int n_loops = loop_structure.n_loops();
    size_t n_coeffs = loop_structure.n_coeffs();
    size_t k_a_idx = n_coeffs - 1;

    /* Diagonal products correspond to Q1*Q1, etc. */
    for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
        bare_dot_prod(i,i) = SQUARE(vars.magnitudes.at(i));
    }
    bare_dot_prod(k_a_idx, k_a_idx) = SQUARE(k_a);

    double value = 0;
    /* k_a is chosen to point in the z-direction, hence products involving k_a
     * and Q_i has the form k_a*Q_i*cos(theta_i) */
    for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
        value = k_a * vars.magnitudes.at(i) * vars.cos_theta.at(i);
        bare_dot_prod(k_a_idx, i) = value;
        bare_dot_prod(i, k_a_idx) = value;
    }

    if (n_loops > 1) {
        /* Compute Q_i * Q_j */
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
                bare_dot_prod(i, j) = value;
                bare_dot_prod(j, i) = value;
            }
        }
    }

    if (loop_structure.spectrum() == BISPECTRUM) {
        size_t k_b_idx = n_coeffs - 2;

        double value = k_a * k_b * cos_ab;
        bare_dot_prod(k_a_idx, k_b_idx) = value;
        bare_dot_prod(k_b_idx, k_a_idx) = value;
        bare_dot_prod(k_b_idx, k_b_idx) = SQUARE(k_b);
        /* Products involving k_b and Q_i (special case since k_b has azimuthal
         * angle 0)*/
        double sin_ab = sqrt(1 - SQUARE(cos_ab));
        for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
            double sin_theta_i = sqrt(1 - SQUARE(vars.cos_theta.at(i)));

            value = sin_ab * cos(vars.phi.at(i)) * sin_theta_i +
                cos_ab * vars.cos_theta.at(i);
            value *= k_b * vars.magnitudes.at(i);
            bare_dot_prod(k_b_idx, i) = value;
            bare_dot_prod(i, k_b_idx) = value;
        }
    }
}



void IntegrandTables::compute_bare_los_proj() {
    int n_loops = loop_structure.n_loops();
    size_t n_coeffs = loop_structure.n_coeffs();

    size_t k_a_idx = n_coeffs - 1;

    double sin_los_angle = sqrt(1 - SQUARE(vars.mu_los));

    /* k_a along L.o.S. is simply k_a * mu_los_ */
    bare_los_projection_.at(k_a_idx) = k_a * vars.mu_los;

    /* Products involving k_a and Q_i has the form k_a*Q_i*cos(theta_i) */
    for (size_t i = 0; i < static_cast<size_t>(n_loops); ++i) {
        double sin_theta_i = sqrt(1 - SQUARE(vars.cos_theta.at(i)));
        bare_los_projection_.at(i) = vars.magnitudes.at(i) * (
                sin_los_angle * sin_theta_i * cos(vars.phi.at(i)) +
                vars.mu_los * vars.cos_theta.at(i)
                );
    }
}



void IntegrandTables::compute_composite_los_proj() {
    size_t n_coeffs = loop_structure.n_coeffs();
    size_t n_configs = loop_structure.n_configs();

    for (size_t a = 0; a < n_configs; ++a) {
        double product_value = 0;

        label2config(static_cast<int>(a), a_coeffs);
        for (size_t i = 0; i < n_coeffs; ++i) {
            product_value += a_coeffs.at(i) * bare_los_projection_.at(i);
        }
        composite_los_projection_.at(a) = product_value;
    }
}



/* Computes table of composite dot products given bare_dot_prod table */
void IntegrandTables::compute_composite_dot_prod()
{
    size_t n_coeffs = loop_structure.n_coeffs();
    size_t n_configs = loop_structure.n_configs();

    // Scalar product matrix is symmetric, hence compute indices [a,b] and
    // [b,a] simultaneously
    for (size_t a = 0; a < n_configs; ++a) {
        for (size_t b = 0; b <= a; ++b) {
            double product_value = 0;

            label2config(static_cast<int>(a), a_coeffs);
            label2config(static_cast<int>(b), b_coeffs);
            for (size_t i = 0; i < n_coeffs; ++i) {
                for (size_t j = 0; j < n_coeffs; ++j) {
                    product_value += a_coeffs.at(i) * b_coeffs.at(j)
                        * bare_dot_prod(i, j);
                }
            }
            if (a == b) {
                composite_dot_prod(a, a) = product_value;
            }
            else {
                composite_dot_prod(a, b) = product_value;
                composite_dot_prod(b, a) = product_value;
            }
        }
    }
}



void IntegrandTables::compute_alpha_beta()
{
    size_t n_configs = loop_structure.n_configs();
    for (size_t a = 0; a < n_configs; ++a) {
        for (size_t b = 0; b < n_configs; ++b) {
            // Special case when a == b
            if (a == b) {
                alpha_(a, b) = 2.0;
                beta_(a, b)  = 2.0;
                continue;
            }

            double alpha_val = 0.0;
            double beta_val  = 0.0;

            // If the first argument is the zero-vector, alpha and beta remains 0
            // If the second argument is the zero-vector, beta remains 0
            if (a != static_cast<size_t>(loop_structure.zero_label())) {
                double product_ab = composite_dot_prod(a, b);
                double product_aa = composite_dot_prod(a, a);

                alpha_val = 1 + product_ab/product_aa;

                if (b != static_cast<size_t>(loop_structure.zero_label())) {
                    double product_bb = composite_dot_prod(b, b);

                    beta_val = product_ab / 2.0
                        * ( 1.0 / product_aa + 1.0 / product_bb
                          + 2.0 * product_ab / (product_aa * product_bb)
                          );
                }
            }
            alpha_(a, b) = alpha_val;
            beta_(a, b) = beta_val;
        }
    }
}



void IntegrandTables::compute_tables()
{
    compute_bare_dot_prod();
    compute_composite_dot_prod();
    compute_alpha_beta();

    if (rsd) {
        compute_bare_los_proj();
        compute_composite_los_proj();
    }
}
