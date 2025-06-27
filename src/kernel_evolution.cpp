#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
    #include <cmath>
    #include <gsl/gsl_odeiv2.h>
    #include <gsl/gsl_errno.h>
}

#include "../include/utilities.hpp"
#include "../include/combinatorics.hpp"
#include "../include/interpolation.hpp"
#include "../include/parameters.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/tables.hpp"
#include "../include/kernel_evolution.hpp"
#include "../include/omega_matrix.hpp"


/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

using std::size_t;


void KernelEvolver::vertex(
    int m_l,
    int m_r,
    const int args_l[],
    const int args_r[],
    int sum_l,
    int sum_r,
    Strided2DVec<double>& partial_rhs_sum
) {
    size_t sl = static_cast<size_t>(sum_l);
    size_t sr = static_cast<size_t>(sum_r);
    double alpha_lr = tables.alpha()(sl, sr);
    double alpha_rl = tables.alpha()(sr, sl);
    double beta     = tables.beta()(sl, sr);
    double kappa    = 0.5 * (alpha_lr + alpha_rl) - beta;

    size_t index_l = static_cast<size_t>(
        compute(args_l, -1, m_l));
    size_t index_r = static_cast<size_t>(
        compute(args_r, -1, m_r));
    const auto& kernel_l =
        tables.kernels.at(index_l).values;
    const auto& kernel_r =
        tables.kernels.at(index_r).values;

    const Interpolation1D& mu2 = tables.ev_params.mu2();
    const EtaGrid& eta_grid = tables.eta_grid;

    auto update_step = [&](size_t i, double mu2_val) {
        if constexpr (COMPONENTS >= 4) {
            partial_rhs_sum(i, 2) = alpha_lr * kernel_l(i, 3) * kernel_r(i, 2);
            partial_rhs_sum(i, 3) = beta     * kernel_l(i, 3) * kernel_r(i, 3);
        }
        if constexpr (COMPONENTS >= 2) {
            partial_rhs_sum(i, 0) = alpha_lr * kernel_l(i, 1) * kernel_r(i, 0);
            partial_rhs_sum(i, 1) = (
                beta * kernel_l(i, 1) * kernel_r(i, 1) +
                mu2_val * kappa * kernel_l(i, 0) * kernel_r(i, 0)
            );
        }
        static_assert(COMPONENTS == 2 || COMPONENTS == 4,
                      "vertex only implemented for COMPONENTS = 2 or 4");

    };

    /* eta_asymp to eta_ini */
    size_t i = 0;
    if (tables.dynamics == EVOLVE_ASYMPTOTIC_ICS) {
        /* Use fixed mu22_val evaluated at eta_ini */
        double mu2_val = mu2(eta_grid.eta_ini());
        for (; i < eta_grid.pre_time_steps(); ++i) {
            update_step(i, mu2_val);
        }
    }

    /* eta_ini to eta_fin */
    for (; i < eta_grid.time_steps(); ++i) {
        double mu2_val = mu2(eta_grid.at(i));
        update_step(i, mu2_val);
    }
}



void KernelEvolver::vertex_cubic(
    int n_a,
    int n_b,
    int n_c,
    const int args_a[],
    const int args_b[],
    const int args_c[],
    int sum_a,
    int sum_b,
    int sum_c,
    int sum_ab,
    Strided2DVec<double>& partial_rhs_sum
) {
    size_t sa = static_cast<size_t>(sum_a);
    size_t sb = static_cast<size_t>(sum_b);
    size_t sc = static_cast<size_t>(sum_c);
    size_t sab = static_cast<size_t>(sum_ab);

    // gamma(a, b)
    double alpha_a_b = tables.alpha()(sa, sb);
    double alpha_b_a = tables.alpha()(sb, sa);
    double beta_a_b   = tables.beta()(sa, sb);
    double kappa_a_b = 0.5 * (alpha_a_b + alpha_b_a) - beta_a_b;

    // gamma(c, ab)
    double alpha_c_ab = tables.alpha()(sc, sab);
    double alpha_ab_c = tables.alpha()(sab, sc);
    double beta_c_ab   = tables.beta()(sc, sab);
    double kappa_c_ab = 0.5 * (alpha_c_ab + alpha_ab_c) - beta_c_ab;

    size_t index_a = static_cast<size_t>(compute(args_a, -1, n_a));
    size_t index_b = static_cast<size_t>(compute(args_b, -1, n_b));
    size_t index_c = static_cast<size_t>(compute(args_c, -1, n_c));

    const auto& kernel_a = tables.kernels.at(index_a).values;
    const auto& kernel_b = tables.kernels.at(index_b).values;
    const auto& kernel_c = tables.kernels.at(index_c).values;

    const Interpolation1D& mu22 = tables.ev_params.mu22();
    const EtaGrid& eta_grid = tables.eta_grid;

    double factor = kappa_a_b * kappa_c_ab;

    /* eta_asymp to eta_ini */
    size_t i = 0;
    if (tables.dynamics == EVOLVE_ASYMPTOTIC_ICS) {
        /* Use fixed mu22_val evaluated at eta_ini */
        double mu22_val = mu22(eta_grid.eta_ini());
        for (; i < eta_grid.pre_time_steps(); ++i) {
            partial_rhs_sum(0, 1) +=
                factor * mu22_val *
                kernel_a(0, 0) * kernel_b(0, 0) * kernel_c(0, 0);
        }
    }

    /* eta_ini to eta_fin */
    for (; i < eta_grid.time_steps(); ++i) {
        // Poisson modification of theta, i.e. component=1
        partial_rhs_sum(i, 1) +=
            factor * mu22(eta_grid.at(i)) *
            kernel_a(i, 0) * kernel_b(i, 0) * kernel_c(i, 0);
    }
}



/* Helper functions for RHS computation */
inline void zero_args(int* args, size_t count, int zero_label) {
    std::fill(args, args + count, zero_label);
}



inline void accumulate_rhs(
    const Strided2DVec<double>& partial,
    std::array<std::vector<double>, COMPONENTS>& rhs,
    int a,
    int b
) {
    size_t a_t = static_cast<size_t>(a);
    size_t b_t = static_cast<size_t>(b);
    for (size_t i = 0; i < COMPONENTS; ++i) {
        for (size_t j = 0; j < partial.rows(); ++j) {
            /* rhs is COMPONENTS x time_steps, partial is time_steps x
             * COMPONENTS because this is faster in vertex, looping over
             * time_steps */
            rhs[i][j] += partial(j, i) / binomial_coeffs[a_t][b_t];
        }
    }
}



inline void accumulate_partial(
    const Strided2DVec<double>& in,
    Strided2DVec<double>& out,
    int a,
    int b
) {
    size_t a_t = static_cast<size_t>(a);
    size_t b_t = static_cast<size_t>(b);
    const double factor = 1.0 / binomial_coeffs[a_t][b_t];
    std::transform(
        in.begin(), in.end(),
        out.begin(),
        out.begin(),
        [factor](double val_in, double val_out) {
            return val_out + val_in * factor;
        }
    );
}



void KernelEvolver::RHS(
    const int arguments[],
    int n,
    std::array<std::vector<double>, COMPONENTS>& rhs /* out, added to */
) {
    const size_t time_steps = tables.eta_grid.time_steps();
    const size_t n_kernel_args = tables.loop_structure.n_kernel_args();
    const int zero_label = tables.loop_structure.zero_label();

    /* partial_rhs_sum is instead a flatted vector for (potential) speedup,
     * use get(partial_rhs_sum, i=time_step, j=component, stride=COMPONENTS)
     * to access */
    Strided2DVec<double> partial_rhs_sum(time_steps, COMPONENTS);

    int args_l[N_KERNEL_ARGS_MAX] = {0};
    int args_r[N_KERNEL_ARGS_MAX] = {0};

    for (int m = 1; m <= n / 2; ++m) {
        zero_args(args_l, n_kernel_args, zero_label);
        zero_args(args_r, n_kernel_args, zero_label);

        std::fill(partial_rhs_sum.begin(), partial_rhs_sum.end(), 0.0);

        /* Go through all ways to pick m (unordered) elements from group of n */
        Combinations comb(n, m);
        do {
            /* Set args_l and args_r from current combination and complement,
             * respectively */
            comb.rearrange_from_current_combination(arguments,
                                                    args_l,
                                                    static_cast<size_t>(m));
            comb.rearrange_from_current_complement(arguments,
                                                   args_r,
                                                   static_cast<size_t>(n-m));

            int sum_l = tables.sum_table(args_l, n_kernel_args);
            int sum_r = tables.sum_table(args_r, n_kernel_args);

            apply_symmetric_vertex(m, n - m, args_l, args_r, sum_l,
                                   sum_r, partial_rhs_sum);
        } while (comb.next());

        accumulate_rhs(partial_rhs_sum, rhs, n, m);
    }
}



void KernelEvolver::RHS_cubic(
    const int arguments[],
    int n,
    std::array<std::vector<double>, COMPONENTS>& rhs /* out, adds to input */
) {
    const size_t time_steps = tables.eta_grid.time_steps();
    const size_t n_kernel_args = tables.loop_structure.n_kernel_args();
    int zero_label = tables.loop_structure.zero_label();

    Strided2DVec<double> outer_partial_rhs_sum(time_steps, COMPONENTS);
    Strided2DVec<double> inner_partial_rhs_sum(time_steps, COMPONENTS);

    int args_a[N_KERNEL_ARGS_MAX] = {0};
    int args_b[N_KERNEL_ARGS_MAX] = {0};
    int args_c[N_KERNEL_ARGS_MAX] = {0};
    int args_l[N_KERNEL_ARGS_MAX] = {0};
    int args_r[N_KERNEL_ARGS_MAX] = {0};

    for (int m = 1; m <= n / 2; ++m) {
        zero_args(args_a, n_kernel_args, zero_label);
        zero_args(args_l, n_kernel_args, zero_label);
        zero_args(args_r, n_kernel_args, zero_label);
        std::fill(outer_partial_rhs_sum.begin(), outer_partial_rhs_sum.end(), 0.0);

        Combinations comb(n, m);
        do {
            comb.rearrange_from_current_combination(arguments, args_l,
                                                    static_cast<size_t>(m));
            comb.rearrange_from_current_complement(arguments, args_r,
                                                   static_cast<size_t>(n - m));

            zero_args(args_a, n_kernel_args, zero_label);
            std::copy(args_l, args_l + m, args_a);
            int sum_a = tables.sum_table(args_a, n_kernel_args);

            for (int p = 1; p <= (n-m)/2; ++p) {
                zero_args(args_b, n_kernel_args, zero_label);
                zero_args(args_c, n_kernel_args, zero_label);
                std::fill(inner_partial_rhs_sum.begin(), inner_partial_rhs_sum.end(), 0.0);

                Combinations comb_r(n-m, p);
                do {
                    comb_r.rearrange_from_current_combination(args_r, args_b,
                                                            static_cast<size_t>(p));
                    comb_r.rearrange_from_current_complement(args_r, args_c,
                                                           static_cast<size_t>(n-m-p));
                    int sum_b = tables.sum_table(args_b, n_kernel_args);
                    int sum_c = tables.sum_table(args_c, n_kernel_args);

                    apply_symmetric_vertex_cubic(m, p, n-m-p, args_a, args_b,
                                                 args_c, sum_a, sum_b, sum_c,
                                                 inner_partial_rhs_sum);
                }
                while (comb_r.next());
                accumulate_partial(inner_partial_rhs_sum, outer_partial_rhs_sum, n-m, p);
            }

            if (m != n - m) {
                zero_args(args_a, n_kernel_args, zero_label);
                std::copy(args_r, args_r + n - m, args_a);
                int sum_a = tables.sum_table(args_a, n_kernel_args);

                for (int p = 1; p <= m/2; ++p) {
                    zero_args(args_b, n_kernel_args, zero_label);
                    zero_args(args_c, n_kernel_args, zero_label);

                    Combinations comb_r(m, p);
                    do {
                        comb_r.rearrange_from_current_combination(args_l, args_b,
                                                                static_cast<size_t>(p));
                        comb_r.rearrange_from_current_complement(args_l, args_c,
                                                               static_cast<size_t>(m-p));
                        int sum_b = tables.sum_table(args_b, n_kernel_args);
                        int sum_c = tables.sum_table(args_c, n_kernel_args);

                        apply_symmetric_vertex_cubic(n-m, p, m-p, args_a, args_b,
                                                     args_c, sum_a, sum_b, sum_c,
                                                     inner_partial_rhs_sum);
                    }
                    while (comb_r.next());
                    accumulate_partial(inner_partial_rhs_sum, outer_partial_rhs_sum, m, p);
                }
            }
        }
        while (comb.next());
        accumulate_rhs(outer_partial_rhs_sum, rhs, n, m);
    }
}



int KernelEvolver::ode_system(double eta, const double y[], double f[], void* ode_input) {
    ODEParameters params = *(static_cast<ODEParameters*>(ode_input));

    int n = params.n;
    double k = params.k;
    Strided2DVec<double>& omega = params.omega;

    update_omega_matrix(eta, k, params.ev_params.kappa(),
        params.ev_params.zeta(), params.ev_params.xi(), omega);

    /* If n == 1, rhs = 0 */
    if (params.n == 1){
        for (size_t i = 0; i < COMPONENTS; ++i) {
            f[i] = 0;
        }
    }
    /* Otherwise, interpolate */
    else {
        for (size_t i = 0; i < COMPONENTS; ++i) {
            f[i] = params.rhs.at(i)(eta);
        }
    }

    for (size_t i = 0; i < COMPONENTS; ++i) {
        /* We rescale the kernels by exp(eta*n) so that the diagonal always
        * includes a term n */
        f[i] -= n * y[i];
        for (std::size_t j = 0; j < COMPONENTS; ++j) {
            f[i] -= omega(i, j) * y[j];
        }
    }
    return GSL_SUCCESS;
}



/* For obtaining the correct initial conditions at eta_ini, evolve system from
 * eta_asymp to eta_ini fixed omega-matrix (i.e. no call to
 * update_omega_matrix) */
int KernelEvolver::ode_system_fixed_eta(double eta, const double y[], double f[], void* ode_input) {
    ODEParameters params = *(static_cast<ODEParameters*>(ode_input));
    int n = params.n;
    Strided2DVec<double>& omega = params.omega;

    /* If n == 1, rhs = 0 */
    if (params.n == 1){
        for (size_t i = 0; i < COMPONENTS; ++i) {
            f[i] = 0;
        }
    }
    /* Otherwise, interpolate */
    else {
        for (size_t i = 0; i < COMPONENTS; ++i) {
            f[i] = params.rhs.at(i)(eta);
        }
    }

    for (size_t i = 0; i < COMPONENTS; ++i) {
        /* We rescale the kernels by exp(eta*n) so that the diagonal always
        * includes a term n */
        f[i] -= n * y[i];
        for (std::size_t j = 0; j < COMPONENTS; ++j) {
            f[i] -= omega(i, j) * y[j];
        }
    }
    return GSL_SUCCESS;
}



void KernelEvolver::set_EdS_ICs(
        const int arguments[],
        int kernel_index,
        int n,
        double k
        )
{
    UNUSED(n);
    UNUSED(k);

    compute_SPT_kernels(arguments, kernel_index, n, tables);

    auto& kernel_vec =
        tables.kernels.at(static_cast<size_t>(kernel_index)).values;
    auto& spt_kernel_vec =
        tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values;

    for (size_t i = 0; i < EDS_SPT_COMPONENTS; ++i) {
        kernel_vec(0, i) = spt_kernel_vec.at(i);
    }
}



void KernelEvolver::set_asymptotic_ICs(
        const int arguments[],
        int kernel_index,
        int n,
        double k
        )
{
    UNUSED(arguments);
    auto& kernel_vec = tables.kernels.at(static_cast<size_t>(kernel_index)).values;
    // Use growing mode IC for F1 at eta_asymp
    if (n == 1) {
        double eigenvalue = tables.omega_eigenspace.eigenvalue()(k);
        double delta_eta = tables.eta_grid.eta_ini() - tables.eta_grid.eta_asymp();
        double F0 = tables.omega_eigenspace.eigenvectors().at(0)(k);
        kernel_vec(0,0) = exp(-eigenvalue * delta_eta);
        for (size_t i = 1; i < COMPONENTS; ++i) {
            /* Normalize by eigenvectors by F0 (PT expansion in first component) */
            double Fi = tables.omega_eigenspace.eigenvectors().at(i)(k) / F0;
            kernel_vec(0, i) = exp(-eigenvalue * delta_eta) * Fi;
        }
    }
    else {
        // Set F(n>1)-kernels to zero at eta_asymp
        std::fill(kernel_vec.begin(), kernel_vec.end(), 0.0);
    }
}



void KernelEvolver::evolve(
        int kernel_index,
        int n,
        double k
    )
{
    const EvolutionParameters& ev_params = tables.ev_params;
    const EtaGrid& eta_grid              = tables.eta_grid;

    const size_t time_steps     = eta_grid.time_steps();
    const size_t pre_time_steps = eta_grid.pre_time_steps();
    const double eta_ini        = eta_grid.eta_ini();
    double eta_current          = eta_grid.at(0);

    auto& k_vec = tables.kernels.at(static_cast<size_t>(kernel_index)).values;

    auto evolve_one_step = [&](size_t i) {
        /* For calculation of kernel_row.at(i), the initial condition is
         * kernel_row.at(i-1). Hence we copy kernel_row.at(i-1) to
         * kernel_row.at(i) */
        if constexpr (COMPONENTS == 2) {
            k_vec(i, 0) = k_vec(i - 1, 0);
            k_vec(i, 1) = k_vec(i - 1, 1);
        } else if constexpr (COMPONENTS == 4) {
            k_vec(i, 0) = k_vec(i - 1, 0);
            k_vec(i, 1) = k_vec(i - 1, 1);
            k_vec(i, 2) = k_vec(i - 1, 2);
            k_vec(i, 3) = k_vec(i - 1, 3);
        } else {
            // Fallback
            std::copy_n(&k_vec(i - 1, 0), COMPONENTS,
                        &k_vec(i, 0));
        }
        /* Then evolve to index i */
        int status = gsl_odeiv2_driver_apply(
            driver, &eta_current, eta_grid.at(i), &k_vec(i, 0)
        );
        if (status != GSL_SUCCESS) {
            throw std::runtime_error("GSL ODE driver gave error value = " +
                                     std::to_string(status));
        }
    };

    size_t i = 1;
    if (tables.dynamics == EVOLVE_ASYMPTOTIC_ICS) {
        if (n == 1) {
            std::array<double, COMPONENTS> k_ini;
            if constexpr (COMPONENTS == 2) {
                k_ini.at(0) = k_vec(0, 0);
                k_ini.at(1) = k_vec(0, 1);
            } else if constexpr (COMPONENTS == 4) {
                k_ini.at(0) = k_vec(0, 0);
                k_ini.at(1) = k_vec(0, 1);
                k_ini.at(2) = k_vec(0, 2);
                k_ini.at(3) = k_vec(0, 3);
            }
            else {
                // fallback
                std::copy_n(&k_vec(0, 0), COMPONENTS,
                            &k_ini[0]);
            }
            for (; i < pre_time_steps + 1; ++i) {
                double factor = std::exp(
                    tables.omega_eigenspace.eigenvalue()(k) *
                    (eta_grid.at(i) - eta_ini));
                for (size_t j = 0; j < COMPONENTS; ++j) {
                    k_vec(i, j) = k_ini[j] * factor;
                }
            }
            eta_current = eta_grid.at(pre_time_steps);
        }
        else {
            /* Evolve from asymptotic eta = eta_asymp to initial eta = eta_ini */
            update_omega_matrix(eta_grid.eta_ini(), k, ev_params.kappa(),
                                ev_params.zeta(), ev_params.xi(), omega);

            sys.function = &ode_system_fixed_eta;
            for (; i < pre_time_steps + 1; ++i) {
                evolve_one_step(i);
            }
        }
    }
    /* eta_ini to eta_fin */
    sys.function = &ode_system;
    for (; i < time_steps; ++i) {
        evolve_one_step(i);
    }
}



KernelEvolver::KernelEvolver(IntegrandTables& tables)
    : tables(tables), omega(COMPONENTS, COMPONENTS)
{
    sys = {ode_system, nullptr, COMPONENTS, nullptr};
    driver = gsl_odeiv2_driver_alloc_y_new( &sys,
            ODE_ROUTINE,
            tables.ev_params.ode_hstart(),
            tables.ev_params.ode_rtol(),
            tables.ev_params.ode_atol()
            );
    if (tables.dynamics == EVOLVE_EDS_ICS) {
        set_ICs = [this](const int arguments[], int kernel_index, int n, double k) {
            return this->set_EdS_ICs(arguments, kernel_index, n, k);
        };
    }
    else if (tables.dynamics == EVOLVE_ASYMPTOTIC_ICS) {
        set_ICs = [this](const int arguments[], int kernel_index, int n, double k) {
            return this->set_asymptotic_ICs(arguments, kernel_index, n, k);
        };
    }
    else {
        throw std::runtime_error(
            "KernelEvolver::KernelEvolver(): KernelEvolver only allows "
            "dynamics = EVOLVE_IC_ASYMP or EVOLVE_ASYMPTOTIC_ICS.");
    }
}



int KernelEvolver::compute(
        const int arguments[],
        int kernel_index,
        int n
        )
{
#if DEBUG >= 1
    kernel_computer_validate_n(arguments, n, tables);
    kernel_computer_validate_kernel_index(arguments, kernel_index, tables);
#endif

    /* If kernel_index is not known, -1 is sent as argument */
    if (kernel_index == -1) {
        kernel_index = tables.loop_structure.args_to_kernel_index(arguments);
    }

    /* Alias reference to kernel we are working with for convenience/readability */
    Kernel& kernel = tables.kernels.at(static_cast<size_t>(kernel_index));

    /* Check if the kernels are already computed */
    if (kernel.computed) return kernel_index;

    /* Compute RHS sum in evolution equation if n > 1. If n == 1, the RHS
     * equals 0 */
    std::array<Interpolation1D, COMPONENTS> rhs;
    if (n > 1) {
        std::array<std::vector<double>, COMPONENTS> rhs_sum;
        for (auto& vec : rhs_sum) {
            vec.resize(tables.eta_grid.time_steps(), 0.0);
        }

        RHS(arguments, n, rhs_sum);
        RHS_cubic(arguments, n, rhs_sum);
        for (size_t i = 0; i < COMPONENTS; ++i) {
            rhs.at(i) =
                Interpolation1D(tables.eta_grid.grid(), rhs_sum.at(i));
        }
    }

    /* Compute k (here: sum of kernel arguments) */
    size_t sum = static_cast<size_t>(tables.sum_table(arguments,
        tables.loop_structure.n_kernel_args()));
    double k = std::sqrt( tables.composite_dot_products()(sum,sum));

    /* Set up ODE input and system */
    ODEParameters params(n, k, tables.eta_grid.eta_ini(),
        tables.ev_params, rhs, omega);
    sys.params = &params;
    /* Set initial conditions */
    set_ICs(arguments, kernel_index, n, k);
    /* Solve ODE */
    evolve(kernel_index, n, k);

    tables.kernels.at(static_cast<size_t>(kernel_index)).computed = true;
    return kernel_index;
}
