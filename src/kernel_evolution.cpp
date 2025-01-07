#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
    #include <gsl/gsl_sf.h>
    #include <gsl/gsl_odeiv2.h>
}

#include "../include/utilities.hpp"
#include "../include/combinatorics.hpp"
#include "../include/interpolation.hpp"
#include "../include/parameters.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/tables.hpp"
#include "../include/kernel_evolution.hpp"


/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

using std::size_t;

void vertex(
        int m_l,
        int m_r,
        const int args_l[],
        const int args_r[],
        int sum_l,
        int sum_r,
        Vec2D<double>& partial_rhs_sum,
        IntegrandTables& tables
        )
{
    double alpha_lr = tables.alpha()
        .at(static_cast<size_t>(sum_l))
        .at(static_cast<size_t>(sum_r));
    double beta = tables.beta()
        .at(static_cast<size_t>(sum_l))
        .at(static_cast<size_t>(sum_r));

    int index_l = compute_gen_kernels(args_l, -1, m_l, tables);
    int index_r = compute_gen_kernels(args_r, -1, m_r, tables);

    size_t a, b, c;

    // Note: all components are zero-indexed
    for (size_t i = 0; i < partial_rhs_sum.size(); ++i) {
        switch (COMPONENTS) {
            case 4:
                /* Component a = 2 */
                a = 2, b = 3, c = 2;
                /* gamma_223 = alpha_lr */
                partial_rhs_sum.at(i).at(a) +=
                    alpha_lr *
                    tables.kernels.at(static_cast<size_t>(index_l))
                        .values.at(i)
                        .at(b) *
                    tables.kernels.at(static_cast<size_t>(index_r))
                        .values.at(i)
                        .at(c);

                /* The term gamma_232 is redundant; due to momentum symmetrization
                 * of the evolution eq. RHS, it will be covered by the term
                 * above. */

                // Component a = 3
                a = 3, b = 3, c = 3;
                // gamma_444 = beta
                partial_rhs_sum.at(i).at(a) +=
                    beta *
                    tables.kernels.at(static_cast<size_t>(index_l))
                        .values.at(i)
                        .at(b) *
                    tables.kernels.at(static_cast<size_t>(index_r))
                        .values.at(i)
                        .at(c);

                /* For case 4 also code in case 2 should be executed */
                /* fall through */
            case 2:
                /* Component a = 0 */
                a = 0, b = 1, c = 0;
                /* gamma_001 = alpha_lr */
                partial_rhs_sum.at(i).at(a) +=
                    alpha_lr *
                    tables.kernels.at(static_cast<size_t>(index_l))
                        .values.at(i)
                        .at(b) *
                    tables.kernels.at(static_cast<size_t>(index_r))
                        .values.at(i)
                        .at(c);

                /* The term gamma_010 is redundant; due to momentum symmetrization
                 * of the evolution eq. RHS, it will be covered by the term
                 * above. */

                /* Component a = 1 */
                a = 1, b = 1, c = 1;
                /* gamma_111 = beta */
                partial_rhs_sum.at(i).at(a) +=
                    beta *
                    tables.kernels.at(static_cast<size_t>(index_l))
                        .values.at(i)
                        .at(b) *
                    tables.kernels.at(static_cast<size_t>(index_r))
                        .values.at(i)
                        .at(c);

                break;
            default:
                throw(std::logic_error("No vertex implemented for components = " +
                            std::to_string(COMPONENTS)));
        }
    }
}



void compute_RHS_sum(
        const int arguments[],
        int n,
        IntegrandTables& tables,
        std::array<Interpolation1D, COMPONENTS>& rhs /* out */
        )
{
    size_t time_steps    = tables.eta_grid.time_steps();
    size_t n_kernel_args = tables.loop_params.n_kernel_args();

    Vec2D<double> rhs_sum;
    Vec2D<double> partial_rhs_sum;

    /* Note different dimensions of rhs_sum and partial_rhs_sum. rhs_sum is
     * interpolated for each component, and needs to be a vector for each
     * component. partial_rhs_sum is opposite due to overhead reduction */
    rhs_sum.assign(COMPONENTS, Vec1D<double>(time_steps, 0.0));
    partial_rhs_sum.assign(time_steps, Vec1D<double>(COMPONENTS));

    int args_l[N_KERNEL_ARGS_MAX] = {0};
    int args_r[N_KERNEL_ARGS_MAX] = {0};

    for (int m = 1; m <= n/2; ++m) {
        /* Initialize partial_rhs_sum to 0 */
        for (auto& el : partial_rhs_sum) {
            std::fill(el.begin(), el.end(), 0.0);
        }

        /* Initialize args_l and args_r */
        std::fill(&args_l[0], &args_l[n_kernel_args],
                tables.loop_params.zero_label());
        std::fill(&args_r[0], &args_r[n_kernel_args],
                tables.loop_params.zero_label());

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

            vertex(m, n-m, args_l, args_r, sum_l, sum_r, partial_rhs_sum,
                    tables);

            // When m != (n - m), we may additionally compute the (n-m)-term by
            // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
            // compute_RHS_sum() only needs to sum up to (including) floor(n/2).
            if (m != n - m) {
                vertex(n - m, m, args_r, args_l, sum_r, sum_l, partial_rhs_sum,
                        tables);
            }
        } while (comb.next());

        // Devide through by symmetrization factor (n choose m)
        int n_choose_m = static_cast<int>(
                gsl_sf_choose(static_cast<unsigned int>(n),
                    static_cast<unsigned int>(m))
                );
        for (size_t i = 0; i < COMPONENTS; ++i) {
            for (size_t j = 0; j < time_steps; ++j) {
                rhs_sum.at(i).at(j) += partial_rhs_sum.at(j).at(i) / n_choose_m;
            }
        }
    }

    for (size_t i = 0; i < COMPONENTS; ++i) {
        rhs.at(i) = Interpolation1D(tables.eta_grid.grid(), rhs_sum.at(i));
    }
}



/* Structure which compiles user input to GSL ODE */
struct ODEInput {
    public:
        const int n;
        const double k;
        const double eta_ini;
        const EvolutionParameters& ev_params;
        const std::array<Interpolation1D, COMPONENTS>& rhs;

        ODEInput(int n,
                double k,
                double eta_ini,
                const EvolutionParameters &ev_params,
                const std::array<Interpolation1D, COMPONENTS>& rhs
                )
            : n(n), k(k), eta_ini(eta_ini), ev_params(ev_params), rhs(rhs) {}
};



namespace EvolveICEdS {
int kernel_gradient(double eta, const double y[], double f[], void *ode_input) {
    ODEInput input = *(static_cast<ODEInput*>(ode_input));
    int n = input.n;

    /* Interpolate RHS */
    double rhs[COMPONENTS] = {0.0};

    // If n == 1, rhs = 0
    if (input.n > 1) {
        for (size_t i = 0; i < COMPONENTS; ++i) {
            rhs[i] = input.rhs.at(i)(eta);
        }
    }

    double zeta = input.ev_params.zeta_at_eta(eta);

    /* etaD parametrization */
    f[0] = rhs[0] - n * y[0] + y[1];
    f[1] = rhs[1] + zeta * y[0] + (-zeta + 1 - n) * y[1];

    return GSL_SUCCESS;
}



void solve_kernel_ODE(
        ODEInput& input,
        const EtaGrid& eta_grid,
        Vec2D<double>& kernels /* time_steps*components table of kernels */
        )
{
    const EvolutionParameters& ev_params = input.ev_params;
    size_t time_steps     = eta_grid.time_steps();

    gsl_odeiv2_system sys = {kernel_gradient, nullptr, COMPONENTS, &input};
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkf45, ev_params.ode_hstart(),
        ev_params.ode_rtol(), ev_params.ode_atol());

    double eta_current = eta_grid.at(0);

    // Evolve system
    for (size_t i = 1; i < time_steps; i++) {
        /* For calculation of kernels.at(i), the initial condition is
         * kernels.at(i-1). Hence we copy kernels.at(i-1) to kernels.at(i) */
        std::copy(kernels.at(i - 1).begin(), kernels.at(i - 1).end(),
                kernels.at(i).begin());
        /* Then evolve to index i */
        int status = gsl_odeiv2_driver_apply(
                driver, &eta_current, eta_grid.at(i), kernels.at(i).data());

        if (status != GSL_SUCCESS) {
            throw(std::runtime_error("GLS ODE driver gave error value = " +
                        std::to_string(status)));
        }
    }

    // Free GSL ODE driver
    gsl_odeiv2_driver_free(driver);
}



static void kernel_initial_conditions(
        int kernel_index,
        IntegrandTables& tables
        )
{
    for (size_t i = 0; i < EDS_SPT_COMPONENTS; ++i) {
        tables.kernels.at(static_cast<size_t>(kernel_index)).values.at(0).at(i) =
            tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[i];
    }
}
}



namespace EvolveICAsymp {
int kernel_gradient(double eta, const double y[], double f[], void *ode_input) {
    ODEInput input = *(static_cast<ODEInput*>(ode_input));

    const EvolutionParameters& ev_params = input.ev_params;
    int n       = input.n;
    double f_nu = ev_params.f_nu();
    double k    = input.k;

    // Between eta_asymp and eta_i, set eta = eta_i
    if (eta < input.eta_ini) {
        eta = input.eta_ini;
    }

    /* Interpolate RHS */
    double rhs[COMPONENTS] = {0.0};

    // If n == 1, rhs = 0
    if (input.n > 1) {
        for (size_t i = 0; i < COMPONENTS; ++i) {
            rhs[i] = input.rhs.at(i)(eta);
        }
    }

    double zeta = ev_params.zeta_at_eta(eta);
    double cs2  = ev_params.cs2(eta, k);

    /* etaD parametrization */
    f[0] = rhs[0] - n * y[0] + y[1];
    f[1] = rhs[1] + zeta * (1 - f_nu) * y[0] + (-zeta + 1 - n) * y[1] +
           zeta * f_nu * y[2];
    f[2] = rhs[2] - n * y[2] + y[3];
    f[3] = rhs[3] + zeta * (1 - f_nu) * y[0] +
           zeta * (f_nu - k * k * cs2) * y[2] + (-zeta + 1 - n) * y[3];

    return GSL_SUCCESS;
}



void solve_kernel_ODE(
        ODEInput& input,
        const EtaGrid& eta_grid,
        Vec2D<double>& kernels /* time_steps*components table of kernels */
        )
{
    const EvolutionParameters& ev_params = input.ev_params;
    size_t time_steps     = eta_grid.time_steps();
    size_t pre_time_steps = eta_grid.pre_time_steps();

    gsl_odeiv2_system sys = {kernel_gradient, nullptr, COMPONENTS, &input};
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkf45, ev_params.ode_hstart(),
        ev_params.ode_rtol(), ev_params.ode_atol());

    // For kernels with n > 1 use ODE SOLVER
    if (input.n > 1) {
        double eta_current = eta_grid.at(0);

        // Evolve system
        for (size_t i = 1; i < time_steps; i++) {
            /* For calculation of kernels.at(i), the initial condition is
             * kernels.at(i-1). Hence we copy kernels.at(i-1) to kernels.at(i) */
            std::copy(kernels.at(i - 1).begin(), kernels.at(i - 1).end(),
                      kernels.at(i).begin());
            /* Then evolve to index i */
            int status = gsl_odeiv2_driver_apply(
                driver, &eta_current, eta_grid.at(i), kernels.at(i).data());

            if (status != GSL_SUCCESS) {
                throw(std::runtime_error("GLS ODE driver gave error value = " +
                                        std::to_string(status)));
            }
        }
    }
    // For n == 1 kernels, we may multipliy by exp("growing mode" eigenvalue)
    // before eta_I
    else {
        for (size_t i = 1; i < pre_time_steps + 1; ++i) {
            for (size_t j = 0; j < COMPONENTS; ++j) {
              kernels.at(i).at(j) =
                  kernels.at(0).at(j) *
                  std::exp(ev_params.omega_eigenvalues_at_k(input.k) *
                           (eta_grid.at(i) - eta_grid.at(0)));
            }
        }

        double eta_current = eta_grid.at(pre_time_steps);

        for (size_t i = pre_time_steps + 1; i < time_steps; i++) {
            /* For calculation of kernels.at(i), the initial condition is
             * kernels.at(i-1). Hence we copy kernels.at(i-1) to kernels.at(i) */
            std::copy(kernels.at(i - 1).begin(), kernels.at(i - 1).end(),
                      kernels.at(i).begin());
            /* Then evolve to index i */
            int status = gsl_odeiv2_driver_apply(
                driver, &eta_current, eta_grid.at(i), kernels.at(i).data());

            if (status != GSL_SUCCESS) {
                throw(std::runtime_error("GLS ODE driver gave error value = " +
                                        std::to_string(status)));
            }
        }
    }

    // Free GSL ODE driver
    gsl_odeiv2_driver_free(driver);
}



static void kernel_initial_conditions(
        int kernel_index,
        int n,
        double k,
        IntegrandTables& tables
        )
{
    // Use growing mode IC for F1 at eta_asymp
    if (n == 1) {
        for (size_t i = 0; i < COMPONENTS; ++i) {
            tables.kernels.at(static_cast<size_t>(kernel_index))
                .values.at(0)
                .at(i) = tables.ev_params.F1_ic_at_k(i, k);
        }
    }
    else {
        // Set F(n>1)-kernels to zero at eta_asymp
        for (size_t i = 0; i < COMPONENTS; ++i) {
            tables.kernels.at(static_cast<size_t>(kernel_index))
                .values.at(0)
                .at(i) = 0;
        }
    }
}
}



int compute_gen_kernels(
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

    /* If kernel_index is not known, -1 is sent as argument */
    if (kernel_index == -1) {
        kernel_index = tables.loop_params.args_2_kernel_index(arguments);
    }

    /* Alias reference to kernel we are working with for convenience/readability */
    Kernel& kernel = tables.kernels.at(static_cast<size_t>(kernel_index));

    /* Check if the kernels are already computed */
    if (kernel.computed) return kernel_index;

    /* Interpolation variables for interpolated RHS */
    std::array<Interpolation1D, COMPONENTS> rhs;

    /* Compute RHS sum in evolution equation if n > 1. If n == 1, the RHS
     * equals 0, which is implemented in solve_kernel_ODE(). */
    if (n > 1) {
        compute_RHS_sum(arguments, n, tables, rhs);
    }

    if (tables.loop_params.dynamics() == EVOLVE_IC_EDS) {
        /* Set up ODE input and system */
        ODEInput input(n, 0, tables.eta_grid.eta_ini(), tables.ev_params, rhs);
        /* Set initial conditions */
        EvolveICEdS::kernel_initial_conditions(kernel_index, tables);
        /* Solve ODE */
        EvolveICEdS::solve_kernel_ODE(
            input, tables.eta_grid,
            tables.kernels.at(static_cast<size_t>(kernel_index)).values);
    }
    else {
        /* tables.loop_params.dynamics() == EVOLVE_ASYMP_IC */
        /* Compute k (sum of kernel arguments) */
        int sum = tables.sum_table.sum_labels(arguments,
                tables.loop_params.n_kernel_args());
        double k = std::sqrt(tables.comp_dot_products()
                                 .at(static_cast<size_t>(sum))
                                 .at(static_cast<size_t>(sum)));
        /* Set up ODE input and system */
        ODEInput input(n, k, tables.eta_grid.eta_ini(), tables.ev_params, rhs);
        /* Set initial conditions */
        EvolveICAsymp::kernel_initial_conditions(kernel_index, n, k, tables);
        /* Solve ODE */
        EvolveICAsymp::solve_kernel_ODE(
            input, tables.eta_grid,
            tables.kernels.at(static_cast<size_t>(kernel_index)).values);
    }

    tables.kernels.at(static_cast<size_t>(kernel_index)).computed = true;
    return kernel_index;
}



/* Only used for EVOLVE_ASYMP_IC */
void compute_F1(
        double k,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid,
        Vec1D<double>& F1_eta_ini, /* out */
        Vec1D<double>& F1_eta_fin  /* out */
        )
{
    size_t time_steps = eta_grid.time_steps();
    size_t pre_time_steps = eta_grid.pre_time_steps();

    Vec2D<double> values;
    values.assign(eta_grid.time_steps(), Vec1D<double>(COMPONENTS, 0.0));

    for (size_t i = 0; i < COMPONENTS; ++i) {
        values.at(0).at(i) = ev_params.F1_ic_at_k(i, k);
    }

    /* Empty rhs */
    std::array<Interpolation1D, COMPONENTS> rhs;

    ODEInput input(1, k, eta_grid.eta_ini(), ev_params, rhs);

    EvolveICAsymp::solve_kernel_ODE(input, eta_grid, values);

    std::copy(values.at(pre_time_steps).begin(),
              values.at(pre_time_steps).end(), F1_eta_ini.begin());
    std::copy(values.at(time_steps - 1).begin(),
              values.at(time_steps - 1).end(), F1_eta_fin.begin());

}
