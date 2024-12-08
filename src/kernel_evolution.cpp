#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

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
        Vec2D<double>& partial_rhs_sum
        )
{
    double alpha_lr = tables.alpha()
        .at(static_cast<size_t>(sum_l))
        .at(static_cast<size_t>(sum_r));
    double beta = tables.beta()
        .at(static_cast<size_t>(sum_l))
        .at(static_cast<size_t>(sum_r));

    int index_l = compute(args_l, -1, m_l);
    int index_r = compute(args_r, -1, m_r);

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



void KernelEvolver::compute_RHS_sum(
        const int arguments[],
        int n,
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

            vertex(m, n-m, args_l, args_r, sum_l, sum_r, partial_rhs_sum);

            // When m != (n - m), we may additionally compute the (n-m)-term by
            // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
            // compute_RHS_sum() only needs to sum up to (including) floor(n/2).
            if (m != n - m) {
                vertex(n - m, m, args_r, args_l, sum_r, sum_l, partial_rhs_sum);
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



int KernelEvolver::ode_system(double eta, const double y[], double f[], void* ode_input) {
    ODEParameters params = *(static_cast<ODEParameters*>(ode_input));

    int n = params.n;
    double k = params.k;
    Vec2D<double>& omega = params.omega;

    double kappa = params.ev_params.get_kappa(0);
    double zeta = params.ev_params.zeta_at_eta(0, eta);
    double xi  = params.ev_params.xi_at_eta_k(0, eta, k);

    UPDATE_OMEGA_MATRIX

    /* If n == 1, rhs = 0 */
    if (params.n == 1){
        for (size_t i = 0; i < COMPONENTS; ++i) {
            f[i] = 0;
            for (std::size_t j = 0; j < COMPONENTS; ++j) {
                f[i] += omega.at(i).at(j) * y[j];
            }
        }
    }
    /* Otherwise, interpolate */
    else {
        for (size_t i = 0; i < COMPONENTS; ++i) {
        /* We rescale the kernels by exp(eta*n) so that the diagonal always
         * includes a term n */
            f[i] = params.rhs.at(i)(eta) - n * y[i];
            for (std::size_t j = 0; j < COMPONENTS; ++j) {
                f[i] += omega.at(i).at(j) * y[j];
            }
        }
    }
    return GSL_SUCCESS;
}



/* For obtaining the correct initial conditions at eta_ini, evolve system from
 * eta_asymp to eta_ini fixed omega-matrix (i.e. no call to
 * UPDATE_OMEGA_MATRIX) */
int KernelEvolver::ode_system_fixed_eta(double eta, const double y[], double f[], void* ode_input) {
    ODEParameters params = *(static_cast<ODEParameters*>(ode_input));
    int n = params.n;
    Vec2D<double>& omega = params.omega;

    /* If n == 1, rhs = 0 */
    if (params.n == 1){
        for (size_t i = 0; i < COMPONENTS; ++i) {
            f[i] = 0;
            for (std::size_t j = 0; j < COMPONENTS; ++j) {
                f[i] += omega.at(i).at(j) * y[j];
            }
        }
    }
    /* Otherwise, interpolate */
    else {
        for (size_t i = 0; i < COMPONENTS; ++i) {
            f[i] = params.rhs.at(i)(eta) - n * y[i];
            for (std::size_t j = 0; j < COMPONENTS; ++j) {
                f[i] += omega.at(i).at(j) * y[j];
            }
        }
    }
    return GSL_SUCCESS;
}



void KernelEvolver::set_EdS_ICs(
        int n,
        double k,
        int kernel_index
        )
{
    UNUSED(n);
    UNUSED(k);
    for (size_t i = 0; i < EDS_SPT_COMPONENTS; ++i) {
        tables.kernels.at(static_cast<size_t>(kernel_index)).values.at(0).at(i) =
            tables.spt_kernels.at(static_cast<size_t>(kernel_index)).values[i];
    }
}



void KernelEvolver::set_asymptotic_ICs(
        int n,
        double k,
        int kernel_index
        )
{
    // Use growing mode IC for F1 at eta_asymp
    if (n == 1) {
        for (size_t i = 0; i < COMPONENTS; ++i) {
            tables.kernels.at(static_cast<size_t>(kernel_index))
                .values.at(0)
                .at(i) = tables.ev_params.omega_eigenvector_at_k(i, k);
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



void KernelEvolver::evolve(
    int n,
    double k,
    int kernel_index
    )
{
    const EvolutionParameters& ev_params = tables.ev_params;
    const EtaGrid& eta_grid              = tables.eta_grid;

    Vec2D<double>& kernel_vec = tables.kernels.at(static_cast<size_t>(kernel_index)).values;

    size_t time_steps     = eta_grid.time_steps();
    size_t pre_time_steps = eta_grid.pre_time_steps();
    double eta_current = eta_grid.at(0);

    size_t i = 0;

    if (tables.loop_params.dynamics() == EVOLVE_ASYMPTOTIC_ICS) {
        /* Evolve from asymptotic eta = eta_asymp to initial eta = eta_ini */
        sys.function = &ode_system_fixed_eta;
        for (i = 1; i < pre_time_steps + 1; ++i) {
            /*For n == 1 kernels, we may multipliy by exp("growing mode" eigenvalue)*/
            /*before eta_I*/
            if (n == 1) {
                for (size_t j = 0; j < COMPONENTS; ++j) {
                    kernel_vec.at(i).at(j) =
                        kernel_vec.at(0).at(j) *
                        std::exp(ev_params.omega_eigenvalues_at_k(k) *
                                (eta_grid.at(i) - eta_grid.at(0)));
                }
            }
            else {
                double eta = eta_grid.eta_ini();
                double kappa = ev_params.get_kappa(0);
                double zeta = ev_params.zeta_at_eta(0, eta);
                double xi  = ev_params.xi_at_eta_k(0, eta, k);

                UPDATE_OMEGA_MATRIX

                /* For calculation of kernel_vec.at(i), the initial condition is
                 * kernel_vec.at(i-1). Hence we copy kernel_vec.at(i-1) to kernel_vec.at(i) */
                std::copy(kernel_vec.at(i - 1).begin(), kernel_vec.at(i - 1).end(),
                        kernel_vec.at(i).begin());
                /* Then evolve to index i */
                int status = gsl_odeiv2_driver_apply( driver, &eta_current, eta_grid.at(i), kernel_vec.at(i).data());

                if (status != GSL_SUCCESS) {
                    throw(std::runtime_error("GLS ODE driver gave error value = " +
                                std::to_string(status)));
                }
            }
        }
    }

    ++i;
    /* eta_ini to eta_fin */
    sys.function = &ode_system;
    for (; i < time_steps; ++i) {
        /* For calculation of kernel_vec.at(i), the initial condition is
            * kernel_vec.at(i-1). Hence we copy kernel_vec.at(i-1) to kernel_vec.at(i) */
        std::copy(kernel_vec.at(i - 1).begin(), kernel_vec.at(i - 1).end(),
                    kernel_vec.at(i).begin());
        /* Then evolve to index i */
        int status = gsl_odeiv2_driver_apply( driver, &eta_current, eta_grid.at(i), kernel_vec.at(i).data());

        if (status != GSL_SUCCESS) {
            throw(std::runtime_error("GLS ODE driver gave error value = " +
                                    std::to_string(status)));
        }
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
        compute_RHS_sum(arguments, n, rhs);
    }

    /* Compute k (sum of kernel arguments) */
    int sum = tables.sum_table.sum_labels(arguments,
        tables.loop_params.n_kernel_args());
    double k = std::sqrt(tables.comp_dot_products()
        .at(static_cast<size_t>(sum))
        .at(static_cast<size_t>(sum)));
    /* Set up ODE input and system */
    initialize_omega(n);
    ODEParameters params(n, k, tables.eta_grid.eta_ini(),
        tables.ev_params, rhs, omega);
    sys.params = &params;
    /* Set initial conditions */
    set_ICs(n, k, kernel_index);
    /* Solve ODE */
    evolve(n, k, kernel_index);

    tables.kernels.at(static_cast<size_t>(kernel_index)).computed = true;
    return kernel_index;
}
