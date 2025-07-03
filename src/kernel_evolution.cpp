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
    double beta     = tables.beta()(sl, sr);

    size_t index_l = static_cast<size_t>(
        compute(args_l, -1, m_l));
    size_t index_r = static_cast<size_t>(
        compute(args_r,
                                                 -1, m_r));
    const auto& kernel_l =
        tables.kernels.at(index_l).values;
    const auto& kernel_r =
        tables.kernels.at(index_r).values;

    for (size_t i = 0; i < tables.eta_grid.time_steps(); ++i) {
        if constexpr (COMPONENTS >= 4) {
            partial_rhs_sum(i, 2) = alpha_lr * kernel_l(i, 3) * kernel_r(i, 2);
            partial_rhs_sum(i, 3) = beta     * kernel_l(i, 3) * kernel_r(i, 3);
        }
        if constexpr (COMPONENTS >= 2) {
            partial_rhs_sum(i, 0) = alpha_lr * kernel_l(i, 1) * kernel_r(i, 0);
            partial_rhs_sum(i, 1) = beta     * kernel_l(i, 1) * kernel_r(i, 1);
        }
        static_assert(COMPONENTS == 2 || COMPONENTS == 4,
            "vertex only implemented for COMPONENTS = 2 or 4");
    }
}



void KernelEvolver::RHS(
    const int arguments[],
    int n,
    std::array<std::vector<double>, COMPONENTS>& rhs /* out, overwritten(!),
    array of size COMPONENTS, containing vectors which will be interpolated
     * later */
) {
    const size_t time_steps = tables.eta_grid.time_steps();
    const size_t n_kernel_args = tables.loop_structure.n_kernel_args();
    int zero_label = tables.loop_structure.zero_label();

    /* Resize and zero rhs */
    for (auto& vec : rhs) {
        vec.resize(time_steps, 0.0);
    }
    /* partial_rhs_sum is instead a flatted vector for (potential) speedup,
     * use get(partial_rhs_sum, i=time_step, j=component, stride=COMPONENTS)
     * to access */
    Strided2DVec<double> partial_rhs_sum(time_steps, COMPONENTS);

    int args_l[N_KERNEL_ARGS_MAX] = {0};
    int args_r[N_KERNEL_ARGS_MAX] = {0};

    for (int m = 1; m <= n / 2; ++m) {
        /* Initialize partial_rhs_sum to 0 */
        std::fill(partial_rhs_sum.begin(), partial_rhs_sum.end(), 0.0);
        /* Initialize args_l and args_r */
        std::fill(&args_l[0], &args_l[n_kernel_args], zero_label);
        std::fill(&args_r[0], &args_r[n_kernel_args], zero_label);

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

            vertex(m, n-m, args_l, args_r, sum_l, sum_r, partial_rhs_sum);

            // When m != (n - m), we may additionally compute the (n-m)-term by
            // swapping args_l, sum_l, m with args_r, sum_r and (n-m). Then
            // compute_RHS_sum() only needs to sum up to (including) floor(n/2).
            if (m != n - m) {
                vertex(n - m, m, args_r, args_l, sum_r, sum_l, partial_rhs_sum);
            }
        } while (comb.next());

        size_t n_t = static_cast<size_t>(n);
        size_t m_t = static_cast<size_t>(m);
        for (size_t i = 0; i < COMPONENTS; ++i) {
            for (size_t j = 0; j < time_steps; ++j) {
                rhs.at(i).at(j) +=
                    partial_rhs_sum(j, i) / binomial_coeffs[n_t][m_t];
            }
        }
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
        RHS(arguments, n, rhs_sum);
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
