#ifndef KERNEL_EVOLUTION_HPP
#define KERNEL_EVOLUTION_HPP

#include <array>
#include <functional>
#include <exception>

#include <gsl/gsl_odeiv2.h>
#include <stdexcept>

#include "utilities.hpp"
#include "parameters.hpp"
#include "tables.hpp"

class Interpolation1D;
class EvolutionParameters;
class EtaGrid;

class KernelEvolver {
    private:
        gsl_odeiv2_system sys;
        gsl_odeiv2_driver* driver;

        IntegrandTables& tables;

        Vec2D<double> omega; /* Used for computing ICs if needed, omega matrix used in ODE is part of ODEParameters */

        struct ODEParameters {
            const int n;
            const double k;
            const double eta_ini;  /* Only used when dynamics = EVOLVE_IC_ASYMP */
            const EvolutionParameters& ev_params;
            std::array<Interpolation1D, COMPONENTS>& rhs;
            Vec2D<double>& omega;
            ODEParameters(
                int n,
                double k,
                double eta_ini,
                const EvolutionParameters &ev_params,
                std::array<Interpolation1D, COMPONENTS>& rhs,
                Vec2D<double>& omega
                )
                : n(n), k(k), eta_ini(eta_ini), ev_params(ev_params), rhs(rhs),
                omega(omega) {}
        };

        void vertex(int m_l, int m_r, const int args_l[], const int args_r[],
                int sum_l, int sum_r, Vec2D<double>& partial_rhs_sum);
        void compute_RHS_sum(const int arguments[], int n, std::array<Interpolation1D, COMPONENTS>& rhs);

        static int ode_system(double eta, const double y[], double f[], void* ode_input);
        static int ode_system_fixed_eta(double eta, const double y[], double f[], void* ode_input);

        std::function<void(int, double, int)> set_ICs;
        void set_EdS_ICs(int n, double k, int kernel_index);
        void set_asymptotic_ICs(int n, double k, int kernel_index);

        void evolve(int n, double k, int kernel_index);
    public:
        KernelEvolver();
        KernelEvolver(IntegrandTables& tables)
            : tables(tables),
            omega(COMPONENTS, Vec1D<double>(COMPONENTS, 0.0))
        {
            sys = {ode_system, nullptr, COMPONENTS, nullptr};
            driver = gsl_odeiv2_driver_alloc_y_new( &sys,
                    gsl_odeiv2_step_rkf45,
                    tables.ev_params.ode_hstart(),
                    tables.ev_params.ode_rtol(),
                    tables.ev_params.ode_atol()
                    );
            if (tables.loop_params.dynamics() == EVOLVE_EDS_ICS) {
                set_ICs = [this](int n, double k, int kernel_index) {
                    return this->set_EdS_ICs(n, k, kernel_index);
                };
            }
            else if (tables.loop_params.dynamics() == EVOLVE_ASYMPTOTIC_ICS) {
                set_ICs = [this](int n, double k, int kernel_index) {
                    return this->set_asymptotic_ICs(n, k, kernel_index);
                };
            }
            else {
              throw std::runtime_error(
                  "include/kernel_evolution.hpp: KernelEvolver only allows "
                  "dynamics = EVOLVE_IC_ASYMP or EVOLVE_ASYMPTOTIC_ICS.");
            }
        }
        ~KernelEvolver() { gsl_odeiv2_driver_free(driver); }

        int compute(
                const int arguments[],
                int kernel_index,
                int n
                );
};

#endif /* ifndef KERNEL_EVOLUTION_HPP */
