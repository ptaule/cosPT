#ifndef KERNEL_EVOLUTION_HPP
#define KERNEL_EVOLUTION_HPP

#include <array>

#include "utilities.hpp"

class Interpolation1D;
class IntegrandTables;
class EvolutionParameters;
class EtaGrid;

class KernelEvolver {
    private:
        const int n;
        const double k;
        const EtaGrid& eta_grid;
        Vec2D<double> omega;
        Vec2D<double> kernels /* time_steps * COMPONENTS table of kernels */

        struct ODEParameters {
            const int n;
            const double k;
            const EvolutionParameters& ev_params;
            const std::array<Interpolation1D, COMPONENTS>& rhs;
            Vec2D<double> omega;
            ODEParameters(
                double k,
                const EvolutionParameters &ev_params,
                const std::array<Interpolation1D, COMPONENTS>& rhs
                )
                : k(k), ev_params(ev_params), rhs(rhs), omega(Vec2D<double>()) : {}
        };
        ODEParameters params;

        void initialize_omega() {
            omega.resize(COMPONENTS);
            /*Initialize linear matrix omega, defined by y' = omega*y +
             * rhs. We rescale the kernels by exp(eta*n) so that the
             * diagonal always includes a term n*/
            for (std::size_t i = 0; i < COMPONENTS; ++i) {
                omega.at(i).resize(COMPONENTS);
                for (std::size_t j = 0; j < COMPONENTS; ++j) {
                    if (i == j) {
                        omega.at(i).at(j) = n;
                    } else {
                        omega.at(i).at(j) = 0;
                    }
                }
            }
        }

        static int ode_system(double eta, const double y[], double f[], void* ode_input);

    public:
        KernelEvolver();
        KernelEvolver(
                int n,
                double k,
                const EtaGrid& eta_grid,
                const EvolutionParameters &ev_params,
                const std::array<Interpolation1D, COMPONENTS>& rhs
                )
            : n(n), k(k), eta_grid(eta_grid), params(k, ev_params, rhs) {
                initialize_omega();
                params.omega = omega;
            }
};


int compute_gen_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );


void compute_F1(
        double k,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid,
        Vec1D<double>& F1_eta_ini, /* out */
        Vec1D<double>& F1_eta_fin  /* out */
        );

#endif /* ifndef KERNEL_EVOLUTION_HPP */
