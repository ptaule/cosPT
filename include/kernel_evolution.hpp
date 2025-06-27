#ifndef KERNEL_EVOLUTION_HPP
#define KERNEL_EVOLUTION_HPP

#include <array>
#include <functional>

#include <gsl/gsl_odeiv2.h>

#include "omega_matrix.hpp"
#include "parameters.hpp"
#include "tables.hpp"
#include "utilities.hpp"

class Interpolation1D;
class EvolutionParameters;
class EtaGrid;

class KernelEvolver {
    private:
        gsl_odeiv2_system sys;
        gsl_odeiv2_driver* driver;

        IntegrandTables& tables;

        /* Used for computing ICs if needed, omega matrix used in ODE is part
         * of ODEParameters */
        Strided2DVec<double> omega;

        struct ODEParameters {
            const int n;
            const double k;
            const double eta_ini;  /* Only used when dynamics = EVOLVE_IC_ASYMP */
            const EvolutionParameters& ev_params;
            std::array<Interpolation1D, COMPONENTS>& rhs;
            Strided2DVec<double>& omega;
            ODEParameters(
                int n,
                double k,
                double eta_ini,
                const EvolutionParameters &ev_params,
                std::array<Interpolation1D, COMPONENTS>& rhs,
                Strided2DVec<double>& omega
                )
                : n(n), k(k), eta_ini(eta_ini), ev_params(ev_params), rhs(rhs),
                omega(omega) {}
        };

        void vertex(int m_l, int m_r, const int args_l[], const int args_r[],
                    int sum_l, int sum_r, Strided2DVec<double>& partial_rhs_sum);
        void vertex_cubic(int n_a, int n_b, int n_c, const int args_a[], const
                          int args_b[], const int args_c[], int sum_a, int
                          sum_b, int sum_c, int sum_ab,
                          Strided2DVec<double>& partial_rhs_sum);

        // Loops splitting kernel arguments are symmetric, because the kernels are
        // symmetric in arguments. Therefore, introduce a helper function to apply the
        // vertex twice, but with "left" and "right" arguments swapped.
        void apply_symmetric_vertex(int n_l, int n_r, int args_l[],
                                    int args_r[], int sum_l, int sum_r,
                                    Strided2DVec<double>& out
        ) {
            vertex(n_l, n_r, args_l, args_r, sum_l, sum_r, out);
            if (n_l != n_r) {
                vertex(n_r, n_l, args_r, args_l,
                       sum_r, sum_l, out);
            }
        }
        inline void apply_symmetric_vertex_cubic(
            int n_a, int n_b, int n_c,
            int a[], int b[], int c[],
            int sum_a, int sum_b, int sum_c,
            Strided2DVec<double>& out
        ) {
            int sum_ab = tables.sum_table(sum_a, sum_b);
            vertex_cubic(n_a, n_b, n_c, a, b,
                         c, sum_a, sum_b, sum_c, sum_ab,
                         out);
            if (n_b != n_c) {
                int sum_ac = tables.sum_table(sum_a, sum_c);
                vertex_cubic(n_a, n_c, n_b, a, c,
                             b, sum_a, sum_c, sum_b,
                             sum_ac, out);
            }
        }

        void RHS(const int arguments[], int n,
                             std::array<std::vector<double>, COMPONENTS>& rhs);
        void RHS_cubic(const int arguments[], int n,
                             std::array<std::vector<double>, COMPONENTS>& rhs);

        static int ode_system(double eta, const double y[], double f[], void* ode_input);
        static int ode_system_fixed_eta(double eta, const double y[], double f[], void* ode_input);

        std::function<void(const int[], int, int, double)> set_ICs;
        void set_EdS_ICs(const int arguments[], int kernel_index , int n,
                         double k);
        void set_asymptotic_ICs(const int arguments[], int kernel_index,
                                int n, double k);

        void evolve(int kernel_index, int n, double k);
    public:
        KernelEvolver();
        KernelEvolver(IntegrandTables& tables);
        ~KernelEvolver() { gsl_odeiv2_driver_free(driver); }

        int compute(
                const int arguments[],
                int kernel_index,
                int n
                );
};


#endif /* ifndef KERNEL_EVOLUTION_HPP */
