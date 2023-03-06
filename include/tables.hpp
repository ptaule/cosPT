/*
   tables.hpp

   Created by Petter Taule on 29.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef TABLES_HPP
#define TABLES_HPP

#include <functional>
#include <vector>

#include "utilities.hpp"

class LoopParameters;
class EvolutionParameters;

struct IntegrationVariables {
    Vec1D<double> magnitudes; /* Loop momenta magnitudes */
    Vec1D<double> cos_theta;  /* Cosine of polar angles of the loop momenta */
    Vec1D<double> phi;        /* Azimutal angles */
    double mu_los;            /* RSD: Cosine of angle between k and L.o.S. */

    IntegrationVariables(std::size_t n_loops) : mu_los(0) {
        magnitudes.assign(n_loops,0);
        cos_theta.assign(n_loops,0);
        phi.assign(n_loops,0);
    }
};


struct SPTKernel {
    double values[EDS_SPT_COMPONENTS] = {0};
    bool computed = false;
};


struct Kernel {
    Vec2D<double> values;
    bool computed = false;
};


struct RSDKernel {
    double value;
    bool computed = false;
};


class SumTable {
    private:
        const int zero_label;
        const std::size_t n_coeffs;

        Vec2D<int> sum_table;

        int sum_two_labels(int a, int b);
    public:
        SumTable(const LoopParameters& loop_params);
        int sum_labels(const int labels[], std::size_t size) const;
};


class EtaGrid {
    private:
        std::size_t pre_time_steps_ = 0;
        std::size_t time_steps_     = 0;
        double eta_ini_     = 0;
        double eta_fin_     = 0;
        double eta_asymp_   = 0;

        Vec1D<double> grid_;
    public:
        EtaGrid() = default;
        EtaGrid(
                std::size_t pre_time_steps,
                std::size_t time_steps,
                double eta_ini,
                double eta_fin,
                double eta_asymp
                );
        EtaGrid(
                std::size_t time_steps,
                double eta_ini,
                double eta_fin
                );

        std::size_t pre_time_steps() const {return pre_time_steps_;}
        std::size_t time_steps() const {return time_steps_;}
        double eta_ini() const {return eta_ini_;}
        double eta_fin() const {return eta_fin_;}
        double eta_asymp() const {return eta_asymp_;}

        const Vec1D<double>& grid() const {return grid_;}

        const double& operator[](std::size_t i) const {return grid_[i];}
        const double& at(size_t i) const {return grid_.at(i);}
};


std::ostream& operator<<(std::ostream& out, const EtaGrid& eta_grid);


/* Tables inside integration, one instance for each integration thread */
class IntegrandTables {
    private:
        double k_a    = 0.0;
        double k_b    = 0.0; // For bispectrum
        double cos_ab = 0.0; // For bispectrum

        double rsd_f; /* Growth factor (at observation redshift) for L.o.S. */

        /* Helper vectors for compute_scalar_products() */
        Vec1D<int> a_coeffs;
        Vec1D<int> b_coeffs;

        Vec2D<double> bare_dot_prod;
        /* N_COEFFS x N_COEFFS, dot products between external/loop wavenumbers */
        Vec2D<double> comp_dot_prod;
        /* N_CONFIGS x N_CONFIGS, dot products between composite wavenumbers */
        Vec2D<double> alpha_;                /* N_CONFIGS x N_CONFIGS */
        Vec2D<double> beta_;                 /* N_CONFIGS x N_CONFIGS */

        Vec1D<double> bare_los_projection_;
        /* N_COEFFS dot products between external/loop wavenumber and L.o.S */
        Vec1D<double> comp_los_projection_;
        /* N_CONFIGS dot products between composite wavenumber and L.o.S */

        void reset_spt_kernels();
        void reset_kernels();
        void reset_rsd_kernels();

        void compute_bare_los_proj();
        void compute_comp_los_proj();

        void compute_bare_dot_prod();

        void compute_comp_dot_prod();
        void compute_alpha_beta();
    public:
        const LoopParameters& loop_params;
        const SumTable& sum_table;

        const EvolutionParameters& ev_params;
        const EtaGrid& eta_grid;

        IntegrationVariables vars;

        Vec1D<SPTKernel> spt_kernels;
        Vec1D<Kernel> kernels;

        Vec1D<RSDKernel> rsd_kernels;
        Vec2D<RSDKernel> vel_power_kernels;

        IntegrandTables(
                double k_a,
                double k_b,
                double cos_ab,
                double rsd_growth_f,
                const LoopParameters& loop_params,
                const SumTable& sum_table,
                const EvolutionParameters& ev_params,
                const EtaGrid& eta_grid
                );

        const Vec2D<double>& bare_dot_products() const { return bare_dot_prod; }
        const Vec2D<double>& comp_dot_products() const { return comp_dot_prod; }
        const Vec2D<double>& alpha() const {return alpha_;}
        const Vec2D<double>& beta() const {return beta_;}

        const Vec1D<double>& bare_los_projection() const { return bare_los_projection_; }
        const Vec1D<double>& comp_los_projection() const { return comp_los_projection_; }

        double rsd_growth_f() const { return rsd_f; }

        void reset();
        void compute_tables();
};


#endif /* ifndef TABLES_HPP */
