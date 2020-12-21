/*
   tables.hpp

   Created by Petter Taule on 29.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef TABLES_HPP
#define TABLES_HPP

#include <vector>

#include "utilities.hpp"
#include "parameters.hpp"


class IntegrationVariables {
    public:
        Vec1D<double> magnitudes; /* Loop momenta magnitudes */
        Vec1D<double> cos_theta;  /* Cosine of polar angles of the loop momenta */
        Vec1D<double> phi;        /* Azimutal angles */

        IntegrationVariables(int n_loops) {
            magnitudes.assign(n_loops,0);
            cos_theta.assign(n_loops,0);
            phi.assign(n_loops,0);
        }
};


class SumTable {
    private:
        const int zero_label;
        const int n_coeffs;

        Vec2D<int> sum_table;

        int sum_two_labels(int a, int b);
    public:
        SumTable(const LoopParameters& loop_params);
        int sum_labels(const int labels[], size_t size) const;
};


class SPTKernel {
    public:
        double values[EDS_SPT_COMPONENTS] = {0};
        bool computed = false;
};


class Kernel {
    public:
        Vec2D<double> values;
        bool computed = false;
};


class EtaGrid {
    private:
        int pre_time_steps_ = 0;
        int time_steps_     = 0;
        double eta_ini_     = 0;
        double eta_fin_     = 0;
        double eta_asymp_   = 0;

        Vec1D<double> grid_;
    public:
        EtaGrid() = default;
        EtaGrid(
                const int pre_time_steps,
                const int time_steps,
                const double eta_ini,
                const double eta_fin,
                const double eta_asymp
                );

        EtaGrid(
                const int time_steps,
                const double eta_ini,
                const double eta_fin
                );

        int pre_time_steps() const {return pre_time_steps_;}
        int time_steps() const {return time_steps_;}
        double eta_ini() const {return eta_ini_;}
        double eta_fin() const {return eta_fin_;}
        double eta_asymp() const {return eta_asymp_;}

        const Vec1D<double>& grid() const {return grid_;}

        const double& operator[](int i) const {return grid_[i];}
        const double& at(int i) const {return grid_.at(i);}
};


std::ostream& operator<<(std::ostream& out, const EtaGrid& eta_grid);


/* Tables inside integration, one instance for each integration thread */
class IntegrandTables {
    private:
        double k_a    = 0.0;
        double k_b    = 0.0; // For bispectrum
        double cos_ab = 0.0; // For bispectrum

        /* Helper vectors for compute_scalar_products() */
        Vec1D<int> a_coeffs;
        Vec1D<int> b_coeffs;

        Vec2D<double> bare_scalar_products_; /* N_COEFFS x N_COEFFS   */
        Vec2D<double> scalar_products_;      /* N_CONFIGS x N_CONFIGS */
        Vec2D<double> alpha_;                /* N_CONFIGS x N_CONFIGS */
        Vec2D<double> beta_;                 /* N_CONFIGS x N_CONFIGS */

        void reset_spt_kernels();
        void reset_kernels();

        void ps_compute_bare_scalar_products(); /* Power spectrum */
        void bs_compute_bare_scalar_products(); /* Bispectrum */

        void compute_scalar_products();
        void compute_alpha_beta();
    public:
        const LoopParameters& loop_params;
        const SumTable& sum_table;

        const EvolutionParameters& ev_params;
        const EtaGrid& eta_grid;

        IntegrationVariables vars;

        Vec1D<SPTKernel> spt_kernels;
        Vec1D<Kernel> kernels;

        IntegrandTables(
                double k_a,
                double k_b,
                double cos_ab,
                const LoopParameters& loop_params,
                const SumTable& sum_table,
                const EvolutionParameters& ev_params,
                const EtaGrid& eta_grid
                );
        IntegrandTables(
                double k_a,
                const LoopParameters& loop_params,
                const SumTable& sum_table,
                const EvolutionParameters& ev_params,
                const EtaGrid& eta_grid
                );

        const Vec2D<double>& bare_scalar_products() const {return
            bare_scalar_products_;}
        const Vec2D<double>& scalar_products() const {return scalar_products_;}
        const Vec2D<double>& alpha() const {return alpha_;}
        const Vec2D<double>& beta() const {return beta_;}

        void reset();
        void compute_tables();
};


#endif /* ifndef TABLES_HPP */
