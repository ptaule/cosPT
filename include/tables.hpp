/*
   tables.hpp

   Created by Petter Taule on 29.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef TABLES_HPP
#define TABLES_HPP

#include <vector>

#include "utilities.hpp"

class Settings {
    public:
        Dynamics dynamics;
        Spectrum spectrum;

        short int n_loops;

        short int n_coeffs;
        short int n_configs;
        short int n_kernels;
        short int n_kernel_args;
        short int zero_label;
        short int single_loop_label_min;
        short int single_loop_label_max;
        short int single_loop_block_size;

        Vec1D<short int> single_loops;

        short int first_composite_block_size = 0; /* Bispectrum */

        short int components     = 0;
        short int time_steps     = 0;
        short int pre_time_steps = 0;
        double eta_i             = 0.0;
        double eta_f             = 0.0;
        double eta_asymp         = 0.0;

        Settings(
                short int n_loops,
                Spectrum spectrum,
                Dynamics dynamics,
                short int time_steps,
                short int pre_time_steps,
                short int components,
                double eta_i,
                double eta_f,
                double eta_asymp
                );
        Settings(
                short int n_loops,
                Spectrum spectrum,
                Dynamics dynamics,
                short int time_steps,
                short int components,
                double eta_i,
                double eta_f
                );
        Settings(
                short int n_loops,
                Spectrum spectrum,
                Dynamics dynamics
                );
};

class IntegrationVariables {
    public:
        Vec1D<double> magnitudes; /* Loop momenta magnitudes */
        Vec1D<double> cos_theta;  /* Cosine of polar angles of the loop momenta */
        Vec1D<double> phi;        /* Azimutal angles */

        IntegrationVariables(short int n_loops) {
            magnitudes.resize(n_loops);
            cos_theta.resize(n_loops);
            phi.resize(n_loops - 1);
        }
};

class SumTable {
    private:
        const Settings& settings;
        Vec2D<short int> sum_table;

        short int sum_two_labels(short int a, short int b);
    public:
        SumTable(const Settings& settings);
        short int sum_labels(const short int labels[], size_t size) const;
};

class SPTKernel {
    public:
        double values[SPT_COMPONENTS] = {0};
        bool computed = false;
};

class Kernel {
    public:
        Vec2D<double> values;
        bool computed = false;
};

/* Tables inside integration, one instance for each integration thread */
class IntegrandTables {
    private:
        void reset_spt_kernels();
        void reset_kernels();

        // Power spectrum
        void ps_compute_bare_scalar_products();
        void ps_compute_scalar_products();
        void ps_compute_alpha_beta();

        // Bi-spectrum
        void bs_compute_bare_scalar_products();
        void bs_compute_scalar_products();
        void bs_compute_alpha_beta();
    public:
        double k_a = 0.0;
        double k_b = 0.0; // For bispectrum

        const Settings& settings;

        const SumTable& sum_table;
        IntegrationVariables vars;

        const Vec1D<double>& eta_grid;

        Vec2D<double> bare_scalar_products; /* N_COEFFS x N_COEFFS   */
        Vec2D<double> scalar_products;      /* N_CONFIGS x N_CONFIGS */
        Vec2D<double> alpha;                /* N_CONFIGS x N_CONFIGS */
        Vec2D<double> beta;                 /* N_CONFIGS x N_CONFIGS */

        Vec1D<SPTKernel> spt_kernels;
        Vec1D<Kernel> kernels;

        IntegrandTables(
                double k_a,
                double k_b,
                const Settings& settings,
                const SumTable& sum_table,
                const Vec1D<double>& eta_grid
                );
        IntegrandTables(
                double k_a,
                const Settings& settings,
                const SumTable& sum_table,
                const Vec1D<double>& eta_grid
                );

        void reset();
        void compute_tables();
};

Vec1D<double> initialize_eta_grid(const Settings& settings);

short int kernel_index_from_arguments(
        const short int arguments[],
        const Settings& settings
        );

#endif /* ifndef TABLES_HPP */
