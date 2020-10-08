/*
   tables.hpp

   Created by Petter Taule on 29.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef TABLES_HPP
#define TABLES_HPP

#include <vector>
#include <functional>

#include "utilities.hpp"
#include "parameters.hpp"


class IntegrationVariables {
    public:
        Vec1D<double> magnitudes; /* Loop momenta magnitudes */
        Vec1D<double> cos_theta;  /* Cosine of polar angles of the loop momenta */
        Vec1D<double> phi;        /* Azimutal angles */

        IntegrationVariables(int n_loops) {
            magnitudes.resize(n_loops);
            cos_theta.resize(n_loops);
            phi.resize(n_loops - 1);
        }
};


class SumTable {
    private:
        const Parameters& params;
        Vec2D<int> sum_table;

        int sum_two_labels(int a, int b);
    public:
        SumTable(const Parameters& params);
        int sum_labels(const int labels[], size_t size) const;
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


class EtaGrid {
    private:
        int pre_time_steps = 0;
        int time_steps     = 0;
        double eta_ini     = 0;
        double eta_fin     = 0;
        double eta_asymp   = 0;

        Vec1D<double> grid;
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

        int get_pre_time_steps() const {return pre_time_steps;}
        int get_time_steps() const {return time_steps;}
        double get_eta_ini() const {return eta_ini;}
        double get_eta_fin() const {return eta_fin;}
        double get_eta_asymp() const {return eta_asymp;}

        const Vec1D<double>& get_grid() const {return grid;}

        const double& operator[](int i) const {return grid[i];}
        const double& at(int i) const {return grid.at(i);}
};


std::ostream& operator<<(std::ostream& out, const EtaGrid& eta_grid);


/* Tables inside integration, one instance for each integration thread */
class IntegrandTables {
    private:
        void reset_spt_kernels();
        void reset_kernels();

        void ps_compute_bare_scalar_products(); /* Power spectrum */
        void bs_compute_bare_scalar_products(); /* Bispectrum */

        void compute_scalar_products();
        void compute_alpha_beta();
    public:
        double k_a = 0.0;
        double k_b = 0.0; // For bispectrum

        const Parameters& params;
        const SumTable& sum_table;

        const EvolutionParameters& ev_params;
        const EtaGrid& eta_grid;

        IntegrationVariables vars;

        Vec2D<double> bare_scalar_products; /* N_COEFFS x N_COEFFS   */
        Vec2D<double> scalar_products;      /* N_CONFIGS x N_CONFIGS */
        Vec2D<double> alpha;                /* N_CONFIGS x N_CONFIGS */
        Vec2D<double> beta;                 /* N_CONFIGS x N_CONFIGS */

        Vec1D<SPTKernel> spt_kernels;
        Vec1D<Kernel> kernels;

        std::function<int(const int[], const Parameters &)>
            kernel_index_from_arguments;

        IntegrandTables(
                double k_a,
                double k_b,
                const Parameters& params,
                const SumTable& sum_table,
                const EvolutionParameters& ev_params,
                const EtaGrid& eta_grid
                );
        IntegrandTables(
                double k_a,
                const Parameters& params,
                const SumTable& sum_table,
                const EvolutionParameters& ev_params,
                const EtaGrid& eta_grid
                );

        void reset();
        void compute_tables();
};

namespace ps {
int kernel_index_from_arguments(const int arguments[],
                                const Parameters &params);
}

namespace bs {
int kernel_index_from_arguments(const int arguments[],
                                const Parameters &params);
}

#endif /* ifndef TABLES_HPP */
