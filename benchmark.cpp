/*
   benchmark.cpp

   Created by Petter Taule on 28.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

#include <benchmark/benchmark.h>

#include "include/utilities.hpp"
#include "include/parameters.hpp"
#include "include/tables.hpp"
#include "include/diagrams.hpp"
#include "include/interpolation.hpp"
#include "include/integrand.hpp"

using std::size_t;

static void BM_PS_EdS_integrand_1loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("ini/benchmark/ps_eds_spt.cfg", -1, -1, -1, 0, 0);

        int n_loops = 1;
        int n_dims = 3 * n_loops - 1;

        LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
        SumTable sum_table(loop_params);
        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        IntegrationInput input(cfg.q_min(), cfg.q_max());

        input.input_ps = Interpolation1D(cfg.input_ps_file(), cfg.input_ps_rescale());

        input.pair_correlations = cfg.pair_correlations();
        input.ps_diagrams = ps::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.k_a(), loop_params, sum_table, ev_params,
                eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with ps::integrand() or bs::integrand() */
            int core = -1;
            ps::integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



static void BM_PS_EdS_integrand_2loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("ini/benchmark/ps_eds_spt.cfg", -1, -1, -1, 0, 0);

        int n_loops = 2;
        int n_dims = 3 * n_loops - 1;

        LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
        SumTable sum_table(loop_params);
        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        IntegrationInput input(cfg.q_min(), cfg.q_max());

        input.input_ps = Interpolation1D(cfg.input_ps_file(), cfg.input_ps_rescale());

        input.pair_correlations = cfg.pair_correlations();
        input.ps_diagrams = ps::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.k_a(), loop_params, sum_table, ev_params,
                eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with ps::integrand() or bs::integrand() */
            int core = -1;
            ps::integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
            return;
    }
}


static void BM_BS_EdS_integrand_1loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("ini/benchmark/bs_eds_spt.cfg", -1, -1, -1, 0, 0);

        int n_loops = 1;
        int n_dims = 3 * n_loops;

        LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
        SumTable sum_table(loop_params);
        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        IntegrationInput input(cfg.q_min(), cfg.q_max());

        input.input_ps = Interpolation1D(cfg.input_ps_file(), cfg.input_ps_rescale());

        input.pair_correlations = cfg.pair_correlations();
        input.bs_diagrams = bs::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.k_a(), cfg.k_b(), cfg.cos_ab(), loop_params,
                sum_table, ev_params, eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with bs::integrand() or bs::integrand() */
            int core = -1;
            bs::integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



static void BM_BS_EdS_integrand_2loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("ini/benchmark/bs_eds_spt.cfg", -1, -1, -1, 0, 0);

        int n_loops = 2;
        int n_dims = 3 * n_loops;

        LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
        SumTable sum_table(loop_params);
        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        IntegrationInput input(cfg.q_min(), cfg.q_max());

        input.input_ps = Interpolation1D(cfg.input_ps_file(), cfg.input_ps_rescale());

        input.pair_correlations = cfg.pair_correlations();
        input.bs_diagrams = bs::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.k_a(), cfg.k_b(), cfg.cos_ab(), loop_params,
                sum_table, ev_params, eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with bs::integrand() or bs::integrand() */
            int core = -1;
            bs::integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



static void BM_PS_2fluid_integrand_1loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("ini/benchmark/m_nu_0.05eV.cfg", -1, -1, -1, 0, 0);

        int n_loops = 1;
        int n_dims = 3 * n_loops - 1;

        LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
        SumTable sum_table(loop_params);
        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        ev_params = EvolutionParameters(cfg.f_nu(), cfg.omega_m_0(),
                cfg.zeta_file(), cfg.redshift_file(),
                cfg.omega_eigenvalues_file(), cfg.F1_ic_files(),
                cfg.effcs2_x_grid(), cfg.effcs2_y_grid(),
                cfg.effcs2_data(), cfg.ode_atol(), cfg.ode_rtol(),
                cfg.ode_hstart());
        eta_grid = EtaGrid(cfg.pre_time_steps(), cfg.time_steps(), cfg.eta_ini(),
                cfg.eta_fin(), cfg.eta_asymp());

        IntegrationInput input(cfg.q_min(), cfg.q_max());

        input.input_ps = Interpolation1D(cfg.input_ps_file(), cfg.input_ps_rescale());

        input.pair_correlations = cfg.pair_correlations();
        input.ps_diagrams = ps::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.k_a(), loop_params, sum_table, ev_params,
                eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with ps::integrand() or bs::integrand() */
            int core = -1;
            ps::integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



static void BM_PS_2fluid_integrand_2loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("ini/benchmark/m_nu_0.05eV.cfg", -1, -1, -1, 0, 0);

        int n_loops = 2;
        int n_dims = 3 * n_loops - 1;

        LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
        SumTable sum_table(loop_params);
        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        ev_params = EvolutionParameters(cfg.f_nu(), cfg.omega_m_0(),
                cfg.zeta_file(), cfg.redshift_file(),
                cfg.omega_eigenvalues_file(), cfg.F1_ic_files(),
                cfg.effcs2_x_grid(), cfg.effcs2_y_grid(),
                cfg.effcs2_data(), cfg.ode_atol(), cfg.ode_rtol(),
                cfg.ode_hstart());
        eta_grid = EtaGrid(cfg.pre_time_steps(), cfg.time_steps(), cfg.eta_ini(),
                cfg.eta_fin(), cfg.eta_asymp());

        IntegrationInput input(cfg.q_min(), cfg.q_max());

        input.input_ps = Interpolation1D(cfg.input_ps_file(), cfg.input_ps_rescale());

        input.pair_correlations = cfg.pair_correlations();
        input.ps_diagrams = ps::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.k_a(), loop_params, sum_table, ev_params,
                eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with ps::integrand() or bs::integrand() */
            int core = -1;
            ps::integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



// Register the functions as benchmarks
BENCHMARK(BM_PS_EdS_integrand_1loop);
BENCHMARK(BM_PS_EdS_integrand_2loop);
BENCHMARK(BM_BS_EdS_integrand_1loop);
BENCHMARK(BM_BS_EdS_integrand_2loop);
BENCHMARK(BM_PS_2fluid_integrand_1loop);
BENCHMARK(BM_PS_2fluid_integrand_2loop);
// Run the benchmark
BENCHMARK_MAIN();
