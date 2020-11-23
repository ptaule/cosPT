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


static void BM_PS_integrand_1loop(benchmark::State& state) {
    // Perform setup here
    Config cfg("ini/benchmark_ps_eds_spt.cfg");

    int n_loops = 1;
    int n_dims = 3 * n_loops - 1;

    LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
    SumTable sum_table(loop_params);
    EvolutionParameters ev_params;
    EtaGrid eta_grid;

    IntegrationInput input(cfg.q_min(), cfg.q_max());

    double twopi_factor = pow(TWOPI, -3);
    input.input_ps = Interpolation1D(cfg.input_ps_file(), twopi_factor, true);

    input.pair_correlations = cfg.pair_correlations();
    input.ps_diagrams = ps::construct_diagrams(loop_params);

    input.tables_vec.emplace_back(cfg.k_a(), loop_params, sum_table, ev_params,
                                  eta_grid);

    int n_correlations = input.pair_correlations.size();

    double* xx = new double[n_dims];
    double* ff = new double[n_correlations];

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



static void BM_PS_integrand_2loop(benchmark::State& state) {
    // Perform setup here
    Config cfg("ini/benchmark_ps_eds_spt.cfg");

    int n_loops = 2;
    int n_dims = 3 * n_loops - 1;

    LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
    SumTable sum_table(loop_params);
    EvolutionParameters ev_params;
    EtaGrid eta_grid;

    IntegrationInput input(cfg.q_min(), cfg.q_max());

    double twopi_factor = pow(TWOPI, -3);
    input.input_ps = Interpolation1D(cfg.input_ps_file(), twopi_factor, true);

    input.pair_correlations = cfg.pair_correlations();
    input.ps_diagrams = ps::construct_diagrams(loop_params);

    input.tables_vec.emplace_back(cfg.k_a(), loop_params, sum_table, ev_params,
                                  eta_grid);

    int n_correlations = input.pair_correlations.size();

    double* xx = new double[n_dims];
    double* ff = new double[n_correlations];

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


static void BM_BS_integrand_1loop(benchmark::State& state) {
    // Perform setup here
    Config cfg("ini/benchmark_bs_eds_spt.cfg");

    int n_loops = 1;
    int n_dims = 3 * n_loops;

    LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
    SumTable sum_table(loop_params);
    EvolutionParameters ev_params;
    EtaGrid eta_grid;

    IntegrationInput input(cfg.q_min(), cfg.q_max());

    double twopi_factor = pow(TWOPI, -3);
    input.input_ps = Interpolation1D(cfg.input_ps_file(), twopi_factor, true);

    input.pair_correlations = cfg.pair_correlations();
    input.bs_diagrams = bs::construct_diagrams(loop_params);

    input.tables_vec.emplace_back(cfg.k_a(), cfg.k_b(), cfg.cos_ab(), loop_params,
            sum_table, ev_params, eta_grid);

    int n_correlations = input.pair_correlations.size();

    double* xx = new double[n_dims];
    double* ff = new double[n_correlations];

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



static void BM_BS_integrand_2loop(benchmark::State& state) {
    // Perform setup here
    Config cfg("ini/benchmark_bs_eds_spt.cfg");

    int n_loops = 2;
    int n_dims = 3 * n_loops;

    LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
    SumTable sum_table(loop_params);
    EvolutionParameters ev_params;
    EtaGrid eta_grid;

    IntegrationInput input(cfg.q_min(), cfg.q_max());

    double twopi_factor = pow(TWOPI, -3);
    input.input_ps = Interpolation1D(cfg.input_ps_file(), twopi_factor, true);

    input.pair_correlations = cfg.pair_correlations();
    input.bs_diagrams = bs::construct_diagrams(loop_params);

    input.tables_vec.emplace_back(cfg.k_a(), cfg.k_b(), cfg.cos_ab(), loop_params,
            sum_table, ev_params, eta_grid);

    int n_correlations = input.pair_correlations.size();

    double* xx = new double[n_dims];
    double* ff = new double[n_correlations];

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



// Register the functions as benchmarks
BENCHMARK(BM_PS_integrand_1loop);
BENCHMARK(BM_PS_integrand_2loop);
BENCHMARK(BM_BS_integrand_1loop);
BENCHMARK(BM_BS_integrand_2loop);
// Run the benchmark
BENCHMARK_MAIN();
