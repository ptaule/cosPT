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
#include "include/tables.hpp"
#include "include/diagrams.hpp"
#include "include/interpolation.hpp"
#include "include/integrand.hpp"


static void BM_kernel_index(benchmark::State& state) {
    // Perform setup here
    Spectrum spectrum = POWERSPECTRUM;
    Dynamics dynamics = EDS_SPT;
    int n_loops = 2;
    LoopParameters loop_params(n_loops, spectrum, dynamics);

    int arguments[] = {22, 16, 10, 14, 12};

    for (auto _ : state) {
        // This code gets timed
        loop_params.arguments_2_kernel_index(arguments);
    }
}



static void BM_ps_integrand(benchmark::State& state) {
    // Perform setup here
    int n_loops = 2;

    double k1 = 0.2;
    double q_min = 1e-4;
    double q_max = 65;

    Interpolation1D input_ps("/home/pettertaule/repos/class_public/output/fiducial/newtonian/z1_pk.dat");

    Vec1D<Pair<int>> pair_correlations = {{0,0}};

    LoopParameters loop_params(n_loops, POWERSPECTRUM, EDS_SPT);
    SumTable sum_table(loop_params);
    EvolutionParameters ev_params;
    EtaGrid eta_grid;

    Vec1D<IntegrandTables> tables_vec;
    tables_vec.emplace_back(k1, loop_params, sum_table, ev_params, eta_grid);

    Vec1D<PowerSpectrumDiagram> diagrams = ps::construct_diagrams(loop_params);

    IntegrationInput input(q_min, q_max, &diagrams, &pair_correlations, input_ps,
            tables_vec);

    double* xx = new double[pair_correlations.size()];
    double* ff = new double[pair_correlations.size()];

    for (size_t i = 0; i < pair_correlations.size(); ++i) {
        xx[i] = 0.5;
        ff[i] = 0;
    }

    for (auto _ : state) {
        // This code gets timed
        int ndim = 3 * n_loops - 1;
        int ncomp = pair_correlations.size();
        int nvec = 1;
        /* core = -1 in accordance with ps::integrand() or bs::integrand() */
        int core = -1;
        ps::integrand(&ndim, xx, &ncomp, ff, &input, &nvec, &core);
    }

    delete[] xx;
    delete[] ff;
}


// Register the functions as benchmarks
/* BENCHMARK(BM_kernel_index); */
BENCHMARK(BM_ps_integrand);
// Run the benchmark
BENCHMARK_MAIN();
