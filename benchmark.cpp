/*
   benchmark.cpp

   Created by Petter Taule on 28.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <vector>

#include <benchmark/benchmark.h>

#include "include/utilities.hpp"
#include "include/tables.hpp"
#include "include/diagrams.hpp"
#include "include/interpolation.hpp"
#include "include/integrand.hpp"

static void BM_kernel_index(benchmark::State& state) {
    // Perform setup here
    Spectrum spectrum = POWERSPECTRUM;
    Dynamics dynamics = SPT;
    int n_loops = 2;
    Settings settings(n_loops, spectrum, dynamics);

    int arguments[] = {13, 7, 1, 5, 3};

    for (auto _ : state) {
        // This code gets timed
        ps::kernel_index_from_arguments(arguments, settings);
    }
}



static void BM_integrand(benchmark::State& state) {
    // Perform setup here
    int n_loops = 2;
    double k1 = 0.2;

    Interpolation1D input_ps("/home/pettertaule/repos/class_public/output/fiducial/newtonian/z1_pk.dat");

    Vec1D<Correlation> correlations = {{0,0}};

    Settings settings(n_loops, POWERSPECTRUM, SPT);
    SumTable sum_table(settings);

    Vec1D<PowerSpectrumDiagram> diagrams = ps::construct_diagrams(settings);

    Vec1D<IntegrandTables> tables_vec;
    tables_vec.push_back(IntegrandTables(k1, settings, sum_table, Vec1D<double>()));

    IntegrationVariables& vars = tables_vec.at(0).vars;
    vars.magnitudes.at(0) = 0.3;
    vars.magnitudes.at(1) = 0.1;
    vars.cos_theta.at(0) = 0.5;
    vars.cos_theta.at(1) = 0.6;
    vars.phi.at(0) = 0.1;

    IntegrationInput input(0, 0, settings, diagrams, input_ps, correlations,
            tables_vec);

    Vec1D<double> results(correlations.size());
    std::fill(results.begin(), results.end(), 0);

    for (auto _ : state) {
        // This code gets timed
        tables_vec.at(0).reset();
        tables_vec.at(0).compute_tables();
        integrand(input, tables_vec.at(0), results);
    }
}


// Register the functions as benchmarks
BENCHMARK(BM_kernel_index);
BENCHMARK(BM_integrand);
// Run the benchmark
BENCHMARK_MAIN();
