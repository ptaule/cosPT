/*
   benchmark.cpp

   Created by Petter Taule on 05.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/


#include <benchmark/benchmark.h>

extern "C" {
#include "include/constants.h"
#include "include/tables.h"
#include "include/integrand.h"
#include "include/io.h"
#include "include/diagrams.h"
#include "include/integrand.h"
}

static void BM_kernel_index(benchmark::State& state) {
    // Perform setup here
    short int arguments[] = {13, 7, 1, 5, 3};

    for (auto _ : state) {
        // This code gets timed
        kernel_index_from_arguments(arguments);
    }
}


static void BM_integrand(benchmark::State& state) {
    // Perform setup here
    const char input_ps_file[] = "/home/pettertaule/repos/class_public/output/fiducial/newtonian/z1_pk.dat";
    gsl_interp_accel* acc;
    gsl_spline* spline;

    read_and_interpolate(input_ps_file,&acc,&spline);

    short int sum_table[N_CONFIGS][N_CONFIGS];
    compute_sum_table(sum_table);

    double k1 = 0.2;

    integration_variables_t vars = {
        .magnitudes = {0.3, 0.1},
        .cos_theta = {0.5, 0.6},
        .phi = {0.1}
    };

    tables_t tables;
    tables.sum_table = sum_table;

    tables.Q_magnitudes = vars.magnitudes;

    // Initialize diagrams to compute at this order in PT
    diagram_t diagrams[N_DIAGRAMS];
    initialize_diagrams(diagrams);

    integration_input_t input = {
        .k = k1,
        .component_a = 0,
        .component_b = 0,
        .acc         = acc,
        .spline      = spline,
        .diagrams    = diagrams,
        .worker_mem  = NULL
    };

    for (auto _ : state) {
        // This code gets timed
        tables_zero_initialize(&tables);
        compute_bare_scalar_products(k1, &vars, tables.bare_scalar_products);
        compute_alpha_beta_tables(tables.bare_scalar_products, tables.alpha,
                tables.beta);
        integrand(&input, &tables);
    }
    diagrams_gc(diagrams);

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}


// Register the functions as benchmarks
BENCHMARK(BM_kernel_index);
BENCHMARK(BM_integrand);
// Run the benchmark
BENCHMARK_MAIN();
