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

using cubareal = double;

int cuba_integrand(
        __attribute__((unused)) const int *ndim,
        const cubareal xx[],
        __attribute__((unused)) const int *ncomp,
        cubareal ff[],
        void *userdata,
        __attribute__((unused)) const int *nvec,
        const int *core
        );

static void BM_kernel_index(benchmark::State& state) {
    // Perform setup here
    Spectrum spectrum = POWERSPECTRUM;
    Dynamics dynamics = SPT;
    int n_loops = 2;
    LoopParameters loop_params(n_loops, spectrum, dynamics);

    int arguments[] = {22, 16, 10, 14, 12};

    for (auto _ : state) {
        // This code gets timed
        loop_params.arguments_2_kernel_index(arguments);
    }
}



static void BM_integrand(benchmark::State& state) {
    // Perform setup here
    int n_loops = 2;

    double k1 = 0.2;
    double q_min = 1e-4;
    double q_max = 65;

    Interpolation1D input_ps("/home/pettertaule/repos/class_public/output/fiducial/newtonian/z1_pk.dat");

    Vec1D<Correlation> correlations = {{0,0}};

    LoopParameters loop_params(n_loops, POWERSPECTRUM, SPT);
    SumTable sum_table(loop_params);
    EvolutionParameters ev_params;
    EtaGrid eta_grid;

    Vec1D<IntegrandTables> tables_vec;
    tables_vec.emplace_back(k1, loop_params, sum_table, ev_params, eta_grid);

    Vec1D<PowerSpectrumDiagram> diagrams = ps::construct_diagrams(loop_params);

    IntegrationInput input(q_min, q_max, &diagrams, correlations, input_ps,
            tables_vec);

    cubareal* xx = new cubareal[correlations.size()];
    cubareal* ff = new cubareal[correlations.size()];

    for (size_t i = 0; i < correlations.size(); ++i) {
        xx[i] = 0.5;
        ff[i] = 0;
    }

    for (auto _ : state) {
        // This code gets timed
        int ndim = 3 * n_loops - 1;
        int ncomp = correlations.size();
        int nvec = 1;
        int core = 0;
        cuba_integrand(&ndim, xx, &ncomp, ff, &input, &nvec, &core);
    }

    delete[] xx;
    delete[] ff;
}


/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif
int cuba_integrand(
        __attribute__((unused)) const int *ndim,
        const cubareal xx[],
        __attribute__((unused)) const int *ncomp,
        cubareal ff[],
        void *userdata,
        __attribute__((unused)) const int *nvec,
        __attribute__((unused)) const int *core
        )
{
    IntegrationInput* input = (IntegrationInput*)userdata;

    /*  For thread <*core + 1> (index 0 is reserved for master), we use the */
    /*  IntegrandTables number *core+1 */
    IntegrandTables& tables = input->tables_vec.at(0);

    int n_loops = tables.loop_params.get_n_loops();
    IntegrationVariables& vars = tables.vars;

    double ratio = input->q_max/input->q_min;
    double log_ratio = log(ratio);
    double jacobian = 0.0;

    switch (n_loops) {
        case 1:
            vars.magnitudes.at(0) = input->q_min * pow(ratio,xx[0]);
            vars.cos_theta.at(0) = xx[1];
            jacobian = log(ratio) * CUBE(vars.magnitudes[0]);
            break;
        case 2:
            vars.magnitudes.at(0) = input->q_min * pow(ratio,xx[0]);
            vars.magnitudes.at(1) = input->q_min * pow(ratio,xx[0] * xx[1]);
            vars.cos_theta.at(0) = xx[2];
            vars.cos_theta.at(1) = xx[3];
            /* We may fix the coordinate system s.t. vars.phi[0] = 0 */
            vars.phi.at(1) = xx[4] * TWOPI;
            jacobian = TWOPI * xx[0]
                * SQUARE(log_ratio)
                * CUBE(vars.magnitudes[0])
                * CUBE(vars.magnitudes[1]);
            break;
        default:
            throw(std::invalid_argument("n_loops is not 1 or 2."));
    }

    Vec1D<double> results(input->correlations.size(), 0.0);
    try {
        /* Zero-initialize kernel tables */
        tables.reset();
        // Compute scalar_products-, alpha- and beta-tables
        tables.compute_tables();

        ps::integrand(*input, tables, results);

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return -999;
    }

    for (size_t i = 0; i < input->correlations.size(); ++i) {
        ff[i] = results.at(i) * jacobian;
    }
    return 0;
}
#undef at


// Register the functions as benchmarks
/* BENCHMARK(BM_kernel_index); */
BENCHMARK(BM_integrand);
// Run the benchmark
BENCHMARK_MAIN();
