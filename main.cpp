/*
   main.cpp

   Created by Petter Taule on 28.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_sf.h>
#include <cuba.h>

#include "include/utilities.hpp"
#include "include/tables.hpp"
#include "include/diagrams.hpp"
#include "include/interpolation.hpp"
#include "include/integrand.hpp"
#include "include/io.hpp"

using std::pow;


int main () {
    double k_a = 1.047129e-01;
    /* double k_b = 1.047129e-02; */
    /* double cos_ab = 1; */
    double q_min = 1e-4;
    double q_max = 65;

    int n_loops = 1;
    Spectrum spectrum = POWERSPECTRUM;
    const std::string input_ps_file =
        "/home/pettertaule/repos/class_public/output/fiducial/newtonian/z1_pk.dat";
    const std::string output_file = "test.dat";

    Interpolation1D input_ps(input_ps_file);

    Vec1D<Pair<int>> correlations = {{0,0}};

    double cuba_epsabs   = 1e-12;
    double cuba_epsrel   = 1e-4;
    double cuba_maxevals = 1e6;
    int cuba_verbose     = 2;
    int n_cores = 4;

    cubacores(n_cores, 10000);

    LoopParameters loop_params(n_loops, spectrum, EDS_SPT);
    SumTable sum_table(loop_params);
    EvolutionParameters ev_params;
    EtaGrid eta_grid;

    Vec1D<IntegrandTables> tables_vec;

    /* (Master + n_cores) instances of IntegrandTables */
    for (int i = 0; i < n_cores + 1; ++i) {
        tables_vec.emplace_back(k_a, loop_params, sum_table, ev_params, eta_grid);
    }

    Vec1D<PowerSpectrumDiagram> diagrams = ps::construct_diagrams(loop_params);

    IntegrationInput input(q_min, q_max, &diagrams, &correlations, input_ps,
            tables_vec);

    /* Non-linear evolution */
    // Overall factors:
    // - Only integrating over cos_theta_i between 0 and 1, multiply by 2 to
    //   obtain [-1,1] (for each loop momenta)
    // - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
    // - Phi integration of first loop momenta gives a factor 2pi
    // - Conventionally divide by ((2pi)^3)^(LOOPS)
    double overall_factor =
        pow(2, n_loops) * gsl_sf_fact(n_loops) * pow(TWOPI, 1 - 3*n_loops);

    int nregions, neval, fail;
    Vec1D<cubareal> integration_results(correlations.size(),0);
    Vec1D<cubareal> integration_errors(correlations.size(),0);
    Vec1D<cubareal> integration_probs(correlations.size(),0);

    /* CUBA settings */
    int n_dims = 0;
    integrand_t integrand;
    if (spectrum == POWERSPECTRUM) {
        n_dims = 3 * n_loops - 1;
        integrand = (integrand_t)ps::integrand;
    }
    else if (spectrum == BISPECTRUM) {
        n_dims = 3 * n_loops;
        integrand = (integrand_t)bs::integrand;
    }
    else {
        std::cerr << "main(): unknown spectrum." << std::endl;
        return EXIT_FAILURE;
    }
#define CUBA_NVEC 1
#define CUBA_LAST 4
#define CUBA_RETAIN_STATEFILE 16
#define CUBA_SEED 0
#define CUBA_MINEVAL 0
#define CUBA_SPIN nullptr
#define CUBA_NNEW 1000
#define CUBA_NMIN 2
#define CUBA_FLATNESS 25.
    Suave(n_dims, correlations.size(), integrand, &input, CUBA_NVEC,
          cuba_epsrel, cuba_epsabs,
          (cuba_verbose | CUBA_LAST | CUBA_RETAIN_STATEFILE), CUBA_SEED,
          CUBA_MINEVAL, cuba_maxevals, CUBA_NNEW, CUBA_NMIN, CUBA_FLATNESS,
          nullptr, CUBA_SPIN, &nregions, &neval, &fail,
          integration_results.data(), integration_errors.data(),
          integration_probs.data());

    Results results(spectrum, correlations);

    for (size_t i = 0; i < correlations.size(); ++i) {
        results.lin_ps.at(i) = input_ps.eval(k_a);
        results.non_lin_ps.at(i) =
            overall_factor * static_cast<double>(integration_results.at(i));
        results.errors.at(i) =
            overall_factor * static_cast<double>(integration_errors.at(i));
    }

    write_results(
            output_file,
            input_ps_file,
            "",
            n_loops,
            cuba_epsabs,
            cuba_epsrel,
            cuba_maxevals,
            k_a,
            input.q_min,
            input.q_max,
            results
            );

    return 0;
}
