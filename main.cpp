/*
   main.cpp

   Created by Petter Taule on 28.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <vector>
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


int cuba_integrand(
        __attribute__((unused)) const int *ndim,
        const cubareal xx[],
        __attribute__((unused)) const int *ncomp,
        cubareal ff[],
        void *userdata,
        __attribute__((unused)) const int *nvec,
        const int *core
        );



int main () {
    double k_a = 1.047129e-01;

    int n_loops = 1;
    const std::string input_ps_file = "/home/pettertaule/repos/class_public/output/fiducial/newtonian/z1_pk.dat";
    Interpolation1D input_ps(input_ps_file);

    Vec1D<Correlation> correlations = {{0,0}};

    double cuba_epsabs   = 1e-12;
    double cuba_epsrel   = 1e-4;
    double cuba_maxevals = 1e6;
    int cuba_verbose     = 2;
    int n_cores = 4;

    cubacores(n_cores, 10000);

    Settings settings(n_loops, POWERSPECTRUM, SPT);
    SumTable sum_table(settings);
    Vec1D<double> eta_grid;
    Vec1D<IntegrandTables> tables_vec;

    /* (Master + n_cores) instances of IntegrandTables */
    for (int i = 0; i < n_cores + 1; ++i) {
        tables_vec.emplace_back(k_a, settings, sum_table, eta_grid);
    }

    Vec1D<PowerSpectrumDiagram> diagrams = ps::construct_diagrams(settings);

    IntegrationInput input(1e-4, 65, settings, diagrams, input_ps, correlations,
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

    // CUBA settings
#define CUBA_NVEC 1
#define CUBA_LAST 4
#define CUBA_RETAIN_STATEFILE 16
#define CUBA_SEED 0
#define CUBA_MINEVAL 0
#define CUBA_SPIN NULL
#define CUBA_NNEW 1000
#define CUBA_NMIN 2
#define CUBA_FLATNESS 25.
    Suave(3 * n_loops - 1, correlations.size(), (integrand_t)cuba_integrand,
            &input, CUBA_NVEC, cuba_epsrel, cuba_epsabs,
            (cuba_verbose | CUBA_LAST | CUBA_RETAIN_STATEFILE), CUBA_SEED,
            CUBA_MINEVAL, cuba_maxevals, CUBA_NNEW, CUBA_NMIN, CUBA_FLATNESS,
            nullptr, CUBA_SPIN, &nregions, &neval, &fail,
            integration_results.data(), integration_errors.data(),
            integration_probs.data());

    Results results(correlations);

    for (size_t i = 0; i < correlations.size(); ++i) {
        results.lin_ps.at(i)     = input_ps.eval(k_a);
        results.non_lin_ps.at(i) = overall_factor * static_cast<double>(integration_results.at(i));
        results.errors.at(i)     = overall_factor * static_cast<double>(integration_errors.at(i));
    }

    write_results(
            "test.dat",
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
        const int *core
        )
{
    IntegrationInput* input = (IntegrationInput*)userdata;

    /*  For thread <*core + 1> (index 0 is reserved for master), we use the */
    /*  IntegrandTables number *core+1 */
    IntegrandTables& tables = input->tables_vec.at(*core + 1);

    int n_loops = input->settings.n_loops;
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
            vars.phi.at(0) = xx[4] * TWOPI;
            jacobian = TWOPI * xx[0]
                * SQUARE(log_ratio)
                * CUBE(vars.magnitudes[0])
                * CUBE(vars.magnitudes[1]);
            break;
        default:
            throw(std::invalid_argument("n_loops is not 1 or 2."));
    }

    /* Zero-initialize kernel tables */
    tables.reset();
    // Compute sum-, bare_scalar_products-, alpha- and beta-tables
    tables.compute_tables();

    Vec1D<double> results(input->correlations.size(), 0);
    integrand(*input, tables, results);

    for (size_t i = 0; i < input->correlations.size(); ++i) {
        ff[i] = results.at(i) * jacobian;
    }
    return 0;
}
#undef at
