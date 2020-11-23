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

#include <getopt.h>
#include <gsl/gsl_sf.h>
#include <cuba.h>

#include "include/utilities.hpp"
#include "include/tables.hpp"
#include "include/diagrams.hpp"
#include "include/interpolation.hpp"
#include "include/kernel_evolution.hpp"
#include "include/integrand.hpp"
#include "include/io.hpp"

using std::pow;


int main(int argc, char* argv[]) {
    /* Default values for k_a_idx etc. indicate that they should be set
     * in ini_file */
    int k_a_idx = -1;
    int k_b_idx = -1;
    int k_c_idx = -1;
    int cuba_maxevals = 0;
    int cuba_cores = -1;
    std::string config_file;

    static struct option long_options[] = {
          {"k_a",     required_argument, 0, 'a'},
          {"k_b",     required_argument, 0, 'b'},
          {"k_c",     required_argument, 0, 'c'},
          {"help",    no_argument,       0, 'h'},
          {"n_evals", required_argument, 0, 'n'},
          {"n_cores", required_argument, 0, 'p'},
          {NULL, 0, NULL, 0}
        };

    int option_index = 0;
    int c = 0;
    while ((c = getopt_long(argc, argv, "a:b:c:hn:p:", long_options, &option_index))
            != -1)
    {
        switch (c) {
            case 'a':
                k_a_idx = atoi(optarg);
                break;
            case 'b':
                k_b_idx = atoi(optarg);
                break;
            case 'c':
                k_c_idx = atoi(optarg);
                break;
            case 'h':
                break;
            case 'n':
                cuba_maxevals = static_cast<int>(atof(optarg));
                break;
            case 'p':
                cuba_cores = atoi(optarg);
                break;
            case '?':
            default:
                std::cerr << "Unknown option given." << std::endl;
                return EXIT_FAILURE;
        }
    }

    if (argc != optind + 1) {
        std::cerr << "Configuration file required as argument." << std::endl;
        return EXIT_FAILURE;
    }
    else {
        config_file = std::string(argv[optind]);
    }

    try {
        std::cout << "Reading configuration file " << config_file << std::endl;
        Config cfg(config_file, k_a_idx, k_b_idx, k_c_idx, cuba_maxevals, cuba_cores);

        int n_loops = cfg.n_loops();
        int n_dims = 0; /* Dimension of integral measure */
        int n_correlations = 0; /* Number of correlations to compute */

        LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics());
        SumTable sum_table(loop_params);

        integrand_t integrand;
        IntegrationInput input(cfg.q_min(), cfg.q_max());
        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        /* Conventionally divide input PS by (2pi)^3 */
        double twopi_factor = pow(TWOPI, -3);
        input.input_ps = Interpolation1D(cfg.input_ps_file(), twopi_factor, true);

        if (cfg.dynamics() == EVOLVE_EDS_IC || cfg.dynamics() == EVOLVE_ASYMP_IC) {
            ev_params = EvolutionParameters(cfg.f_nu(), cfg.omega_m_0(),
                    cfg.zeta_file(), cfg.redshift_file(),
                    cfg.omega_eigenvalues_file(), cfg.F1_ic_files(),
                    cfg.effcs2_x_grid(), cfg.effcs2_y_grid(),
                    cfg.effcs2_data(), cfg.ode_atol(), cfg.ode_rtol(),
                    cfg.ode_hstart());
            eta_grid = EtaGrid(cfg.pre_time_steps(), cfg.time_steps(), cfg.eta_ini(),
                    cfg.eta_fin(), cfg.eta_asymp());
        }

        if (cfg.spectrum() == POWERSPECTRUM) {
            n_dims = 3 * n_loops - 1;
            integrand = (integrand_t)ps::integrand;

            input.pair_correlations = cfg.pair_correlations();
            n_correlations = input.pair_correlations.size();
            input.ps_diagrams = ps::construct_diagrams(loop_params);

            /* (Master + n_cores) instances of IntegrandTables */
            for (int i = 0; i < cfg.cuba_cores() + 1; ++i) {
                input.tables_vec.emplace_back(cfg.k_a(), loop_params,
                        sum_table, ev_params, eta_grid);
            }
        }
        else if (cfg.spectrum() == BISPECTRUM) {
            n_dims = 3 * n_loops;
            integrand = (integrand_t)bs::integrand;

            input.triple_correlations = cfg.triple_correlations();
            n_correlations = input.triple_correlations.size();
            input.bs_diagrams = bs::construct_diagrams(loop_params);

            /* (Master + n_cores) instances of IntegrandTables */
            for (int i = 0; i < cfg.cuba_cores() + 1; ++i) {
                input.tables_vec.emplace_back(cfg.k_a(), cfg.k_b(),
                        cfg.cos_ab(), loop_params, sum_table, ev_params,
                        eta_grid);
            }
        }
        else {
            throw ConfigException("Unknown spectrum.");
        }

        Vec1D<double> lin_ps(n_correlations, 0);
        Vec1D<double> non_lin_ps(n_correlations, 0);
        Vec1D<double> errors(n_correlations, 0);

        if (cfg.spectrum() == POWERSPECTRUM) {
            /* Linear power spectrum */
            for (auto& el : lin_ps) {
                el = input.input_ps.eval(cfg.k_a());
            }
            if (cfg.dynamics() == EVOLVE_ASYMP_IC) {
                Vec1D<double> F1_eta_ini(COMPONENTS, 0);
                Vec1D<double> F1_eta_fin(COMPONENTS, 0);
                compute_F1(cfg.k_a(), ev_params, eta_grid, F1_eta_ini, F1_eta_fin);

                for (int i = 0; i < n_correlations; ++i) {
                  lin_ps.at(i) *=
                      F1_eta_fin.at(input.pair_correlations.at(i).first()) *
                      F1_eta_fin.at(input.pair_correlations.at(i).second());
                }
            }
        }


        int cuba_cores = cfg.cuba_cores();
        int cuba_points = 10000;
        cubacores(&cuba_cores, &cuba_points);
        int cuba_retain_statefile = 0;
        std::string cuba_statefile = cfg.cuba_statefile();
        if (cfg.cuba_retain_statefile()) {
            cuba_retain_statefile = 16;
        }

        Vec1D<cubareal> integration_results(n_correlations, 0);
        Vec1D<cubareal> integration_errors(n_correlations, 0);
        Vec1D<cubareal> integration_probs(n_correlations, 0);
#define CUBA_NVEC 1
#define CUBA_LAST 4
#define CUBA_SEED 0
#define CUBA_MINEVAL 0
#define CUBA_SPIN nullptr
#define CUBA_NNEW 1000
#define CUBA_NMIN 2
#define CUBA_FLATNESS 25.
        Suave(n_dims, n_correlations, integrand, &input, CUBA_NVEC,
                cfg.cuba_rtol(), cfg.cuba_atol(),
                (cfg.cuba_verbose() | CUBA_LAST | cuba_retain_statefile), CUBA_SEED,
                CUBA_MINEVAL, cfg.cuba_maxevals(), CUBA_NNEW, CUBA_NMIN, CUBA_FLATNESS,
                (cuba_statefile.empty() ? nullptr : cuba_statefile.c_str()),
                CUBA_SPIN, &cfg.cuba_subregions(), &cfg.cuba_evals(),
                &cfg.cuba_fail(), integration_results.data(),
                integration_errors.data(), integration_probs.data());

        /* Overall factors:
         * - Only integrating over cos_theta_i between 0 and
         *   1, multiply by 2 to obtain [-1,1] (for each loop momenta)
         * - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
         * - Phi integration of first loop momenta gives a factor 2pi
         *   (powerspectrum) */
        double overall_factor = pow(2, n_loops) * gsl_sf_fact(n_loops) *
                                (cfg.spectrum() == POWERSPECTRUM ? TWOPI : 1);

        for (size_t i = 0; i < static_cast<size_t>(n_correlations); ++i) {
            non_lin_ps.at(i) =
                overall_factor * static_cast<double>(integration_results.at(i));
            errors.at(i) =
                overall_factor * static_cast<double>(integration_errors.at(i));
        }

        std::cout << "Integration probability/probabilities: ";
        for (auto& el : integration_probs) {
            std::cout << el << ", ";
        }
        std::cout << std::endl;

        write_results(cfg, lin_ps, non_lin_ps, errors);
        std::cout << "Results written to " << cfg.output_file() << "." << std::endl;
    }
    catch (const ConfigException& e) {
        std::cerr << e.what() << "\nExiting." << std::endl;
        return EXIT_FAILURE;
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\nExiting." << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}
