#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include <getopt.h>
#include <gsl/gsl_sf.h>
#include <cuba.h>

#include "include/diagrams.hpp"
#include "include/integrand.hpp"
#include "include/io.hpp"
#include "include/ir_resum.hpp"
#include "include/parameters.hpp"
#include "include/tables.hpp"
#include "include/tree_level.hpp"
#include "include/utilities.hpp"

using std::pow;

void print_help();

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
          {"k_a_idx", required_argument, 0, 'a'},
          {"k_b_idx", required_argument, 0, 'b'},
          {"k_c_idx", required_argument, 0, 'c'},
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
                print_help();
                return EXIT_SUCCESS;
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
        std::cerr << "Configuration file required as argument. "
            << "Use --help to see instructions.\nExiting." << std::endl;
        return EXIT_FAILURE;
    }
    config_file = std::string(argv[optind]);

    try {
        std::cout << "Reading configuration file \"" << config_file << "\"." << std::endl;
        Config cfg(config_file, k_a_idx, k_b_idx, k_c_idx, cuba_maxevals, cuba_cores);

        int n_loops = cfg.n_loops();
        int n_dims = 0; /* Dimension of integral measure */
        size_t n_comp = 0; /* Integrand dimension */

        LoopParameters loop_params(n_loops, cfg.spectrum(), cfg.dynamics(), cfg.rsd());
        SumTable sum_table(loop_params);

        IRresumSettings ir_settings(cfg.k_s(), cfg.k_osc());

        InputPowerSpectrum ps(cfg.input_ps_file(), cfg.input_ps_rescale(),
                cfg.ir_resum(), ir_settings, n_loops, cfg.pt_order(),
                cfg.rsd(), cfg.rsd_growth_f());

        IntegrationInput input(ps, cfg.q_min(), cfg.q_max(), cfg.single_hard_limit());

        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        if (cfg.dynamics() == EVOLVE_IC_EDS) {
            if (COMPONENTS != 2) {
                std::cerr <<
                    "Dynamics = EVOLVE_IC_EDS not implemented for COMPONENTS != 2"
                    << std::endl;
                return EXIT_FAILURE;
            }
            ev_params = EvolutionParameters(cfg.zeta_file(), cfg.ode_atol(),
                    cfg.ode_rtol(), cfg.ode_hstart());
            eta_grid = EtaGrid(cfg.time_steps(), cfg.eta_ini(), cfg.eta_fin());
        }
        else if (cfg.dynamics() == EVOLVE_IC_ASYMP) {
            if (COMPONENTS != 4) {
                std::cerr <<
                    "Dynamics = EVOLVE_IC_EDS not implemented for COMPONENTS != 2"
                    << std::endl;
                return EXIT_FAILURE;
            }
            ev_params = EvolutionParameters(cfg.f_nu(), cfg.omega_m_0(),
                    cfg.zeta_file(), cfg.redshift_file(),
                    cfg.omega_eigenvalues_file(), cfg.F1_ic_files(),
                    cfg.effcs2_x_grid(), cfg.effcs2_y_grid(),
                    cfg.effcs2_data(), cfg.ode_atol(), cfg.ode_rtol(),
                    cfg.ode_hstart());
            eta_grid = EtaGrid(cfg.pre_time_steps(), cfg.time_steps(),
                    cfg.eta_ini(), cfg.eta_fin(), cfg.eta_asymp());
        }

        if (cfg.spectrum() == POWERSPECTRUM) {
            input.pair_correlations = cfg.pair_correlations();

            if (loop_params.rsd()) {
                n_comp = 3; /* Monopole, quadrupole, hexadecapole */
            }
            else {
                n_comp = input.pair_correlations.size();
            }

            if (n_loops > 0) {
                if (loop_params.rsd()) {
                    n_dims = 3 * n_loops + 1;
                }
                else {
                    n_dims = 3 * n_loops - 1;
                }

                input.ps_diagrams = ps::construct_diagrams(loop_params);

                /* (Master + n_cores) instances of IntegrandTables */
                for (int i = 0; i < cfg.cuba_cores() + 1; ++i) {
                    input.tables_vec.emplace_back(cfg.k_a(), 0, 0,
                            cfg.rsd_growth_f(), loop_params, sum_table,
                            ev_params, eta_grid);
                }
            }
        }
        else if (cfg.spectrum() == BISPECTRUM) {
            input.triple_correlations = cfg.triple_correlations();
            n_comp = input.triple_correlations.size();

            if (cfg.rsd()) {
                std::cerr << "RSD is not implemented for bispectrum. Exiting."
                    << std::endl;
                return EXIT_FAILURE;
            }

            if (n_loops > 0) {
                n_dims = 3 * n_loops;
                input.bs_diagrams = bs::construct_diagrams(loop_params);

                /* (Master + n_cores) instances of IntegrandTables */
                for (int i = 0; i < cfg.cuba_cores() + 1; ++i) {
                    input.tables_vec.emplace_back(cfg.k_a(), cfg.k_b(),
                            cfg.cos_ab(), cfg.rsd_growth_f(), loop_params,
                            sum_table, ev_params,
                            eta_grid);
                }
            }
        }
        else {
            throw ConfigException("Unknown spectrum.");
        }

        if (cfg.compute_eft_displacement_dispersion()) {
            /* Compute small scale displacement dispersion,
             * sigma_d^2 = \int_{q_max}^{infty} dq P_{input}(q)
             * useful for EFT parameter RGEs */
            cfg.eft_displacement_dispersion() = ps.integral(cfg.q_max(), 10);
        }

        Vec1D<double> tree_level_result(n_comp, 0);
        Vec1D<double> loop_result(n_comp, 0);
        Vec1D<double> errors(n_comp, 0);

        /* Tree-level results */
        if (cfg.spectrum() == POWERSPECTRUM) {
            ps::tree_level(cfg.k_a(), cfg.dynamics(), ps, eta_grid, ev_params,
                           input.pair_correlations, tree_level_result);
        }
        else if (cfg.spectrum() == BISPECTRUM) {
            /* Tree level bispectrum */
            IntegrandTables tables(cfg.k_a(), cfg.k_b(), cfg.cos_ab(),
                    0, loop_params, sum_table, ev_params, eta_grid);
            bs::tree_level(tables, ps, input.triple_correlations,
                                  tree_level_result);
        }
        else {
            throw ConfigException("Unknown spectrum.");
        }

        /* Single hard limit */
        if (cfg.single_hard_limit()) {
            n_dims -= 1;
            for (auto& el : input.tables_vec) {
                el.vars.magnitudes.at(0) = cfg.sh_Q1();
            }
        }

        if (n_loops > 0) {
            cuba_cores = cfg.cuba_cores();
            int cuba_points = 10000;
            cubacores(&cuba_cores, &cuba_points);
            int cuba_retain_statefile = 0;
            std::string cuba_statefile = cfg.cuba_statefile();
            if (cfg.cuba_retain_statefile()) {
                cuba_retain_statefile = 16;
            }

            Vec1D<cubareal> integration_results(n_comp, 0);
            Vec1D<cubareal> integration_errors(n_comp, 0);
            Vec1D<cubareal> integration_probs(n_comp, 0);
#define CUBA_NVEC 1
#define CUBA_LAST 4
#define CUBA_SEED 0
#define CUBA_MINEVAL 0
#define CUBA_SPIN nullptr
#define CUBA_NNEW 1000
#define CUBA_NMIN 2
#define CUBA_FLATNESS 25.
            Suave(n_dims, static_cast<int>(n_comp), (integrand_t)integrand, &input,
                  CUBA_NVEC, cfg.cuba_rtol(), cfg.cuba_atol(),
                  (cfg.cuba_verbose() | CUBA_LAST | cuba_retain_statefile),
                  CUBA_SEED, CUBA_MINEVAL, cfg.cuba_maxevals(), CUBA_NNEW,
                  CUBA_NMIN, CUBA_FLATNESS,
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
            double overall_factor =
                pow(2, n_loops) * gsl_sf_fact(static_cast<unsigned>(n_loops));
            if (cfg.spectrum() == POWERSPECTRUM && !cfg.rsd()) {
                overall_factor *= TWOPI;
            }

            for (size_t i = 0; i < n_comp; ++i) {
                loop_result.at(i) = overall_factor *
                    static_cast<double>(integration_results.at(i));
                errors.at(i) = overall_factor *
                    static_cast<double>(integration_errors.at(i));

                if (input.single_hard_limit) {
                    loop_result.at(i) *= SQUARE(cfg.sh_Q1());
                    errors.at(i)      *= SQUARE(cfg.sh_Q1());
                }

                if (cfg.rsd()) {
                    loop_result.at(i) *= 4 * static_cast<double>(i) + 1;
                    errors.at(i)      *= 4 * static_cast<double>(i) + 1;
                }
            }

            std::cout << "Integration probability/probabilities: ";
            for (auto& el : integration_probs) {
                std::cout << el << ", ";
            }
            std::cout << std::endl;
        }

        write_results(cfg, tree_level_result, loop_result, errors);
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
    return EXIT_SUCCESS;
}



void print_help() {
    const char* help = R"(CosPT

CosPT computes loop corrections to the power spectrum or bispectrum in
cosmological perturbation theory. The program takes a configuration file as
argument, see README.md for how to enter settings in this file. Some settings
may also be set as a command line options, see below.

Usage:
  cosPT config_file
  cosPT config_file [--k_a_idx=IDX] [--n_evals=NUM] [--n_cores=NUM]
  cosPT config_file [--k_a_idx=IDX] [--k_b_idx=IDX] [--k_c_idx=IDX]
                    [--n_evals=NUM] [--n_cores=NUM]
  cosPT -h | --help

Options:
  -h --help     Show this screen.
  -a --k_a_idx  External wavenumber to compute, index in k_a_grid
                file. Overrides index given in configuration file.
  -b --k_b_idx  Similarly for second external wavenumber (bispectrum).
  -c --k_c_idx  Similarly for third external wavenumber (bispectrum).
  -n --n_evals  Number of evaluations in Monte Carlo integration.
  -p --n_cores  Number of cores/threads for CUBA to spawn, in addition to master.
            )";
    std::cout << help;
    return;
}
