#include <cmath>
#include <fstream>
#include <iostream>

extern "C" {
    #include <getopt.h>
    #include <gsl/gsl_sf.h>
    #include <cuba.h>
}

#include "include/diagrams.hpp"
#include "include/integrand.hpp"
#include "include/io.hpp"
#include "include/ir_resum.hpp"
#include "include/parameters.hpp"
#include "include/tables.hpp"
#include "include/tree_level.hpp"

using std::pow;
using std::string;

void print_help();

int main(int argc, char* argv[]) {
    /* Default values -1 for k_a_idx etc. indicate that they should be set
     * in ini_file */
    int k_a_idx = -1;
    int k_b_idx = -1;
    int k_c_idx = -1;
    int cuba_maxevals = 0;
    int cuba_cores = -1;
    string config_file;
    int stdout_mode = -1; /* Not set (-1), write to file (0), only print to stdout (1) */
    int verbosity_level = -1;

    static struct option long_options[] = {
          {"k_a_idx", required_argument, 0, 'a'},
          {"k_b_idx", required_argument, 0, 'b'},
          {"k_c_idx", required_argument, 0, 'c'},
          {"help",    no_argument,       0, 'h'},
          {"n_evals", required_argument, 0, 'n'},
          {"n_cores", required_argument, 0, 'p'},
          {"stdout_mode", no_argument, 0, 'P'},
          {"verbosity", required_argument, 0, 'v'},
          {NULL, 0, NULL, 0}
        };

    int option_index = 0;
    int c = 0;
    while ((c = getopt_long(argc, argv, "a:b:c:hn:p:Pv:", long_options, &option_index))
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
            case 'P':
                stdout_mode = 1;
                break;
            case 'v':
                verbosity_level = atoi(optarg);
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
    config_file = string(argv[optind]);

    try {
        if (verbosity_level > 0) {
            std::cout << "Reading configuration file \"" << config_file << "\"." << std::endl;
        }
        Config cfg(config_file, k_a_idx, k_b_idx, k_c_idx, cuba_maxevals, cuba_cores);

        /* If stdout_mode/verbosity_level is not set with command line options,
         * use the one from the config file */
        if (stdout_mode == -1) {
            stdout_mode = cfg.get<bool>("stdout_mode");
        }
        if (verbosity_level == -1) {
            verbosity_level = cfg.get<int>("cuba_verbosity_level");
        }

        int n_loops = cfg.get<int>("n_loops");
        int n_dims = 0; /* Dimension of integral measure */
        size_t n_comp = 0; /* Integrand dimension */

        LoopParameters loop_params(n_loops, cfg.get<Spectrum>("spectrum"),
                                   cfg.get<Dynamics>("dynamics"),
                                   cfg.get<bool>("rsd"));
        SumTable sum_table(loop_params);

        IRresumSettings ir_settings(
            n_loops,
            cfg.get<int>("pt_order"),
            cfg.get<double>("k_s"),
            cfg.get<double>("k_osc"),
        );

        InputPowerSpectrum ps(cfg.get<string>("input_ps_file"),
                              cfg.get<double>("input_ps_rescale_num"),
                              cfg.get<bool>("ir_resum"),
                              ir_settings, cfg.get<bool>("rsd"),
                              cfg.get<double>("rsd_growth_f"));

        IntegrationInput input(ps, cfg.get<double>("q_min"),
                               cfg.get<double>("q_max"),
                               cfg.get<bool>("single_hard_limit"));

        EtaGrid eta_grid(cfg.get<double>("eta_ini"),
                         cfg.get<double>("eta_fin"),
                         cfg.get<size_t>("time_steps"),
                         cfg.get<size_t>("pre_time_steps"),
                         cfg.get<double>("eta_asymp"));
        EvolutionParameters ev_params(cfg.kappa(),
                                      cfg.zeta_files(),
                                      cfg.xi_files(),
                                      cfg.get<double>("ode_abs_tolerance"),
                                      cfg.get<double>("ode_rel_tolerance"),
                                      cfg.get<double>("ode_start_step"));
        OmegaEigenspace omega_eigenspace(
            cfg.get<Dynamics>("dynamics"),
            eta_grid.eta_ini(),
            ev_params,
            cfg.get<int>("omega_eigenmode"),
            cfg.get<double>("omega_k_min"),
            cfg.get<double>("omega_k_max"),
            cfg.get<int>("omega_N"),
            cfg.get<double>("omega_imag_threshold")
            );

        if (cfg.get<Spectrum>("spectrum") == POWERSPECTRUM) {
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
                for (int i = 0; i < cfg.get<int>("cuba_n_cores") + 1; ++i) {
                    input.tables_vec.emplace_back(cfg.get<double>("k_a"), 0, 0,
                                                  cfg.get<double>("rsd_growth_f"),
                                                  loop_params, sum_table,
                                                  ev_params, eta_grid,
                                                  omega_eigenspace);
                }
            }
        }
        else if (cfg.get<Spectrum>("spectrum") == BISPECTRUM) {
            input.triple_correlations = cfg.triple_correlations();
            n_comp = input.triple_correlations.size();

            if (cfg.get<bool>("rsd")) {
                std::cerr << "RSD is not implemented for bispectrum. Exiting."
                    << std::endl;
                return EXIT_FAILURE;
            }

            if (n_loops > 0) {
                n_dims = 3 * n_loops;
                input.bs_diagrams = bs::construct_diagrams(loop_params);

                /* (Master + n_cores) instances of IntegrandTables */
                for (int i = 0; i < cfg.get<int>("cuba_n_cores") + 1; ++i) {
                    input.tables_vec.emplace_back(cfg.get<double>("k_a"),
                                                  cfg.get<double>("k_b"),
                                                  cfg.get<double>("cos_ab"),
                                                  cfg.get<double>("rsd_growth_f"),
                                                  loop_params, sum_table,
                                                  ev_params, eta_grid,
                                                  omega_eigenspace);
                }
            }
        }
        else {
            throw ConfigException("Unknown spectrum.");
        }

        if (cfg.get<bool>("compute_eft_displacement_dispersion")) {
            /* Compute small scale displacement dispersion,
             * sigma_d^2 = \int_{q_max}^{infty} dq P_{input}(q)
             * useful for EFT parameter RGEs */
            cfg.set("eft_displacement_dispersion", ps.integral(cfg.get<double>("q_max"), 10));
        }

        Vec1D<double> tree_level_result(n_comp, 0);
        Vec1D<double> loop_result(n_comp, 0);
        Vec1D<double> errors(n_comp, 0);

        /* Tree-level results */
        if (cfg.get<Spectrum>("spectrum") == POWERSPECTRUM) {
            IntegrandTables tables(cfg.get<double>("k_a"), 0, 0,
                    cfg.get<double>("rsd_growth_f"), loop_params, sum_table,
                    ev_params, eta_grid, omega_eigenspace);
            ps::tree_level(tables, ps, input.pair_correlations, tree_level_result);
        }
        else if (cfg.get<Spectrum>("spectrum") == BISPECTRUM) {
            /* Tree level bispectrum */
            IntegrandTables tables(cfg.get<double>("k_a"),
                                   cfg.get<double>("k_b"),
                                   cfg.get<double>("cos_ab"), 0, loop_params,
                                   sum_table, ev_params, eta_grid, omega_eigenspace);
            bs::tree_level(tables, ps, input.triple_correlations,
                           tree_level_result);
        }
        else {
            throw ConfigException("Unknown spectrum.");
        }

        /* Single hard limit */
        if (cfg.get<bool>("single_hard_limit")) {
            n_dims -= 1;
            for (auto& el : input.tables_vec) {
                el.vars.magnitudes.at(0) = cfg.get<double>("sh_Q1");
            }
        }

        if (n_loops > 0) {
            cuba_cores = cfg.get<int>("cuba_n_cores");
            int cuba_points = 10000;
            cubacores(&cuba_cores, &cuba_points);
            int cuba_retain_statefile = 0;
            string cuba_statefile = cfg.get<string>("cuba_statefile");
            if (cfg.get<bool>("cuba_retain_statefile")) {
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
            Suave(n_dims,
                  static_cast<int>(n_comp),
                  (integrand_t)integrand,
                  &input,
                  CUBA_NVEC,
                  cfg.get<double>("cuba_rel_tolerance"),
                  cfg.get<double>("cuba_abs_tolerance"),
                  (verbosity_level | CUBA_LAST | cuba_retain_statefile),
                  CUBA_SEED,
                  CUBA_MINEVAL,
                  cfg.get<int>("cuba_max_evaluations"),
                  CUBA_NNEW,
                  CUBA_NMIN,
                  CUBA_FLATNESS,
                  (cuba_statefile.empty() ? nullptr : cuba_statefile.c_str()),
                  CUBA_SPIN,
                  &cfg.cuba_subregions(),
                  &cfg.cuba_evals(),
                  &cfg.cuba_fail(),
                  integration_results.data(),
                  integration_errors.data(),
                  integration_probs.data()
                  );

            /* Overall factors:
             * - Only integrating over cos_theta_i between 0 and
             *   1, multiply by 2 to obtain [-1,1] (for each loop momenta)
             * - Assuming Q1 > Q2 > ..., hence multiply result by LOOPS factorial
             * - Phi integration of first loop momenta gives a factor 2pi
             *   (powerspectrum) */
            double overall_factor =
                pow(2, n_loops) * gsl_sf_fact(static_cast<unsigned>(n_loops));
            if (cfg.get<Spectrum>("spectrum") == POWERSPECTRUM && !cfg.get<bool>("rsd")) {
                overall_factor *= TWOPI;
            }

            for (size_t i = 0; i < n_comp; ++i) {
                loop_result.at(i) = overall_factor *
                    static_cast<double>(integration_results.at(i));
                errors.at(i) = overall_factor *
                    static_cast<double>(integration_errors.at(i));

                if (input.single_hard_limit) {
                    loop_result.at(i) *= SQUARE(cfg.get<double>("sh_Q1"));
                    errors.at(i)      *= SQUARE(cfg.get<double>("sh_Q1"));
                }

                if (cfg.get<bool>("rsd")) {
                    loop_result.at(i) *= 4 * static_cast<double>(i) + 1;
                    errors.at(i)      *= 4 * static_cast<double>(i) + 1;
                }
            }

            if (verbosity_level > 0) {
                std::cout << "Integration probability/probabilities: ";
                for (auto& el : integration_probs) {
                    std::cout << el << ", ";
                }
                std::cout << std::endl;
            }
        }

        if (stdout_mode == 1) {
            write_results(cfg, tree_level_result, loop_result, errors,
                    std::cout, false);
        } else {
            const std::string filename = cfg.get<std::string>("output_file");
            std::ofstream out(filename);
            if (!out) {
                throw std::runtime_error("Failed to open " + filename);
            }

            write_results(cfg, tree_level_result, loop_result, errors,
                    out, cfg.get<bool>("write_header"));

            if (verbosity_level > 0) {
                std::cout << "Results written to " << filename << "." << std::endl;
            }
        }
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
argument. For details on how to set up the configuration file, see README.md.
Some settings may also be modified via command-line options, as described below.

Usage:
  cosPT config_file
  cosPT config_file [--k_a_idx=IDX] [--n_evals=NUM] [--n_cores=NUM]
  cosPT config_file [--k_a_idx=IDX] [--k_b_idx=IDX] [--k_c_idx=IDX]
                    [--n_evals=NUM] [--n_cores=NUM]
  cosPT -h | --help

Options:
  -h --help        Show this help message.
  -a --k_a_idx     External wavenumber to compute, index in k_a_grid
                   file. Overrides the index specified in the configuration file.
  -b --k_b_idx     Similarly for second external wavenumber (for bispectrum).
  -c --k_c_idx     Similarly for third external wavenumber (for bispectrum).
  -n --n_evals     Number of evaluations for Monte Carlo integration.
  -p --n_cores     Number of cores/threads to spawn for CUBA, in addition to the master.
  -P --stdout_mode Output results to stdout instead of to a file.
  -v --verbosity   Set verbosity level for standard output. Higher values provide
                   more detailed output.

)";
    std::cout << help;
    return;
}
