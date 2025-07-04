#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

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
using std::exp;
using std::string;

struct CLIOptions {
    /* Default values -1 for k_a_idx etc. indicate that they should be set in ini_file */
    int k_a_idx = -1;
    int k_b_idx = -1;
    int k_c_idx = -1;
    int cuba_maxevals = -1;
    int cuba_n_cores = -1;
    int stdout_mode = -1;
    int verbosity = -1;
    std::string config_file;
};

struct CalculationSetup {
    Spectrum spectrum;
    Dynamics dynamics;
    bool rsd;

    int n_loops;
    int n_dims;
    std::size_t n_components;

    LoopStructure loop_structure;
    SumTable sum_table;

    EtaGrid eta_grid;
    EvolutionParameters evolution_params;
    OmegaEigenspace omega_eigenspace;

    double growth_factor_multiplier;
    IRresumSettings ir_settings;
    InputPowerSpectrum ps;

    CalculationSetup(const Config& cfg);
};

void print_help();
CLIOptions parse_cli_args(int argc, char* argv[]);
IntegrationInput initialize_integration_input(const Config& cfg,
                                              const CalculationSetup& setup);
void run_loop_integration(Config& cfg, const CalculationSetup& setup,
                          IntegrationInput& input, Vec1D<double>& loop_result,
                          Vec1D<double>& errors);
void write_output(const Config& cfg, const Vec1D<double>& tree_level_result,
                  const Vec1D<double>& loop_result, const Vec1D<double>& errors,
                  int stdout_mode, int verbosity);
void compute_tree_level(const CalculationSetup& setup, const Config& cfg, const
                        IntegrationInput& input, std::vector<double>&
                        tree_level_result);
void small_scale_dispersion(Config& cfg, const CalculationSetup& setup);



int main(int argc, char* argv[]) {
    Config cfg;
    CLIOptions opts;
    try {
        opts = parse_cli_args(argc, argv);
        cfg = Config(opts.config_file, opts.k_a_idx,
                     opts.k_b_idx, opts.k_c_idx);

        // If read from command line, overwrite config values
        if (opts.stdout_mode != -1)   cfg.set("stdout_mode", opts.stdout_mode);
        if (opts.verbosity != -1)     cfg.set("cuba_verbosity_level", opts.verbosity);
        if (opts.cuba_n_cores != -1)  cfg.set("cuba_n_cores", opts.cuba_n_cores);
        if (opts.cuba_maxevals != -1) cfg.set("cuba_max_evaluations", opts.cuba_maxevals);

        CalculationSetup setup(cfg);
        IntegrationInput input = initialize_integration_input(cfg, setup);

        /* Tree-level result */
        std::vector<double> tree_level_result(setup.n_components, 0);
        compute_tree_level(setup, cfg, input, tree_level_result);

        /* Loop-level result */
        std::vector<double> loop_result(setup.n_components, 0);
        std::vector<double> errors(setup.n_components, 0);
        run_loop_integration(cfg, setup, input, loop_result, errors);

        /* Compute small scale displacement dispersion useful for EFT parameter RGE,
         * sigma_d^2 = \int_{q_max}^{infty} dq P_{input}(q) */
        small_scale_dispersion(cfg, setup);

        write_output(cfg, tree_level_result, loop_result, errors,
                     opts.stdout_mode, opts.verbosity);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\nExiting." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}



CLIOptions parse_cli_args(int argc, char* argv[]) {
    CLIOptions opts;

    static struct option long_options[] = {
        {"k_a_idx",      required_argument, 0, 'a'},
        {"k_b_idx",      required_argument, 0, 'b'},
        {"k_c_idx",      required_argument, 0, 'c'},
        {"help",         no_argument,       0, 'h'},
        {"n_evals",      required_argument, 0, 'n'},
        {"n_cores",      required_argument, 0, 'p'},
        {"stdout_mode",  no_argument,       0, 'P'},
        {"verbosity",    required_argument, 0, 'v'},
        {nullptr, 0, nullptr, 0}
    };

    int option_index = 0;
    int c = 0;
    while ((
        c = getopt_long(argc, argv, "a:b:c:hn:p:Pv:",
                            long_options, &option_index))
        != -1)
    {
        switch (c) {
            case 'a': opts.k_a_idx = std::stoi(optarg); break;
            case 'b': opts.k_b_idx = std::stoi(optarg); break;
            case 'c': opts.k_c_idx = std::stoi(optarg); break;
            case 'h':
                print_help();
                std::exit(EXIT_SUCCESS);
            case 'n': opts.cuba_maxevals = static_cast<int>(std::stod(optarg)); break;
            case 'p': opts.cuba_n_cores = std::stoi(optarg); break;
            case 'P': opts.stdout_mode = 1; break;
            case 'v': opts.verbosity = std::stoi(optarg); break;
            case '?':
            default:
                throw std::invalid_argument("Unknown or malformed CLI option");
        }
    }

    if (optind >= argc) {
        throw std::invalid_argument("Configuration file is required.\nUse --help to see usage.");
    }

    opts.config_file = argv[optind];
    return opts;
}



CalculationSetup::CalculationSetup(const Config& cfg) :
    spectrum(cfg.get<Spectrum>("spectrum")),
    dynamics(cfg.get<Dynamics>("dynamics")),
    rsd(cfg.get<bool>("rsd")),
    n_loops(cfg.get<int>("n_loops")),
    loop_structure(n_loops, spectrum),
    sum_table(loop_structure),
    eta_grid(cfg.get<double>("eta_ini"),
             cfg.get<double>("eta_fin"),
             cfg.get<size_t>("time_steps"),
             cfg.get<size_t>("pre_time_steps"),
             cfg.get<double>("eta_asymp")),
    evolution_params(cfg.kappa(), cfg.zeta_files(),
                     cfg.xi_files(),
                     cfg.get<double>("ode_abs_tolerance"),
                     cfg.get<double>("ode_rel_tolerance"),
                     cfg.get<double>("ode_start_step")),
    omega_eigenspace( dynamics, eta_grid.eta_ini(),
                     evolution_params,
                     cfg.get<int>("omega_eigenmode"),
                     cfg.get<double>("omega_k_min"),
                     cfg.get<double>("omega_k_max"),
                     cfg.get<int>("omega_N"),
                     cfg.get<double>("omega_imag_threshold")),
    growth_factor_multiplier(
        (dynamics == EVOLVE_EDS_ICS || dynamics == EVOLVE_ASYMPTOTIC_ICS) ?
            std::exp(eta_grid.eta_fin() - eta_grid.eta_ini()) : 1
    ),
    ir_settings( n_loops, cfg.get<int>("pt_order"),
                cfg.get<double>("k_s"),
                cfg.get<double>("k_osc"),
                growth_factor_multiplier),
    ps( cfg.get<std::string>("input_ps_file"),
       cfg.get<double>("input_ps_rescale_num"),
       cfg.get<bool>("ir_resum"), ir_settings, rsd,
       cfg.get<double>("rsd_growth_f"))
{
    // Determine output dimensionality
    if (spectrum == POWERSPECTRUM) {
        n_components = rsd ? 3 : cfg.pair_correlations().size();
        if (n_loops > 0) {
            n_dims = rsd ? 3 * n_loops + 1 : 3 * n_loops - 1;
        }
    } else if (spectrum == BISPECTRUM) {
        if (rsd) {
            throw std::runtime_error("RSD is not implemented for bispectrum.");
        }
        if (cfg.get<bool>("biased_tracers")) {
            throw std::runtime_error("Biased_tracers are not implemented for bispectrum.");
        }
        n_components = cfg.triple_correlations().size();
        n_dims = 3 * n_loops;
    } else {
        throw std::runtime_error("Unknown spectrum type.");
    }
}



IntegrationInput initialize_integration_input(const Config& cfg, const CalculationSetup& setup) {
    IntegrationInput input(
        setup.ps,
        cfg.get<double>("q_min"),
        cfg.get<double>("q_max"),
        cfg.get<bool>("single_hard_limit")
    );

    if (setup.spectrum == POWERSPECTRUM) {
        input.pair_correlations = cfg.pair_correlations();
        if (setup.n_loops > 0) {
            input.ps_diagrams = ps::construct_diagrams(setup.loop_structure);
        }
    } else if (setup.spectrum == BISPECTRUM) {
        input.triple_correlations = cfg.triple_correlations();
        if (setup.n_loops > 0) {
            input.bs_diagrams = bs::construct_diagrams(setup.loop_structure);
        }
    }

    // Create (master + worker cores) IntegrandTables
    for (int i = 0; i < cfg.get<int>("cuba_n_cores") + 1; ++i) {
        input.tables_vec.emplace_back(
            cfg.get<double>("k_a"),
            (setup.spectrum == BISPECTRUM ? cfg.get<double>("k_b") : 0.0),
            (setup.spectrum == BISPECTRUM ? cfg.get<double>("cos_ab") : 0.0),
            setup.rsd,
            cfg.get<double>("rsd_growth_f"),
            setup.dynamics,
            setup.loop_structure,
            setup.sum_table,
            setup.evolution_params,
            setup.eta_grid,
            setup.omega_eigenspace
        );
    }

    return input;
}



IntegrandTables make_tables_for_tree_level(const Config& cfg, const CalculationSetup& setup) {
    return IntegrandTables(
        cfg.get<double>("k_a"),
        (setup.spectrum == BISPECTRUM ? cfg.get<double>("k_b") : 0.0),
        (setup.spectrum == BISPECTRUM ? cfg.get<double>("cos_ab") : 0.0),
        setup.rsd,
        cfg.get<double>("rsd_growth_f"),
        setup.dynamics,
        setup.loop_structure,
        setup.sum_table,
        setup.evolution_params,
        setup.eta_grid,
        setup.omega_eigenspace
    );
}



void compute_tree_level(const CalculationSetup& setup, const Config& cfg,
                        const IntegrationInput& input, std::vector<double>& tree_level_result) {
    if (setup.spectrum == POWERSPECTRUM) {
        IntegrandTables tables(cfg.get<double>("k_a"), 0, 0, setup.rsd,
                               cfg.get<double>("rsd_growth_f"),
                               setup.dynamics, setup.loop_structure,
                               setup.sum_table, setup.evolution_params,
                               setup.eta_grid, setup.omega_eigenspace);
        ps::tree_level(tables, setup.ps, input.pair_correlations, tree_level_result);
    } else if (setup.spectrum == BISPECTRUM) {
        IntegrandTables tables(cfg.get<double>("k_a"),
                               cfg.get<double>("k_b"),
                               cfg.get<double>("cos_ab"), false, 0,
                               setup.dynamics,
                               setup.loop_structure, setup.sum_table,
                               setup.evolution_params, setup.eta_grid,
                               setup.omega_eigenspace);
        bs::tree_level(tables, setup.ps, input.triple_correlations, tree_level_result);
    }
    for (auto& el : tree_level_result) {
        el *= SQUARE(setup.growth_factor_multiplier);
    }
}



void run_loop_integration(
    Config& cfg,
    const CalculationSetup& setup,
    IntegrationInput& input,
    Vec1D<double>& loop_result,
    Vec1D<double>& errors
) {
    int cuba_n_cores = cfg.get<int>("cuba_n_cores");
    int cuba_maxevals = cfg.get<int>("cuba_max_evaluations");
    int cuba_points = 10000;

    cubacores(&cuba_n_cores, &cuba_points);

    int cuba_retain_statefile = cfg.get<bool>("cuba_retain_statefile") ? 16 : 0;
    std::string cuba_statefile = cfg.get<std::string>("cuba_statefile");
    int verbosity_level = cfg.get<int>("cuba_verbosity_level");

    std::size_t n_comp = setup.n_components;
    int n_dims = setup.n_dims;

    // Single hard limit
    if (cfg.get<bool>("single_hard_limit")) {
        --n_dims;
        for (auto& el : input.tables_vec) {
            el.vars.magnitudes.at(0) = cfg.get<double>("sh_Q1");
        }
    }

    Vec1D<cubareal> integration_results(n_comp, 0);
    Vec1D<cubareal> integration_errors(n_comp, 0);
    Vec1D<cubareal> integration_probs(n_comp, 0);

    /* Silence gcc warning converting integrand to integrand_t. integrand has
     * nvec and core as additional arguments compared to integrand_t (defined
     * in cuba.h), but this is fine according to the documentation of
     * (parallel) CUBA. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
    Suave(n_dims,
          static_cast<int>(n_comp),
          (integrand_t)integrand,
          &input,
          1, // CUBA_NVEC
          cfg.get<double>("cuba_rel_tolerance"),
          cfg.get<double>("cuba_abs_tolerance"),
          (verbosity_level | 4 | cuba_retain_statefile), // CUBA_LAST = 4
          0, // CUBA_SEED
          0, // CUBA_MINEVAL
          cuba_maxevals,
          1000, 2, 25., // CUBA_NNEW, NMIN, FLATNESS
          (cuba_statefile.empty() ? nullptr : cuba_statefile.c_str()),
          nullptr,
          &cfg.cuba_subregions(),
          &cfg.cuba_evals(),
          &cfg.cuba_fail(),
          integration_results.data(),
          integration_errors.data(),
          integration_probs.data());
#pragma GCC diagnostic pop

    // Compute overall prefactor
    double factor = pow(2, setup.n_loops) *
        gsl_sf_fact(static_cast<unsigned int>(setup.n_loops));
    if (setup.spectrum == POWERSPECTRUM && !setup.rsd) {
        factor *= TWOPI;
    }

    // Convert results
    loop_result.resize(n_comp);
    errors.resize(n_comp);
    for (size_t i = 0; i < n_comp; ++i) {
        double result = factor * static_cast<double>(integration_results.at(i));
        double error  = factor * static_cast<double>(integration_errors.at(i));

        /* Relevant for dynamical results, otherwise growth_multiplier is 1 */
        double gpow = pow(setup.growth_factor_multiplier, 2 * (setup.n_loops + 1));
        result *= gpow;
        error  *= gpow;

        if (input.single_hard_limit) {
            double shQ1_sq = SQUARE(cfg.get<double>("sh_Q1"));
            result *= shQ1_sq;
            error  *= shQ1_sq;
        }

        if (setup.rsd) {
            double l_factor = 4 * static_cast<double>(i) + 1;
            result *= l_factor;
            error  *= l_factor;
        }

        loop_result[i] = result;
        errors[i] = error;
    }

    if (verbosity_level > 0) {
        std::cout << "Integration probability/probabilities: ";
        for (auto& p : integration_probs) std::cout << p << ", ";
        std::cout << "\n";
    }
}



void small_scale_dispersion(Config& cfg, const CalculationSetup& setup) {
    if (cfg.get<bool>("compute_eft_displacement_dispersion")) {
        constexpr int eft_integral_limit = 10;
        cfg.set("eft_displacement_dispersion",
                setup.ps.integral(cfg.get<double>("q_max"), eft_integral_limit));
    }
}



void write_output(
    const Config& cfg,
    const Vec1D<double>& tree_result,
    const Vec1D<double>& loop_result,
    const Vec1D<double>& errors,
    int stdout_mode,
    int verbosity_level
) {
    if (stdout_mode == 1) {
        write_results(cfg, tree_result,
                      loop_result, errors,
                      std::cout, false);
    } else {
        const std::string filename = cfg.get<std::string>("output_file");
        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Failed to open " + filename);
        }

        write_results(cfg, tree_result,
                      loop_result, errors, out,
                      cfg.get<bool>("write_header"));

        if (verbosity_level > 0) {
            std::cout << "Results written to " << filename << "." << std::endl;
        }
    }
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
