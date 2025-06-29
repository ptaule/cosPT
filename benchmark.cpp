#include <iostream>

#include <benchmark/benchmark.h>

#include "include/diagrams.hpp"
#include "include/integrand.hpp"
#include "include/ir_resum.hpp"
#include "include/parameters.hpp"
#include "include/tables.hpp"
#include "include/utilities.hpp"

using std::size_t;
using std::string;

void run_integrand_benchmark(
    benchmark::State& state,
    int n_loops,
    Spectrum spectrum,
    const std::string& config_path
) {
    try {
        Config cfg(config_path, -1, -1, -1);

        bool rsd = cfg.get<bool>("rsd");

        LoopParameters loop_params(n_loops, cfg.get<Spectrum>("spectrum"),
                                   cfg.get<Dynamics>("dynamics"),
                                   rsd);
        SumTable sum_table(loop_params);

        IRresumSettings ir_settings(
            n_loops,
            cfg.get<int>("pt_order"),
            cfg.get<double>("k_s"),
            cfg.get<double>("k_osc"),
            1
        );

        InputPowerSpectrum ps(cfg.get<std::string>("input_ps_file"),
                              cfg.get<double>("input_ps_rescale_num"),
                              cfg.get<bool>("ir_resum"),
                              ir_settings, rsd,
                              rsd);

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

        int n_dims = 0; /* Dimension of integral measure */
        int n_correlations = 0;
        switch (spectrum) {
            case POWERSPECTRUM:
                input.pair_correlations = cfg.pair_correlations();
                input.ps_diagrams = ps::construct_diagrams(loop_params);
                input.tables_vec.emplace_back(
                    cfg.get<double>("k_a"), 0, 0,
                    cfg.get<double>("rsd_growth_f"),
                    loop_params, sum_table,
                    ev_params, eta_grid, omega_eigenspace);

                if (rsd) {
                    n_dims = 3 * n_loops + 1;
                    /* (not correlations, but monopole, quadrupole, hexadecapole */
                    n_correlations = 3;
                }
                else {
                    n_dims = 3 * n_loops - 1;
                    n_correlations = static_cast<int>(input.pair_correlations.size());
                }

                break;
            case BISPECTRUM:
                n_dims = 3 * n_loops;
                input.triple_correlations = cfg.triple_correlations();
                n_correlations = static_cast<int>(input.triple_correlations.size());
                input.bs_diagrams = bs::construct_diagrams(loop_params);
                input.tables_vec.emplace_back(
                    cfg.get<double>("k_a"),
                    cfg.get<double>("k_b"),
                    cfg.get<double>("cos_ab"),
                    cfg.get<double>("rsd_growth_f"),
                    loop_params, sum_table,
                    ev_params, eta_grid, omega_eigenspace);
                break;
        }

        std::vector<double> xx(static_cast<size_t>(n_dims), 0.5);
        std::vector<double> ff(static_cast<size_t>(n_correlations));

        for (auto _ : state) {
            int nvec = 1;
            int core = -1;
            integrand(&n_dims, xx.data(), &n_correlations,
                      ff.data(), &input, &nvec,
                      &core);
        }
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
    }
}



#define REGISTER_BENCHMARK(NAME, N_LOOPS, spectrum, CONFIG_FILE) \
    static void NAME(benchmark::State& state) { \
        run_integrand_benchmark(state, N_LOOPS, spectrum, CONFIG_FILE); \
    } \
    BENCHMARK(NAME);

REGISTER_BENCHMARK(BM_PS_EdS_integrand_1loop, 1, POWERSPECTRUM, "benchmark/ini/ps_eds_spt.cfg")
REGISTER_BENCHMARK(BM_PS_EdS_integrand_2loop, 2, POWERSPECTRUM, "benchmark/ini/ps_eds_spt.cfg")
REGISTER_BENCHMARK(BM_PS_RSD_EdS_integrand_1loop, 1, POWERSPECTRUM, "benchmark/ini/ps_eds_spt_rsd.cfg")
REGISTER_BENCHMARK(BM_PS_RSD_EdS_integrand_2loop, 2, POWERSPECTRUM, "benchmark/ini/ps_eds_spt_rsd.cfg")
REGISTER_BENCHMARK(BM_BS_EdS_integrand_1loop, 1, BISPECTRUM,   "benchmark/ini/bs_eds_spt.cfg")
REGISTER_BENCHMARK(BM_BS_EdS_integrand_2loop, 2, BISPECTRUM,   "benchmark/ini/bs_eds_spt.cfg")
REGISTER_BENCHMARK(BM_PS_2fluid_integrand_1loop, 1, POWERSPECTRUM, "benchmark/ini/quijote_Mnu_0.1eV.cfg")
REGISTER_BENCHMARK(BM_PS_2fluid_integrand_2loop, 2, POWERSPECTRUM, "benchmark/ini/quijote_Mnu_0.1eV.cfg")

// Run the benchmark
BENCHMARK_MAIN();
