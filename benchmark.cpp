#include <iostream>
#include <cmath>

#include <benchmark/benchmark.h>

#include "include/diagrams.hpp"
#include "include/integrand.hpp"
#include "include/ir_resum.hpp"
#include "include/parameters.hpp"
#include "include/tables.hpp"
#include "include/utilities.hpp"

using std::size_t;
using std::string;

static void BM_PS_EdS_integrand_1loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("benchmark/ini/ps_eds_spt.cfg", -1, -1, -1, 0, 0);

        int n_loops = 1;
        int n_dims = 3 * n_loops - 1;

        LoopParameters loop_params(n_loops, cfg.get<Spectrum>("spectrum"),
                                   cfg.get<Dynamics>("dynamics"),
                                   cfg.get<bool>("rsd"));
        SumTable sum_table(loop_params);

        IRresumSettings ir_settings(cfg.get<double>("k_s"), cfg.get<double>("k_osc"));

        InputPowerSpectrum ps(cfg.get<string>("input_ps_file"),
                              cfg.get<double>("input_ps_rescale_num"),
                              cfg.get<bool>("ir_resum"), ir_settings, n_loops,
                              cfg.get<int>("pt_order"), cfg.get<bool>("rsd"),
                              cfg.get<double>("rsd_growth_f"));

        IntegrationInput input(ps, cfg.get<double>("q_min"),
                               cfg.get<double>("q_max"),
                               cfg.get<bool>("single_hard_limit"));

        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        input.pair_correlations = cfg.pair_correlations();
        input.ps_diagrams = ps::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.get<double>("k_a"), 0, 0,
                                        cfg.get<double>("rsd_growth_f"),
                                        loop_params, sum_table,
                                        ev_params, eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with ps::integrand() or bs::integrand() */
            int core = -1;
            integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



static void BM_PS_EdS_integrand_2loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("benchmark/ini/ps_eds_spt.cfg", -1, -1, -1, 0, 0);

        int n_loops = 2;
        int n_dims = 3 * n_loops - 1;

        LoopParameters loop_params(n_loops, cfg.get<Spectrum>("spectrum"),
                                   cfg.get<Dynamics>("dynamics"),
                                   cfg.get<bool>("rsd"));
        SumTable sum_table(loop_params);

        IRresumSettings ir_settings(cfg.get<double>("k_s"), cfg.get<double>("k_osc"));

        InputPowerSpectrum ps(cfg.get<string>("input_ps_file"),
                              cfg.get<double>("input_ps_rescale_num"),
                              cfg.get<bool>("ir_resum"), ir_settings, n_loops,
                              cfg.get<int>("pt_order"), cfg.get<bool>("rsd"),
                              cfg.get<double>("rsd_growth_f"));

        IntegrationInput input(ps, cfg.get<double>("q_min"),
                               cfg.get<double>("q_max"),
                               cfg.get<bool>("single_hard_limit"));

        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        input.pair_correlations = cfg.pair_correlations();
        input.ps_diagrams = ps::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.get<double>("k_a"), 0, 0,
                                        cfg.get<double>("rsd_growth_f"),
                                        loop_params, sum_table,
                                        ev_params, eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with ps::integrand() or bs::integrand() */
            int core = -1;
            integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
            return;
    }
}


static void BM_BS_EdS_integrand_1loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("benchmark/ini/bs_eds_spt.cfg", -1, -1, -1, 0, 0);

        int n_loops = 1;
        int n_dims = 3 * n_loops;

        LoopParameters loop_params(n_loops, cfg.get<Spectrum>("spectrum"),
                                   cfg.get<Dynamics>("dynamics"),
                                   cfg.get<bool>("rsd"));
        SumTable sum_table(loop_params);

        IRresumSettings ir_settings(cfg.get<double>("k_s"), cfg.get<double>("k_osc"));

        InputPowerSpectrum ps(cfg.get<string>("input_ps_file"),
                              cfg.get<double>("input_ps_rescale_num"),
                              cfg.get<bool>("ir_resum"), ir_settings, n_loops,
                              cfg.get<int>("pt_order"), cfg.get<bool>("rsd"),
                              cfg.get<double>("rsd_growth_f"));

        IntegrationInput input(ps, cfg.get<double>("q_min"),
                               cfg.get<double>("q_max"),
                               cfg.get<bool>("single_hard_limit"));

        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        input.pair_correlations = cfg.pair_correlations();
        input.bs_diagrams = bs::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.get<double>("k_a"),
                                      cfg.get<double>("k_b"),
                                      cfg.get<double>("cos_ab"),
                                      cfg.get<double>("rsd_growth_f"),
                                      loop_params, sum_table, ev_params,
                                      eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with bs::integrand() or bs::integrand() */
            int core = -1;
            integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



static void BM_BS_EdS_integrand_2loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("benchmark/ini/bs_eds_spt.cfg", -1, -1, -1, 0, 0);

        int n_loops = 2;
        int n_dims = 3 * n_loops;

        LoopParameters loop_params(n_loops, cfg.get<Spectrum>("spectrum"),
                                   cfg.get<Dynamics>("dynamics"),
                                   cfg.get<bool>("rsd"));
        SumTable sum_table(loop_params);

        IRresumSettings ir_settings(cfg.get<double>("k_s"), cfg.get<double>("k_osc"));

        InputPowerSpectrum ps(cfg.get<string>("input_ps_file"),
                              cfg.get<double>("input_ps_rescale_num"),
                              cfg.get<bool>("ir_resum"), ir_settings, n_loops,
                              cfg.get<int>("pt_order"), cfg.get<bool>("rsd"),
                              cfg.get<double>("rsd_growth_f"));

        IntegrationInput input(ps, cfg.get<double>("q_min"),
                               cfg.get<double>("q_max"),
                               cfg.get<bool>("single_hard_limit"));

        EvolutionParameters ev_params;
        EtaGrid eta_grid;

        input.pair_correlations = cfg.pair_correlations();
        input.bs_diagrams = bs::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.get<double>("k_a"),
                                      cfg.get<double>("k_b"),
                                      cfg.get<double>("cos_ab"),
                                      cfg.get<double>("rsd_growth_f"),
                                      loop_params, sum_table, ev_params,
                                      eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with bs::integrand() or bs::integrand() */
            int core = -1;
            integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



static void BM_PS_2fluid_integrand_1loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("benchmark/ini/quijote_Mnu_0.1eV.cfg",
                   -1, -1, -1, 0,
                   0);

        int n_loops = 1;
        int n_dims = 3 * n_loops - 1;

        LoopParameters loop_params(n_loops, cfg.get<Spectrum>("spectrum"),
                                   cfg.get<Dynamics>("dynamics"),
                                   cfg.get<bool>("rsd"));
        SumTable sum_table(loop_params);

        IRresumSettings ir_settings(cfg.get<double>("k_s"), cfg.get<double>("k_osc"));

        InputPowerSpectrum ps(cfg.get<string>("input_ps_file"),
                              cfg.get<double>("input_ps_rescale_num"),
                              cfg.get<bool>("ir_resum"), ir_settings, n_loops,
                              cfg.get<int>("pt_order"), cfg.get<bool>("rsd"),
                              cfg.get<double>("rsd_growth_f"));

        IntegrationInput input(ps, cfg.get<double>("q_min"),
                               cfg.get<double>("q_max"),
                               cfg.get<bool>("single_hard_limit"));

        EvolutionParameters ev_params(cfg.get<double>("f_nu"),
                                            cfg.get<double>("omega_m_0"),
                                            cfg.get<string>("zeta_file"),
                                            cfg.get<string>("redshift_file"),
                                            cfg.get<string>("omega_eigenvalues_file"),
                                            cfg.F1_ic_files(),
                                            cfg.get<string>("effcs2_x_grid"),
                                            cfg.get<string>("effcs2_y_grid"),
                                            cfg.get<string>("effcs2_data"),
                                            cfg.get<double>("ode_atol"),
                                            cfg.get<double>("ode_rtol"),
                                            cfg.get<double>("ode_hstart"));
        EtaGrid eta_grid(cfg.get<size_t>("pre_time_steps"),
                               cfg.get<size_t>("time_steps"),
                               cfg.get<double>("eta_ini"),
                               cfg.get<double>("eta_fin"),
                               cfg.get<double>("eta_asymp"));

        input.pair_correlations = cfg.pair_correlations();
        input.ps_diagrams = ps::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.get<double>("k_a"), 0, 0,
                                        cfg.get<double>("rsd_growth_f"),
                                        loop_params, sum_table,
                                        ev_params, eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with ps::integrand() or bs::integrand() */
            int core = -1;
            integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



static void BM_PS_2fluid_integrand_2loop(benchmark::State& state) {
    // Perform setup here
    try {
        Config cfg("benchmark/ini/quijote_Mnu_0.1eV.cfg",
                   -1, -1, -1, 0,
                   0);

        int n_loops = 2;
        int n_dims = 3 * n_loops - 1;

        LoopParameters loop_params(n_loops, cfg.get<Spectrum>("spectrum"),
                                   cfg.get<Dynamics>("dynamics"),
                                   cfg.get<bool>("rsd"));
        SumTable sum_table(loop_params);

        IRresumSettings ir_settings(cfg.get<double>("k_s"), cfg.get<double>("k_osc"));

        InputPowerSpectrum ps(cfg.get<string>("input_ps_file"),
                              cfg.get<double>("input_ps_rescale_num"),
                              cfg.get<bool>("ir_resum"), ir_settings, n_loops,
                              cfg.get<int>("pt_order"), cfg.get<bool>("rsd"),
                              cfg.get<double>("rsd_growth_f"));

        IntegrationInput input(ps, cfg.get<double>("q_min"),
                               cfg.get<double>("q_max"),
                               cfg.get<bool>("single_hard_limit"));

        EvolutionParameters ev_params(cfg.get<double>("f_nu"),
                                            cfg.get<double>("omega_m_0"),
                                            cfg.get<string>("zeta_file"),
                                            cfg.get<string>("redshift_file"),
                                            cfg.get<string>("omega_eigenvalues_file"),
                                            cfg.F1_ic_files(),
                                            cfg.get<string>("effcs2_x_grid"),
                                            cfg.get<string>("effcs2_y_grid"),
                                            cfg.get<string>("effcs2_data"),
                                            cfg.get<double>("ode_atol"),
                                            cfg.get<double>("ode_rtol"),
                                            cfg.get<double>("ode_hstart"));
        EtaGrid eta_grid(cfg.get<size_t>("pre_time_steps"),
                               cfg.get<size_t>("time_steps"),
                               cfg.get<double>("eta_ini"),
                               cfg.get<double>("eta_fin"),
                               cfg.get<double>("eta_asymp"));

        input.pair_correlations = cfg.pair_correlations();
        input.ps_diagrams = ps::construct_diagrams(loop_params);

        input.tables_vec.emplace_back(cfg.get<double>("k_a"), 0, 0,
                                        cfg.get<double>("rsd_growth_f"),
                                        loop_params, sum_table,
                                        ev_params, eta_grid);

        int n_correlations = static_cast<int>(input.pair_correlations.size());

        double* xx = new double[static_cast<size_t>(n_dims)];
        double* ff = new double[static_cast<size_t>(n_correlations)];

        for (int i = 0; i < n_dims; ++i) {
            xx[i] = 0.5;
        }

        for (auto _ : state) {
            // This code gets timed
            int nvec = 1;
            /* core = -1 in accordance with ps::integrand() or bs::integrand() */
            int core = -1;
            integrand(&n_dims, xx, &n_correlations, ff, &input, &nvec, &core);
        }

        delete[] xx;
        delete[] ff;
    }
    catch (const ConfigException& cfgex) {
        std::cerr << cfgex.what() << std::endl;
            return;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return;
    }
}



// Register the functions as benchmarks
BENCHMARK(BM_PS_EdS_integrand_1loop);
BENCHMARK(BM_PS_EdS_integrand_2loop);
BENCHMARK(BM_BS_EdS_integrand_1loop);
BENCHMARK(BM_BS_EdS_integrand_2loop);
BENCHMARK(BM_PS_2fluid_integrand_1loop);
BENCHMARK(BM_PS_2fluid_integrand_2loop);
// Run the benchmark
BENCHMARK_MAIN();
