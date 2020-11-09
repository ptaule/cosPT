/*
   parameters.cpp

   Created by Petter Taule on 04.10.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <utility>

#include <libconfig.h++>

#include "../include/utilities.hpp"
#include "../include/version.hpp"
#include "../include/wavenumbers.hpp"
#include "../include/parameters.hpp"

using std::size_t;


Config::Config(const std::string& ini_file, int k_a_idx, int k_b_idx)
    : k_a_idx(k_a_idx), k_b_idx(k_b_idx)
{
    libconfig::Config cfg;

    /* Read config file */
    try {
        cfg.readFile(ini_file.c_str());
    }
    catch (const libconfig::FileIOException& ioex) {
        throw ConfigException("Unable to read " + ini_file + ".");
    }
    catch (const libconfig::ParseException& pex) {
        throw ConfigException("Parse error at " + std::string(pex.getFile()) +
                              ":" + std::to_string(pex.getLine()) + " - " +
                              std::string(pex.getError()));
    }

    /* Number of loops */
    try {
        n_loops_ = cfg.lookup("loops");
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        throw ConfigException("No loops setting in configuration file.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for loops setting.");
    }
    /* Dynamics */
    try {
        std::string dynamics_str = cfg.lookup("dynamics");
        std::transform(dynamics_str.begin(), dynamics_str.end(),
                dynamics_str.begin(), [](unsigned char c) {return
                tolower(c);});

        if (dynamics_str == "eds-spt") {
            dynamics_ = EDS_SPT;
        }
        else if (dynamics_str == "evolve-asymp-ic") {
            dynamics_ = EVOLVE_ASYMP_IC;
        }
        else if (dynamics_str == "eolve-eds-ic") {
            throw ConfigException("eolve-eds-ic not implemented.");
        }
        else {
            throw ConfigException("Unknown dynamics in configuration file.");
        }
    }
    catch (const libconfig::SettingNotFoundException& nfex){
        throw ConfigException("No dynamics setting in configuration file.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for dynamics setting.");
    }
    catch (const ConfigException& e) {
        throw e;
    }

    /* Spectrum */
    try {
        std::string spetrcum_str = cfg.lookup("spectrum");
        std::transform(spetrcum_str.begin(), spetrcum_str.end(),
                spetrcum_str.begin(), [](unsigned char c) {return
                tolower(c);});

        if (spetrcum_str == "powerspectrum") {
            spectrum_ = POWERSPECTRUM;
        }
        else if (spetrcum_str == "bispectrum") {
            spectrum_ = BISPECTRUM;
        }
        else {
            throw ConfigException("Unknown spectrum in configuration file.");
        }
    }
    catch (const libconfig::SettingNotFoundException& nfex){
        throw ConfigException("No spectrum setting in configuration file.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for spectrum setting.");
    }
    catch (const ConfigException& e) {
        throw e;
    }

    /* First wavenumber */
    if (k_a_idx != -1) {
        k_a_ = Wavenumbers::grid[k_a_idx];
    }
    else if (cfg.lookupValue("k_a_idx", k_a_idx)) {
        k_a_ = Wavenumbers::grid[k_a_idx];
    }
    else if (cfg.lookupValue("k_a", k_a_)) {}
    else {
        throw ConfigException("Neither k_a nor k_a_idx in configuration file, "
                              "and no k_a index given as command line argument.");
    }

    /* For power spectrum, two-point correlations */
    if (spectrum_ == POWERSPECTRUM) {
        try {
            const libconfig::Setting& correlation = cfg.lookup("correlations");
            int count = correlation.getLength();
            if (count == 0) {
                throw ConfigException("No correlations in configuration file.");
            }
            for (int i = 0; i < count; ++i) {
                if (correlation[i].getLength() != 2) {
                    throw ConfigException(
                        "Correlation must be two indices (powerspectrum)");
                }
                int a = correlation[i][0];
                int b = correlation[i][1];
                if (a < 0 || a >= COMPONENTS || b < 0 || b >= COMPONENTS) {
                    throw ConfigException("a and b must be element in [0," +
                            std::to_string(COMPONENTS) + "]");
                }
                pair_correlations_.push_back({a,b});
            }
        }
        catch (const libconfig::SettingNotFoundException& nfex) {
            throw ConfigException("No correlations in configuration file.");
        }
        catch (const ConfigException& e) {
            throw e;
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception for correlations setting.");
        }

        if (k_b_idx != -1 || cfg.exists("k_b_idx") || cfg.exists("k_b") ||
            cfg.exists("cos_ab")) {
            std::cerr << "For powerspectrum, k_b/cos_ab settings are ignored."
                    << std::endl;
        }
    }
    /* For bispectrum: k_b, cos_ab and three-point correlations */
    else {
        try {
            const libconfig::Setting& correlation = cfg.lookup("correlations");
            int count = correlation.getLength();
            if (count == 0) {
                throw ConfigException("No correlations in configuration file.");
            }
            for (int i = 0; i < count; ++i) {
                if (correlation[i].getLength() != 3) {
                    throw ConfigException(
                        "Correlation must be three indices (bispectrum)");
                }
                int a = correlation[i][0];
                int b = correlation[i][1];
                int c = correlation[i][2];
                if (a < 0 || a >= COMPONENTS || b < 0 || b >= COMPONENTS ||
                        c < 0 || c >= COMPONENTS) {
                    throw ConfigException("a, b and c must be element in [0," +
                            std::to_string(COMPONENTS) + "]");
                }
                triple_correlations_.push_back({a,b,c});
            }
        }
        catch (const ConfigException& e) {
            throw e;
        }
        catch (const libconfig::SettingNotFoundException& nfex) {
            throw ConfigException("No correlations in configuration file.");
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception for correlations setting.");
        }

        if (k_b_idx != -1) {
            k_b_ = Wavenumbers::grid[k_b_idx];
        }
        else if (cfg.lookupValue("k_b_idx", k_b_idx)) {
            k_b_ = Wavenumbers::grid[k_b_idx];
        }
        else if (cfg.lookupValue("k_b", k_b_)) {}
        else {
            throw ConfigException(
                "Neither k_b nor k_b_idx in configuration file, "
                "and no k_b index given as command line argument.");
        }
        try {
            cos_ab_ = cfg.lookup("cos_ab");
        }
        catch (const libconfig::SettingNotFoundException& nfex) {
            throw ConfigException(
                "No value set for cos_ab (spectrum = BISPECTRUM).");
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception for cos_ab setting.");
        }
    }

    /* Integration ranges */
    try {
        q_min_ = static_cast<double>(cfg.lookup("q_min"));
        q_max_ = static_cast<double>(cfg.lookup("q_max"));
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        throw ConfigException("Missing q_min and/or q_max in configuration.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for q_min/q_max setting.");
    }

    /* Output file */
    try {
        output_file_ = static_cast<std::string>(cfg.lookup("output_file"));
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        throw ConfigException("No setting for output file in configuration.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for output_file setting.");
    }
    /* CUBA settings */
    if (cfg.exists("cuba_settings")) {
        const libconfig::Setting& cuba_settings = cfg.lookup("cuba_settings");
        cuba_settings.lookupValue("abs_tolerance", cuba_atol_);
        cuba_settings.lookupValue("rel_tolerance", cuba_rtol_);
        cuba_settings.lookupValue("max_evaluations", cuba_maxevals_);
        cuba_settings.lookupValue("verbosity_level", cuba_verbose_);
        cuba_settings.lookupValue("threads", cuba_cores_);
        cuba_settings.lookupValue("statefile", cuba_statefile_);
        cuba_settings.lookupValue("retain_statefile", cuba_retain_statefile_);
    }

    /* ODE settings */
    if (cfg.exists("ODE settings")) {
        const libconfig::Setting& ode_settings = cfg.lookup("cuba_settings");
        ode_settings.lookupValue("abs_tolerance", ode_atol_);
        ode_settings.lookupValue("rel_tolerance", ode_rtol_);
        ode_settings.lookupValue("start_step", ode_hstart_);
    }

    /* Input power spectrum */
    try {
        input_ps_file_ = static_cast<std::string>(cfg.lookup("input_ps_file"));
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        throw ConfigException("No input PS file in configuration.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for input_ps_file setting.");
    }

    /* Settings for evolution dynamics */
    if (dynamics_ == EVOLVE_ASYMP_IC || dynamics_ == EVOLVE_EDS_IC) {
        try {
            time_steps_ = cfg.lookup("time_steps");
            eta_ini_ = static_cast<double>(cfg.lookup("eta_ini"));
            eta_fin_ = static_cast<double>(cfg.lookup("eta_fin"));
        }
        catch (const libconfig::SettingNotFoundException& nfex) {
            throw ConfigException(
                "Missing time grid settings in configuration file.");
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception in time grid settings.");
        }
        if (dynamics_ == EVOLVE_ASYMP_IC) {
            try {
                pre_time_steps_ = cfg.lookup("pre_time_steps");
                eta_asymp_ = static_cast<double>(cfg.lookup("eta_asymp"));
            }
            catch (const libconfig::SettingNotFoundException& nfex) {
                throw ConfigException(
                    "Missing time grid settings in configuration file.");
            }
            catch (const libconfig::SettingTypeException& tex) {
                throw ConfigException("Encountered type exception in time grid settings.");
            }
        }
        try {
            zeta_file_ = static_cast<std::string>(cfg.lookup("zeta_file"));
            redshift_file_ = static_cast<std::string>(cfg.lookup("redshift_file"));
            omega_eigenvalues_file_ =
                static_cast<std::string>(cfg.lookup("omega_eigenvalues_file"));

            const libconfig::Setting& F1_ic_files_list = cfg.lookup("F1_ic_files");
            int count = F1_ic_files_list.getLength();
            if (count != COMPONENTS) {
                throw ConfigException("There should be " +
                                      std::to_string(COMPONENTS) +
                                      " F1 ic files in the configuration");
            }
            for (int i = 0; i < count; ++i) {
                F1_ic_files_.push_back(
                    static_cast<std::string>(F1_ic_files_list[i]));
            }

            const libconfig::Setting& effcs2_files =
                cfg.lookup("effective_cs2_files");
            effcs2_x_grid_ = static_cast<std::string>(effcs2_files.lookup("x_grid"));
            effcs2_y_grid_ = static_cast<std::string>(effcs2_files.lookup("y_grid"));
            effcs2_data_   = static_cast<std::string>(effcs2_files.lookup("data"));
        }
        catch (const libconfig::SettingNotFoundException& nfex) {
            throw ConfigException(
                "Missing file for interpolation in configuration.");
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception in settings for "
                                  "interpolation files.");
        }
        try {
            f_nu_ = static_cast<double>(cfg.lookup("f_nu"));
            omega_m_0_ = static_cast<double>(cfg.lookup("omega_m_0"));
        }
        catch (const libconfig::SettingNotFoundException& nfex) {
            throw ConfigException("Missing f_nu or omega_m_0 in configuration.");
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception for f_nu or omega_m_0.");
        }
    }

    /* Store potential description */
    cfg.lookupValue("description", description_);
}



std::ostream& operator<<(std::ostream& out, const Config& cfg) {
    if (cfg.spectrum() == POWERSPECTRUM) {
        out << "# Matter power spectrum P(k) at " << cfg.n_loops()
            << "-loop for k = " << cfg.k_a() << " (h/Mpc)" << "\n";
    }
    else {
      out << "# Matter bispectrum B(k) at " << cfg.n_loops() << "-loop for\n"
          << "#\n# k_a = " << cfg.k_a() << " (h/Mpc)"
          << "\n# k_b = " << cfg.k_b() << " (h/Mpc)"
          << "\n# cos_ab = " << cfg.cos_ab() << "\n";
    }
    out << "#\n# Description: " << cfg.description() << "\n";
    out <<    "# Git hash:    " << build_git_sha << "\n";
    out <<    "# Build time:  " << build_git_time << "\n";

    out << "#\n# Correlations (zero-indexed components):\n# ";
    if (cfg.spectrum() == POWERSPECTRUM) {
        for (auto& el : cfg.pair_correlations()) {
            out << " " << el <<  " ,";
        }
    }
    else {
        for (auto& el : cfg.triple_correlations()) {
            out << " " << el <<  " ,";
        }
    }

    out << "\n#\n# Input power spectrum read from " << cfg.input_ps_file() << "\n";
    out << std::scientific;
    out << "# Integration limits:\n";
    out << "#\t q_min = " << cfg.q_min() << "\n";
    out << "#\t q_max = " << cfg.q_max() << "\n";

    out << "#\n# Cuba settings:\n";
    out << "#\t abs tolerance = " << cfg.cuba_atol() << "\n";
    out << "#\t rel tolerance = " << cfg.cuba_rtol() << "\n";
    out << "#\t max. evals    = " << cfg.cuba_maxevals() << "\n";

    if (!cfg.cuba_statefile().empty()) {
        out << "#\t statefile     = " << cfg.cuba_statefile() << "\n";
    }

    if (cfg.dynamics() == EVOLVE_EDS_IC || cfg.dynamics() == EVOLVE_ASYMP_IC) {
        out << "#\n# ODE settings:\n";
        out << std::scientific;
        out << "#\t abs tolerance = " << cfg.ode_atol() << "\n";
        out << "#\t rel tolerance = " << cfg.ode_rtol() << "\n";
        out << "#\t start step    = " << cfg.ode_hstart() << "\n";

        out << "#\n# Time grid settings:\n";
        out << "#\t time steps     = " << cfg.time_steps() << "\n";
        if (cfg.dynamics() == EVOLVE_ASYMP_IC) {
            out << "#\t pre time steps = " << cfg.pre_time_steps() << "\n";
            out << "#\t eta asymp      = " << cfg.eta_asymp() << "\n";
        }
        out << "#\t eta ini        = " << cfg.eta_ini() << "\n";
        out << "#\t eta fin        = " << cfg.eta_fin() << "\n";

        out << "#\n# Dynamics settings:\n";
        out << "# f_nu      = " << cfg.f_nu() << "\n";
        out << "# omega_m_0 = " << cfg.omega_m_0() << "\n";

        out << "#\n# zeta file              = " << cfg.zeta_file() << "\n";
        out << "# redshift file          = " << cfg.redshift_file() << "\n";
        out << "# omega eigenvalues file = " << cfg.omega_eigenvalues_file() << "\n";
        for (int i = 0; i < COMPONENTS; ++i) {
            out << "# F1 ic files[" << i << "]\t\t = " << cfg.F1_ic_files().at(i)
                << "\n";
        }
        out << "# effective cs2 x grid   = " << cfg.effcs2_x_grid() << "\n";
        out << "# effective cs2 y grid   = " << cfg.effcs2_y_grid() << "\n";
        out << "# effective cs2 data     = " << cfg.effcs2_data() << "\n";
    }

    return out;
}



LoopParameters::LoopParameters(int n_loops, Spectrum spectrum, Dynamics dynamics)
    : dynamics_(dynamics), spectrum_(spectrum), n_loops_(n_loops)
{
    if (spectrum_ == POWERSPECTRUM && (n_loops_ < 1 || n_loops_ > 2)) {
        throw(std::invalid_argument(
            "LoopParameters::LoopParameters(): POWERSPECTRUM only "
            "implemented for n_loops = 1,2."));
    }
    if (spectrum_ == BISPECTRUM && (n_loops_ != 1)) {
        throw(std::invalid_argument(
            "LoopParameters::LoopParameters(): BISPECTRUM only "
            "implemented for n_loops = 1."));
    }
    if (spectrum_ == POWERSPECTRUM) {
        n_coeffs_      = n_loops_ + 1;
        n_configs_     = pow(3, n_coeffs_);
        n_kernels_     = (n_configs_/3 + 1) * pow(4,n_loops_);
        n_kernel_args_ = 2 * n_loops_ + 1;
        zero_label_    = get_zero_label(n_coeffs_);

        /* Define block size for kernel indexing. A block consists of all */
        /* single loop label argument combinations. */
        single_loop_block_size = pow(4, n_loops_);
        /* Max/min single_loops labels */
        single_loop_label_min = n_configs_ / 3;
        single_loop_label_max = 2 * n_configs_ / 3 - 1;
    }
    else if (spectrum_ == BISPECTRUM) {
        n_coeffs_      = n_loops_ + 2;
        n_configs_     = pow(3, n_coeffs_);
        n_kernels_     = pow(n_configs_ + 1, 2) * pow(4, n_loops_);
        n_kernel_args_ = 2 * n_loops_ + 2;
        zero_label_    = get_zero_label(n_coeffs_);

        single_loop_block_size = pow(4, n_loops_);
        /* Max/min single_loops labels */
        single_loop_label_min = 4 * n_configs_ / 9;
        single_loop_label_max = 5 * n_configs_ / 9 - 1;
        /* Composite label has either k_a, k_b, or multiple loop momenta
         * (applicable for overall loop diagram) */
        first_composite_block_size = (n_configs_ + 1) * single_loop_block_size;
    }
    else {
        throw(std::invalid_argument(
            "LoopParameters::LoopParameters(): invalid spectrum."));
    }

    /* List of single loop labels */
    Vec1D<int> config(n_coeffs_, 0);
    int label;

    for (int i = 0; i < n_loops_; ++i) {
        config.at(i) = -1;
        label = config2label(config);
        single_loops.push_back(label);

        config.at(i) = 1;
        label = config2label(config);
        single_loops.push_back(label);

        config.at(i) = 0;
    }
}



int LoopParameters::ps_arguments_2_kernel_index(const int arguments[]) const
{
   /* Precompute powers of two for speedup */
    int pow2[] = {1,2,4,8,16,32,64,128};

    // In DEBUG-mode, check that non-zero arguments (zero_label) are unique
#if DEBUG >= 1
    if (!unique_elements(arguments, n_kernel_args_, zero_label_))
        throw(std::logic_error(
            "LoopParameters::ps_arguments_2_kernel_index(): duplicate "
            "vector arguments passed."));
    int n_k_labels = 0;
#endif

    int index = 0;

    for (int i = 0; i < n_kernel_args_; ++i) {
        /* First, check if argument is a zero vector */
        if (arguments[i] == zero_label_) continue;

        /* Argument is a k-type vector (i.e. on the form k + c_i Q_i) if k is
         * present. In our vector-label convention, k is the last coefficient,
         * hence +k is present if label > single_loop_label_max */
        if (arguments[i] > single_loop_label_max) {
            index += (arguments[i] - single_loop_label_max) * single_loop_block_size;
#if DEBUG >= 1
            /* Count k-type labels */
            ++n_k_labels;
#endif
        }
#if DEBUG >= 1
        /* We should not get -k in power spectrum computation */
        else if (arguments[i] < single_loop_label_min) {
            throw(std::logic_error("LoopParameters::ps_arguments_2_kernel_index()"
                                   ": got argument with -k."));
        }
#endif
        else {
            /* Single loop */
            for (size_t j = 0; j < single_loops.size(); ++j) {
                if (arguments[i] == single_loops[j]) {
                    index += pow2[j];
                    break;
                }
            }
#if DEBUG >= 1
            /* Check that this is in fact a single loop vector */
            if(!single_loop_label(arguments[i], n_coeffs_, spectrum_))
                throw(std::logic_error(
                    "LoopParameters::ps_arguments_2_kernel_index(): argument is "
                    "neither 0, composite type, or single loop."));
#endif
        }
    }
#if DEBUG >= 1
    if (n_k_labels > 1)
        throw(std::logic_error("LoopParameters::ps_arguments_2_kernel_index(): "
                               "more than one argument is of composite type."));
#endif

    return index;
}



int LoopParameters::bs_arguments_2_kernel_index(const int arguments[]) const
{
   /* Precompute powers of two for speedup */
    int pow2[] = {1,2,4,8,16,32,64,128};

    /* In DEBUG-mode, check that non-zero_label arguments are unique */
#if DEBUG >= 1
    if (!unique_elements(arguments, n_kernel_args_, zero_label_))
        throw(std::logic_error(
            "LoopParameters::bs_arguments_2_kernel_index(): duplicate "
            "vector arguments passed."));
#endif

    int index = 0;

    /* Counter of composite arguments, i.e. not single loop momenta */
    int n_composite = 0;

    for (int i = 0; i < n_kernel_args_; ++i) {
        // First, check if argument is a zero vector
        if (arguments[i] == zero_label_) continue;

        /* Argument is not single loop if label < single_loop_label_min or
         * label > single_loop_label_max */
        if (arguments[i] < single_loop_label_min) {
          index += (arguments[i] + 1) * (n_composite > 0
                                             ? single_loop_block_size
                                             : first_composite_block_size);
          ++n_composite;
        }
        else if (arguments[i] > single_loop_label_max) {
            int block = n_composite > 0 ? single_loop_block_size
                                        : first_composite_block_size;
            /* #  */
            index += (arguments[i] + 1) * block;

            ++n_composite;
        }
        else {
            /* Single loop */
            for (size_t j = 0; j < single_loops.size(); ++j) {
                if (arguments[i] == single_loops[j]) {
                    index += pow2[j];
                    goto found_single_loop;
                }
            }
            /* Went through for loop, which means that argument is composite,
             * hence this corresponds to a connecting line with overall loop */
            index += (arguments[i] + 1) * (n_composite > 0
                                               ? single_loop_block_size
                                               : first_composite_block_size);
            ++n_composite;
found_single_loop: ;
        }
    }
#if DEBUG >= 1
    if (n_composite > 2)
        throw(std::logic_error(
            "LoopParameters::bs_arguments_2_kernel_index(): more than two "
            "arguments is of composite type."));
#endif

    return index;
}



EvolutionParameters::EvolutionParameters(EvolutionParameters&& other) noexcept
    : f_nu_(other.f_nu_), cs2_factor(other.cs2_factor),
    cg2_factor(other.cg2_factor), ode_atol_(other.ode_atol_),
    ode_rtol_(other.ode_rtol_), ode_hstart_(other.ode_hstart_),
    zeta(std::move(other.zeta)), redshift(std::move(other.redshift)),
    omega_eigenvalues(std::move(other.omega_eigenvalues)),
    effcs2(std::move(other.effcs2))
{
    for (auto&& el : other.F1_ic) {
        F1_ic.emplace_back(std::move(el));
    }
}



EvolutionParameters& EvolutionParameters::operator=(EvolutionParameters&& other)
{
    if (this != &other) {
        f_nu_       = other.f_nu_;
        cs2_factor  = other.cs2_factor;
        cg2_factor  = other.cg2_factor;
        ode_atol_   = other.ode_atol_;
        ode_rtol_   = other.ode_rtol_;
        ode_hstart_ = other.ode_hstart_;
        sound_speed = other.sound_speed;

        zeta              = std::move(other.zeta);
        redshift          = std::move(other.redshift);
        omega_eigenvalues = std::move(other.omega_eigenvalues);
        effcs2            = std::move(other.effcs2);

        for (auto&& el : other.F1_ic) {
            F1_ic.emplace_back(std::move(el));
        }
    }
    return *this;
}



EvolutionParameters::EvolutionParameters(
        double f_nu,
        double omega_m_0,
        double ode_atol,
        double ode_rtol,
        double ode_hstart,
        const std::string& zeta_file,
        const std::string& redshift_file,
        const std::string& omega_eigenvalues_file,
        const Vec1D<std::string>& F1_ic_files,
        const std::string& effcs2_x_file,
        const std::string& effcs2_y_file,
        const std::string& effcs2_data_file
        ) :
    f_nu_(f_nu), ode_atol_(ode_atol), ode_rtol_(ode_rtol),
    ode_hstart_(ode_hstart), zeta(zeta_file), redshift(redshift_file),
    omega_eigenvalues(omega_eigenvalues_file), effcs2(effcs2_x_file,
            effcs2_y_file, effcs2_data_file)
{
    for (auto& F1_ic_file : F1_ic_files) {
        F1_ic.emplace_back(F1_ic_file);
    }

    sound_speed = EXACT;

    /* cs2_factor = 2/3 * 1/(omegaM a^2 H^2) */
    cs2_factor = 2.0/3.0 * SQUARE(3e3) / omega_m_0;
}



EvolutionParameters::EvolutionParameters(
        double f_nu,
        double omega_m_0,
        const std::string& zeta_file,
        const std::string& redshift_file,
        const std::string& omega_eigenvalues_file,
        const Vec1D<std::string>& F1_ic_files,
        const std::string& effcs2_x_file,
        const std::string& effcs2_y_file,
        const std::string& effcs2_data_file
        ) :
    EvolutionParameters(f_nu, omega_m_0, 1e-6, 1e-4, 1e-3, zeta_file,
            redshift_file, omega_eigenvalues_file, F1_ic_files, effcs2_x_file,
            effcs2_y_file, effcs2_data_file)
{}



EvolutionParameters::EvolutionParameters(
        double m_nu,
        double f_nu,
        double omega_m_0,
        double ode_atol,
        double ode_rtol,
        double ode_hstart,
        const std::string& zeta_file,
        const std::string& redshift_file,
        const std::string& omega_eigenvalues_file,
        const Vec1D<std::string>& F1_ic_files
        ) :
    f_nu_(f_nu), ode_atol_(ode_atol), ode_rtol_(ode_rtol),
    ode_hstart_(ode_hstart), zeta(zeta_file), redshift(redshift_file),
    omega_eigenvalues(omega_eigenvalues_file)
{
    for (auto& F1_ic_file : F1_ic_files) {
        F1_ic.emplace_back(F1_ic_file);
    }

    sound_speed = ADIABATIC;

    /* cs2_factor = 2/3 * 1/(omegaM a^2 H^2) */
    cg2_factor = 2.0/3.0 * SQUARE(3e3) / omega_m_0;

    /* Neutrino temperature today [eV] */
    double T_nu0 = 1.67734976e-4;
    /* 5/9 Zeta(5)/Zeta(3) = 7.188565369 */
    double a = 7.188565369;
    cg2_factor *= a * SQUARE(T_nu0) / SQUARE(m_nu);
}



EvolutionParameters::EvolutionParameters(
        double m_nu,
        double f_nu,
        double omega_m_0,
        const std::string& zeta_file,
        const std::string& redshift_file,
        const std::string& omega_eigenvalues_file,
        const Vec1D<std::string>& F1_ic_files
        ) :
    EvolutionParameters(m_nu, f_nu, omega_m_0, 1e-6, 1e-4, 1e-3, zeta_file,
            redshift_file, omega_eigenvalues_file, F1_ic_files)
{}



EvolutionParameters::EvolutionParameters(
        double ode_atol,
        double ode_rtol,
        double ode_hstart,
        const std::string& zeta_file
        ) :
    ode_atol_(ode_atol), ode_rtol_(ode_rtol), ode_hstart_(ode_hstart),
    zeta(zeta_file)
{}



EvolutionParameters::EvolutionParameters(const std::string& zeta_file)
    : EvolutionParameters(1e-6, 1e-4, 1e-3, zeta_file)
{}
