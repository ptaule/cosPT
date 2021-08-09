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
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <filesystem>

#include <libconfig.h++>

#include "../include/utilities.hpp"
#include "../include/io.hpp"
#include "../include/version.hpp"
#include "../include/parameters.hpp"

using std::size_t;
using std::pow;
namespace fs = std::filesystem;


Config::Config(const std::string& ini_file,
        int k_a_idx,
        int k_b_idx,
        int k_c_idx,
        int cuba_maxevals,
        int cuba_cores
        )
    : k_a_idx(k_a_idx), k_b_idx(k_b_idx), k_c_idx(k_c_idx),
    cuba_maxevals_(cuba_maxevals), cuba_cores_(cuba_cores)
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

    /* First wavenumber */
    std::string k_a_grid_file;
    Vec2D<double> k_a_grid;
    if (cfg.lookupValue("k_a_grid", k_a_grid_file)) {
        read_columns_from_file(k_a_grid_file, 1, k_a_grid);
    }
    if (!k_a_grid.empty() && k_a_idx != -1) {
        try {
            k_a_ = k_a_grid.at(0).at(static_cast<size_t>(k_a_idx));
        }
        catch (const std::out_of_range& e) {
            throw ConfigException("k_a_idx out of range (of k_a_grid).");
        }
    }
    else if (!k_a_grid.empty() && cfg.lookupValue("k_a_idx", k_a_idx)) {
        try {
            k_a_ = k_a_grid.at(0).at(static_cast<size_t>(k_a_idx));
        }
        catch (const std::out_of_range& e) {
            throw ConfigException("k_a_idx out of range (of k_a_grid).");
        }
        k_a_ = k_a_grid.at(0).at(static_cast<size_t>(k_a_idx));
    }
    else if (cfg.lookupValue("k_a", k_a_)) {}
    else {
        throw ConfigException("Did not obtain any value for k_a. Either provide "
                              "k_a, or k_a_idx and k_a_grid.");
    }

    /* Integration ranges */
    try {
        int q;
        if (cfg.lookupValue("q_min", q_min_)) {}
        /* If not found as double, try int */
        else if (cfg.lookupValue("q_min", q)) {
            q_min_ = static_cast<double>(q);
        }
        else {
            throw ConfigException("Missing q_min in configuration.");
        }
        if (cfg.lookupValue("q_max", q_max_)) {}
        /* If not found as double, try int */
        else if (cfg.lookupValue("q_max", q)) {
            q_max_ = static_cast<double>(q);
        }
        else {
            throw ConfigException("Missing q_max in configuration.");
        }
    }
    catch (const ConfigException& ex) {
        throw ex;
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for q_min/q_max setting.");
    }

    set_spectrum(cfg);
    set_dynamics(cfg);

    /* CUBA settings */
    try {
        if (cfg.exists("cuba_settings")) {
            const libconfig::Setting& cuba_settings = cfg.lookup("cuba_settings");
            cuba_settings.lookupValue("abs_tolerance", cuba_atol_);
            cuba_settings.lookupValue("rel_tolerance", cuba_rtol_);
            cuba_settings.lookupValue("verbosity_level", cuba_verbose_);
            cuba_settings.lookupValue("retain_statefile", cuba_retain_statefile_);

            /* If cuba_maxevals is not already set, look up value */
            if (cuba_maxevals_ == 0) {
                if (!cuba_settings.exists("max_evaluations")) {
                    /* Default value */
                    std::cerr << "No cuba max. evaluations given. Using "
                                 "default value: 1e6." << std::endl;
                    cuba_maxevals_ = 1e6;
                }
                else {
                    /* Casting first to double, so that libconfig recognizes
                     * input such as 1e6, then cast to int */
                    cuba_maxevals_ = static_cast<int>(static_cast<double>(
                        cuba_settings.lookup("max_evaluations")));
                }
            }

            /* If cuba_cores_ is not already set, look up value */
            if (cuba_cores_ == -1) {
                if (!cuba_settings.lookupValue("n_cores", cuba_cores_)) {
                    std::cerr << "No n_cores value given. Using default value: 0"
                              << std::endl;
                    cuba_cores_ = 0;
                }
            }
            set_cuba_statefile(cuba_settings);
        }
        else {
            if (cuba_cores_ == -1) {
                std::cerr << "No n_cores value given. Using default value: 0"
                    << std::endl;
                cuba_cores_ = 0;
            }
        }
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for cuba_settings.");
    }

    /* Input power spectrum */
    try {
        input_ps_file_ = static_cast<std::string>(cfg.lookup("input_ps_file"));
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        throw ConfigException("No input_ps_file setting found in configuration.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException(
            "Encountered type exception for input_ps_file setting.");
    }
    try {
        /* input_ps_rescale setting. Rescaling factor set to 1 by default. */
        /* First, check if given as string "2pi^-3" or "2pi^3" */
        int input_ps_rescale_int;
        if (cfg.lookupValue("input_ps_rescale", input_ps_rescale_str_)) {
            if (input_ps_rescale_str_.compare("2pi^-3") == 0) {
                input_ps_rescale_num = pow(TWOPI,-3);
            }
            else if (input_ps_rescale_str_.compare("2pi^3") == 0) {
                input_ps_rescale_num = pow(TWOPI,3);
            }
            else {
                throw ConfigException("Got input_ps_rescale string \"" + input_ps_rescale_str_ +
                        "\" which is neither \"2pi^-3\" nor \"2pi^3\".");
            }
        }
        /* If not found as string, try double */
        else if (cfg.lookupValue("input_ps_rescale", input_ps_rescale_num)) {}
        /* Last, try integer */
        else if (cfg.lookupValue("input_ps_rescale", input_ps_rescale_int)) {
            input_ps_rescale_num = static_cast<double>(input_ps_rescale_int);
        }
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException(
            "Encountered type exception for input_ps_rescale setting.");
    }
    catch (const ConfigException& e) {
        throw e;
    }

    /* Store potential description */
    cfg.lookupValue("description", description_);

    set_output_file(cfg);
}



void Config::set_spectrum(const libconfig::Config& cfg)
{
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
            k_c_idx != -1 || cfg.exists("k_c_idx") || cfg.exists("k_c") ||
            cfg.exists("cos_ab"))
        {
            std::cerr << "For powerspectrum, k_b/k_c/cos_ab settings are ignored."
                    << std::endl;
        }
    }
    /* For bispectrum: k_b, k_c/cos_ab and three-point correlations */
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
            throw ConfigException(
                "Encountered type exception for correlations setting.");
        }

        std::string k_b_grid_file;
        Vec2D<double> k_b_grid;
        if (cfg.lookupValue("k_b_grid", k_b_grid_file)) {
            read_columns_from_file(k_b_grid_file, 1, k_b_grid);
        }
        if (!k_b_grid.empty() && k_b_idx != -1) {
            try {
                k_b_ = k_b_grid.at(0).at(static_cast<size_t>(k_b_idx));
            }
            catch (const std::out_of_range& e) {
                throw ConfigException("k_b_idx out of range (of k_b_grid).");
            }
        }
        else if (!k_b_grid.empty() && cfg.lookupValue("k_b_idx", k_b_idx)) {
            try {
                k_b_ = k_b_grid.at(0).at(static_cast<size_t>(k_b_idx));
            }
            catch (const std::out_of_range& e) {
                throw ConfigException("k_b_idx out of range (of k_b_grid).");
            }
            k_b_ = k_b_grid.at(0).at(static_cast<size_t>(k_b_idx));
        }
        else if (!k_b_grid.empty() && cfg.lookupValue("k_b_idx", k_b_idx)) {
            k_b_ = k_b_grid.at(0).at(static_cast<size_t>(k_b_idx));
        }
        else if (cfg.lookupValue("k_b", k_b_)) {}
        else {
            throw ConfigException(
                "Did not obtain any value for k_b. Either provide "
                "k_b, or k_b_idx and k_b_grid.");
        }

        bool k_c_given = false;
        try {
            std::string k_c_grid_file;
            Vec2D<double> k_c_grid;
            if (cfg.lookupValue("k_c_grid", k_c_grid_file)) {
                read_columns_from_file(k_c_grid_file, 1, k_c_grid);
            }
            if (!k_c_grid.empty() && k_c_idx != -1) {
                try {
                    k_c_ = k_c_grid.at(0).at(static_cast<size_t>(k_c_idx));
                }
                catch (const std::out_of_range& e) {
                    throw ConfigException("k_c_idx out of range (of k_c_grid).");
                }
                k_c_given = true;
            }
            else if (!k_c_grid.empty() && cfg.lookupValue("k_c_idx", k_c_idx)) {
                try {
                    k_c_ = k_c_grid.at(0).at(static_cast<size_t>(k_c_idx));
                }
                catch (const std::out_of_range& e) {
                    throw ConfigException("k_c_idx out of range (of k_c_grid).");
                }
                k_c_ = k_c_grid.at(0).at(static_cast<size_t>(k_c_idx));
                k_c_given = true;
            }
            else if (!k_c_grid.empty() && cfg.lookupValue("k_c_idx", k_c_idx)) {
                k_c_ = k_c_grid.at(0).at(static_cast<size_t>(k_c_idx));
                k_c_given = true;
            }
            else if (cfg.lookupValue("k_c", k_c_)) {
                k_c_given = true;
            }
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception for k_c setting.");
        }

        bool cos_ab_given = false;
        try {
            if (cfg.lookupValue("cos_ab", cos_ab_)) {
                cos_ab_given = true;

                if (cos_ab_ < -1 || cos_ab_ > 1) {
                    throw ConfigException(
                        "Got cos_ab = " + std::to_string(cos_ab_) +
                        ", which is not between -1 and 1.");
                }
            }
        }
        catch (const ConfigException& e) {
            throw e;
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException(
                "Encountered type exception for cos_ab setting.");
        }

        /* If neither/both k_c and cos_ab was given, throw exception */
        if ((!k_c_given) && (!cos_ab_given)) {
            throw ConfigException(
                "Neither k_c nor cos_ab given. Please provide one of them.");
        }
        if (k_c_given && cos_ab_given) {
            throw ConfigException(
                "Both k_c and cos_ab given. Please choose one.");
        }

        /* k_c = - k_a - k_b */
        if (k_c_given) {
            cos_ab_ = 0.5 * (SQUARE(k_c_)/(k_a_*k_b_) - k_a_/k_b_ - k_b_/k_a_);

            if (cos_ab_ < -1 || cos_ab_ > 1) {
                throw ConfigException(
                    "The k_a, k_b, k_c values given does not constitute a valid "
                    "configuration (cos_ab out of range).");
            }
        }
        else {
            /* cos_ab given */
            k_c_ = std::sqrt(SQUARE(k_a_) + SQUARE(k_b_) + 2*k_a_*k_b_*cos_ab_);
        }
    }
}



void Config::set_dynamics(const libconfig::Config& cfg)
{
    try {
        std::string dynamics_str = cfg.lookup("dynamics");
        std::transform(dynamics_str.begin(), dynamics_str.end(),
                dynamics_str.begin(), [](unsigned char c) {return
                tolower(c);});

        if (dynamics_str == "eds-spt") {
            dynamics_ = EDS_SPT;
        }
        else if (dynamics_str == "evolve-ic-asymp") {
            dynamics_ = EVOLVE_IC_ASYMP;
        }
        else if (dynamics_str == "evolve-ic-eds") {
            dynamics_ = EVOLVE_IC_EDS;
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

    /* Settings for evolution dynamics */
    if (dynamics_ == EVOLVE_IC_ASYMP || dynamics_ == EVOLVE_IC_EDS) {

        /* ODE settings */
        try {
            if (cfg.exists("ode_settings")) {
                const libconfig::Setting& ode_settings = cfg.lookup("ode_settings");
                ode_settings.lookupValue("abs_tolerance", ode_atol_);
                ode_settings.lookupValue("rel_tolerance", ode_rtol_);
                ode_settings.lookupValue("start_step", ode_hstart_);
            }
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception in ODE settings.");
        }

        try {
            /* Need to read setting as int, then convert so size_t */
            int ts = cfg.lookup("time_steps");
            time_steps_ = static_cast<size_t>(ts);
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
        if (dynamics_ == EVOLVE_IC_ASYMP) {
            try {
                /* Need to read setting as int, then convert so size_t */
                int pts = cfg.lookup("pre_time_steps");
                pre_time_steps_ = static_cast<size_t>(pts);
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

            if (dynamics_ == EVOLVE_IC_ASYMP) {
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
        }
        catch (const libconfig::SettingNotFoundException& nfex) {
            throw ConfigException(
                "Missing file for interpolation in configuration.");
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception in settings for "
                                  "interpolation files.");
        }
        if (dynamics_ == EVOLVE_IC_ASYMP) {
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
    }
}



void Config::set_output_file(const libconfig::Config& cfg)
{
    try {
        if (cfg.lookupValue("output_file", output_file_)) {
            /* Check that directory exists */
            fs::path p(output_file_);
            if (p.has_parent_path() && !fs::exists(p.parent_path())) {
                throw ConfigException("Output file directory " +
                                      std::string(p.parent_path()) +
                                      " does not exist.");
            }
        }
        else if (cfg.lookupValue("output_path", output_path)) {
            /* Check that directory exists */
            if (!fs::exists(fs::path(output_path))) {
                throw ConfigException("Output directory " + output_path +
                                      " does not exist.");
            }

            std::stringstream ss;
            ss << output_path;
            ss << "/";
            /* Add k_a_idx (& k_b_idx) to end of file */
            if (k_a_idx != -1) {
                ss << std::setfill('0') << std::setw(3) << k_a_idx;
            }
            else {
                ss << std::scientific << std::setprecision(6) << k_a_;
            }
            if (spectrum_ == BISPECTRUM) {
                ss << "_";
                if (k_b_idx != -1) {
                    ss << std::setfill('0') << std::setw(3) << k_b_idx;
                }
                else {
                    ss << std::scientific << std::setprecision(6) << k_b_;
                }
                ss << "_";
                if (k_c_idx != -1) {
                    ss << std::setfill('0') << std::setw(3) << k_c_idx;
                }
                else {
                    ss << std::scientific << std::setprecision(6) << k_c_;
                }
            }
            ss << ".dat";
            output_file_ = ss.str();
        }
        else {
            throw ConfigException("No output path/file given in configuration.");
        }
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        throw ConfigException("No setting for output file or path in configuration.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for output_file setting.");
    }
    catch (const ConfigException& ex) {
        throw ex;
    }
}



void Config::set_cuba_statefile(const libconfig::Setting& cuba_settings)
{
    try {
        if (cuba_settings.lookupValue("statefile", cuba_statefile_)) {
            /* Check that directory exists */
            fs::path p(cuba_statefile_);
            if (!fs::exists(p.parent_path())) {
                throw ConfigException("CUBA statefile file directory " +
                                      std::string(p.parent_path()) +
                                      " does not exist.");
            }
        }
        else if (cuba_settings.lookupValue("statefile_path",
                                           cuba_statefile_path)) {
            /* Check that directory exists */
            if (!fs::exists(fs::path(cuba_statefile_path))) {
                throw ConfigException("CUBA statefile directory " +
                                      cuba_statefile_path + " does not exist.");
            }
            cuba_statefile_ = cuba_statefile_path;
            cuba_statefile_ += "/";
            /* Add k_a_idx (& k_b_idx) to end of file */
            if (k_a_idx != -1) {
                cuba_statefile_ += std::to_string(k_a_idx);
            }
            else {
                cuba_statefile_ += std::to_string(k_a_);
            }
            if (spectrum_ == BISPECTRUM) {
                if (k_b_idx != -1) {
                    cuba_statefile_ += "_" + std::to_string(k_b_idx);
                }
                else {
                    cuba_statefile_ += "_" + std::to_string(k_b_);
                }
                if (k_c_idx != -1) {
                    cuba_statefile_ += "_" + std::to_string(k_c_idx);
                }
                else {
                    cuba_statefile_ += "_" + std::to_string(k_c_);
                }
            }
            cuba_statefile_ += ".state";
        }
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        throw ConfigException(
            "No setting for CUBA statefile or path in configuration.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException(
            "Encountered type exception for CUBA statefile setting.");
    }
    catch (const ConfigException& ex) {
        throw ex;
    }
}



std::ostream& operator<<(std::ostream& out, const Config& cfg) {
    if (cfg.spectrum() == POWERSPECTRUM) {
        out << "# Matter power spectrum P(k) at " << cfg.n_loops()
            << "-loop for k = " << cfg.k_a() << "\n";
    }
    else {
      out << "# Matter bispectrum B(k) at " << cfg.n_loops() << "-loop for\n"
          << "#\n# k_a    = " << cfg.k_a() << "\n"
          << "# k_b    = " << cfg.k_b() << "\n"
          << "# k_c    = " << cfg.k_c() << "\n"
          << "# cos_ab = " << cfg.cos_ab() << "\n";
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
    if (cfg.input_ps_rescale() != 1) {
        out << "# Input power spectrum rescaled by a factor ";
        if (!cfg.input_ps_rescale_str().empty()) {
            out << cfg.input_ps_rescale_str();
        }
        else {
            out << cfg.input_ps_rescale();
        }
        out << " before interpolation.";
    }

    out << std::scientific;
    out << "\n#\n# Integration limits:\n";
    out << "#\t q_min = " << cfg.q_min() << "\n";
    out << "#\t q_max = " << cfg.q_max() << "\n";

    out << "#\n# Cuba settings:\n";
    out << "#\t abs tolerance = " << cfg.cuba_atol() << "\n";
    out << "#\t rel tolerance = " << cfg.cuba_rtol() << "\n";
    out << "#\t max. evals    = " << cfg.cuba_maxevals() << "\n";

    if (!cfg.cuba_statefile().empty()) {
        out << "#\t statefile     = " << cfg.cuba_statefile() << "\n";
    }

    if (cfg.dynamics() == EVOLVE_IC_ASYMP || cfg.dynamics() == EVOLVE_IC_EDS) {
        out << "#\n# ODE settings:\n";
        out << std::scientific;
        out << "#\t abs tolerance = " << cfg.ode_atol() << "\n";
        out << "#\t rel tolerance = " << cfg.ode_rtol() << "\n";
        out << "#\t start step    = " << cfg.ode_hstart() << "\n";

        out << "#\n# Time grid settings:\n";
        out << "#\t time steps     = " << cfg.time_steps() << "\n";
        if (cfg.dynamics() == EVOLVE_IC_ASYMP) {
            out << "#\t pre time steps = " << cfg.pre_time_steps() << "\n";
            out << "#\t eta asymp      = " << cfg.eta_asymp() << "\n";
        }
        out << "#\t eta ini        = " << cfg.eta_ini() << "\n";
        out << "#\t eta fin        = " << cfg.eta_fin() << "\n";

        out << "#\n# Dynamics settings:\n";
        if (cfg.dynamics() == EVOLVE_IC_ASYMP) {
            out << "# f_nu      = " << cfg.f_nu() << "\n";
            out << "# omega_m_0 = " << cfg.omega_m_0() << "\n";
        }

        out << "#\n# zeta file              = " << cfg.zeta_file() << "\n";
        if (cfg.dynamics() == EVOLVE_IC_ASYMP) {
            out << "# redshift file          = " << cfg.redshift_file() << "\n";
            out << "# omega eigenvalues file = " << cfg.omega_eigenvalues_file()
                << "\n";
            for (size_t i = 0; i < COMPONENTS; ++i) {
                out << "# F1 ic files[" << i
                    << "]\t\t = " << cfg.F1_ic_files().at(i) << "\n";
            }
            out << "# effective cs2 x grid   = " << cfg.effcs2_x_grid() << "\n";
            out << "# effective cs2 y grid   = " << cfg.effcs2_y_grid() << "\n";
            out << "# effective cs2 data     = " << cfg.effcs2_data() << "\n";
        }
    }
    out << "#\n# Information from Cuba integration:\n";
    if (cfg.cuba_fail() == 0) {
        out << "#\t Accuracy reached.\n";
    }
    else if (cfg.cuba_fail() > 0) {
        out << "#\t Accuracy not reached.\n";
    }
    else if (cfg.cuba_fail() == -1) {
        out << "#\t Error: dimension out of range.\n";
    }
    else {
        out << "#\t Unknown error status from CUBA.\n";
    }
    out << "#\t Num. evaluations = " << cfg.cuba_evals() << "\n";
    out << "#\t Num. subregions  = " << cfg.cuba_subregions() << "\n";

    return out;
}



/* Helper conversions function */
static inline int intpow(int a, int b) {
    return static_cast<int>(pow(a,b));
}
static inline size_t uintpow(int a, int b) {
    return static_cast<size_t>(pow(a,b));
}



LoopParameters::LoopParameters(int n_loops, Spectrum spectrum, Dynamics dynamics)
    : dynamics_(dynamics), spectrum_(spectrum), n_loops_(n_loops)
{
    if (n_loops_ < 0 || n_loops_ > 2) {
        throw(std::invalid_argument("LoopParameters::LoopParameters(): "
                                    "implementation for n_loops = 0,1,2 only."));
    }

    if (spectrum_ == POWERSPECTRUM) {
        n_coeffs_      = static_cast<size_t>(n_loops_) + 1;
        n_configs_     = uintpow(3, static_cast<int>(n_coeffs_));
        n_kernels_     = (n_configs_/3 + 1) * uintpow(4,n_loops_);
        n_kernel_args_ = 2 * static_cast<size_t>(n_loops_) + 1;
        zero_label_    = get_zero_label(n_coeffs_);

        /* Define block size for kernel indexing. A block consists of all */
        /* single loop label argument combinations. */
        single_loop_block_size = intpow(4, n_loops_);
        /* Max/min single_loops labels */
        single_loop_label_min = static_cast<int>(n_configs_) / 3;
        single_loop_label_max = 2 * static_cast<int>(n_configs_) / 3 - 1;
    }
    else if (spectrum_ == BISPECTRUM) {
        n_coeffs_  = static_cast<size_t>(n_loops_) + 2;
        n_configs_ = uintpow(3, static_cast<int>(n_coeffs_));
        n_kernels_ = uintpow(static_cast<int>(n_configs_) + 1, 2) *
                     uintpow(4, n_loops_);
        n_kernel_args_ = 2 * static_cast<size_t>(n_loops_) + 2;
        zero_label_    = get_zero_label(n_coeffs_);

        single_loop_block_size = intpow(4, n_loops_);
        /* Max/min single_loops labels */
        single_loop_label_min = 4 * static_cast<int>(n_configs_) / 9;
        single_loop_label_max = 5 * static_cast<int>(n_configs_) / 9 - 1;
        /* Composite label has either k_a, k_b, or multiple loop momenta
         * (applicable for overall loop diagram) */
        first_composite_block_size = (static_cast<int>(n_configs_) + 1) *
            single_loop_block_size;
    }
    else {
        throw(std::invalid_argument(
            "LoopParameters::LoopParameters(): invalid spectrum."));
    }

    /* List of single loop labels */
    Vec1D<int> config(n_coeffs_, 0);
    int label;

    for (size_t i = 0; i < static_cast<size_t>(n_loops_); ++i) {
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

    for (size_t i = 0; i < n_kernel_args_; ++i) {
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

    for (size_t i = 0; i < n_kernel_args_; ++i) {
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
        const std::string& zeta_file,
        const std::string& redshift_file,
        const std::string& omega_eigenvalues_file,
        const Vec1D<std::string>& F1_ic_files,
        const std::string& effcs2_x_file,
        const std::string& effcs2_y_file,
        const std::string& effcs2_data_file,
        double ode_atol,
        double ode_rtol,
        double ode_hstart
        ) :
    f_nu_(f_nu), ode_atol_(ode_atol), ode_rtol_(ode_rtol),
    ode_hstart_(ode_hstart), redshift(redshift_file),
    omega_eigenvalues(omega_eigenvalues_file), effcs2(effcs2_x_file,
            effcs2_y_file, effcs2_data_file)
{
    /* zeta always enters with prefactor 1.5, hence we redefine and multiply
     * here once */
    zeta = Interpolation1D(zeta_file, 1.5);

    for (auto& F1_ic_file : F1_ic_files) {
        F1_ic.emplace_back(F1_ic_file);
    }

    sound_speed = EXACT;

    /* cs2_factor = 2/3 * 1/(omegaM a^2 H^2) */
    cs2_factor = 2.0/3.0 * SQUARE(3e3) / omega_m_0;
}



EvolutionParameters::EvolutionParameters(
        double m_nu,
        double f_nu,
        double omega_m_0,
        const std::string& zeta_file,
        const std::string& redshift_file,
        const std::string& omega_eigenvalues_file,
        const Vec1D<std::string>& F1_ic_files,
        double ode_atol,
        double ode_rtol,
        double ode_hstart
        ) :
    f_nu_(f_nu), ode_atol_(ode_atol), ode_rtol_(ode_rtol),
    ode_hstart_(ode_hstart), redshift(redshift_file),
    omega_eigenvalues(omega_eigenvalues_file)
{
    /* zeta always enters with prefactor 1.5, hence we redefine and multiply
     * here once */
    zeta = Interpolation1D(zeta_file, 1.5);

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
        const std::string& zeta_file,
        double ode_atol,
        double ode_rtol,
        double ode_hstart
        ) :
    ode_atol_(ode_atol), ode_rtol_(ode_rtol), ode_hstart_(ode_hstart),
    zeta(zeta_file)
{
    /* zeta always enters with prefactor 1.5, hence we redefine and multiply
     * here once */
    zeta = Interpolation1D(zeta_file, 1.5);
}
