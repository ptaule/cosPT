#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <sstream>
#include <stdexcept>
#include <utility>

#if (__cplusplus > 201703L)
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include <libconfig.h++>

#include "../include/utilities.hpp"
#include "../include/io.hpp"
#include "../include/version.hpp"
#include "../include/parameters.hpp"

using std::size_t;
using std::pow;
using std::string;


/* Special definition for double, allowing also int and converting */
template<>
bool Config::set_param_value<double>(
    const libconfig::Config& cfg,
    const std::string& param,
    bool required
    )
{
    /* Try double */
    try {
        set<double>(param, cfg.lookup(param));
        return true;
    }
    catch (const libconfig::SettingTypeException& tex) {
        /* If not recognized, try int */
        try {
            int value = cfg.lookup(param);
            set<double>(param, static_cast<int>(value));
            return true;
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception parsing " + param +
                    " setting. ");
        }
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        if (required) {
            throw ConfigException("No " + param + " (required) setting found in "
                    "configuration file.");
        }
    }
    return true;
}



/* Special definition for size_t, converting from int */
template<>
bool Config::set_param_value<size_t>(
    const libconfig::Config& cfg,
    const std::string& param,
    bool required
    )
{
    try {
        int value = cfg.lookup(param);
        if (value < 0) {
            throw ConfigException("Negative value for " + param + " setting "
                    "not accepted.");
        }
        set<size_t>(param, static_cast<size_t>(value));
    }
    catch (const ConfigException& e) {
        throw e;
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        if (required) {
            throw ConfigException("No " + param + " (required) setting found in "
                                  "configuration file.");
        }
        return false;
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception parsing " + param +
                              " setting.");
    }
    return true;
}



template<typename T>
bool Config::set_param_value(
    const libconfig::Config& cfg,
    const string& param,
    bool required
    )
{
    try {
        set<T>(param, cfg.lookup(param));
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        if (required) {
            throw ConfigException("No " + param + " (required) setting found in "
                                  "configuration file.");
        }
        return false;
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception parsing " + param +
                              " setting.");
    }
    return true;
}



template<typename T>
bool Config::set_param_value(
    const libconfig::Setting& cfg,
    const string& prefix,
    const string& param,
    bool required
    )
{
    try {
        set<T>(prefix + "_" + param, cfg.lookup(param));
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        if (required) {
            throw ConfigException("No " + param + " (required) setting found in "
                                  "configuration file.");
        }
        return false;
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception parsing " + param +
                              " setting.");
    }
    return true;
}



template<typename T>
T Config::get_param_value(
    const libconfig::Config& cfg,
    const string& param,
    bool required
    )
{
    try {
        return cfg.lookup(param);
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        if (required) {
            throw ConfigException("No " + param + " (required) setting found in "
                                  "configuration file.");
        }
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception parsing " + param +
                              " setting.");
    }
    return T();
}



Vec1D<string> Config::keys_not_recognized(const libconfig::Config& cfg) const {
    const libconfig::Setting& root = cfg.getRoot();
    Vec1D<string> keys;
    /* Keys that are converted into other options in Config, these are allowed
     * as well */
    Vec1D<string> additional_keys = {
        "loops",
        "k_a_grid",
        "k_a_grid_file",
        "k_b_grid",
        "k_b_grid_file",
        "k_c_grid",
        "k_c_grid_file",
        "correlations",
        "input_ps_rescale",
        "cuba_settings",
        "ode_settings",
        "f_nu",
        "omega_m_0",
        "redshift_file",
        "omega_eigenvalues_file",
        "zeta_file",
        "effective_cs2_files",
        "F1_ic_files",
    };

    for (int i = 0; i < root.getLength(); ++i) {
        string key = root[i].getName();
        if (params.find(key) == params.end() &&
                std::find(additional_keys.begin(), additional_keys.end(), key)
                == additional_keys.end())
        {
            keys.push_back(key);
        }
    }
    return keys;
}



void Config::set_spectrum(const libconfig::Config& cfg)
{
    string spectrum_str = get_param_value<string>(cfg, "spectrum", true);

    std::transform(spectrum_str.begin(), spectrum_str.end(),
            spectrum_str.begin(),
            [](unsigned char c) {return tolower(c);}
            );

    if (spectrum_str == "powerspectrum") {
        set("spectrum", POWERSPECTRUM);
    }
    else if (spectrum_str == "bispectrum") {
        set("spectrum", BISPECTRUM);
    }
    else {
        throw ConfigException("Unknown spectrum parsed from configuration file.");
    }

    try {
        const libconfig::Setting& correlation = cfg.lookup("correlations");
        int count = correlation.getLength();
        if (count == 0 && !get<bool>("rsd")) {
            std::cout << "Info: correlations setting empty, computing";
            if (get<Spectrum>("spectrum") == POWERSPECTRUM) {
                std::cout << " [0,0] ";
                pair_correlations_.push_back({0,0});
            }
            else if (get<Spectrum>("spectrum") == BISPECTRUM) {
                std::cout << " [0,0,0] ";
                triple_correlations_.push_back({0,0,0});
            }
            std::cout << std::endl;
        }
        else if (count > 0 && get<bool>("rsd")) {
            std::cout << "Warning: correlations setting ignored for RSD."
                << std::endl;
        }
        else {
            for (int i = 0; i < count; ++i) {
                int a = 0;
                int b = 0;
                int c = 0;

                if (get<Spectrum>("spectrum") == POWERSPECTRUM) {
                    if (correlation[i].getLength() != 2) {
                        throw ConfigException(
                            "Correlation must be two indices (powerspectrum)");
                    }
                    a = correlation[i][0];
                    b = correlation[i][1];
                    pair_correlations_.push_back({a,b});
                }
                else if (get<Spectrum>("spectrum") == BISPECTRUM) {
                    if (correlation[i].getLength() != 3) {
                        throw ConfigException(
                            "Correlation must be three indices (bispectrum)");
                    }
                    a = correlation[i][0];
                    b = correlation[i][1];
                    c = correlation[i][1];
                    triple_correlations_.push_back({a,b,c});
                }
                if (a < 0 || a >= COMPONENTS ||
                    b < 0 || b >= COMPONENTS ||
                    c < 0 || c >= COMPONENTS
                    ) {
                    throw ConfigException("a, b and c must be element in [0," +
                            std::to_string(COMPONENTS) + "]");
                }
            }
        }
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        if (!get<bool>("rsd")) {
            std::cout << "Info: correlations setting empty, computing";
            if (get<Spectrum>("spectrum") == POWERSPECTRUM) {
                std::cout << " [0,0] ";
                pair_correlations_.push_back({0,0});
            }
            else if (get<Spectrum>("spectrum") == BISPECTRUM) {
                std::cout << " [0,0,0] ";
                triple_correlations_.push_back({0,0,0});
            }
            std::cout << std::endl;
        }
    }
    catch (const ConfigException& e) {
        throw e;
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception for correlations setting.");
    }
}



void Config::set_bispectrum_ext_momenta(const libconfig::Config& cfg)
{
    double k_a = get<double>("k_a");

    string k_b_grid_file;
    Vec2D<double> k_b_grid;
    if (cfg.lookupValue("k_b_grid", k_b_grid_file)) {
        read_delimited_file(k_b_grid_file, k_b_grid);
    }
    /* If k_b_idx != -1, i.e. given as command line option, use this value */
    if (!k_b_grid.empty() && (k_b_idx != -1 || cfg.lookupValue("k_b_idx", k_b_idx))) {
        try {
            set("k_b", k_b_grid.at(static_cast<size_t>(k_b_idx)).at(0));
        }
        catch (const std::out_of_range& e) {
            throw ConfigException("k_b_idx out of range (of k_b_grid).");
        }
    }
    else if (cfg.exists("k_b")) {
        set_param_value<double>(cfg, "k_b");
    }
    else {
        set("k_b", k_a);
        std::cout << "Info: no k_b read, setting k_b = k_a." << std::endl;
    }

    double k_b = get<double>("k_b");
    double k_c = 0;
    double cos_ab = 0;
    bool k_c_given = false;
    bool cos_ab_given = false;

    string k_c_grid_file;
    Vec2D<double> k_c_grid;
    if (cfg.lookupValue("k_c_grid", k_c_grid_file)) {
        read_delimited_file(k_c_grid_file, k_c_grid);
    }
    /* If k_c_idx != -1, i.e. given as command line option, use this value */
    if (!k_c_grid.empty() && (k_c_idx != -1 || cfg.lookupValue("k_c_idx", k_c_idx))) {
        try {
            k_c = k_c_grid.at(static_cast<size_t>(k_c_idx)).at(0);
        }
        catch (const std::out_of_range& e) {
            throw ConfigException("k_c_idx out of range (of k_c_grid).");
        }
        k_c_given = true;
    }
    else if (cfg.exists("k_c")) {
        set_param_value<double>(cfg, "k_c");
        k_c = get<double>("k_c");
        k_c_given = true;
    }

    if (cfg.exists("cos_ab")) {
        set_param_value<double>(cfg, "cos_ab");
        cos_ab = get<double>("cos_ab");
        cos_ab_given = true;

        if (cos_ab < -1 || cos_ab > 1) {
            throw ConfigException(
                "Got cos_ab = " + std::to_string(cos_ab) + ", which is not"
                "between -1 and 1.");
        }
    }

    /* If neither/both k_c and cos_ab was given, set k_c = k_a */
    if ((!k_c_given) && (!cos_ab_given)) {
        set("k_c", k_a);
        k_c = k_a;
        std::cout << "Info: no k_c read, setting k_c = k_a." << std::endl;
    }
    if (k_c_given && cos_ab_given) {
        throw ConfigException(
            "Both k_c and cos_ab given. Please choose one.");
    }

    /* k_c = - k_a - k_b */
    if (k_c_given) {
        cos_ab = 0.5 * ((k_c*k_c)/(k_a*k_b) - k_a/k_b - k_b/k_a);
        set<double>("cos_ab", cos_ab);

        if (cos_ab < -1 || cos_ab > 1) {
            throw ConfigException(
                "The k_a, k_b, k_c values given does not constitute a valid "
                "configuration (cos_ab out of range).");
        }
    }
    else {
        /* cos_ab given */
        k_c = std::sqrt(SQUARE(k_a) + SQUARE(k_b) + 2*k_a*k_b*cos_ab);
        set<double>("k_c", k_c);
    }
}



void Config::set_input_ps(const libconfig::Config& cfg)
{
    /* Input power spectrum */
    set_param_value<string>(cfg, "input_ps_file", true);
    try {
        /* input_ps_rescale setting. Rescaling factor set to 1 by default. */
        /* First, check if given as string "2pi^-3" or "2pi^3" */
        string r_str;
        double r_double;
        int r_int;
        if (cfg.lookupValue("input_ps_rescale", r_str)) {
            set<string>("input_ps_rescale_str", r_str);
            if (r_str.compare("2pi^-3") == 0) {
                set("input_ps_rescale_num", pow(TWOPI,-3));
            }
            else if (r_str.compare("2pi^3") == 0) {
                set("input_ps_rescale_num", pow(TWOPI,3));
            }
            else {
                throw ConfigException("Got input_ps_rescale string \"" + r_str +
                                      "\" which is neither \"2pi^-3\" nor \"2pi^3\".");
            }
        }
        /* If not found as string, try double */
        else if (cfg.lookupValue("input_ps_rescale", r_double)) {
            set("input_ps_rescale_num", r_double);
        }
        /* Last, try integer */
        else if (cfg.lookupValue("input_ps_rescale", r_int)) {
            set("input_ps_rescale_num", r_int);
        }
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException(
            "Encountered type exception parsing input_ps_rescale setting.");
    }
    catch (const ConfigException& e) {
        throw e;
    }
}



void Config::set_output_file(const libconfig::Config& cfg)
{
    try {
        /* Instead of running the entire integration before checking whether
         * the output file/path can be written to, perform certain checks
         * immediately */
        string output_path;
        if (set_param_value<string>(cfg, "output_file")) {
            /* For c++ version > 2017 use routines in <filesystem> to check
             * that directory exists */
#if (__cplusplus > 201703L)
            fs::path p(get<string>("output_file"));
            if (p.has_parent_path() && !fs::exists(p.parent_path())) {
                throw ConfigException("Output file directory " +
                                      string(p.parent_path()) +
                                      " does not exist.");
            }
#endif
        }
        else if (cfg.lookupValue("output_path", output_path)) {
            /* For c++ version > 2017 use routines in <filesystem> to check
             * that directory exists */
#if (__cplusplus > 201703L)
            if (!fs::exists(fs::path(output_path))) {
                throw ConfigException("Output directory \"" + output_path +
                                      "\" does not exist.");
            }
            if (!fs::is_directory(output_path)) {
                throw ConfigException("Output directory \"" + output_path +
                                      "\" is not a directory.");
            }
#endif

            set<string>("output_file",
                        create_filename_from_wavenumbers(output_path, ".dat"));
        }
        else {
            throw ConfigException("No output path/file given in configuration.");
        }
    }
    catch (const libconfig::SettingNotFoundException& nfex) {
        throw ConfigException("No output_file or output_path (required) found "
                              "in configuration file.");
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException("Encountered type exception parsing output_file setting.");
    }
    catch (const ConfigException& ex) {
        throw ex;
    }
}



void Config::set_dynamics(const libconfig::Config& cfg)
{
    Dynamics dynamics;
    string dynamics_str = get_param_value<string>(cfg, "dynamics", true);

    std::transform(dynamics_str.begin(), dynamics_str.end(),
            dynamics_str.begin(), [](unsigned char c) {return
            tolower(c);});

    if (dynamics_str == "eds-spt") {
        dynamics = EDS_SPT;
    }
    else if (dynamics_str == "evolve-ic-asymp") {
        dynamics = EVOLVE_IC_ASYMP;
    }
    else if (dynamics_str == "evolve-ic-eds") {
        dynamics = EVOLVE_IC_EDS;
    }
    else {
        throw ConfigException("Unknown dynamics in configuration file.");
    }
    set("dynamics", dynamics);

    /* Settings for evolution dynamics */
    if (dynamics == EVOLVE_IC_ASYMP || dynamics == EVOLVE_IC_EDS) {
        /* ODE settings */
        try {
            if (cfg.exists("ode_settings")) {
                const libconfig::Setting& ode_settings = cfg.lookup("ode_settings");
                set<double>("ode_atol", ode_settings.lookup("abs_tolerance"));
                set<double>("ode_rtol", ode_settings.lookup("rel_tolerance"));
                set<double>("ode_hstart", ode_settings.lookup("start_step"));
            }
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception parsing ODE settings.");
        }

        set_param_value<size_t>(cfg, "time_steps", true);
        set_param_value<double>(cfg, "eta_ini", true);
        set_param_value<double>(cfg, "eta_fin", true);

        /* Need additional information about time steps before eta_ini to set ICs */
        if (dynamics == EVOLVE_IC_ASYMP) {
            set_param_value<size_t>(cfg, "pre_time_steps", true);
            set_param_value<double>(cfg, "eta_asymp", true);
        }

        /* Read files containing scale/time-dependent functions depending on cosmology.
         * zeta(eta) depends only on time, effcs2(eta, k) depends also on scale */
        set_param_value<string>(cfg, "zeta_file", true);

        try {
            if (dynamics == EVOLVE_IC_ASYMP) {
                set_param_value<string>(cfg, "redshift_file", true);
                set_param_value<string>(cfg, "omega_eigenvalues_file");

                const libconfig::Setting& F1_ic_files_list = cfg.lookup("F1_ic_files");
                int count = F1_ic_files_list.getLength();
                if (count != COMPONENTS) {
                    throw ConfigException("There should be " +
                                          std::to_string(COMPONENTS) + " F1 ic files in the configuration");
                }
                for (int i = 0; i < count; ++i) {
                    F1_ic_files_.push_back(F1_ic_files_list[i].c_str());
                }

                if (cfg.exists("effective_cs2_files")) {
                    const libconfig::Setting& effcs2_files = cfg.lookup("effective_cs2_files");
                    set<string>("effcs2_x_grid", effcs2_files.lookup("x_grid"));
                    set<string>("effcs2_y_grid", effcs2_files.lookup("y_grid"));
                    set<string>("effcs2_data", effcs2_files.lookup("data"));
                }
            }
        }
        catch (const libconfig::SettingNotFoundException& nfex) {
            throw ConfigException(
                "Missing file for interpolation in configuration.");
        }
        catch (const libconfig::SettingTypeException& tex) {
            throw ConfigException("Encountered type exception parsing interpolation files.");
        }
        if (dynamics == EVOLVE_IC_ASYMP) {
            set_param_value<double>(cfg, "f_nu", true);
            set_param_value<double>(cfg, "omega_m_0", true);
        }
    }
}



void Config::set_cuba_config(
    const libconfig::Setting& cuba_settings,
    int cuba_max_evaluations,
    int cuba_cores
    )
{
    set_param_value<double>(cuba_settings, "cuba", "abs_tolerance");
    set_param_value<double>(cuba_settings, "cuba", "rel_tolerance");
    set_param_value<int>(cuba_settings, "cuba", "verbosity_level");
    set_param_value<bool>(cuba_settings, "cuba", "retain_statefile");

    /* If cuba_max_evaluations is not already set, look up value */
    if (cuba_max_evaluations == 0) {
        /* Try int */
        try {
            set<int>("cuba_max_evaluations", cuba_settings.lookup("max_evaluations"));
        }
        catch (const libconfig::SettingTypeException& tex) {
            /* If not recognized, try double */
            try {
                double value = cuba_settings.lookup("max_evaluations");
                if (value > std::numeric_limits<int>::max() || value < std::numeric_limits<int>::min()) {
                    throw std::overflow_error("Value is out of range for int");
                }
                set<int>("cuba_max_evaluations", static_cast<int>(std::round(value)));
            }
            catch (std::overflow_error& oe) {
                throw ConfigException("Parsing cuba_settings: "
                        "max_evaluations: value is out of range for int");
            }
            catch (const libconfig::SettingTypeException& tex) {
                throw ConfigException("Encountered type exception parsing "
                        "cuba_settings: max_evaluations setting.");
            }
        }
        catch (const libconfig::SettingNotFoundException& nfex) {
            std::cout << "No cuba max. evaluations given. Using default value: "
                << get<double>("cuba_maxeval") << std::endl;
        }
    }
    /* If cuba_cores is not already set, look up value */
    if (cuba_cores == 0 && !set_param_value<int>(cuba_settings, "cuba", "cores")) {
        std::cout << "No CUBA number of cores option given. Using default value: "
            << get<double>("cuba_cores")
            << " (0 defaults to the number of cores available)"<< std::endl;
    }
}



void Config::set_cuba_statefile(const libconfig::Setting& cuba_settings)
{
    try {
        string statefile_path;
        if (set_param_value<string>(cuba_settings, "cuba", "statefile")) {
            /* For c++ version > 2017 use routines in <filesystem> to check
             * that directory exists */
#if (__cplusplus > 201703L)
            fs::path p(get<string>("statefile"));
            if (!fs::exists(p.parent_path())) {
                throw ConfigException("CUBA statefile file directory " +
                                      string(p.parent_path()) + " does not exist.");
            }
#endif
        }
        else if (cuba_settings.lookupValue("statefile_path", statefile_path)) {
            /* For c++ version > 2017 use routines in <filesystem> to check
             * that directory exists */
#if (__cplusplus > 201703L)
            if (!fs::exists(fs::path(statefile_path))) {
                throw ConfigException("CUBA statefile directory \"" +
                                      statefile_path + "\" does not exist.");
            }
            if (!fs::is_directory(statefile_path)) {
                throw ConfigException("CUBA statefile directory \"" +
                                      statefile_path + "\" is not a directory.");
            }
#endif

            set<string>("output_file",
                        create_filename_from_wavenumbers(statefile_path, ".state"));
        }
    }
    catch (const libconfig::SettingTypeException& tex) {
        throw ConfigException(
            "Encountered type exception for CUBA statefile setting.");
    }
    catch (const ConfigException& ex) {
        throw ex;
    }
}



string Config::create_filename_from_wavenumbers(
    const string& base_path,
    const string& file_extension
    )
{
    std::stringstream ss;
    ss << base_path;
    ss << "/";
    /* Add k_a_idx (& k_b_idx) to end of file */
    if (k_a_idx != -1) {
        ss << std::setfill('0') << std::setw(3) << k_a_idx;
    }
    else {
        ss << std::scientific << std::setprecision(6) << get<double>("k_a");
    }
    if (get<Spectrum>("spectrum") == BISPECTRUM) {
        ss << "_";
        if (k_b_idx != -1) {
            ss << std::setfill('0') << std::setw(3) << k_b_idx;
        }
        else {
            ss << std::scientific << std::setprecision(6) << get<double>("k_b");
        }
        ss << "_";
        if (k_c_idx != -1) {
            ss << std::setfill('0') << std::setw(3) << k_c_idx;
        }
        else {
            ss << std::scientific << std::setprecision(6) << get<double>("k_c");
        }
    }
    ss << file_extension;
    return ss.str();
}



Config::Config()
{
    /* Set default values */
    set("n_loops", 0);
    set("dynamics", EDS_SPT);
    set("spectrum", POWERSPECTRUM);

    set("single_hard_limit", false);
    set("sh_Q1", 10.0);

    set("k_a", 0.0);
    set("k_b", 0.0);
    set("k_c", 0.0);
    set("cos_ab", 0.0);

    /* Integration limits */
    set("q_min", 1e-4);
    set("q_max", 1.0);

    set("input_ps_file", string());
    set("input_ps_rescale_num", 1.0);
    set("input_ps_rescale_str", string());

    set("output_path", string());
    set("output_file", string());

    set("rsd", false);
    set("rsd_growth_f", 1.0);

    set("ir_resum", false);
    /* Values from 1605.02149 */
    set("k_s", 0.2);
    set("k_osc", 1.0/110.0);
    /* Which PT order are we working at? Relevant for avoiding overcounting
     * of IR contributions */
    set("pt_order", 0);

    set("compute_eft_displacement_dispersion", false);
    set("eft_displacement_dispersion", 0);
    /* Upper limit for EFT displacement dispersion integral. Lower limit = cutoff */
    set("eft_displacement_dispersion_infty_", 10);

    set("cuba_abs_tolerance", 1e-12);
    set("cuba_rel_tolerance", 1e-4);
    set<int>("cuba_max_evaluations", 1e6);
    set("cuba_verbosity_level", 0);
    set("cuba_n_cores", 0);
    set("cuba_retain_statefile", false);
    set("cuba_statefile", string());
    set("cuba_statefile_path", string());

    set("description", string());

    set("ode_abs_tolerance", 1e-6);
    set("ode_rel_tolerance", 1e-4);
    set("ode_start_step", 1e-3);

    set("eta_ini", 0.0);
    set("eta_fin", 0.0);
    set<size_t>("time_steps", 0);
    set<size_t>("pre_time_steps", 0);
    set("eta_asymp", 0.0);

    set("omega_eigenmode", 0);
    set("omega_k_min", get<double>("q_min"));
    set("omega_k_max", get<double>("q_max"));
    set("omega_N", 100);
    set("omega_imag_threshold", 1e-3);
}



Config::Config(const string& ini_file,
        int k_a_idx,
        int k_b_idx,
        int k_c_idx,
        int cuba_max_evaluations,
        int cuba_cores
        )
    : Config()
{
    this->k_a_idx = k_a_idx;
    this->k_b_idx = k_b_idx;
    this->k_c_idx = k_c_idx;

    libconfig::Config cfg;
    /* Read config file */
    try {
        cfg.readFile(ini_file.c_str());
    }
    catch (const libconfig::FileIOException& ioex) {
        throw ConfigException("Unable to read \"" + ini_file + "\".");
    }
    catch (const libconfig::ParseException& pex) {
        throw ConfigException("Parse error at " + string(pex.getFile()) +
                              ":" + std::to_string(pex.getLine()) + " - " +
                              string(pex.getError()));
    }

    Vec1D<string> keys = keys_not_recognized(cfg);
    if (!keys.empty()) {
        string error_msg = "Keys in configuration file not recognized: ";
        for (auto key : keys) {
            error_msg += key + ", ";
        }
        throw ConfigException(error_msg);
    }

    /* Number of loops */
    int n_loops = get_param_value<int>(cfg, "loops", true);
    set("n_loops", n_loops);

    /* First wavenumber */
    string k_a_grid_file;
    Vec2D<double> k_a_grid;
    if (cfg.lookupValue("k_a_grid", k_a_grid_file)) {
        read_delimited_file(k_a_grid_file, k_a_grid);
    }
    /* If k_a_idx != -1, i.e. given as command line option, use this value */
    if (!k_a_grid.empty() && (k_a_idx != -1 || cfg.lookupValue("k_a_idx", k_a_idx))) {
        try {
            set("k_a", k_a_grid.at(static_cast<size_t>(k_a_idx)).at(0));
        }
        catch (const std::out_of_range& e) {
            throw ConfigException("k_a_idx out of range (of k_a_grid).");
        }
    }
    else if (cfg.exists("k_a")) {
        set_param_value<double>(cfg, "k_a");
    }
    else {
        throw ConfigException("Did not obtain any value for k_a. Either provide "
                              "k_a, or k_a_idx and k_a_grid.");
    }

    /* Integration ranges */
    set_param_value<double>(cfg, "q_min", true);
    set_param_value<double>(cfg, "q_max", true);

    /* RSD. Read before spectrum settings, so that warning messages can be
     * given if rsd AND correlations are set */
    if(set_param_value<bool>(cfg, "rsd")) {
        if(!set_param_value<double>(cfg, "rsd_growth_f")) {
            std::cout <<
                "Info: No value for rsd_growth_f read, using default value = "
                << get<double>("rsd_growth_f") << "."<< std::endl;
        }
    }

    /* IR resummation */
    if(set_param_value<bool>(cfg, "ir_resum")) {
        set_param_value<double>(cfg, "k_s");
        set_param_value<double>(cfg, "k_osc");
        set_param_value<int>(cfg, "pt_order", true);
    }

    /* Compute eft_displacement_dispersion? */
    if(set_param_value<bool>(cfg, "compute_eft_displacement_dispersion")) {
        if(!set_param_value<double>(cfg, "eft_displacement_dispersion_infty")) {
            std::cout <<
                "Info: No value for eft_displacement_dispersion_infty read, "
                "using default value = " <<
                get<double>("eft_displacement_dispersion_infty") << "."<<
                std::endl;
        }
    }

    /* Compute single hard limit? */
    if(set_param_value<bool>(cfg, "single_hard_limit")) {
        if(!set_param_value<double>(cfg, "single_hard_limit_q")) {
            set<double>("sh_Q1", get<double>("q_max") * 10);
            std::cout << "Info: Single hard limit: Setting q = 10 * q_max = "
                << get<double>("sh_Q1") << "."<< std::endl;
        }
        else if (get<double>("sh_Q1") < get<double>("q_max")) {
            std::cout << "Warning: Single hard limit: Fixed q = " <<
                get<double>("sh_Q1") << " is less than q_max = " <<
                get<double>("q_max") << "."<< std::endl;
        }
    }

    /* Store potential description */
    set_param_value<string>(cfg, "description");

    if (cfg.exists("cuba_settings")) {
        const libconfig::Setting& cuba_settings = cfg.lookup("cuba_settings");
        set_cuba_config(cuba_settings, cuba_max_evaluations, cuba_cores);
        set_cuba_statefile(cuba_settings);
    }

    set_spectrum(cfg);
    if (get<Spectrum>("spectrum") == BISPECTRUM) {
        set_bispectrum_ext_momenta(cfg);
    }
    set_input_ps(cfg);
    set_dynamics(cfg);
    set_output_file(cfg);
}



std::ostream& operator<<(std::ostream& out, const Config& c) {
    if (c.get<Spectrum>("spectrum") == POWERSPECTRUM) {
        if (c.get<bool>("rsd")) {
            out << "# Power spectrum multipoles Pl(k) in redshift space at "
                << c.get<int>("n_loops")
                << "-loop for k = " << c.get<double>("k_a") << "\n";
        }
        else {
            out << "# Matter power spectrum P(k) at " << c.get<int>("n_loops")
                << "-loop for k = " << c.get<double>("k_a") << "\n";
        }
    }
    else {
        out << "# Matter bispectrum B(k) at " << c.get<int>("n_loops") << "-loop for\n"
            << "#\n# k_a    = " << c.get<double>("k_a") << "\n"
            << "# k_b    = " << c.get<double>("k_b") << "\n"
            << "# k_c    = " << c.get<double>("k_c") << "\n"
            << "# cos_ab = " << c.get<double>("cos_ab") << "\n";
    }
    out << "#\n";
    out << "# Description: " << c.get<string>("description") << "\n";
    out << "# Git hash:    " << build_git_sha     << "\n";
    out << "# Build time:  " << build_git_time    << "\n";
    out << "#\n";

    if (!c.get<bool>("rsd")) {
        out << "# Correlations (zero-indexed components):\n# ";
        if (c.get<Spectrum>("spectrum") == POWERSPECTRUM) {
            for (auto& el : c.pair_correlations()) {
                out << " " << el <<  " ,";
            }
        }
        else {
            for (auto& el : c.triple_correlations()) {
                out << " " << el <<  " ,";
            }
        }
        out << "\n#\n";
    }

    if (c.get<bool>("single_hard_limit")) {
        out << "# Computing single-hard limit with fixed Q1 = "
            << c.get<double>("sh_Q1")
            << "\n#\n";
    }

    out << "# Input power spectrum read from " << c.get<string>("input_ps_file") << "\n";
    if (c.get<double>("input_ps_rescale") != 1) {
        out << "# Input power spectrum rescaled by a factor ";
        if (!c.get<string>("input_ps_rescale_str").empty()) {
            out << c.get<string>("input_ps_rescale_str");
        }
        else {
            out << c.get<double>("input_ps_rescale");
        }
        out << " before interpolation." << "\n";
    }
    out << "#\n";

    out << std::scientific;
    out << "# Integration limits:\n";
    out << "#\t q_min = " << c.get<double>("q_min") << "\n";
    out << "#\t q_max = " << c.get<double>("q_max") << "\n";

    out << "#\n# Cuba settings:\n";
    out << "#\t abs tolerance = " << c.get<double>("cuba_abs_tolerance") << "\n";
    out << "#\t rel tolerance = " << c.get<double>("cuba_rel_tolerance") << "\n";
    out << "#\t max. evals    = " << c.get<int>("cuba_max_evaluations") << "\n";

    if (!c.get<string>("cuba_statefile").empty()) {
        out << "#\t statefile     = " << c.get<string>("cuba_statefile") << "\n";
    }
    out << "#\n";

    if (c.get<bool>("rsd")) {
        out << "# (RSD) growth_f = " << c.get<double>("rsd_growth_f") << "\n#\n";
    }
    if (c.get<bool>("ir_resum")) {
        out << "# IR resummation on, working at PT order N = " <<
            c.get<int>("pt_order") << "\n";
        out << "#\tk_s   = " << c.get<double>("k_s") << "\n";
        out << "#\tk_osc = " << c.get<double>("k_osc") << "\n";
        out << "#\n";
    }

    if (c.get<Dynamics>("dynamics") == EVOLVE_IC_ASYMP
        || c.get<Dynamics>("dynamics") == EVOLVE_IC_EDS) {
        out << "# ODE settings:\n";
        out << std::scientific;
        out << "#\t abs tolerance = " << c.get<double>("ode_abs_tolerance") << "\n";
        out << "#\t rel tolerance = " << c.get<double>("ode_rel_tolerance") << "\n";
        out << "#\t start step    = " << c.get<double>("ode_start_step") << "\n#\n";

        out << "# Time grid settings:\n";
        out << "#\t time steps     = " << c.get<size_t>("time_steps") << "\n";
        if (c.get<Dynamics>("dynamics") == EVOLVE_IC_ASYMP) {
            out << "#\t pre time steps = " << c.get<size_t>("pre_time_steps") << "\n";
            out << "#\t eta asymp      = " << c.get<double>("eta_asymp") << "\n";
        }
        out << "#\t eta ini        = " << c.get<double>("eta_ini") << "\n";
        out << "#\t eta fin        = " << c.get<double>("eta_fin") << "\n#\n";

        out << "# Dynamics settings:\n";
        if (c.get<Dynamics>("dynamics") == EVOLVE_IC_ASYMP) {
            out << "# f_nu      = " << c.get<double>("f_nu") << "\n";
            out << "# omega_m_0 = " << c.get<double>("omega_m_0") << "\n";
        }
        out << "#\n";

        out << "# zeta file              = " << c.get<string>("zeta_file") << "\n";
        if (c.get<Dynamics>("dynamics") == EVOLVE_IC_ASYMP) {
            out << "# redshift file          = " << c.get<string>("redshift_file") << "\n";
            out << "# omega eigenvalues file = " << c.get<string>("omega_eigenvalues_file")
                << "\n";
            for (size_t i = 0; i < COMPONENTS; ++i) {
                out << "# F1 ic files[" << i
                    << "]\t\t = " << c.F1_ic_files().at(i) << "\n";
            }
            out << "# effective cs2 x grid   = " << c.get<string>("effcs2_x_grid") << "\n";
            out << "# effective cs2 y grid   = " << c.get<string>("effcs2_y_grid") << "\n";
            out << "# effective cs2 data     = " << c.get<string>("effcs2_data") << "\n";
        }
        out << "#\n";
    }

    out << "# Information from Cuba integration:\n";
    if (c.cuba_fail() == 0) {
        out << "#\t Accuracy reached.\n";
    }
    else if (c.cuba_fail() > 0) {
        out << "#\t Accuracy not reached.\n";
    }
    else if (c.cuba_fail() == -1) {
        out << "#\t Error: dimension out of range.\n";
    }
    else {
        out << "#\t Unknown error status from CUBA.\n";
    }
    out << "#\t Num. evaluations = " << c.cuba_evals() << "\n";
    out << "#\t Num. subregions  = " << c.cuba_subregions() << "\n";

    if (c.get<bool>("compute_eft_displacement_dispersion")) {
        out << "#\n# Small scale displacement dispersion (relevant for EFT):\n";
        out << "# sigma_d^2 = \\int_{q_max}^{infty} dq P_{input}(q) = " <<
            c.get<double>("eft_displacement_dispersion") << "\twith infty = " <<
            c.get<double>("eft_displacement_dispersion_infty") << "\n";
    }
    out << "#\n";

    return out;
}



/* Helper conversions function */
static inline int intpow(int a, int b) {
    return static_cast<int>(pow(a,b));
}
static inline size_t uintpow(int a, int b) {
    return static_cast<size_t>(pow(a,b));
}



LoopParameters::LoopParameters(int n_loops, Spectrum spectrum,
        Dynamics dynamics, bool rsd)
    : dynamics_(dynamics), spectrum_(spectrum), rsd_(rsd), n_loops_(n_loops)
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

        /* Set argument-to-kernel index converter to a lambda function for the PS case */
        args_2_kernel_index = [this](const int args[]) {
            return this->ps_args_2_kernel_index(args);
        };
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

        /* Set argument-to-kernel index converter to a lambda function for the BS case */
        args_2_kernel_index = [this](const int args[]) {
            return this->bs_args_2_kernel_index(args);
        };
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



int LoopParameters::ps_args_2_kernel_index(const int arguments[]) const
{
   /* Precompute powers of two for speedup */
    int pow2[] = {1,2,4,8,16,32,64,128};

    // In DEBUG-mode, check that non-zero arguments (zero_label) are unique
#if DEBUG >= 1
    if (!unique_elements(arguments, n_kernel_args_, zero_label_))
        throw(std::logic_error(
            "LoopParameters::ps_args_2_kernel_index(): duplicate "
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
            throw(std::logic_error("LoopParameters::ps_args_2_kernel_index()"
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
                    "LoopParameters::ps_args_2_kernel_index(): argument is "
                    "neither 0, composite type, or single loop."));
#endif
        }
    }
#if DEBUG >= 1
    if (n_k_labels > 1)
        throw(std::logic_error("LoopParameters::ps_args_2_kernel_index(): "
                               "more than one argument is of composite type."));
#endif

    return index;
}



int LoopParameters::bs_args_2_kernel_index(const int arguments[]) const
{
   /* Precompute powers of two for speedup */
    int pow2[] = {1,2,4,8,16,32,64,128};

    /* In DEBUG-mode, check that non-zero_label arguments are unique */
#if DEBUG >= 1
    if (!unique_elements(arguments, n_kernel_args_, zero_label_))
        throw(std::logic_error(
            "LoopParameters::bs_args_2_kernel_index(): duplicate "
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
            "LoopParameters::bs_args_2_kernel_index(): more than two "
            "arguments is of composite type."));
#endif

    return index;
}



EvolutionParameters::EvolutionParameters(EvolutionParameters&& other) noexcept
    : f_nu_(other.f_nu_), cs2_factor(other.cs2_factor),
    ode_atol_(other.ode_atol_), ode_rtol_(other.ode_rtol_),
    ode_hstart_(other.ode_hstart_),
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
        ode_atol_   = other.ode_atol_;
        ode_rtol_   = other.ode_rtol_;
        ode_hstart_ = other.ode_hstart_;

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
        const string& zeta_file,
        const string& redshift_file,
        const string& omega_eigenvalues_file,
        const Vec1D<string>& F1_ic_files,
        const string& effcs2_x_file,
        const string& effcs2_y_file,
        const string& effcs2_data_file,
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

    /* cs2_factor = 2/3 * 1/(omegaM a^2 H^2) */
    cs2_factor = 2.0/3.0 * SQUARE(3e3) / omega_m_0;
}



EvolutionParameters::EvolutionParameters(
        const string& zeta_file,
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
