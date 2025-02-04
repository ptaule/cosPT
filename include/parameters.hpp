#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <functional>
#include <string>
#include <map>
#include <any>

#include "utilities.hpp"
#include "interpolation.hpp"

namespace libconfig {
class Config;
class Setting;
}


class ConfigException : std::exception {
    private:
        std::string ex;
    public:
        ConfigException (const std::string& ex) : ex(ex) {}
        const char* what() const noexcept {return ex.c_str();}
};


class Config {
    private:
        std::map<std::string, std::any> params;

        /* External wavenumber indices */
        int k_a_idx = -1;
        int k_b_idx = -1;
        int k_c_idx = -1;

        /* Information from CUBA integration */
        int cuba_evals_      = 0;
        int cuba_fail_       = 0;
        int cuba_subregions_ = 0;

        /* Correlations to compute */
        Vec1D<Pair<int>> pair_correlations_;     /* Power spectrum */
        Vec1D<Triple<int>> triple_correlations_; /* Bispectrum */

        Vec1D<std::string> F1_ic_files_;

        /* T is the type of the parameter */
        template<typename T>
        bool set_param_value(
            const libconfig::Config& cfg,
            const std::string& param,
            bool required = false
            );
        template<typename T>
        bool set_param_value(
            const libconfig::Setting& cfg,
            const std::string& prefix,
            const std::string& param,
            bool required = false
            );

        template<typename T>
        T get_param_value(
            const libconfig::Config& cfg,
            const std::string& param,
            bool required = false
            );

        Vec1D<std::string> keys_not_recognized(const libconfig::Config& cfg) const;

        void set_spectrum(const libconfig::Config& cfg);
        void set_bispectrum_ext_momenta(const libconfig::Config& cfg);

        void set_input_ps(const libconfig::Config& cfg);
        void set_output_file(const libconfig::Config& cfg);

        void set_dynamics(const libconfig::Config& cfg);

        void set_cuba_config(
            const libconfig::Setting& cuba_settings,
            int cuba_maxevals = 0,
            int cuba_cores = -1
            );
        void set_cuba_statefile(const libconfig::Setting& cuba_settings);
        std::string create_filename_from_wavenumbers(
            const std::string& base_path,
            const std::string& file_extension
        );
   public:
        Config();
        Config(const std::string& ini_file,
                int k_a_idx = -1,
                int k_b_idx = -1,
                int k_c_idx = -1,
                int cuba_maxevals = 0,
                int cuba_cores = -1
                );

        Vec1D<Pair<int>> pair_correlations() const {return pair_correlations_;}
        Vec1D<Triple<int>> triple_correlations() const {return triple_correlations_;}

        Vec1D<std::string> F1_ic_files() const {return F1_ic_files_;}

        /* Read CUBA info */
        int cuba_evals() const {return cuba_evals_;}
        int cuba_fail() const {return cuba_fail_;}
        int cuba_subregions() const {return cuba_subregions_;}

        /* Write CUBA info */
        int& cuba_evals() {return cuba_evals_;}
        int& cuba_fail() {return cuba_fail_;}
        int& cuba_subregions() {return cuba_subregions_;}

        /* Setters/getters for params map */
        template<typename T>
        void set(const std::string& key, const T& value) {
            params[key] = value;
        }

        template<typename T>
        T get(const std::string& key) const {
            auto it = params.find(key);
            if (it != params.end()) {
                try {
                    return std::any_cast<T>(it->second);
                } catch (const std::bad_any_cast&) {
                    throw ConfigException(
                    "Config::get(): bad cast with key = " + key);
                }
            }
            return T(); // Return default value if key not found or casting fails
        }
};


std::ostream& operator<<(std::ostream& out, const Config& cfg);


class LoopParameters {
    private:
        const Dynamics dynamics_;
        const Spectrum spectrum_;
        const bool rsd_;

        const int n_loops_;

        std::size_t n_coeffs_;
        std::size_t n_configs_;
        std::size_t n_kernels_;
        std::size_t n_kernel_args_;
        int zero_label_;
        int single_loop_label_min;
        int single_loop_label_max;
        int single_loop_block_size;

        Vec1D<int> single_loops;

        int first_composite_block_size = 0; /* Bispectrum */

        int ps_args_2_kernel_index(const int arguments[]) const;
        int bs_args_2_kernel_index(const int arguments[]) const;

    public:
        LoopParameters(int n_loops, Spectrum spectrum, Dynamics dynamics,
                bool rsd);

        Dynamics dynamics() const { return dynamics_; }
        Spectrum spectrum() const { return spectrum_; }
        bool rsd() const { return rsd_; }

        int n_loops() const { return n_loops_; }
        std::size_t n_coeffs() const { return n_coeffs_; }
        std::size_t n_configs() const { return n_configs_; }
        std::size_t n_kernels() const { return n_kernels_; }
        std::size_t n_kernel_args() const { return n_kernel_args_; }
        int zero_label() const { return zero_label_; }

        std::function<int(const int[])> args_2_kernel_index;
};


class EvolutionParameters {
    private:
        double f_nu_       = 0.0;
        double cs2_factor  = 0.0;

        double ode_atol_   = 0.0;
        double ode_rtol_   = 0.0;
        double ode_hstart_ = 0.0;

        Interpolation1D zeta;
        Interpolation1D redshift;
        Interpolation1D omega_eigenvalues;
        Vec1D<Interpolation1D> F1_ic;
        Interpolation2D effcs2;
    public:
        EvolutionParameters() = default;
        EvolutionParameters(const EvolutionParameters&) = delete;
        EvolutionParameters& operator=(const EvolutionParameters&) = delete;
        EvolutionParameters(EvolutionParameters&& other) noexcept;
        EvolutionParameters& operator=(EvolutionParameters&& other);

        /* Effective sound speed constructors */
        EvolutionParameters(
                double f_nu,
                double omega_m_0,
                const std::string& zeta_file,
                const std::string& redshift_file,
                const std::string& omega_eigenvalues_file,
                const Vec1D<std::string>& F1_ic_files,
                const std::string& effcs2_x_file,
                const std::string& effcs2_y_file,
                const std::string& effcs2_data_file,
                double ode_atol = 1e-6,
                double ode_rtol = 1e-4,
                double ode_hstart = 1e-3
                );

        /* No sound speed constructors */
        EvolutionParameters(
                const std::string& zeta_file,
                double ode_atol = 1e-6,
                double ode_rtol = 1e-4,
                double ode_hstart = 1e-3
                );

        double f_nu() const {return f_nu_;}
        double ode_atol() const {return ode_atol_;}
        double ode_rtol() const {return ode_rtol_;}
        double ode_hstart() const {return ode_hstart_;}

        double zeta_at_eta(double eta) const {return zeta(eta);}
        double omega_eigenvalues_at_k(double k) const {
            return omega_eigenvalues(k);
        }

        double F1_ic_at_k(std::size_t i, double k) const {
            return F1_ic[i](k);
        }

        double cs2(double eta, double k) const {
            return cs2_factor * effcs2(eta, k) / (1 + redshift(eta));
        }
};


#endif /* ifndef PARAMETERS_HPP */
