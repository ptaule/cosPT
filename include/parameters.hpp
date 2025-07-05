#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <any>
#include <functional>
#include <stdexcept>
#include <string>
#include <map>

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

        Vec1D<double> kappa_;
        Vec1D<std::string> zeta_files_;
        Vec1D<std::string> xi_files_;

        Vec1D<double> bias_parameters_;

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
        bool set_output_file(const libconfig::Config& cfg);

        void set_dynamics(const libconfig::Config& cfg);

        void set_cuba_config(const libconfig::Setting& cuba_settings);
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
                int k_c_idx = -1
                );

        Vec1D<Pair<int>> pair_correlations() const {return pair_correlations_;}
        Vec1D<Triple<int>> triple_correlations() const {return triple_correlations_;}

        Vec1D<double> kappa() const {return kappa_;}
        Vec1D<std::string> zeta_files() const {return zeta_files_;}
        Vec1D<std::string> xi_files() const {return xi_files_;}
        Vec1D<double> bias_parameters() const {return bias_parameters_;}

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
        void set(const std::string& key, T value) {
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
        bool exists(const std::string& key) const {
            return params.find(key) != params.end();
        }
};


std::ostream& operator<<(std::ostream& out, const Config& cfg);


// Represents the structural and indexing logic for n-loop computations.
// Provides configuration metadata and converts argument lists into
// hashed kernel indices (args_to_kernel_index).
class LoopStructure {
    private:
        const Spectrum spectrum_;

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

        int ps_args_to_kernel_index(const int arguments[]) const;
        int bs_args_to_kernel_index(const int arguments[]) const;

    public:
        LoopStructure(int n_loops, Spectrum spectrum);

        Spectrum spectrum() const { return spectrum_; }

        int n_loops() const { return n_loops_; }
        std::size_t n_coeffs() const { return n_coeffs_; }
        std::size_t n_configs() const { return n_configs_; }
        std::size_t n_kernels() const { return n_kernels_; }
        std::size_t n_kernel_args() const { return n_kernel_args_; }
        int zero_label() const { return zero_label_; }

        std::function<int(const int[])> args_to_kernel_index;
};


class EvolutionParameters {
    private:
        Vec1D<double> kappa_ = {0.0};
        Vec1D<Interpolation1D> zeta_; /* Time-dependent functions */
        Vec1D<Interpolation2D> xi_;   /* Space- and time-dependent functions */

        double ode_atol_   = 0.0;
        double ode_rtol_   = 0.0;
        double ode_hstart_ = 0.0;
    public:
        EvolutionParameters() = default;
        EvolutionParameters(const EvolutionParameters&) = delete;
        EvolutionParameters& operator=(const EvolutionParameters&) = delete;
        EvolutionParameters(EvolutionParameters&& other) noexcept;
        EvolutionParameters& operator=(EvolutionParameters&& other);

        EvolutionParameters(
                const Vec1D<double>& kappa,
                const Vec1D<std::string>& zeta_files,
                const Vec1D<std::string>& xi_files,
                double ode_atol = 1e-6,
                double ode_rtol = 1e-4,
                double ode_hstart = 1e-3
                );

        double ode_atol() const { return ode_atol_; }
        double ode_rtol() const { return ode_rtol_; }
        double ode_hstart() const { return ode_hstart_; }

        const Vec1D<double>& kappa() const { return kappa_; }
        const Vec1D<Interpolation1D>& zeta() const { return zeta_; }
        const Vec1D<Interpolation2D>& xi() const { return xi_; }
};


class OmegaEigenspace {
    private:
        Interpolation1D eigenvalue_;
        Vec1D<Interpolation1D> eigenvectors_;

        double eta_ini = 0;           /* Time eta at which to compute eigenspace */
        int eigenmode = 0;            /* Which eigenmode to compute [0,COMPONENTS-1] */
        double imag_threshold = 1e-3; /* Threshold for imaginary part of eigenvalues/eigenvectors */
        const Dynamics dynamics;
        const EvolutionParameters& ev_params;

        void omega_eigenspace_at_k(
            double k,
            double& eigenvalue,        /* out */
            Vec1D<double>& eigenvector /* out */
            );
    public:
        OmegaEigenspace() = delete;
        OmegaEigenspace(const OmegaEigenspace&) = delete;
        OmegaEigenspace& operator=(const OmegaEigenspace&) = delete;
        OmegaEigenspace(OmegaEigenspace&& other) noexcept;
        OmegaEigenspace& operator=(OmegaEigenspace&& other) = delete;

        OmegaEigenspace(
                const Dynamics dynamics,
                double eta_ini,                       /* Time eta at which to compute eigenspace */
                const EvolutionParameters& ev_params,
                int eigenmode,                        /* Which eigenmode to compute [0,COMPONENTS-1] */
                double k_min,
                double k_max,                         /* Over which k-range to compute */
                int N,                                /* How many grid points to use in interpolation */
                double imag_threshold = 1e-3          /* Threshold for imaginary part of eigenvalues/eigenvectors */
                );

        const Interpolation1D& eigenvalue() const {
#if DEBUG > 0
            if (dynamics != EVOLVE_ASYMPTOTIC_ICS) {
                throw std::runtime_error("OmegaEigenspace::eigenvalue(): "
                        "Only for dynamics = EVOLVE_ASYMPTOTIC_ICS is Omega "
                        "eigenspace computed.");
            }
#endif
            return eigenvalue_;
        };

        const Vec1D<Interpolation1D>& eigenvectors() const {
#if DEBUG > 0
            if (dynamics != EVOLVE_ASYMPTOTIC_ICS) {
                throw std::runtime_error("OmegaEigenspace::eigenvalue(): "
                        "Only for dynamics = EVOLVE_ASYMPTOTIC_ICS is Omega "
                        "eigenspace computed.");
            }
#endif
            return eigenvectors_;
        };
};


#endif /* ifndef PARAMETERS_HPP */
