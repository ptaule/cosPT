/*
   parameters.hpp

   Created by Petter Taule on 04.10.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <stdexcept>
#include <string>

#include "utilities.hpp"
#include "interpolation.hpp"

namespace libconfig {
class Config;
class Setting;
}

class Config {
    private:
        int n_loops_;
        Dynamics dynamics_;
        Spectrum spectrum_;

        bool single_hard_limit_;
        double sh_Q1_ = 0;

        int k_a_idx = -1;
        int k_b_idx = -1;
        int k_c_idx = -1;

        double k_a_;
        double k_b_ = 0;
        double k_c_ = 0;
        double cos_ab_ = 0;

        double q_min_ = 1e-4;
        double q_max_ = 1;

        Vec1D<Pair<int>> pair_correlations_;     /* Power spectrum */
        Vec1D<Triple<int>> triple_correlations_; /* Bispectrum */

        double cuba_atol_           = 1e-12;
        double cuba_rtol_           = 1e-4;
        int cuba_maxevals_          = 0;
        int cuba_verbose_           = 1;
        int cuba_cores_             = 0;
        bool cuba_retain_statefile_ = false;
        std::string cuba_statefile_;
        std::string cuba_statefile_path;

        /* Information from CUBA integration */
        int cuba_evals_      = 0;
        int cuba_fail_       = 0;
        int cuba_subregions_ = 0;

        std::string description_;

        double ode_atol_   = 1e-6;
        double ode_rtol_   = 1e-4;
        double ode_hstart_ = 1e-3;

        std::string input_ps_file_;
        double input_ps_rescale_num = 1;
        std::string input_ps_rescale_str_;

        std::string output_path;
        std::string output_file_;

        std::size_t time_steps_     = 0;
        std::size_t pre_time_steps_ = 0;
        double eta_asymp_   = 0;
        double eta_ini_     = 0;
        double eta_fin_     = 0;

        double f_nu_ = 0;
        double omega_m_0_ = 0;

        std::string zeta_file_;
        std::string redshift_file_;
        std::string omega_eigenvalues_file_;
        Vec1D<std::string> F1_ic_files_;
        std::string effcs2_x_grid_;
        std::string effcs2_y_grid_;
        std::string effcs2_data_;

        void set_spectrum(const libconfig::Config& cfg);
        void set_dynamics(const libconfig::Config& cfg);
        void set_output_file(const libconfig::Config& cfg);
        void set_cuba_statefile(const libconfig::Setting& cuba_settings);
    public:
        Config(const std::string& ini_file,
                int k_a_idx = -1,
                int k_b_idx = -1,
                int k_c_idx = -1,
                int cuba_maxevals = 0,
                int cuba_cores = -1
                );
        int n_loops() const {return n_loops_;}
        Dynamics dynamics() const {return dynamics_;}
        Spectrum spectrum() const {return spectrum_;}
        bool single_hard_limit() const {return single_hard_limit_;}
        double sh_Q1() const {return sh_Q1_;}

        double k_a() const {return k_a_;}
        double k_b() const {return k_b_;}
        double k_c() const {return k_c_;}
        double cos_ab() const {return cos_ab_;};

        double q_min() const {return q_min_;}
        double q_max() const {return q_max_;}

        Vec1D<Pair<int>> pair_correlations() const {return pair_correlations_;}
        Vec1D<Triple<int>> triple_correlations() const {return triple_correlations_;}

        double cuba_atol() const {return cuba_atol_ ;}
        double cuba_rtol() const {return cuba_rtol_;}
        int cuba_maxevals() const {return cuba_maxevals_;}
        int cuba_verbose() const {return cuba_verbose_;}
        int cuba_cores() const {return cuba_cores_;}
        bool cuba_retain_statefile() const {return cuba_retain_statefile_;}
        std::string cuba_statefile() const {return cuba_statefile_;}

        int cuba_subregions() const {return cuba_subregions_;}
        int cuba_evals() const {return cuba_evals_;}
        int cuba_fail() const {return cuba_fail_;}

        /* Functions to write MC integration info */
        int& cuba_subregions() {return cuba_subregions_;}
        int& cuba_evals() {return cuba_evals_;}
        int& cuba_fail() {return cuba_fail_;}

        std::string description() const {return description_;}
        std::string& description() {return description_;}

        double ode_atol() const {return ode_atol_ ;}
        double ode_rtol() const {return ode_rtol_;}
        double ode_hstart() const {return ode_hstart_;}

        std::string input_ps_file() const {return input_ps_file_;}
        double input_ps_rescale() const {return input_ps_rescale_num;}
        std::string input_ps_rescale_str() const {return input_ps_rescale_str_;}

        std::string output_file() const {return output_file_;}

        std::size_t time_steps() const {return time_steps_;}
        std::size_t pre_time_steps() const {return pre_time_steps_;}
        double eta_asymp() const {return eta_asymp_;}
        double eta_ini() const {return eta_ini_;}
        double eta_fin() const {return eta_fin_;}

        double f_nu() const {return f_nu_;}
        double omega_m_0() const {return omega_m_0_;}

        std::string zeta_file() const {return zeta_file_;}
        std::string redshift_file() const {return redshift_file_;}
        std::string omega_eigenvalues_file() const {return omega_eigenvalues_file_;}
        Vec1D<std::string> F1_ic_files() const {return F1_ic_files_;}
        std::string effcs2_x_grid() const {return effcs2_x_grid_;}
        std::string effcs2_y_grid() const {return effcs2_y_grid_;}
        std::string effcs2_data() const {return effcs2_data_;}
};

std::ostream& operator<<(std::ostream& out, const Config& cfg);


class ConfigException : std::exception {
    private:
        std::string ex;
    public:
        ConfigException (const std::string& ex) : ex(ex) {}
        const char* what() const noexcept {return ex.c_str();}
};


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

        int ps_arguments_2_kernel_index(const int arguments[]) const;
        int bs_arguments_2_kernel_index(const int arguments[]) const;

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

        int arguments_2_kernel_index(const int arguments[]) const {
            if (spectrum_ == POWERSPECTRUM)
                return ps_arguments_2_kernel_index(arguments);
            else if (spectrum_ == BISPECTRUM)
                return bs_arguments_2_kernel_index(arguments);
            else
                throw(std::logic_error(
                    "Parameters::arguments_2_kernel_index(): invalid spectrum."));
        }

        /* std::vector argument wrapper function */
        int arguments_2_kernel_index(const std::vector<int> arguments) const
        {
            return arguments_2_kernel_index(arguments.data());
        }
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
