/*
   parameters.hpp

   Created by Petter Taule on 04.10.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <functional>
#include <stdexcept>

#include "utilities.hpp"
#include "interpolation.hpp"

class LoopParameters {
    private:
        const Dynamics dynamics_;
        const Spectrum spectrum_;

        const int n_loops_;

        int n_coeffs_;
        int n_configs_;
        int n_kernels_;
        int n_kernel_args_;
        int zero_label_;
        int single_loop_label_min;
        int single_loop_label_max;
        int single_loop_block_size;

        Vec1D<int> single_loops;

        int first_composite_block_size = 0; /* Bispectrum */

        int ps_arguments_2_kernel_index(const int arguments[]) const;
        int bs_arguments_2_kernel_index(const int arguments[]) const;

    public:
        LoopParameters(int n_loops, Spectrum spectrum, Dynamics dynamics);

        Dynamics dynamics() const { return dynamics_; }
        Spectrum spectrum() const { return spectrum_; }

        int n_loops() const { return n_loops_; }
        int n_coeffs() const { return n_coeffs_; }
        int n_configs() const { return n_configs_; }
        int n_kernels() const { return n_kernels_; }
        int n_kernel_args() const { return n_kernel_args_; }
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
};



class EvolutionParameters {
    private:
        double f_nu_       = 0.0;
        double cs2_factor  = 0.0;
        double cg2_factor  = 0.0;

        double ode_atol_   = 0.0;
        double ode_rtol_   = 0.0;
        double ode_hstart_ = 0.0;

        Interpolation1D zeta;
        Interpolation1D redshift;
        Interpolation1D omega_eigenvalues;
        Vec1D<Interpolation1D> F1_ic;
        Interpolation2D effcs2;

        double eff_sound_speed(double eta, double k) const {
            return cs2_factor * effcs2.eval(eta, k) / (1 + redshift.eval(eta));
        }
        double ad_sound_speed(double eta,
                              __attribute__((unused)) double k) const {
            return cg2_factor * (1 + redshift.eval(eta));
        }

    public:
        EvolutionParameters() = default;

        /* Effective sound speed constructors */
        EvolutionParameters(
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
                );
        EvolutionParameters(
                double f_nu,
                double omega_m_0,
                const std::string& zeta_file,
                const std::string& redshift_file,
                const std::string& omega_eigenvalues_file,
                const Vec1D<std::string>& F1_ic_files,
                const std::string& effcs2_x_file,
                const std::string& effcs2_y_file,
                const std::string& effcs2_data_file
                );

        /* Adiabatic sound speed constructors */
        EvolutionParameters(
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
                );
        EvolutionParameters(
                double m_nu,
                double f_nu,
                double omega_m_0,
                const std::string& zeta_file,
                const std::string& redshift_file,
                const std::string& omega_eigenvalues_file,
                const Vec1D<std::string>& F1_ic_files
                );

        /* No sound speed constructors */
        EvolutionParameters(
                double ode_atol,
                double ode_rtol,
                double ode_hstart,
                const std::string& zeta_file
                );
        EvolutionParameters(const std::string& zeta_file);

        double f_nu() const {return f_nu_;}
        double ode_atol() const {return ode_atol_;}
        double ode_rtol() const {return ode_rtol_;}
        double ode_hstart() const {return ode_hstart_;}

        double zeta_at_eta(double eta) const {return zeta.eval(eta);}
        double omega_eigenvalues_at_k(double k) const {
            return omega_eigenvalues.eval(k);
        }

        double F1_ic_at_k(int i, double k) const {
            return F1_ic[i].eval(k);
        }

        std::function<double(double, double)> cs2;
};


#endif /* ifndef PARAMETERS_HPP */
