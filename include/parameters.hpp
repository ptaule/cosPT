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
        const Dynamics dynamics;
        const Spectrum spectrum;

        const int n_loops;

        int n_coeffs;
        int n_configs;
        int n_kernels;
        int n_kernel_args;
        int zero_label;
        int single_loop_label_min;
        int single_loop_label_max;
        int single_loop_block_size;

        Vec1D<int> single_loops;

        int first_composite_block_size = 0; /* Bispectrum */

        int ps_arguments_2_kernel_index(const int arguments[]) const;
        int bs_arguments_2_kernel_index(const int arguments[]) const;

    public:
        LoopParameters(int n_loops, Spectrum spectrum, Dynamics dynamics);

        Dynamics get_dynamics() const { return dynamics; }
        Spectrum get_spectrum() const { return spectrum; }

        int get_n_loops() const { return n_loops; }
        int get_n_coeffs() const { return n_coeffs; }
        int get_n_configs() const { return n_configs; }
        int get_n_kernels() const { return n_kernels; }
        int get_n_kernel_args() const { return n_kernel_args; }
        int get_zero_label() const { return zero_label; }

        int arguments_2_kernel_index(const int arguments[]) const {
            if (spectrum == POWERSPECTRUM)
                return ps_arguments_2_kernel_index(arguments);
            else if (spectrum == BISPECTRUM)
                return bs_arguments_2_kernel_index(arguments);
            else
                throw(std::logic_error(
                            "Parameters::arguments_2_kernel_index(): invalid spectrum."));
        }
};



class EvolutionParameters {
    private:
        double f_nu       = 0.0;
        double cs2_factor = 0.0;
        double cg2_factor = 0.0;

        double ode_atol   = 0.0;
        double ode_rtol   = 0.0;
        double ode_hstart = 0.0;

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

        double get_f_nu() const {return f_nu;}
        double get_ode_atol() const {return ode_atol;}
        double get_ode_rtol() const {return ode_rtol;}
        double get_ode_hstart() const {return ode_hstart;}

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
