/*
   parameters.cpp

   Created by Petter Taule on 04.10.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>

#include "../include/utilities.hpp"
#include "../include/parameters.hpp"

using std::size_t;


LoopParameters::LoopParameters(int n_loops, Spectrum spectrum, Dynamics dynamics)
    : dynamics(dynamics), spectrum(spectrum), n_loops(n_loops)
{
    if (spectrum == POWERSPECTRUM && (n_loops < 1 || n_loops > 2)) {
        throw(std::invalid_argument(
            "LoopParameters::LoopParameters(): POWERSPECTRUM only "
            "implemented for n_loops = 1,2."));
    }
    if (spectrum == BISPECTRUM && (n_loops != 1)) {
        throw(std::invalid_argument(
            "LoopParameters::LoopParameters(): BISPECTRUM only "
            "implemented for n_loops = 1."));
    }
    if (spectrum == POWERSPECTRUM) {
        n_coeffs = n_loops + 1;
        n_configs = pow(3, n_coeffs);
        n_kernels = (n_configs/3 + 1) * pow(4,n_loops);
        n_kernel_args = 2 * n_loops + 1;
        zero_label = ::get_zero_label(n_coeffs);

        /* Define block size for kernel indexing. A block consists of all */
        /* single loop label argument combinations. */
        single_loop_block_size = pow(4, n_loops);
        /* Max/min single_loops labels */
        single_loop_label_min = n_configs / 3;
        single_loop_label_max = 2 * n_configs / 3 - 1;
    }
    else if (spectrum == BISPECTRUM) {
        n_coeffs = n_loops + 2;
        n_configs = pow(3, n_coeffs);
        n_kernels = pow(n_configs - pow(3, n_loops), 2) * pow(4, n_loops);
        n_kernel_args = 2 * n_loops + 2;
        zero_label = ::get_zero_label(n_coeffs);

        single_loop_block_size = pow(4, n_loops);
        /* Max/min single_loops labels */
        single_loop_label_min = 4 * n_configs / 9;
        single_loop_label_max = 5 * n_configs / 9 - 1;
        /* Composite label has either k_a, k_b, or multiple loop momenta
         * (applicable for overall loop diagram) */
        first_composite_block_size = 26 * n_configs/27 * single_loop_block_size;
    }
    else {
        throw(std::invalid_argument(
            "LoopParameters::LoopParameters(): invalid spectrum."));
    }

    /* List of single loop labels */
    int coeffs[N_COEFFS_MAX] = {zero_label};
    int label;

    for (int i = 0; i < n_loops; ++i) {
        coeffs[i] = -1;
        label = config2label(coeffs, n_coeffs);
        single_loops.push_back(label);

        coeffs[i] = 1;
        label = config2label(coeffs, n_coeffs);
        single_loops.push_back(label);

        coeffs[i] = 0;
    }
}



int LoopParameters::ps_arguments_2_kernel_index(const int arguments[]) const
{
   /* Precompute powers of two for speedup */
    int pow2[] = {1,2,4,8,16,32,64,128};

    // In DEBUG-mode, check that non-zero arguments (zero_label) are unique
#if DEBUG >= 1
    if (!unique_elements(arguments, n_kernel_args, zero_label))
        throw(std::logic_error(
            "LoopParameters::ps_arguments_2_kernel_index(): duplicate "
            "vector arguments passed."));
    int n_k_labels = 0;
#endif

    int index = 0;

    for (int i = 0; i < n_kernel_args; ++i) {
        // First, check if argument is a zero vector
        if (arguments[i] == zero_label) continue;

        // Argument is a k-type vector (i.e. on the form k + c_i Q_i) if k is
        // present. In our vector-label convention, k is the last coefficient,
        // hence +k is present if label > single_loop_label_max
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
            if(!single_loop_label(arguments[i], n_coeffs, spectrum))
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

    // In DEBUG-mode, check that non-zero arguments (zero_label) are unique
#if DEBUG >= 1
    if (!unique_elements(arguments, n_kernel_args, zero_label))
        throw(std::logic_error(
            "LoopParameters::bs_arguments_2_kernel_index(): duplicate "
            "vector arguments passed."));
#endif

    int index = 0;

    /* Counter of composite arguments, i.e. not single loop momenta */
    int n_composite = 0;

    for (int i = 0; i < n_kernel_args; ++i) {
        // First, check if argument is a zero vector
        if (arguments[i] == zero_label) continue;

        // Argument is not single loop if label < single_loop_label_min or
        // label > single_loop_label_max */
        if (arguments[i] < single_loop_label_min) {
            index += (arguments[i] + 1) *
                (n_composite > 0 ? single_loop_block_size :
                 first_composite_block_size);

            ++n_composite;
        }
        else if (arguments[i] > single_loop_label_max) {
            int block = n_composite > 0 ? single_loop_block_size : first_composite_block_size;
            index += (arguments[i] - n_configs/9 + 1) * block;

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
            /* Went through loop, which means that argument is composite, hence
             * connecting line with overall loop */
            index += (arguments[i] - single_loop_label_min + 1) *
                (n_composite > 0 ? single_loop_block_size :
                 first_composite_block_size);
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
    f_nu(f_nu), ode_atol(ode_atol), ode_rtol(ode_rtol), ode_hstart(ode_hstart),
    zeta(zeta_file), redshift(redshift_file),
    omega_eigenvalues(omega_eigenvalues_file), effcs2(effcs2_x_file,
            effcs2_y_file, effcs2_data_file)
{
    for (auto& F1_ic_file : F1_ic_files) {
        F1_ic.emplace_back(F1_ic_file);
    }
    /* cs2_factor = 2/3 * 1/(omegaM a^2 H^2) */
    cs2_factor = 2.0/3.0 * SQUARE(3e3) / omega_m_0;

    using namespace std::placeholders;
    cs2 = std::bind(&EvolutionParameters::eff_sound_speed, this, _1, _2);
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
    f_nu(f_nu), ode_atol(ode_atol), ode_rtol(ode_rtol), ode_hstart(ode_hstart),
    zeta(zeta_file), redshift(redshift_file),
    omega_eigenvalues(omega_eigenvalues_file)
{
    for (auto& F1_ic_file : F1_ic_files) {
        F1_ic.emplace_back(F1_ic_file);
    }
    /* cs2_factor = 2/3 * 1/(omegaM a^2 H^2) */
    cg2_factor = 2.0/3.0 * SQUARE(3e3) / omega_m_0;

    /* Neutrino temperature today [eV] */
    double T_nu0 = 1.67734976e-4;
    /* 5/9 Zeta(5)/Zeta(3) = 7.188565369 */
    double a = 7.188565369;
    cg2_factor *= a * SQUARE(T_nu0) / SQUARE(m_nu);

    using namespace std::placeholders;
    cs2 = std::bind(&EvolutionParameters::ad_sound_speed, this, _1, _2);
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
    ode_atol(ode_atol), ode_rtol(ode_rtol), ode_hstart(ode_hstart), zeta(zeta_file)
{}



EvolutionParameters::EvolutionParameters(const std::string& zeta_file)
    : EvolutionParameters(1e-6, 1e-4, 1e-3, zeta_file)
{}
