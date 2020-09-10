/*
   diagrams.cpp

   Created by Petter Taule on 30.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
   */

#include <iostream>
#include <numeric>
#include <cmath>
#include <stdexcept>

#include <gsl/gsl_sf.h>

#include "../include/utilities.hpp"
#include "../include/combinatorics.hpp"
#include "../include/tables.hpp"
#include "../include/diagrams.hpp"


void PowerSpectrumDiagram::kernel_arguments(
        short int n_coeffs,
        short int a,
        short int b
        )
{
    Vec1D<short int>& rearrangement = rearrangements[a];
    Vec1D<bool>& signs = sign_configs[b];

    Vec1D<short int>& arguments_l = arg_configs_l[a][b].args;
    Vec1D<short int>& arguments_r = arg_configs_r[a][b].args;

    /* First argument : k1 - (±k2) - (±k3) - ... */
    short int config[N_COEFFS_MAX] = {0};
    config[n_coeffs - 1] = 1;

    for (short int i = 2; i <= m; ++i) {
        // Note extra minus sign from kernel expression
        config[rearrangement[i-2]] = (signs[i-2] ? - 1 : 1);
    }

    arguments_l[0] = config2label(config, n_coeffs);
    arguments_r[0] = config2label(config, n_coeffs);

    // Reset config
    std::fill(config, config + n_coeffs, 0);

    // Argument indices
    size_t index_l = 1;
    size_t index_r = 1;

    // k2,k3,...,km arguments, the ordering of which is stored by the first
    // (m-1) entries of rearrangement[]. Note that the signs are opposite of
    // corresponding terms in the first argument
    for (short int i = 2; i <= m; ++i) {
        config[rearrangement[i-2]] = (signs[i-2] ? 1 : -1);
        arguments_l[index_l++] = config2label(config, n_coeffs);
        arguments_r[index_r++] = config2label(config, n_coeffs);
        std::fill(config, config + n_coeffs, 0);
    }

    // l-loop arguments
    for (short int i = 0; i < l; ++i) {
        short int loop_momentum_index = rearrangement[i + m - 1];
        config[loop_momentum_index] = 1;
        arguments_l[index_l++] = config2label(config, n_coeffs);
        config[loop_momentum_index] = -1;
        arguments_l[index_l++] = config2label(config, n_coeffs);

        std::fill(config, config + n_coeffs, 0);
    }
    // r-loop arguments
    for (short int i = 0; i < r; ++i) {
        short int loop_momentum_index = rearrangement[i + m - 1 + l];
        config[loop_momentum_index] = 1;
        arguments_r[index_r++] = config2label(config, n_coeffs);
        config[loop_momentum_index] = -1;
        arguments_r[index_r++] = config2label(config, n_coeffs);

        std::fill(config, config + n_coeffs, 0);
    }

    // Fill remaining spots with zero-label
    short int zero_label = get_zero_label(n_coeffs);
    size_t n_kernel_args = static_cast<size_t>(settings.n_kernel_args);
    while (index_l < n_kernel_args) {
        arguments_l[index_l++] = zero_label;
    }
    while (index_r < n_kernel_args) {
        arguments_r[index_r++] = zero_label;
    }

#if DEBUG >= 1
    /* kernel_index_from_arguments() assumes that the length of arguments[] is
     * n_kernel_args. Checking this explicitly: */
    if (arguments_l.size() != n_kernel_args) {
        throw(std::logic_error("PowerSpectrumDiagram::kernel_arguments(): Size of left argument vector does not equal n_kernel_args."));
    }
    if (arguments_r.size() != n_kernel_args) {
        throw(std::logic_error("PowerSpectrumDiagram::kernel_arguments(): Size of right argument vector does not equal n_kernel_args."));
    }
#endif

    arg_configs_l[a][b].kernel_index =
        kernel_index_from_arguments(arguments_l.data(), settings);
    arg_configs_r[a][b].kernel_index =
        kernel_index_from_arguments(arguments_r.data(), settings);
}



void PowerSpectrumDiagram::compute_rearrangements(short int n_loops)
{
    // Allocate memory for rearrangements
    rearrangements.resize(n_rearrangements);
    for (short int i = 0; i < n_rearrangements; ++i) {
        rearrangements[i].resize(n_loops);
    }

    Vec1D<short int> group_sizes;

    /* Connecting loops? */
    if (m > 1) {
        group_sizes.push_back(m-1);
    }
    if (l > 0) {
        group_sizes.push_back(l);
    }
    if (r > 0) {
        group_sizes.push_back(r);
    }

    Orderings orderings(n_loops, group_sizes);

    int counter = 0;
    do {
        rearrangements.at(counter) = orderings.get_current();
        ++counter;
    } while (orderings.next());

    if (counter != n_rearrangements) {
        throw(std::logic_error("PowerSpectrumDiagram::compute_rearrangements(): Created "
                    + std::to_string(counter) +
                    " rearrangements, which does not equal n_rearrangements = "
                    + std::to_string(n_rearrangements) + "."));
    }
}



void PowerSpectrumDiagram::compute_sign_flips() {
    // Allocate memory for sign configurations
    sign_configs.resize(n_sign_configs);
    for (short int i = 0; i < n_sign_configs; ++i) {
        // Signs of k2,...,km can be flipped, hence allocate (m-1)-length array
        sign_configs[i].resize(m-1);
    }

    // Loop over possible sign configurations and store in sign table
    for (short int i = 0; i < n_sign_configs; ++i) {
        for (short int j = 0; j < (m-1); ++j) {
            // Add 1 so that the first sign configuration is +1,+1,...,+1
            sign_configs[i][j] = (i/(int)pow(2,j) + 1) % 2;
        }
    }
}



PowerSpectrumDiagram::PowerSpectrumDiagram(
        const Settings& settings,
        short int m,
        short int l,
        short int r
        ) : settings(settings), m(m), l(l), r(r)
{
    short int n_loops = settings.n_loops;
    short int n_kernel_args = settings.n_kernel_args;

    diagram_factor = (gsl_sf_fact(2*l + m) * gsl_sf_fact(2*r + m)) /
        (pow(2,l+r) * gsl_sf_fact(l) * gsl_sf_fact(r) * gsl_sf_fact(m));

    n_rearrangements = gsl_sf_fact(n_loops) /
        (gsl_sf_fact(m-1) * gsl_sf_fact(l) * gsl_sf_fact(r));
    n_sign_configs = pow(2,m-1);

    compute_rearrangements(n_loops);
    compute_sign_flips();

    // Allocate memory
    arg_configs_l.resize(n_rearrangements);
    arg_configs_r.resize(n_rearrangements);

    for (short int a = 0; a < n_rearrangements; ++a) {
        arg_configs_l[a].resize(n_sign_configs);
        arg_configs_r[a].resize(n_sign_configs);

        for (short int b = 0; b < n_sign_configs; ++b) {
            arg_configs_l[a][b].args.resize(n_kernel_args);
            arg_configs_r[a][b].args.resize(n_kernel_args);

            // Initialize arguments and kernel indices for this configuration
            kernel_arguments(settings.n_coeffs, a, b);
        }
    }
}



void PowerSpectrumDiagram::print_diagram_tags(std::ostream& out) const
{
    out << COLOR_MAGENTA 
        << "(m,l,r) = (" << m << "," << l << "," << r << ")" 
        << COLOR_RESET;
}



void PowerSpectrumDiagram::print_argument_configuration(
        std::ostream& out,
        short int a,
        short int b
        ) const
{
    out << COLOR_BLUE;
    print_labels(arg_configs_l[a][b].args.data(),
            arg_configs_l[a][b].args.size(), settings.n_coeffs,
            settings.zero_label, out);
    out << " ";
    print_labels(arg_configs_r[a][b].args.data(),
            arg_configs_r[a][b].args.size(), settings.n_coeffs,
            settings.zero_label, out);
    out << COLOR_RESET;
}



std::ostream& operator<<(std::ostream& out, const PowerSpectrumDiagram& diagram) {
    diagram.print_diagram_tags(out);
    out << std::endl;

    for (short int a = 0; a < diagram.n_rearrangements; ++a) {
        for (short int b = 0; b < diagram.n_sign_configs; ++b) {
            out << "\t";
            diagram.print_argument_configuration(out, a, b);
            out << std::endl;
        }
    }
    return out;
}



Vec1D<PowerSpectrumDiagram> construct_diagrams(const Settings& settings) {
    Vec1D<PowerSpectrumDiagram> diagrams;

    short int n_loops = settings.n_loops;
    short int m = 0;
    short int index = 0;

    // Find (distinct) power spectrum diagrams at L-loop
    // They satisfy: m >= 1; l,r > 0; l + r + m = L + 1
    for (m = 1; m <= n_loops + 1; ++m) {
        short int l = n_loops + 1 - m;
        short int r = 0;
        while (l >= r) {
            if (index >= 2 * n_loops) {
                throw(std::logic_error("construct_diagrams(): Index larger than 2 * n_loops."));
            }
            diagrams.push_back(PowerSpectrumDiagram(settings, m, l, r));

            l = n_loops + 1 - m - (++r);
            index++;
        };
    }
    return diagrams;
}
