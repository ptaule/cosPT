/*
   diagrams.cpp

   Created by Petter Taule on 30.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
   */

#include <iostream>
#include <algorithm>
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
    Vec1D<short int>& rearrangement = rearrangements.at(a);
    Vec1D<bool>& signs = sign_configs.at(b);

    Vec1D<short int>& arguments_l = arg_configs_l.at(a).at(b).args;
    Vec1D<short int>& arguments_r = arg_configs_r.at(a).at(b).args;

    /* First argument : k1 - (±k2) - (±k3) - ... */
    short int config[N_COEFFS_MAX] = {0};
    config[n_coeffs - 1] = 1;

    for (short int i = 2; i <= m; ++i) {
        // Note extra minus sign from kernel expression
        config[rearrangement.at(i-2)] = (signs.at(i-2) ? - 1 : 1);
    }

    arguments_l.at(0) = config2label(config, n_coeffs);
    arguments_r.at(0) = config2label(config, n_coeffs);

    // Reset config
    std::fill(config, config + n_coeffs, 0);

    // Argument indices
    size_t index_l = 1;
    size_t index_r = 1;

    // k2,k3,...,km arguments, the ordering of which is stored by the first
    // (m-1) entries of rearrangement[]. Note that the signs are opposite of
    // corresponding terms in the first argument
    for (short int i = 2; i <= m; ++i) {
        config[rearrangement.at(i-2)] = (signs.at(i-2) ? 1 : -1);
        arguments_l.at(index_l++) = config2label(config, n_coeffs);
        arguments_r.at(index_r++) = config2label(config, n_coeffs);
        std::fill(config, config + n_coeffs, 0);
    }

    // l-loop arguments
    for (short int i = 0; i < l; ++i) {
        short int loop_momentum_index = rearrangement.at(i + m - 1);
        config[loop_momentum_index] = 1;
        arguments_l.at(index_l++) = config2label(config, n_coeffs);
        config[loop_momentum_index] = -1;
        arguments_l.at(index_l++) = config2label(config, n_coeffs);

        std::fill(config, config + n_coeffs, 0);
    }
    // r-loop arguments
    for (short int i = 0; i < r; ++i) {
        short int loop_momentum_index = rearrangement.at(i + m - 1 + l);
        config[loop_momentum_index] = 1;
        arguments_r.at(index_r++) = config2label(config, n_coeffs);
        config[loop_momentum_index] = -1;
        arguments_r.at(index_r++) = config2label(config, n_coeffs);

        std::fill(config, config + n_coeffs, 0);
    }

    // Fill remaining spots with zero-label
    short int zero_label = get_zero_label(n_coeffs);
    size_t n_kernel_args = static_cast<size_t>(settings.n_kernel_args);
    while (index_l < n_kernel_args) {
        arguments_l.at(index_l++) = zero_label;
    }
    while (index_r < n_kernel_args) {
        arguments_r.at(index_r++) = zero_label;
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

    arg_configs_l.at(a).at(b).kernel_index =
        kernel_index_from_arguments(arguments_l.data(), settings);
    arg_configs_r.at(a).at(b).kernel_index =
        kernel_index_from_arguments(arguments_r.data(), settings);
}



void PowerSpectrumDiagram::compute_rearrangements(short int n_loops)
{
    // Allocate memory for rearrangements
    rearrangements.resize(n_rearrangements);
    for (short int i = 0; i < n_rearrangements; ++i) {
        rearrangements.at(i).resize(n_loops);
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
        sign_configs.at(i).resize(m-1);
    }

    // Loop over possible sign configurations and store in sign table
    for (short int i = 0; i < n_sign_configs; ++i) {
        for (short int j = 0; j < (m-1); ++j) {
            // Add 1 so that the first sign configuration is +1,+1,...,+1
            sign_configs.at(i).at(j) = (i/(int)pow(2,j) + 1) % 2;
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

    if (m + l + r != n_loops + 1) {
        throw(std::invalid_argument("m + l + r != n_loops + 1"));
    }

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
        arg_configs_l.at(a).resize(n_sign_configs);
        arg_configs_r.at(a).resize(n_sign_configs);

        for (short int b = 0; b < n_sign_configs; ++b) {
            arg_configs_l.at(a).at(b).args.resize(n_kernel_args);
            arg_configs_r.at(a).at(b).args.resize(n_kernel_args);

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
    print_labels(arg_configs_l.at(a).at(b).args.data(),
            arg_configs_l.at(a).at(b).args.size(), settings.n_coeffs,
            settings.zero_label, settings.spectrum, out);
    out << " ";
    print_labels(arg_configs_r.at(a).at(b).args.data(),
            arg_configs_r.at(a).at(b).args.size(), settings.n_coeffs,
            settings.zero_label, settings.spectrum, out);
    out << COLOR_RESET;
}



void BiSpectrumDiagram::compute_rearrangements(short int n_loops) {
    // Allocate memory for rearrangements
    rearrangements.resize(n_rearrangements);
    for (short int i = 0; i < n_rearrangements; ++i) {
        rearrangements.at(i).resize(n_loops);
    }

    Vec1D<short int> group_sizes;

    /* Overall loops? */
    if (overall_loop) {
        group_sizes.push_back(1);
    }
    /* Connecting loops? */
    if (n_ab > 1) {
        group_sizes.push_back(n_ab - 1);
    }
    if (n_bc > 1) {
        group_sizes.push_back(n_bc - 1);
    }
    if (n_ca > 1) {
        group_sizes.push_back(n_ca - 1);
    }
    /* Self loops? */
    if (n_a > 0) {
        group_sizes.push_back(n_a);
    }
    if (n_b > 0) {
        group_sizes.push_back(n_b);
    }
    if (n_c > 0) {
        group_sizes.push_back(n_c);
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



void BiSpectrumDiagram::compute_sign_flips(
        short int connecting_lines,
        short int n_sign_configs,
        Vec2D<bool> sign_configs
        )
{
    if (connecting_lines == 0) {
        return;
    }
    /* Allocate memory for sign configurations */
    sign_configs.resize(n_sign_configs);
    for (short int i = 0; i < n_sign_configs; ++i) {
        sign_configs.at(i).resize(connecting_lines - 1);
    }

    /* Loop over possible sign configurations and store */
    for (short int i = 0; i < n_sign_configs; ++i) {
        for (short int j = 0; j < (connecting_lines - 1); ++j) {
            // Add 1 so that the first sign configuration is +1,+1,...,+1
            sign_configs.at(i).at(j) = (i/(int)pow(2,j) + 1) % 2;
        }
    }
}



BiSpectrumDiagram::BiSpectrumDiagram(
                const Settings& settings,
                short int n_ab,
                short int n_bc,
                short int n_ca,
                short int n_a,
                short int n_b,
                short int n_c
        ) : settings(settings),
    n_ab(n_ab), n_bc(n_bc), n_ca(n_ca),
    n_a(n_a), n_b(n_b), n_c(n_c)
{
    short int n_loops = settings.n_loops;
    short int n_kernel_args = settings.n_kernel_args;

    if (n_ab + n_bc + n_ca + n_a + n_b + n_c != n_loops + 2) {
        throw(std::invalid_argument(
                    "n_ab + n_bc + n_ca + n_a + n_b + n_c != n_loops + 2"));
    }
    if ((n_ab > 0) + (n_bc > 0) + (n_ca > 0) < 2) {
        throw(std::invalid_argument("More than one number out of {n_ab, n_bc, n_ca} is zero."));
    }
    overall_loop = true;
    if (n_ab == 0 || n_bc == 0 || n_ca == 0) {
        overall_loop = false;
    }

    diagram_factor =
        gsl_sf_fact(2 * n_a + n_ab + n_ca) *
        gsl_sf_fact(2 * n_b + n_bc + n_ab) *
        gsl_sf_fact(2 * n_c + n_ca + n_bc) /
        (
         pow(2,n_a + n_b + n_c) *
         gsl_sf_fact(n_a) * gsl_sf_fact(n_b) * gsl_sf_fact(n_c) *
         gsl_sf_fact(n_ab) * gsl_sf_fact(n_bc) * gsl_sf_fact(n_ca)
         );

    n_rearrangements = gsl_sf_fact(n_loops) /
        (
         gsl_sf_fact(n_ab > 0 ? n_ab - 1 : 0) *
         gsl_sf_fact(n_bc > 0 ? n_bc - 1 : 0) *
         gsl_sf_fact(n_ca > 0 ? n_ca - 1 : 0) *
         gsl_sf_fact(n_a) * gsl_sf_fact(n_b) * gsl_sf_fact(n_c)
        );

    n_sign_configs_ab = pow(2, n_ab > 0 ? n_ab - 1 : 0);
    n_sign_configs_bc = pow(2, n_bc > 0 ? n_bc - 1 : 0);
    n_sign_configs_ca = pow(2, n_ca > 0 ? n_ca - 1 : 0);

    compute_rearrangements(n_loops);
    compute_sign_flips(n_ab, n_sign_configs_ab, sign_configs_ab);
    compute_sign_flips(n_bc, n_sign_configs_bc, sign_configs_bc);
    compute_sign_flips(n_ab, n_sign_configs_ab, sign_configs_ab);

    // Allocate memory
    arg_configs_a.resize(n_rearrangements);
    arg_configs_b.resize(n_rearrangements);
    arg_configs_c.resize(n_rearrangements);

    for (short int a = 0; a < n_rearrangements; ++a) {
        arg_configs_a.at(a).resize(n_sign_configs_ab * n_sign_configs_ca);
        arg_configs_b.at(a).resize(n_sign_configs_bc * n_sign_configs_ab);
        arg_configs_c.at(a).resize(n_sign_configs_ca * n_sign_configs_bc);

        /* for (short int b = 0; b < n_sign_configs; ++b) { */
        /*     arg_configs_a[a][b].args.resize(n_kernel_args); */
        /*     arg_configs_b[a][b].args.resize(n_kernel_args); */
        /*     arg_configs_c[a][b].args.resize(n_kernel_args); */

        /*     // Initialize arguments and kernel indices for this configuration */
        /*     kernel_arguments(n_loops, a, b); */
        /* } */
    }
}



void BiSpectrumDiagram::print_diagram_tags(std::ostream& out) const
{
    out << COLOR_MAGENTA
        << "(n_a, n_b, n_c, n_ab, n_bc, n_ca) = ("
        << n_a << "," << n_b << "," << n_c << "," << n_ab << "," << n_bc << ","
        << n_ca << ")" << COLOR_RESET;
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



Vec1D<PowerSpectrumDiagram> construct_ps_diagrams(const Settings& settings) {
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



Vec1D<BiSpectrumDiagram> construct_bs_diagrams(const Settings& settings) {
    Vec1D<BiSpectrumDiagram> diagrams;

    short int n_loops = settings.n_loops;

    /* Need to place L + 2 lines */
    /* First, consider diagram where a connecting line is zero, e.g. n_ab = 0 */
    for (short int n_bc = 1; n_bc <= n_loops + 1; ++n_bc) {
        for (short int n_ca = 1; n_ca <= n_loops - n_bc + 2; ++n_ca) {
            /* Any remaining loops are self-loops */
            short int n_self = n_loops + 2 - n_bc - n_ca;
            /* Assining loops to n_a */
            for (short int n_a = 0; n_a <= n_self; ++n_a) {
                /* Split remaining loops between n_b and n_c */
                short int n_b = n_self - n_a;
                short int n_c = 0;
                while (n_b >= n_c) {
                    /* Three diagrams, picking n_ab, n_bc, or n_ca as line being zero */
                    diagrams.push_back(BiSpectrumDiagram(settings, 0,    n_bc, n_ca, n_a, n_b, n_c));
                    diagrams.push_back(BiSpectrumDiagram(settings, n_bc, 0,    n_ca, n_a, n_b, n_c));
                    diagrams.push_back(BiSpectrumDiagram(settings, n_ca, n_bc, 0,    n_a, n_b, n_c));
                    /* If n_b != n_c, there are diagrams with n_b <-> n_c */
                    if (n_b != n_c) {
                        diagrams.push_back(BiSpectrumDiagram(settings, 0,    n_bc, n_ca, n_a, n_c, n_b));
                        diagrams.push_back(BiSpectrumDiagram(settings, n_bc, 0,    n_ca, n_a, n_c, n_b));
                        diagrams.push_back(BiSpectrumDiagram(settings, n_ca, n_bc, 0,    n_a, n_c, n_b));
                    }

                    n_b = n_self - n_a - (++n_c);
                }
            }
        }
    }

    /* Diagrams with n_ab, n_bc, n_ca > 0 */
    for (short int n_ab = 1; n_ab <= n_loops + 1; ++n_ab) {
        for (short int n_bc = 1; n_bc <= n_loops - n_ab + 2; ++n_bc) {
            for (short int n_ca = 1; n_ca <= n_loops - n_ab - n_bc + 2; ++n_ca) {
                /* Any remaining loops are self-loops */
                short int n_self = n_loops + 2 - n_ab - n_bc - n_ca;
                /* Assining loops to n_a */
                for (short int n_a = 0; n_a <= n_self; ++n_a) {
                    /* Split remaining loops between n_b and n_c */
                    short int n_b = n_self - n_a;
                    short int n_c = 0;
                    while (n_b >= n_c) {
                        diagrams.push_back(BiSpectrumDiagram(settings, n_ab, n_bc, n_ca, n_a, n_b, n_c));
                        /* If n_b != n_c, there is a diagram with n_b <-> n_c */
                        if (n_b != n_c) {
                            diagrams.push_back(BiSpectrumDiagram(settings, n_ab, n_bc, n_ca, n_a, n_c, n_b));
                        }

                        n_b = n_self - n_a - (++n_c);
                    }
                }
            }
        }
    }

    return diagrams;
}
