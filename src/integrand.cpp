/*
   integrand.cpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <stdexcept>
#include <cmath>
#include <iostream>

#include "../include/integrand.hpp"

#include "../include/utilities.hpp"
#include "../include/tables.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/diagrams.hpp"
#include "../include/interpolation.hpp"
#include "../include/integrand.hpp"


using std::size_t;

/* // For debuggin purposes */
/* __attribute__((unused)) */
/* void print_evolved_kernel( */
/*         const short int arguments[], */
/*         short int index, */
/*         short int n, */
/*         const tables_t* tables */
/*         ) */
/* { */
/*     printf("F%d",n); */
/*     print_labels(arguments, N_KERNEL_ARGS); */
/*     printf("\n"); */
/*     for (int i = 0; i < TIME_STEPS; ++i) { */
/*         for (int j = 0; j < COMPONENTS; ++j) { */
/*             printf("%.5e\t",tables->kernels[index].values[i][j]); */
/*         } */
/*         printf("\n"); */
/*     } */
/* } */



double compute_k1(
        short int m,
        short int n_coeffs,
        const Vec1D<short int>& rearrangement,
        const Vec1D<bool>& signs,
        const Vec2D<double>& bare_scalar_products
        )
{
    double k1 = bare_scalar_products[n_coeffs - 1][n_coeffs - 1];
    for (int i = 2; i <= m; ++i) {
        int index = rearrangement[i-2];
        k1 += bare_scalar_products[index][index];
        // Note that elements in the signs-array correspond to rearranged loop
        // momenta, thus we use 'i-2', not 'index' as index here
        k1 -= 2 * (signs[i-2] ? 1 : -1) * bare_scalar_products[n_coeffs - 1][index];
    }

    for (int i = 3; i <= m; ++i) {
        for (int j = 2; j < i; ++j) {
            k1 += 2 * (signs[i-2] ? 1 : -1) * (signs[j-2] ? 1 : -1) *
                bare_scalar_products[rearrangement[i-2]][rearrangement[j-2]];
        }
    }

    return std::sqrt(k1);
}



inline int heaviside_theta(
        short int m,
        double k1,
        const Vec1D<short int>& rearrangement,
        const Vec1D<double>& Q_magnitudes
        )
{
    if (m == 1) return 1;

    // Heaviside-theta (k1 - k2)
    if (k1 <= Q_magnitudes[rearrangement[0]]) return 0;

#if DEBUG >= 1
    // Check that the heaviside-theta (k2 - k3) etc. are satisfied by
    // (reparametrized) momenta from CUBA
    for (int i = 3; i <= m; ++i) {
        if ( Q_magnitudes[rearrangement[i-3]]
                < Q_magnitudes[rearrangement[i-2]])
            throw(std::logic_error("Heaviside theta: Q" +
                        std::to_string(rearrangement[i-3] + 1) + " < Q" +
                        std::to_string(rearrangement[i-2] + 1) + "."));
    }
#endif
    return m;
}



void integrand_term(
        const PowerSpectrumDiagram& diagram,
        short int a,
        short int b,
        const Vec1D<Correlation>& correlations,
        IntegrandTables& tables,
        Vec1D<double>& term_results
        )
{
    short int kernel_index_l = diagram.arg_configs_l[a][b].kernel_index;
    short int kernel_index_r = diagram.arg_configs_r[a][b].kernel_index;

    const short int* arguments_l = diagram.arg_configs_l[a][b].args.data();
    const short int* arguments_r = diagram.arg_configs_r[a][b].args.data();

    /* Compute kernels */
    compute_SPT_kernels(arguments_l, kernel_index_l, 2*diagram.l + diagram.m, tables);
    compute_SPT_kernels(arguments_r, kernel_index_r, 2*diagram.r + diagram.m, tables);

    // Get values for two-point correlators
    // If l != r, there are two diagrams corresponding to l <-> r
    if (diagram.l == diagram.r) {
        for (size_t i = 0; i < correlations.size(); ++i) {
            term_results[i] = 
                tables.spt_kernels[kernel_index_l].values[correlations[i].first] *
                tables.spt_kernels[kernel_index_r].values[correlations[i].second];
        }
    } else {
        for (size_t i = 0; i < correlations.size(); ++i) {
            term_results[i] = 
                tables.spt_kernels[kernel_index_l].values[correlations[i].first] *
                tables.spt_kernels[kernel_index_r].values[correlations[i].second]
                +
                tables.spt_kernels[kernel_index_l].values[correlations[i].second] *
                tables.spt_kernels[kernel_index_r].values[correlations[i].first];
        }
    }
}



void integrand(
        const IntegrationInput& input, 
        IntegrandTables& tables,
        Vec1D<double>& results
        )
{
    size_t n_correlations = input.correlations.size();
    // Loop over all diagrams
    for (auto& dg : input.diagrams) {
#if DEBUG >= 2
        dg.print_diagram_tags(std::cout);
        std::cout << std::endl;
#endif
        Vec1D<double> diagram_results(n_correlations, 0);
        // Loop over momentum rearrangement and sign flips 
        for (short int a = 0; a < dg.n_rearrangements; ++a) {
            for (short int b = 0; b < dg.n_sign_configs; ++b) {
#if DEBUG >= 2
                dg.print_argument_configuration(std::cout, a, b);
#endif
                Vec1D<double> term_results(n_correlations);

                double k1 = compute_k1(dg.m, input.settings.n_coeffs,
                        dg.rearrangements[a], dg.sign_configs[b],
                        tables.bare_scalar_products);
                int h_theta = heaviside_theta(dg.m, k1, dg.rearrangements[a],
                        tables.vars.magnitudes);
                if (h_theta == 0) {
#if DEBUG >= 2
                    std::cout << "\t" << 0 << std::endl;
#endif
                    continue;
                }
                integrand_term(dg, a, b, input.correlations, tables, term_results);
                for (auto& el : term_results) {
                    el *= h_theta * input.input_ps.eval(k1);
                }
                for (size_t j = 0; j < n_correlations; ++j) {
                    diagram_results[j] += term_results[j];
                }
#if DEBUG >= 2
                for (auto& el : term_results) {
                    std::cout << "\t" << el;
                }
                std::cout << std::endl;
#endif
            }
        }

        for (size_t j = 0; j < n_correlations; ++j) {
            diagram_results[j] *= dg.diagram_factor;
            diagram_results[j] /= dg.n_rearrangements * dg.n_sign_configs;
            results[j] += diagram_results[j];
        }
    }
    for (short int i = 0; i < input.settings.n_loops; ++i) {
        for (auto& el : results) {
            el *= input.input_ps.eval(tables.vars.magnitudes[i]);
        }
    }
}



std::ostream& operator<<(std::ostream& out, const Correlation& correlation) {
    out << "<" << correlation.first << "," << correlation.second << ">";
    return out;
}
