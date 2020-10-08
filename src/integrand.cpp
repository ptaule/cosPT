/*
   integrand.cpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <stdexcept>
#include <utility>
#include <cmath>

#include "../include/integrand.hpp"

#include "../include/utilities.hpp"
#include "../include/tables.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/diagrams.hpp"
#include "../include/interpolation.hpp"
#include "../include/integrand.hpp"

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

using std::size_t;


double compute_k1(
        int m,
        int n_coeffs,
        const Vec1D<int>& rearrangement,
        const Vec1D<bool>& signs,
        const Vec2D<double>& bare_scalar_products
        )
{
    double k1 = bare_scalar_products.at(n_coeffs - 1).at(n_coeffs - 1);
    for (int i = 2; i <= m; ++i) {
        int index = rearrangement.at(i-2);
        k1 += bare_scalar_products.at(index).at(index);
        // Note that elements in the signs-array correspond to rearranged loop
        // momenta, thus we use 'i-2', not 'index' as index here
        k1 -= 2 * (signs.at(i-2) ? 1 : -1)
            * bare_scalar_products.at(n_coeffs - 1).at(index);
    }

    for (int i = 3; i <= m; ++i) {
        for (int j = 2; j < i; ++j) {
            k1 += 2 * (signs.at(i-2) ? 1 : -1) * (signs.at(j-2) ? 1 : -1) *
                bare_scalar_products.at(rearrangement.at(i-2)).at(rearrangement.at(j-2));
        }
    }

    return std::sqrt(k1);
}



inline int heaviside_theta(
        int m,
        double k1,
        const Vec1D<int>& rearrangement,
        const Vec1D<double>& Q_magnitudes
        )
{
    if (m == 1) return 1;

    // Heaviside-theta (k1 - k2)
    if (k1 <= Q_magnitudes.at(rearrangement.at(0))) return 0;

#if DEBUG >= 1
    // Check that the heaviside-theta (k2 - k3) etc. are satisfied by
    // (reparametrized) momenta from CUBA
    for (int i = 3; i <= m; ++i) {
        if ( Q_magnitudes.at(rearrangement.at(i-3))
                < Q_magnitudes.at(rearrangement.at(i-2)))
            throw(std::logic_error(
                "Heaviside theta: Q" +
                std::to_string(rearrangement.at(i - 3) + 1) + " < Q" +
                std::to_string(rearrangement.at(i - 2) + 1) + "."));
    }
#endif
    return m;
}



void integrand_term(
        const PowerSpectrumDiagram& diagram,
        int a,
        int b,
        const Vec1D<Correlation>& correlations,
        IntegrandTables& tables,
        Vec1D<double>& term_results
        )
{
    int kernel_index_l = diagram.arg_configs_l.at(a).at(b).kernel_index;
    int kernel_index_r = diagram.arg_configs_r.at(a).at(b).kernel_index;

    const int* arguments_l = diagram.arg_configs_l.at(a).at(b).args.data();
    const int* arguments_r = diagram.arg_configs_r.at(a).at(b).args.data();

    /* Compute kernels */
    compute_SPT_kernels(arguments_l, kernel_index_l, 2*diagram.l + diagram.m, tables);
    compute_SPT_kernels(arguments_r, kernel_index_r, 2*diagram.r + diagram.m, tables);

    // Get values for two-point correlators
    // If l != r, there are two diagrams corresponding to l <-> r
    if (diagram.l == diagram.r) {
        for (size_t i = 0; i < correlations.size(); ++i) {
            term_results.at(i) =
                tables.spt_kernels.at(kernel_index_l).values[correlations.at(i).first] *
                tables.spt_kernels.at(kernel_index_r).values[correlations.at(i).second];
        }
    } else {
        for (size_t i = 0; i < correlations.size(); ++i) {
            term_results[i] =
                tables.spt_kernels.at(kernel_index_l).values[correlations.at(i).first] *
                tables.spt_kernels.at(kernel_index_r).values[correlations.at(i).second]
                +
                tables.spt_kernels.at(kernel_index_l).values[correlations.at(i).second] *
                tables.spt_kernels.at(kernel_index_r).values[correlations.at(i).first];
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
        Vec1D<double> diagram_results(n_correlations, 0.0);

        // Loop over momentum rearrangement and sign flips
        for (size_t a = 0; a < dg.rearrangements.size(); ++a) {
            for (size_t b = 0; b < dg.sign_configs.size(); ++b) {
#if DEBUG >= 2
                dg.print_argument_configuration(std::cout, a, b);
#endif

                double k1 = compute_k1(dg.m, tables.params.n_coeffs,
                        dg.rearrangements.at(a), dg.sign_configs.at(b),
                        tables.bare_scalar_products);
                int h_theta = heaviside_theta(dg.m, k1, dg.rearrangements.at(a),
                        tables.vars.magnitudes);
                if (h_theta == 0) {
#if DEBUG >= 2
                    std::cout << "\t" << 0 << std::endl;
#endif
                    continue;
                }

                Vec1D<double> term_results(n_correlations, 0.0);
                integrand_term(dg, a, b, input.correlations, tables, term_results);
                for (auto& el : term_results) {
                    el *= h_theta * input.input_ps.eval(k1);
                }
                for (size_t j = 0; j < n_correlations; ++j) {
                    diagram_results.at(j) += term_results.at(j);
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
            diagram_results.at(j) *= dg.diagram_factor;
            diagram_results.at(j) /= dg.rearrangements.size() * dg.sign_configs.size();
            results.at(j) += diagram_results.at(j);
        }
    }
    for (int i = 0; i < tables.params.n_loops; ++i) {
        for (auto& el : results) {
            el *= input.input_ps.eval(tables.vars.magnitudes[i]);
        }
    }
}



std::ostream& operator<<(std::ostream& out, const Correlation& correlation) {
    out << "<" << correlation.first << "," << correlation.second << ">";
    return out;
}
