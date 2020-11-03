/*
   integrand.cpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <stdexcept>
#include <utility>

#include "../include/integrand.hpp"

#include "../include/utilities.hpp"
#include "../include/parameters.hpp"
#include "../include/tables.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/kernel_evolution.hpp"
#include "../include/diagrams.hpp"
#include "../include/interpolation.hpp"
#include "../include/integrand.hpp"

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

using std::size_t;


namespace ps {
void integrand_term(
        const PowerSpectrumDiagram& diagram,
        int rearr_idx,
        int sign_idx,
        const Vec1D<Correlation>& correlations,
        IntegrandTables& tables,
        Vec1D<double>& term_results
        )
{
    const ArgumentConfiguration& arg_config_l =
        diagram.get_arg_config_l(rearr_idx, sign_idx);
    const ArgumentConfiguration& arg_config_r =
        diagram.get_arg_config_r(rearr_idx, sign_idx);

    /* Pointers to SPTKernel vector or last time step of Kernel vector */
    double* values_l = nullptr;
    double* values_r = nullptr;

    if (tables.loop_params.get_dynamics() == EDS_SPT) {
        compute_SPT_kernels(arg_config_l.args.data(),
                arg_config_l.kernel_index, 2 * diagram.l + diagram.m, tables);
        compute_SPT_kernels(arg_config_r.args.data(),
                arg_config_r.kernel_index, 2 * diagram.r + diagram.m, tables);

        values_l = tables.spt_kernels.at(arg_config_l.kernel_index).values;
        values_r = tables.spt_kernels.at(arg_config_r.kernel_index).values;
    }
    else if (tables.loop_params.get_dynamics() == EVOLVE_ASYMP_IC ||
             tables.loop_params.get_dynamics() == EVOLVE_EDS_IC) {
        kernel_evolution(arg_config_l.args.data(), arg_config_l.kernel_index,
                2 * diagram.l + diagram.m, tables);
        kernel_evolution(arg_config_r.args.data(), arg_config_r.kernel_index,
                2 * diagram.r + diagram.m, tables);

        int time_steps = tables.eta_grid.get_time_steps();
        values_l = tables.kernels.at(arg_config_l.kernel_index)
                       .values.at(time_steps - 1)
                       .data();
        values_r = tables.kernels.at(arg_config_r.kernel_index)
                       .values.at(time_steps - 1)
                       .data();
    }
    else {
        throw std::runtime_error("integrand_term(): Unknown dynamics.");
    }

    // Get values for two-point correlators
    // If l != r, there are two diagrams corresponding to l <-> r
    if (diagram.l == diagram.r) {
        for (size_t i = 0; i < correlations.size(); ++i) {
            term_results.at(i) = values_l[correlations.at(i).first] *
                                 values_r[correlations.at(i).second];
        }
    } else {
        for (size_t i = 0; i < correlations.size(); ++i) {
            term_results[i] = values_l[correlations.at(i).first] *
                              values_r[correlations.at(i).second]
                                    +
                              values_l[correlations.at(i).second] *
                              values_r[correlations.at(i).first];
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
    for (auto& dg : *input.ps_diagrams) {
#if DEBUG >= 2
        dg.print_diagram_tags(std::cout);
        std::cout << std::endl;
#endif
        Vec1D<double> diagram_results(n_correlations, 0.0);

        // Loop over momentum rearrangement and sign flips
        for (size_t rearr_idx = 0; rearr_idx < dg.n_rearrangements(); ++rearr_idx) {
            for (size_t sign_idx = 0; sign_idx < dg.n_sign_configs(); ++sign_idx) {
#if DEBUG >= 2
                dg.print_argument_configuration(std::cout, rearr_idx, sign_idx);
#endif

                double q_m1 = dg.q_m1(rearr_idx, sign_idx,
                        tables.scalar_products);
                int heaviside_theta = dg.heaviside_theta(q_m1, rearr_idx,
                        tables.vars.magnitudes);
                if (heaviside_theta == 0) {
#if DEBUG >= 2
                    std::cout << "\t" << 0 << std::endl;
#endif
                    continue;
                }

                Vec1D<double> term_results(n_correlations, 0.0);
                integrand_term(dg, rearr_idx, sign_idx, input.correlations,
                        tables, term_results);

                for (auto& el : term_results) {
                    el *= heaviside_theta * input.input_ps.eval(q_m1);
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
            diagram_results.at(j) *= dg.diagram_factor();
            diagram_results.at(j) /= dg.n_rearrangements() * dg.n_sign_configs();
            results.at(j) += diagram_results.at(j);
        }
    }
    for (int i = 0; i < tables.loop_params.get_n_loops(); ++i) {
        for (auto& el : results) {
            el *= input.input_ps.eval(tables.vars.magnitudes[i]);
        }
    }
}
}



std::ostream& operator<<(std::ostream& out, const Correlation& correlation) {
    out << "<" << correlation.first << "," << correlation.second << ">";
    return out;
}
