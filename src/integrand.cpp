/*
   integrand.cpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cuba.h>

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
void configuration_term(
        const PowerSpectrumDiagram& diagram,
        size_t rearr_idx,
        size_t sign_idx,
        const Vec1D<Pair<int>>& pair_correlations,
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

    if (tables.loop_params.dynamics() == EDS_SPT) {
        compute_SPT_kernels(arg_config_l.args.data(),
                arg_config_l.kernel_index, 2 * diagram.l + diagram.m, tables);
        compute_SPT_kernels(arg_config_r.args.data(),
                arg_config_r.kernel_index, 2 * diagram.r + diagram.m, tables);

        values_l = tables.spt_kernels
                       .at(static_cast<size_t>(arg_config_l.kernel_index))
                       .values;
        values_r = tables.spt_kernels
                       .at(static_cast<size_t>(arg_config_r.kernel_index))
                       .values;
    }
    else if (tables.loop_params.dynamics() == EVOLVE_IC_ASYMP ||
             tables.loop_params.dynamics() == EVOLVE_IC_EDS) {
        /* If EdS-SPT initial conditions, compute EdS-kernels */
        if (tables.loop_params.dynamics() == EVOLVE_IC_EDS) {
            compute_SPT_kernels(arg_config_l.args.data(),
                    arg_config_l.kernel_index, 2 * diagram.l + diagram.m, tables);
            compute_SPT_kernels(arg_config_r.args.data(),
                    arg_config_r.kernel_index, 2 * diagram.r + diagram.m, tables);
        }

        compute_gen_kernels(arg_config_l.args.data(), arg_config_l.kernel_index,
                2 * diagram.l + diagram.m, tables);
        compute_gen_kernels(arg_config_r.args.data(), arg_config_r.kernel_index,
                2 * diagram.r + diagram.m, tables);

        size_t time_steps = tables.eta_grid.time_steps();
        values_l = tables.kernels.at(static_cast<size_t>(arg_config_l.kernel_index))
                       .values.at(time_steps - 1)
                       .data();
        values_r = tables.kernels.at(static_cast<size_t>(arg_config_r.kernel_index))
                       .values.at(time_steps - 1)
                       .data();
    }
    else {
        throw std::runtime_error("integrand_term(): Unknown dynamics.");
    }

    // Get values for two-point correlators
    // If l != r, there are two diagrams corresponding to l <-> r
    if (diagram.l == diagram.r) {
        for (size_t i = 0; i < pair_correlations.size(); ++i) {
            term_results.at(i) = values_l[pair_correlations.at(i).first()] *
                                 values_r[pair_correlations.at(i).second()];
        }
    } else {
        for (size_t i = 0; i < pair_correlations.size(); ++i) {
            term_results[i] = values_l[pair_correlations.at(i).first()] *
                              values_r[pair_correlations.at(i).second()]
                                    +
                              values_l[pair_correlations.at(i).second()] *
                              values_r[pair_correlations.at(i).first()];
        }
    }
}



void diagram_term(
        const PowerSpectrumDiagram& diagram,
        const IntegrationInput& input,
        IntegrandTables& tables,
        Vec1D<double>& diagram_results
        )
{
#if DEBUG >= 2
    diagram.print_diagram_tags(std::cout);
    std::cout << std::endl;
#endif
    size_t n_correlations = input.pair_correlations.size();

    // Loop over momentum rearrangement and sign flips
    for (size_t i = 0; i < diagram.n_rearrangements(); ++i) {
        for (size_t j = 0; j < diagram.n_sign_configs(); ++j) {
#if DEBUG >= 2
            diagram.print_argument_configuration(std::cout, i, j);
#endif

            double q_m1 = diagram.q_m1(i, j, tables.comp_dot_products());
            int heaviside_theta = diagram.heaviside_theta(q_m1, i,
                    tables.vars.magnitudes);

            /* We set P(k < q_min) and P(k > q_max) to zero, which means we can
             * skip this configuration if q_m1 < qmin or q_m1 > q_max.
             * Furthermore if heaviside_theta == 0, we can also skip. */
            if (heaviside_theta == 0 || q_m1 < input.q_min ||
                q_m1 > input.q_max) {
#if DEBUG >= 2
                std::cout << "\t" << 0 << std::endl;
#endif
                continue;
            }

            Vec1D<double> term_results(n_correlations, 0.0);
            configuration_term(diagram, i, j, input.pair_correlations, tables,
                    term_results);

            for (auto& el : term_results) {
                el *= heaviside_theta * input.input_ps(q_m1);
            }
            for (size_t a = 0; a < n_correlations; ++a) {
                diagram_results.at(a) += term_results.at(a);
            }
#if DEBUG >= 2
            for (auto& el : term_results) {
                std::cout << "\t" << el;
            }
            std::cout << std::endl;
#endif
        }
    }

    for (size_t i = 0; i < n_correlations; ++i) {
        diagram_results.at(i) *= diagram.diagram_factor();
        diagram_results.at(i) /= static_cast<double>(
                diagram.n_rearrangements() * diagram.n_sign_configs());
    }
}
} /* namespace ps */



namespace bs {
void configuration_term(
        const BiSpectrumDiagram& diagram,
        size_t rearr_idx,
        size_t sign_idx,
        size_t overall_loop_idx,
        const Vec1D<Triple<int>>& triple_correlations,
        IntegrandTables& tables,
        Vec1D<double>& term_results
        )
{
    const Triple<ArgumentConfiguration>& arg_config =
        diagram.get_arg_config(rearr_idx, sign_idx, overall_loop_idx);

    /* Pointers to SPTKernel vector or last time step of Kernel vector */
    double* values_a = nullptr;
    double* values_b = nullptr;
    double* values_c = nullptr;

    if (tables.loop_params.dynamics() == EDS_SPT) {
        compute_SPT_kernels(arg_config.a().args.data(),
                arg_config.a().kernel_index, 2 * diagram.n_a + diagram.n_ab +
                diagram.n_ca, tables);
        compute_SPT_kernels(arg_config.b().args.data(),
                arg_config.b().kernel_index, 2 * diagram.n_b + diagram.n_ab +
                diagram.n_bc, tables);
        compute_SPT_kernels(arg_config.c().args.data(),
                arg_config.c().kernel_index, 2 * diagram.n_c + diagram.n_bc +
                diagram.n_ca, tables);

        values_a = tables.spt_kernels
                       .at(static_cast<size_t>(arg_config.a().kernel_index))
                       .values;
        values_b = tables.spt_kernels
                       .at(static_cast<size_t>(arg_config.b().kernel_index))
                       .values;
        values_c = tables.spt_kernels
                       .at(static_cast<size_t>(arg_config.c().kernel_index))
                       .values;
    }
    else if (tables.loop_params.dynamics() == EVOLVE_IC_ASYMP ||
             tables.loop_params.dynamics() == EVOLVE_IC_EDS) {
        /* If EdS-SPT initial conditions, compute EdS-kernels */
        if (tables.loop_params.dynamics() == EVOLVE_IC_EDS) {
            compute_SPT_kernels(arg_config.a().args.data(),
                    arg_config.a().kernel_index, 2 * diagram.n_a + diagram.n_ab +
                    diagram.n_ca, tables);
            compute_SPT_kernels(arg_config.b().args.data(),
                    arg_config.b().kernel_index, 2 * diagram.n_b + diagram.n_ab +
                    diagram.n_bc, tables);
            compute_SPT_kernels(arg_config.c().args.data(),
                    arg_config.c().kernel_index, 2 * diagram.n_c + diagram.n_bc +
                    diagram.n_ca, tables);
        }

        compute_gen_kernels(arg_config.a().args.data(),
                arg_config.a().kernel_index, 2 * diagram.n_a + diagram.n_ab +
                diagram.n_ca, tables);
        compute_gen_kernels(arg_config.b().args.data(),
                arg_config.b().kernel_index, 2 * diagram.n_b + diagram.n_ab +
                diagram.n_bc, tables);
        compute_gen_kernels(arg_config.c().args.data(),
                arg_config.c().kernel_index, 2 * diagram.n_c + diagram.n_bc +
                diagram.n_ca, tables);

        size_t time_steps = tables.eta_grid.time_steps();
        values_a =
            tables.kernels.at(static_cast<size_t>(arg_config.a().kernel_index))
                .values.at(time_steps - 1)
                .data();
        values_b =
            tables.kernels.at(static_cast<size_t>(arg_config.b().kernel_index))
                .values.at(time_steps - 1)
                .data();
        values_c =
            tables.kernels.at(static_cast<size_t>(arg_config.c().kernel_index))
                .values.at(time_steps - 1)
                .data();
    }
    else {
        throw std::runtime_error("integrand_term(): Unknown dynamics.");
    }

    for (size_t i = 0; i < triple_correlations.size(); ++i) {
        term_results[i] = values_a[triple_correlations.at(i).first()] *
                          values_b[triple_correlations.at(i).second()] *
                          values_c[triple_correlations.at(i).third()];
    }
}



void diagram_term(
        const BiSpectrumDiagram& diagram,
        const IntegrationInput& input,
        IntegrandTables& tables,
        Vec1D<double>& diagram_results
        )
{
#if DEBUG >= 2
    diagram.print_diagram_tags(std::cout);
    std::cout << std::endl;
#endif
    size_t n_correlations = input.triple_correlations.size();

    /* Loop over momentum rearrangement, sign flips and overall loop assosiations */
    size_t overall_loop_assosiations = diagram.overall_loop() ? 6 : 1;
    for (size_t i = 0; i < diagram.n_rearrangements(); ++i) {
        for (size_t j = 0; j < diagram.n_sign_configs(); ++j) {
            for (size_t k = 0; k < overall_loop_assosiations; ++k) {
#if DEBUG >= 2
                diagram.print_argument_configuration(std::cout, i, j, k);
#endif
                Triple<double> q_xy1 = {0,0,0};
                int heaviside_theta = 1;
                diagram.connecting_lines_factors(i, j, k,
                        tables.vars.magnitudes, tables.comp_dot_products(),
                        q_xy1, heaviside_theta);

                if (heaviside_theta == 0) {
#if DEBUG >= 2
                    std::cout << "\t" << 0 << std::endl;
#endif
                    continue;
                }

                Vec1D<double> term_results(n_correlations, 0.0);
                configuration_term(diagram, i, j, k,
                        input.triple_correlations, tables, term_results);

                /* Multiply by P_lin(q_xy1) if xy connecting line is present.
                 * Use Interpolation1D::eval(x, min, max) which returns 0 when
                 * x is not between min and max */
                if (diagram.n_ab > 0) {
                    for (auto& el : term_results) {
                        el *= input.input_ps(q_xy1.a(), input.q_min,
                                input.q_max);
                    }
                }
                if (diagram.n_bc > 0) {
                    for (auto& el : term_results) {
                        el *= input.input_ps(q_xy1.b(), input.q_min,
                                input.q_max);
                    }
                }
                if (diagram.n_ca > 0) {
                    for (auto& el : term_results) {
                        el *= input.input_ps(q_xy1.c(), input.q_min,
                                input.q_max);
                    }
                }
                /* For overall-loop diagrams: divide by pLin(rearr(Q)) which is
                 * multiplied also in bs::integrand() */
                if (diagram.overall_loop()) {
                    for (auto& el : term_results) {
                        el /= input.input_ps(
                            diagram.q1_magnitude(i, tables.vars.magnitudes));
                    }
                }
                for (auto& el : term_results) {
                    el *= heaviside_theta;
                }
                for (size_t a = 0; a < n_correlations; ++a) {
                    diagram_results.at(a) += term_results.at(a);
                }
#if DEBUG >= 2
                for (auto& el : term_results) {
                    std::cout << "\t" << el;
                }
                std::cout << std::endl;
#endif
            }
        }
    }

    for (size_t i = 0; i < n_correlations; ++i) {
        diagram_results.at(i) *= diagram.diagram_factor();
        /* Divide by # rearrangments, sign configurations for connecting lines
         * and divide by 2 for switching sign of overall loop (if present) */
        diagram_results.at(i) /= static_cast<double>(
            diagram.n_rearrangements() * diagram.n_sign_configs() *
            (diagram.overall_loop() ? 2 : 1));
    }
}
} /* namespace bs */



int integrand(
        __attribute__((unused)) const int *ndim,
        const cubareal xx[],
        __attribute__((unused)) const int *ncomp,
        cubareal ff[],
        void *userdata,
        __attribute__((unused)) const int *nvec,
        const int *core
        )
{
    IntegrationInput& input = *static_cast<IntegrationInput*>(userdata);

    /*  For thread <*core + 1> (index 0 is reserved for master), we use the */
    /*  IntegrandTables number *core+1 */
    IntegrandTables& tables = input.tables_vec.at(static_cast<size_t>(*core + 1));

    int n_loops = tables.loop_params.n_loops();
    IntegrationVariables& vars = tables.vars;


#if DEBUG > 1
    /* Check that n_loops > 0 */
    if (n_loops < 1) {
        throw std::runtime_error("n_loops is not larger than or equal to 1.");
    }
#endif

    double ratio = input.q_max/input.q_min;
    double log_ratio = std::log(ratio);
    double jacobian = 1.0;

    /* 1-loop */
    if (!input.single_hard_limit) {
        vars.magnitudes.at(0) = input.q_min * pow(ratio,xx[0]);
        jacobian *= log_ratio * CUBE(vars.magnitudes.at(0));
    }

    int idx = 1;
    vars.cos_theta.at(0) = xx[idx++];
    if (tables.loop_params.rsd()) {
        vars.phi.at(0) = xx[idx++] * TWOPI;
        jacobian *= TWOPI;
    }

    /* 2-loop */
    if (n_loops > 1){
        if (!input.single_hard_limit) {
            vars.magnitudes.at(1) = pow(vars.magnitudes.at(0), xx[idx++]);
            jacobian *= xx[0];
        }
        else {
            vars.magnitudes.at(1) = input.q_min * pow(ratio,xx[idx]);
        }
        jacobian *= log_ratio * CUBE(vars.magnitudes.at(1));
        vars.cos_theta.at(1) = xx[idx++];
        vars.phi.at(1) = xx[idx++] * TWOPI;
        jacobian *= TWOPI;
    }

    size_t n_correlations = static_cast<size_t>(*ncomp);

    Vec1D<double> results(n_correlations, 0.0);
    Vec1D<double> diagram_results(n_correlations, 0.0);

    try {
        /* Zero-initialize kernel tables */
        tables.reset();
        // Compute dot_products-, alpha- and beta-tables
        tables.compute_tables();

        // Loop over all diagrams
        if (tables.loop_params.spectrum() == POWERSPECTRUM) {
            for (auto& diagram : input.ps_diagrams) {
                ps::diagram_term(diagram, input, tables, diagram_results);

                /* Add diagram results to results and reset diagram results */
                for (size_t j = 0; j < n_correlations; ++j) {
                    results.at(j) += diagram_results.at(j);
                }
                std::fill(diagram_results.begin(), diagram_results.end(), 0.0);
            }
        }
        else {
            /* Bispectrum */
            for (auto& diagram : input.bs_diagrams) {
                bs::diagram_term(diagram, input, tables, diagram_results);

                /* Add diagram results to results and reset diagram results */
                for (size_t j = 0; j < n_correlations; ++j) {
                    results.at(j) += diagram_results.at(j);
                }
                std::fill(diagram_results.begin(), diagram_results.end(), 0.0);
            }
        }

        int i = 0;
        /* Skip P_lin(Q1) if single_hard_limit = true */
        if (input.single_hard_limit) ++i;
        for (; i < tables.loop_params.n_loops(); ++i) {
            for (auto& el : results) {
                el *= input.input_ps(
                    tables.vars.magnitudes.at(static_cast<size_t>(i)));
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        /* Tell CUBA an error occured */
        return -999;
    }

    for (size_t i = 0; i < n_correlations; ++i) {
        ff[i] = results.at(i) * jacobian;
    }

    /* Return success */
    return 0;
}
