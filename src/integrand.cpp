/*
   integrand.cpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <stdexcept>
#include <utility>
#include <array>
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
        int rearr_idx,
        int sign_idx,
        const Vec1D<PairCorrelation>& pair_correlations,
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
        for (size_t i = 0; i < pair_correlations.size(); ++i) {
            term_results.at(i) = values_l[pair_correlations.at(i).first] *
                                 values_r[pair_correlations.at(i).second];
        }
    } else {
        for (size_t i = 0; i < pair_correlations.size(); ++i) {
            term_results[i] = values_l[pair_correlations.at(i).first] *
                              values_r[pair_correlations.at(i).second]
                                    +
                              values_l[pair_correlations.at(i).second] *
                              values_r[pair_correlations.at(i).first];
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
    size_t n_correlations = input.pair_correlations->size();

    // Loop over momentum rearrangement and sign flips
    for (size_t i = 0; i < diagram.n_rearrangements(); ++i) {
        for (size_t j = 0; j < diagram.n_sign_configs(); ++j) {
#if DEBUG >= 2
            diagram.print_argument_configuration(std::cout, i, j);
#endif

            double q_m1 = diagram.q_m1(i, j, tables.scalar_products);
            int heaviside_theta = diagram.heaviside_theta(q_m1, i,
                    tables.vars.magnitudes);
            if (heaviside_theta == 0) {
#if DEBUG >= 2
                std::cout << "\t" << 0 << std::endl;
#endif
                continue;
            }

            Vec1D<double> term_results(n_correlations, 0.0);
            configuration_term(diagram, i, j, *input.pair_correlations, tables,
                    term_results);

            for (auto& el : term_results) {
                el *= heaviside_theta * input.input_ps.eval(q_m1);
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
        diagram_results.at(i) /= diagram.n_rearrangements() *
            diagram.n_sign_configs();
    }
}



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
    IntegrationInput* input = (IntegrationInput*)userdata;

    /*  For thread <*core + 1> (index 0 is reserved for master), we use the */
    /*  IntegrandTables number *core+1 */
    IntegrandTables& tables = input->tables_vec.at(*core + 1);

    int n_loops = tables.loop_params.get_n_loops();
    IntegrationVariables& vars = tables.vars;

    double ratio = input->q_max/input->q_min;
    double log_ratio = std::log(ratio);
    double jacobian = 0.0;

    switch (n_loops) {
        case 1:
            vars.magnitudes.at(0) = input->q_min * pow(ratio,xx[0]);
            vars.cos_theta.at(0) = xx[1];
            jacobian = log(ratio) * CUBE(vars.magnitudes[0]);
            break;
        case 2:
            vars.magnitudes.at(0) = input->q_min * pow(ratio,xx[0]);
            vars.magnitudes.at(1) = input->q_min * pow(ratio,xx[0] * xx[1]);
            vars.cos_theta.at(0) = xx[2];
            vars.cos_theta.at(1) = xx[3];
            /* We may fix the coordinate system s.t. vars.phi[0] = 0 */
            vars.phi.at(1) = xx[4] * TWOPI;
            jacobian = TWOPI * xx[0]
                * SQUARE(log_ratio)
                * CUBE(vars.magnitudes[0])
                * CUBE(vars.magnitudes[1]);
            break;
        default:
            throw(std::invalid_argument("ps::integrand(): n_loops is not 1 or 2."));
    }

    size_t n_correlations = input->pair_correlations->size();
    Vec1D<double> results(n_correlations, 0.0);
    Vec1D<double> diagram_results(n_correlations, 0.0);
    try {
        /* Zero-initialize kernel tables */
        tables.reset();
        // Compute scalar_products-, alpha- and beta-tables
        tables.compute_tables();

        // Loop over all diagrams
        for (auto& diagram : *input->ps_diagrams) {
            diagram_term(diagram, *input, tables, diagram_results);

            /* Add diagram results to results and reset diagram results */
            for (size_t j = 0; j < n_correlations; ++j) {
                results.at(j) += diagram_results.at(j);
            }
            std::fill(diagram_results.begin(), diagram_results.end(), 0.0);
        }
        for (int i = 0; i < tables.loop_params.get_n_loops(); ++i) {
            for (auto& el : results) {
                el *= input->input_ps.eval(tables.vars.magnitudes[i]);
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
}



namespace bs {
void configuration_term(
        const BiSpectrumDiagram& diagram,
        int rearr_idx,
        int sign_idx,
        int overall_loop_idx,
        const Vec1D<TripleCorrelation>& triple_correlations,
        IntegrandTables& tables,
        Vec1D<double>& term_results
        )
{
    const ArgumentConfiguration& arg_config_a =
        diagram.get_arg_config_a(rearr_idx, sign_idx, overall_loop_idx);
    const ArgumentConfiguration& arg_config_b =
        diagram.get_arg_config_b(rearr_idx, sign_idx, overall_loop_idx);
    const ArgumentConfiguration& arg_config_c =
        diagram.get_arg_config_c(rearr_idx, sign_idx, overall_loop_idx);

    /* Pointers to SPTKernel vector or last time step of Kernel vector */
    double* values_a = nullptr;
    double* values_b = nullptr;
    double* values_c = nullptr;

    if (tables.loop_params.get_dynamics() == EDS_SPT) {
        compute_SPT_kernels(arg_config_a.args.data(),
                arg_config_a.kernel_index, 2 * diagram.n_a + diagram.n_ab +
                diagram.n_ca, tables);
        compute_SPT_kernels(arg_config_b.args.data(),
                arg_config_b.kernel_index, 2 * diagram.n_b + diagram.n_ab +
                diagram.n_bc, tables);
        compute_SPT_kernels(arg_config_c.args.data(),
                arg_config_c.kernel_index, 2 * diagram.n_c + diagram.n_bc +
                diagram.n_ca, tables);

        values_a = tables.spt_kernels.at(arg_config_a.kernel_index).values;
        values_b = tables.spt_kernels.at(arg_config_b.kernel_index).values;
        values_c = tables.spt_kernels.at(arg_config_c.kernel_index).values;
    }
    else if (tables.loop_params.get_dynamics() == EVOLVE_ASYMP_IC ||
             tables.loop_params.get_dynamics() == EVOLVE_EDS_IC) {
        kernel_evolution(arg_config_a.args.data(),
                arg_config_a.kernel_index, 2 * diagram.n_a + diagram.n_ab +
                diagram.n_ca, tables);
        kernel_evolution(arg_config_b.args.data(),
                arg_config_b.kernel_index, 2 * diagram.n_b + diagram.n_ab +
                diagram.n_bc, tables);
        kernel_evolution(arg_config_c.args.data(),
                arg_config_c.kernel_index, 2 * diagram.n_c + diagram.n_bc +
                diagram.n_ca, tables);

        int time_steps = tables.eta_grid.get_time_steps();
        values_a = tables.kernels.at(arg_config_a.kernel_index)
                       .values.at(time_steps - 1).data();
        values_b = tables.kernels.at(arg_config_b.kernel_index)
                       .values.at(time_steps - 1).data();
        values_c = tables.kernels.at(arg_config_c.kernel_index)
                       .values.at(time_steps - 1).data();
    }
    else {
        throw std::runtime_error("integrand_term(): Unknown dynamics.");
    }

    for (size_t i = 0; i < triple_correlations.size(); ++i) {
        term_results[i] = values_a[triple_correlations.at(i).at(0)] *
                          values_b[triple_correlations.at(i).at(1)] *
                          values_c[triple_correlations.at(i).at(2)];
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
    size_t n_correlations = input.triple_correlations->size();

    /* Loop over momentum rearrangement, sign flips and overall loop assosiations */
    int overall_loop_assosiations = diagram.overall_loop() ? 3 : 1;
    for (size_t i = 0; i < diagram.n_rearrangements(); ++i) {
        for (size_t j = 0; j < diagram.n_sign_configs(); ++j) {
            for (int k = 0; k < overall_loop_assosiations; ++k) {
#if DEBUG >= 2
                diagram.print_argument_configuration(std::cout, i, j, k);
#endif
                double q_ab1, q_bc1, q_ca1;
                int heaviside_theta = 1;
                diagram.connecting_lines_factors(i, j, k,
                        tables.scalar_products, q_ab1, q_bc1, q_ca1,
                        heaviside_theta);

                if (heaviside_theta == 0) {
#if DEBUG >= 2
                    std::cout << "\t" << 0 << std::endl;
#endif
                    continue;
                }

                Vec1D<double> term_results(n_correlations, 0.0);
                configuration_term(diagram, i, j, k,
                        *input.triple_correlations, tables, term_results);

                for (auto& el : term_results) {
                    el *= heaviside_theta;
                    el *= input.input_ps.eval(q_ab1);
                    el *= input.input_ps.eval(q_bc1);
                    el *= input.input_ps.eval(q_ca1);
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
        diagram_results.at(i) /= diagram.n_rearrangements() *
            diagram.n_sign_configs();
    }
}



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
    IntegrationInput* input = (IntegrationInput*)userdata;

    /*  For thread <*core + 1> (index 0 is reserved for master), we use the */
    /*  IntegrandTables number *core+1 */
    IntegrandTables& tables = input->tables_vec.at(*core + 1);

    int n_loops = tables.loop_params.get_n_loops();
    IntegrationVariables& vars = tables.vars;

    double ratio = input->q_max/input->q_min;
    double jacobian = 0.0;

    switch (n_loops) {
        case 1:
            vars.magnitudes.at(0) = input->q_min * pow(ratio,xx[0]);
            vars.cos_theta.at(0) = xx[1];
            vars.phi.at(0) = xx[2];
            jacobian = TWOPI * log(ratio) * CUBE(vars.magnitudes[0]);
            break;
        default:
            throw(std::invalid_argument("bs::integrand(): n_loops is not 1."));
    }

    size_t n_correlations = input->triple_correlations->size();
    Vec1D<double> results(n_correlations, 0.0);
    Vec1D<double> diagram_results(n_correlations, 0.0);
    try {
        /* Zero-initialize kernel tables */
        tables.reset();
        // Compute scalar_products-, alpha- and beta-tables
        tables.compute_tables();

        // Loop over all diagrams
        for (auto& diagram : *input->bs_diagrams) {
            diagram_term(diagram, *input, tables, diagram_results);

            /* Add diagram results to results and reset diagram results */
            for (size_t j = 0; j < n_correlations; ++j) {
                results.at(j) += diagram_results.at(j);
            }
            std::fill(diagram_results.begin(), diagram_results.end(), 0.0);
        }
        for (int i = 0; i < tables.loop_params.get_n_loops(); ++i) {
            for (auto& el : results) {
                el *= input->input_ps.eval(tables.vars.magnitudes[i]);
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
}



std::ostream& operator<<(std::ostream& out, const PairCorrelation& pair_correlation) {
    out << "<" << pair_correlation.first << "," << pair_correlation.second << ">";
    return out;
}



std::ostream& operator<<(std::ostream& out, const TripleCorrelation&
        triple_correlation)
{
    out << "<" << triple_correlation.at(0) << "," << triple_correlation.at(1) << ","
        << triple_correlation.at(2) << ">";
    return out;
}
