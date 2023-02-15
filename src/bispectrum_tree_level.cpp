/*
   bispectrum_tree_level.cpp

   Created by Petter Taule on 19.12.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <stdexcept>

#include "../include/interpolation.hpp"
#include "../include/diagrams.hpp"
#include "../include/kernel_evolution.hpp"
#include "../include/parameters.hpp"
#include "../include/tables.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/bispectrum_tree_level.hpp"

using std::size_t;

Triple<ArgumentConfiguration> kernel_arguments(
        int diagram_idx,
        const LoopParameters& loop_params
        )
{
    size_t n_coeffs      = loop_params.n_coeffs();
    size_t n_kernel_args = loop_params.n_kernel_args();
    int zero_label    = loop_params.zero_label();

    size_t k_a_idx = n_coeffs - 1;
    size_t k_b_idx = n_coeffs - 2;

    /* Configs for k_a, k_b and k_c = -k_a - k_b */
    Triple<Vec1D<int>> configs(
        Vec1D<int>(n_coeffs, 0),
        Vec1D<int>(n_coeffs, 0),
        Vec1D<int>(n_coeffs, 0)
    );
    configs.a().at(k_a_idx) = 1;
    configs.b().at(k_b_idx) = 1;
    configs.c().at(k_a_idx) = -1;
    configs.c().at(k_b_idx) = -1;

    Triple<int> labels(
        config2label(configs.a()),
        config2label(configs.b()),
        config2label(configs.c())
    );

    Triple<ArgumentConfiguration> arg_config;
    /* Allocate memory */
    arg_config.a().args.resize(n_kernel_args);
    arg_config.b().args.resize(n_kernel_args);
    arg_config.c().args.resize(n_kernel_args);

    /* Arg config a has always 2 arguments, i.e. corresponds to two incoming
     * connecting lines */
    switch (diagram_idx) {
        case 0:
            /* B000110 */
            arg_config.a().args.at(0) = flip_signs(labels.a(), n_coeffs);
            arg_config.a().args.at(1) = flip_signs(labels.c(), n_coeffs);
            arg_config.b().args.at(0) = labels.a();
            arg_config.c().args.at(0) = labels.c();
            break;
        case 1:
            /* B000101 */
            arg_config.a().args.at(0) = flip_signs(labels.b(), n_coeffs);
            arg_config.a().args.at(1) = flip_signs(labels.c(), n_coeffs);
            arg_config.b().args.at(0) = labels.b();
            arg_config.c().args.at(0) = labels.c();
            break;
        case 2:
            /* B000011 */
            arg_config.a().args.at(0) = flip_signs(labels.a(), n_coeffs);
            arg_config.a().args.at(1) = flip_signs(labels.b(), n_coeffs);
            arg_config.b().args.at(0) = labels.a();
            arg_config.c().args.at(0) = labels.b();
            break;
        default:
            throw(std::logic_error(
                "kernel_arguments(): got argument i which is not 0,1,2."));
    }

    /* Set remaining arguments to zero (zero_label) */
    arg_config.b().args.at(1) = zero_label;
    arg_config.c().args.at(1) = zero_label;
    for (size_t j = 2; j < n_kernel_args; ++j) {
        arg_config.a().args.at(j) = zero_label;
        arg_config.b().args.at(j) = zero_label;
        arg_config.c().args.at(j) = zero_label;
    }

    arg_config.a().kernel_index =
        loop_params.arguments_2_kernel_index(arg_config.a().args);
    arg_config.b().kernel_index =
        loop_params.arguments_2_kernel_index(arg_config.b().args);
    arg_config.c().kernel_index =
        loop_params.arguments_2_kernel_index(arg_config.c().args);

    return arg_config;
}



void diagram_term(
        int diagram_idx,
        IntegrandTables& tables,
        const Interpolation1D& input_ps,
        const Vec1D<Triple<int>>& triple_correlations,
        Vec1D<double>& diagram_results /* out */
        )
{
    size_t n_coeffs = tables.loop_params.n_coeffs();
    size_t k_a_idx = n_coeffs - 1;
    size_t k_b_idx = n_coeffs - 2;

    Triple<ArgumentConfiguration> arg_config =
        kernel_arguments(diagram_idx, tables.loop_params);

    /* Pointers to SPTKernel vector or last time step of Kernel vector */
    double* values_a = nullptr;
    double* values_b = nullptr;
    double* values_c = nullptr;

    if (tables.loop_params.dynamics() == EDS_SPT) {
        compute_SPT_kernels(arg_config.a().args.data(),
                arg_config.a().kernel_index, 2, tables);
        compute_SPT_kernels(arg_config.b().args.data(),
                arg_config.b().kernel_index, 1, tables);
        compute_SPT_kernels(arg_config.c().args.data(),
                arg_config.c().kernel_index, 1, tables);

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
        if (tables.loop_params.dynamics() == EVOLVE_IC_EDS) {
            compute_SPT_kernels(arg_config.a().args.data(),
                    arg_config.a().kernel_index, 2, tables);
            compute_SPT_kernels(arg_config.b().args.data(),
                    arg_config.b().kernel_index, 1, tables);
            compute_SPT_kernels(arg_config.c().args.data(),
                    arg_config.c().kernel_index, 1, tables);
        }

        compute_gen_kernels(arg_config.a().args.data(),
                arg_config.a().kernel_index, 2, tables);
        compute_gen_kernels(arg_config.b().args.data(),
                arg_config.b().kernel_index, 1, tables);
        compute_gen_kernels(arg_config.c().args.data(),
                arg_config.c().kernel_index, 1, tables);

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

    /* External momentum values */
    double k_a = std::sqrt(tables.bare_dot_products().at(k_a_idx).at(k_a_idx));
    double k_b = std::sqrt(tables.bare_dot_products().at(k_b_idx).at(k_b_idx));
    double k_c = std::sqrt(SQUARE(k_a) + SQUARE(k_b) + 2 *
            tables.bare_dot_products().at(k_a_idx).at(k_b_idx));


    for (size_t i = 0; i < triple_correlations.size(); ++i) {
        diagram_results.at(i) = 2 *
            values_a[triple_correlations.at(i).first()] *
            values_b[triple_correlations.at(i).second()] *
            values_c[triple_correlations.at(i).third()];

        switch (diagram_idx) {
            case 0:
                diagram_results.at(i) *= input_ps(k_a);
                diagram_results.at(i) *= input_ps(k_c);
              break;
            case 1:
                diagram_results.at(i) *= input_ps(k_b);
                diagram_results.at(i) *= input_ps(k_c);
                break;
            case 2:
                diagram_results.at(i) *= input_ps(k_a);
                diagram_results.at(i) *= input_ps(k_b);
                break;
            default:
                throw(std::logic_error(
                            "diagram_term(): got argument i which is not 0,1,2."));
        }
    }
}



void tree_level_bispectrum(
        IntegrandTables& tables,
        const Interpolation1D& input_ps,
        const Vec1D<Triple<int>>& triple_correlations,
        Vec1D<double>& results /* out */
        )
{
    /* Set some irrelevant values for non-existent loops */
    for (size_t i = 0; i < static_cast<size_t>(tables.loop_params.n_loops()); ++i) {
        tables.vars.magnitudes.at(i) = 0;
        tables.vars.cos_theta.at(i) = 0;
        tables.vars.phi.at(i) = 0;
    }
    /* Zero-initialize kernel tables */
    tables.reset();
    /* Compute dot_products-, alpha- and beta-tables */
    tables.compute_tables();

    /* Three diagrams: 000110, 000101, 000011 */
    for (int i = 0; i < 3; ++i) {
        Vec1D<double> diagram_results(triple_correlations.size(), 0.0);
        diagram_term(i, tables, input_ps, triple_correlations, diagram_results);

        for (size_t j = 0; j < triple_correlations.size(); ++j) {
            results.at(j) += diagram_results.at(j);
        }
    }

    /* Zero-initialize kernel tables again, to ensure that values do not
     * corrupt futher computations */
    tables.reset();
}
