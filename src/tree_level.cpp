/*
   tree_level.cpp

   Created by Petter Taule on 19.12.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <stdexcept>

#include <gsl/gsl_integration.h>

#include "../include/diagrams.hpp"
#include "../include/ir_resum.hpp"
#include "../include/kernel_evolution.hpp"
#include "../include/parameters.hpp"
#include "../include/tables.hpp"
#include "../include/tree_level.hpp"
#include "../include/spt_kernels.hpp"

using std::size_t;


namespace ps {

void rsd_tree_level(
    double k,
    const InputPowerSpectrum& ps,
    Vec1D<double>& results /* out */
    )
{
    double rsd_growth_f = ps.rsd_growth_f();
    for (auto& el : results) el = ps.tree_level(k, 0);
    /* Linear monopole */
    results.at(0) *= (
            1 + 2.0/3.0 * rsd_growth_f +
            1.0/5.0 * SQUARE(rsd_growth_f)
            );
    /* Linear quadrupole */
    results.at(1) *= (
            4.0/3.0 * rsd_growth_f
            + 4.0/7.0 * SQUARE(rsd_growth_f)
            );
    /* Linear hexadecapole */
    results.at(2) *= (8.0/35.0 * SQUARE(rsd_growth_f));
}



void rsd_tree_level_ir_resum(
    double k,
    const InputPowerSpectrum& ps,
    Vec1D<double>& results, /* out */
    std::size_t integration_sub_regions = 10000,
    double integration_atol = 0,
    double integration_rtol = 1e-6,
    int integration_key = GSL_INTEG_GAUSS61
    )
{
    gsl_integration_workspace* workspace =
        gsl_integration_workspace_alloc(integration_sub_regions);

    /* l=0 integration */
    auto integral_l0 = [&k, &ps](double mu) {
        return ps.tree_level(k, mu);
    };

    gsl_function F;
    F.function = [] (double x, void* p) {
        return (*(decltype(integral_l0)*)p)(x);
    };
    F.params = &integral_l0;

    double abserr;
    int status = gsl_integration_qag(&F, 0, 1, integration_atol,
                                     integration_rtol,
                                     integration_sub_regions, integration_key,
                                     workspace, &results.at(0), &abserr);

    if (status != 0) {
        throw std::runtime_error("RSD tree-level integration l=0 failed \
                with error code" + std::to_string(status));
    }

    /* l=2 integration */
    auto integral_l2 = [&k, &ps](double mu) {
        return 0.5 * (3*mu*mu - 1) * ps.tree_level(k, mu);
    };

    F.function = [] (double x, void* p) {
        return (*(decltype(integral_l2)*)p)(x);
    };
    F.params = &integral_l2;

    status = gsl_integration_qag(&F, 0, 1, integration_atol,
                                 integration_rtol,
                                 integration_sub_regions, integration_key,
                                 workspace, &results.at(1), &abserr);

    if (status != 0) {
        throw std::runtime_error("RSD tree-level integration l=2 failed \
                with error code" + std::to_string(status));
    }

    /* l=4 integration */
    auto integral_l4 = [&k, &ps](double mu) {
        return 0.125 * (35 * POW4(mu) - 30 * mu * mu + 3)
            * ps.tree_level(k, mu);
    };

    F.function = [] (double x, void* p) {
        return (*(decltype(integral_l4)*)p)(x);
    };
    F.params = &integral_l4;

    status = gsl_integration_qag(&F, 0, 1, integration_atol,
                                 integration_rtol,
                                 integration_sub_regions, integration_key,
                                 workspace, &results.at(2), &abserr);

    if (status != 0) {
        throw std::runtime_error("RSD tree-level integration l=4 failed \
                with error code" + std::to_string(status));
    }

    /* Multiply by two for integration mu=[-1,0] (integrands are even) */
    for (auto& el : results) el *= 2;
}



void tree_level(
    double k_a,
    Dynamics dynamics,
    const InputPowerSpectrum& ps,
    const EtaGrid& eta_grid,
    const EvolutionParameters& ev_params,
    const Vec1D<Pair<int>>& pair_correlations,
    Vec1D<double>& results /* out */
)
{
    if (ps.rsd()) {
        if (ps.ir_resum()) {
            size_t sub_regions = 10000;
            double atol = std::min(1e-10, ps(k_a,0));
            double rtol = 1e-6;
            int integration_key = GSL_INTEG_GAUSS61;

            ps::rsd_tree_level_ir_resum(k_a, ps, results, sub_regions, atol,
                    rtol, integration_key);
        }
        else {
            ps::rsd_tree_level(k_a, ps, results);
        }
    }
    else {
        for (auto& el : results) el = ps.tree_level(k_a, 0);
    }

    if (dynamics == EVOLVE_IC_ASYMP) {
        Vec1D<double> F1_eta_ini(COMPONENTS, 0);
        Vec1D<double> F1_eta_fin(COMPONENTS, 0);
        compute_F1(k_a, ev_params, eta_grid, F1_eta_ini,
                   F1_eta_fin);

        for (size_t i = 0; i < results.size(); ++i) {
            results.at(i) *=
                F1_eta_fin.at(static_cast<size_t>(
                    pair_correlations.at(i).first())) *
                F1_eta_fin.at(static_cast<size_t>(
                    pair_correlations.at(i).second()));
        }
    }
}
} /* namespace ps */


namespace bs {
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
        loop_params.args_2_kernel_index(arg_config.a().args.data());
    arg_config.b().kernel_index =
        loop_params.args_2_kernel_index(arg_config.b().args.data());
    arg_config.c().kernel_index =
        loop_params.args_2_kernel_index(arg_config.c().args.data());

    return arg_config;
}



void diagram_term(
        int diagram_idx,
        IntegrandTables& tables,
        const InputPowerSpectrum& ps,
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
                diagram_results.at(i) *= ps(k_a, 0);
                diagram_results.at(i) *= ps(k_c, 0);
              break;
            case 1:
                diagram_results.at(i) *= ps(k_b, 0);
                diagram_results.at(i) *= ps(k_c, 0);
                break;
            case 2:
                diagram_results.at(i) *= ps(k_a, 0);
                diagram_results.at(i) *= ps(k_b, 0);
                break;
            default:
                throw(std::logic_error(
                            "diagram_term(): got argument i which is not 0,1,2."));
        }
    }
}



void tree_level(
        IntegrandTables& tables,
        const InputPowerSpectrum& ps,
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
        diagram_term(i, tables, ps, triple_correlations, diagram_results);

        for (size_t j = 0; j < triple_correlations.size(); ++j) {
            results.at(j) += diagram_results.at(j);
        }
    }

    /* Zero-initialize kernel tables again, to ensure that values do not
     * corrupt futher computations */
    tables.reset();
}
} /* namespace bs */
