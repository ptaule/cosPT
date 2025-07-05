#include <cmath>
#include <iostream>
#include <stdexcept>

extern "C" {
    #include <cuba.h>
}

#include "../include/biased_tracers.hpp"
#include "../include/diagrams.hpp"
#include "../include/kernel_evolution.hpp"
#include "../include/ir_resum.hpp"
#include "../include/parameters.hpp"
#include "../include/rsd.hpp"
#include "../include/spt_kernels.hpp"
#include "../include/tables.hpp"
#include "../include/utilities.hpp"

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
    const Dynamics dynamics = tables.dynamics;

    const ArgumentConfiguration& arg_config_l =
        diagram.get_arg_config_l(rearr_idx, sign_idx);
    const ArgumentConfiguration& arg_config_r =
        diagram.get_arg_config_r(rearr_idx, sign_idx);

    int k_idx_l = arg_config_l.kernel_index;
    int k_idx_r = arg_config_r.kernel_index;
    const auto args_l = arg_config_l.args.data();
    const auto args_r = arg_config_r.args.data();

    int n_l = 2 * diagram.l + diagram.m;
    int n_r = 2 * diagram.r + diagram.m;

    /* If dynamics is EdS-SPT or the corresponding kernels are used for initial
     * conditions, compute the EdS-SPT kernels */
    if (dynamics == EDS_SPT || dynamics == EVOLVE_EDS_ICS) {
        compute_SPT_kernels(args_l, k_idx_l, n_l, tables);
        compute_SPT_kernels(args_r, k_idx_r, n_r, tables);
    }
    /* If dynamics is not EdS-SPT, solve general ODE system for kernels */
    if (dynamics != EDS_SPT) {
        KernelEvolver kernel_evolver(tables);
        kernel_evolver.compute(args_l, k_idx_l, n_l);
        kernel_evolver.compute(args_r, k_idx_r, n_r);
    }

    if (tables.biased_tracers) {
        compute_rsd_biased_kernels(args_l, k_idx_l, n_l, tables);
        compute_rsd_biased_kernels(args_r, k_idx_r, n_r, tables);
        compute_rsd_biased_kernels(arg_config_l.args.data(), arg_config_l.kernel_index,
                2 * diagram.l + diagram.m, tables);
        compute_rsd_biased_kernels(arg_config_r.args.data(), arg_config_r.kernel_index,
                2 * diagram.r + diagram.m, tables);

        /* Read values for two-point correlators.
         * If l != r, there are two diagrams corresponding to l <-> r */
        term_results.at(0) = tables.rsd_kernels
            .at(static_cast<size_t>(k_idx_l)).value
            * tables.rsd_kernels
            .at(static_cast<size_t>(k_idx_r)).value
            * (diagram.l == diagram.r ? 1 : 2);
        return;
    }
    /* Generic RSD computation (general loop order) */
    else if (tables.rsd) {
        compute_rsd_kernels(args_l, k_idx_l, n_l, tables);
        compute_rsd_kernels(args_r, k_idx_r, n_r, tables);

        /* Read values for two-point correlators.
         * If l != r, there are two diagrams corresponding to l <-> r */
        term_results.at(0) = tables.rsd_kernels
            .at(static_cast<size_t>(k_idx_l)).value
            * tables.rsd_kernels
            .at(static_cast<size_t>(k_idx_r)).value
            * (diagram.l == diagram.r ? 1 : 2);
        return;
    }

    auto get_kernel_values_ptr = [&](size_t k_idx) -> double* {
        if (dynamics == EDS_SPT) {
            return tables.spt_kernels.at(k_idx).values.data();  // std::array<double, N>
        } else {
            size_t last = tables.eta_grid.time_steps() - 1;
            auto& kernel_vec = tables.kernels.at(k_idx).values;
            return &kernel_vec(last, 0);
        }
    };

    /* Pointers to SPTKernel vector or last time step of Kernel vector */
    double* vals_l = get_kernel_values_ptr(static_cast<size_t>(k_idx_l));
    double* vals_r = get_kernel_values_ptr(static_cast<size_t>(k_idx_r));

    bool symmetric = diagram.l == diagram.r;
    for (size_t i = 0; i < pair_correlations.size(); ++i) {
        int fst = pair_correlations.at(i).first();
        int snd = pair_correlations.at(i).second();

        double result = vals_l[fst] * vals_r[snd];
        if (!symmetric)
            result += vals_l[snd] * vals_r[fst];

        term_results.at(i) = result;
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
    std::cout << Colors::MAGENTA << diagram.tags()
        << Colors::RESET << std::endl;
#endif
    size_t n_comp = diagram_results.size();

    // Loop over momentum rearrangement and sign flips
    for (size_t i = 0; i < diagram.n_rearrangements(); ++i) {
        for (size_t j = 0; j < diagram.n_sign_configs(); ++j) {
#if DEBUG >= 2
        std::cout << Colors::BLUE << diagram.argument_configuration(i, j)
            << Colors::RESET << std::endl;
#endif

            double q_m1 = diagram.q_m1(i, j,
                    tables.composite_dot_products());
            int heaviside_theta = diagram.heaviside_theta(q_m1, i,
                    tables.vars.magnitudes);

            /* We set P(k < q_min) and P(k > q_max) to zero, hence we
             * skip this configuration if q_m1 < qmin or q_m1 > q_max.
             * Furthermore if heaviside_theta == 0, we can also skip. */
            if (heaviside_theta == 0 || q_m1 < input.q_min ||
                q_m1 > input.q_max) {
#if DEBUG >= 2
                std::cout << "\t" << 0 << std::endl;
#endif
                continue;
            }

            Vec1D<double> term_results(n_comp, 0.0);
            configuration_term(diagram, i, j, input.pair_correlations, tables,
                    term_results);

            for (auto& el : term_results) {
                el *= heaviside_theta * input.ps(q_m1, tables.vars.mu_los);
            }
            for (size_t a = 0; a < n_comp; ++a) {
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

    for (size_t i = 0; i < n_comp; ++i) {
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
    const Dynamics dynamics = tables.dynamics;
    const Triple<ArgumentConfiguration>& arg_config =
        diagram.get_arg_config(rearr_idx, sign_idx, overall_loop_idx);

    int k_idx_a = arg_config.a().kernel_index;
    int k_idx_b = arg_config.b().kernel_index;
    int k_idx_c = arg_config.c().kernel_index;
    const auto args_a = arg_config.a().args.data();
    const auto args_b = arg_config.b().args.data();
    const auto args_c = arg_config.c().args.data();

    int n_a = 2 * diagram.n_a + diagram.n_ab + diagram.n_ca;
    int n_b = 2 * diagram.n_b + diagram.n_ab + diagram.n_bc;
    int n_c = 2 * diagram.n_c + diagram.n_bc + diagram.n_ca;

    if (dynamics == EDS_SPT || dynamics == EVOLVE_EDS_ICS) {
        compute_SPT_kernels(args_a, k_idx_a, n_a, tables);
        compute_SPT_kernels(args_b, k_idx_b, n_b, tables);
        compute_SPT_kernels(args_c, k_idx_c, n_c, tables);
    }
    if (dynamics != EDS_SPT) {
        KernelEvolver kernel_evolver(tables);

        kernel_evolver.compute(args_a, k_idx_a, n_a);
        kernel_evolver.compute(args_b, k_idx_b, n_b);
        kernel_evolver.compute(args_c, k_idx_c, n_c);
    }

    auto get_kernel_values_ptr = [&](size_t k_idx) -> double* {
        if (dynamics == EDS_SPT) {
            return tables.spt_kernels.at(k_idx).values.data();  // std::array<double, N>
        } else {
            size_t last = tables.eta_grid.time_steps() - 1;
            auto& kernel_vec = tables.kernels.at(k_idx).values;
            return &kernel_vec(last, 0);
        }
    };

    /* Pointers to SPTKernel vector or last time step of Kernel vector */
    double* vals_a = get_kernel_values_ptr(static_cast<size_t>(k_idx_a));
    double* vals_b = get_kernel_values_ptr(static_cast<size_t>(k_idx_b));
    double* vals_c = get_kernel_values_ptr(static_cast<size_t>(k_idx_c));

    for (size_t i = 0; i < triple_correlations.size(); ++i) {
        term_results[i] = vals_a[triple_correlations.at(i).first()] *
                          vals_b[triple_correlations.at(i).second()] *
                          vals_c[triple_correlations.at(i).third()];
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
    std::cout << Colors::MAGENTA << diagram.tags()
        << Colors::RESET << std::endl;
#endif
    size_t n_comp = diagram_results.size();

    /* Loop over momentum rearrangement, sign flips and overall loop assosiations */
    size_t overall_loop_assosiations = diagram.overall_loop() ? 6 : 1;
    for (size_t i = 0; i < diagram.n_rearrangements(); ++i) {
        for (size_t j = 0; j < diagram.n_sign_configs(); ++j) {
            for (size_t k = 0; k < overall_loop_assosiations; ++k) {
#if DEBUG >= 2
                std::cout << Colors::BLUE <<
                    diagram.argument_configuration(i, j,
                                                   k)
                    << Colors::RESET << std::endl;
#endif
                Triple<double> q_xy1(0,0,0);
                int heaviside_theta = 1;
                diagram.connecting_lines_factors(i, j, k,
                        tables.vars.magnitudes,
                        tables.composite_dot_products(),
                        q_xy1, heaviside_theta);

                if (heaviside_theta == 0) {
#if DEBUG >= 2
                    std::cout << "\t" << 0 << std::endl;
#endif
                    continue;
                }

                Vec1D<double> term_results(n_comp, 0.0);
                configuration_term(diagram, i, j, k,
                        input.triple_correlations, tables, term_results);

                /* Multiply by P_lin(q_xy1) if xy connecting line is present.
                 * In addition, if q_xy1 is outside integration range [q_min,
                 * q_max], we set P_lin to zero, i.e. no contribution from this
                 * term. */
                if (diagram.n_ab > 0) {
                    if (q_xy1.a() < input.q_min || q_xy1.a() > input.q_max) {
#if DEBUG >= 2
                        std::cout << "\t" << 0 << std::endl;
#endif
                        continue;
                    }
                    for (auto& el : term_results) {
                        el *= input.ps(q_xy1.a(), tables.vars.mu_los);
                    }
                }
                if (diagram.n_bc > 0) {
                    if (q_xy1.b() < input.q_min || q_xy1.b() > input.q_max) {
#if DEBUG >= 2
                        std::cout << "\t" << 0 << std::endl;
#endif
                        continue;
                    }
                    for (auto& el : term_results) {
                        el *= input.ps(q_xy1.b(), tables.vars.mu_los);
                    }
                }
                if (diagram.n_ca > 0) {
                    if (q_xy1.c() < input.q_min || q_xy1.c() > input.q_max) {
#if DEBUG >= 2
                        std::cout << "\t" << 0 << std::endl;
#endif
                        continue;
                    }
                    for (auto& el : term_results) {
                        el *= input.ps(q_xy1.c(), tables.vars.mu_los);
                    }
                }
                /* For overall-loop diagrams: divide by pLin(rearr(Q)) which is
                 * multiplied also in bs::integrand() */
                if (diagram.overall_loop()) {
                    for (auto& el : term_results) {
                        el /= input.ps(
                                diagram.q1_magnitude(i, tables.vars.magnitudes),
                                tables.vars.mu_los
                                );
                    }
                }
                for (auto& el : term_results) {
                    el *= heaviside_theta;
                }
                for (size_t a = 0; a < n_comp; ++a) {
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

    for (size_t i = 0; i < n_comp; ++i) {
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
        const int *ndim,
        const cubareal xx[],
        const int *ncomp,
        cubareal ff[],
        void *userdata,
        const int *nvec,
        const int *core
        )
{
    UNUSED(ndim);
    UNUSED(ncomp);
    UNUSED(nvec);

    IntegrationInput& input = *static_cast<IntegrationInput*>(userdata);

    /*  For thread <*core + 1> (index 0 is reserved for master), we use the */
    /*  IntegrandTables number *core+1 */
    IntegrandTables& tables = input.tables_vec.at(static_cast<size_t>(*core + 1));

    size_t n_comp = static_cast<size_t>(*ncomp);
    int n_loops = tables.loop_structure.n_loops();
    IntegrationVariables& vars = tables.vars;

    Spectrum spectrum = tables.loop_structure.spectrum();
    bool rsd = tables.rsd;

#if DEBUG > 1
    /* Check that n_loops > 0 */
    if (n_loops < 1) {
        throw std::runtime_error("n_loops is not larger than or equal to 1.");
    }
#endif

    double ratio = input.q_max/input.q_min;
    double log_ratio = std::log(ratio);
    double jacobian = 1.0;

    int idx = 0; /* Index for reading xx-array */
    double Q1_xx = 0; /* Temp storage of xx-element for Q1 if not single_hard_limit */

    if (rsd) {
        vars.mu_los = xx[idx++];
    }
    /* 1-loop */
    if (!input.single_hard_limit) {
        Q1_xx = xx[idx++];
        vars.magnitudes.at(0) = input.q_min * pow(ratio, Q1_xx);
        jacobian *= log_ratio * CUBE(vars.magnitudes.at(0));
    }

    vars.cos_theta.at(0) = xx[idx++];
    /* For bispectrum or RSD we cannot fix phi[0] = 0 */
    if (spectrum == BISPECTRUM || rsd) {
        vars.phi.at(0) = xx[idx++] * TWOPI;
        jacobian *= TWOPI;
    }

    /* 2-loop */
    if (n_loops > 1){
        if (!input.single_hard_limit) {
            vars.magnitudes.at(1) = input.q_min *
                pow(vars.magnitudes.at(0) / input.q_min, xx[idx++]);
            jacobian *= Q1_xx;
        }
        else {
            vars.magnitudes.at(1) = input.q_min * pow(ratio,xx[idx++]);
        }
        jacobian *= log_ratio * CUBE(vars.magnitudes.at(1));
        vars.cos_theta.at(1) = xx[idx++];
        vars.phi.at(1) = xx[idx++] * TWOPI;
        jacobian *= TWOPI;
    }

    Vec1D<double> results(n_comp, 0.0);
    Vec1D<double> diagram_results;

    /* For RSD, we only compute the [0,0] (or [0,0,0] for bispectrum)
     * correlation. Therefore we only need diagram_results of size 1 */
    if (rsd) {
        diagram_results.resize(1);
    }
    else {
        diagram_results.resize(n_comp);
    }

    try {
        /* Zero-initialize kernel tables */
        tables.reset();
        // Compute dot_products-, alpha- and beta-tables
        tables.compute_tables();

        // Loop over all diagrams
        if (spectrum == POWERSPECTRUM) {
            for (auto& diagram : input.ps_diagrams) {
                ps::diagram_term(diagram, input, tables, diagram_results);

                /* Add diagram results to results and reset diagram results */
                for (size_t j = 0; j < n_comp; ++j) {
                    /* If RSD, only one element of diagram_results (index 0),
                     * if not use j */
                    size_t idx = rsd ? 0 : j;
                    results.at(j) += diagram_results.at(idx);
                }
                std::fill(diagram_results.begin(), diagram_results.end(), 0.0);
            }
        }
        else {
            /* Bispectrum */
            for (auto& diagram : input.bs_diagrams) {
                bs::diagram_term(diagram, input, tables, diagram_results);

                /* Add diagram results to results and reset diagram results */
                for (size_t j = 0; j < n_comp; ++j) {
                    results.at(j) += diagram_results.at(j);
                }
                std::fill(diagram_results.begin(), diagram_results.end(), 0.0);
            }
        }

        int i = 0;
        /* Skip P_lin(Q1) if single_hard_limit = true */
        if (input.single_hard_limit) ++i;
        for (; i < tables.loop_structure.n_loops(); ++i) {
            for (auto& el : results) {
                el *= input.ps(
                        vars.magnitudes.at(static_cast<size_t>(i)),
                        vars.mu_los
                        );
            }
        }

        if (tables.biased_tracers) {
            /* We define b2 such that the constant contribution to Id2d2 as k->
             * 0 is removed */
            /* TODO: atm only implemented for n_loops = 1 */
            /* Compute F1(k) (which is 1 for EdS-SPT) */
            double F1 = 1;
            if (tables.dynamics != EDS_SPT) {
                KernelEvolver kernel_evolver(tables);

                Vec1D<int> config(tables.loop_structure.n_coeffs(), 0);
                if (config.empty()) {
                    throw std::runtime_error("integrand(): n_coeffs = 0");
                }
                config.at(0) = 1; /* Config for q1 */

                Vec1D<int> arguments(tables.loop_structure.n_kernel_args(),
                    tables.loop_structure.zero_label());
                arguments.at(0) = config2label(config);

                int kernel_index =
                    kernel_evolver.compute(arguments.data(),
                                           -1, 1);
                size_t k_idx_t = static_cast<size_t>(kernel_index);
                size_t last = tables.eta_grid.time_steps() - 1;
                F1 =
                    tables.kernels.at(k_idx_t).values(last, 0);
            }

            double b2_subtract = 0.5 *
                SQUARE(F1) *
                SQUARE(tables.bias_parameters.at(1)) *
                SQUARE(input.ps(vars.magnitudes.at(0), vars.mu_los));
            results.at(0) -= b2_subtract;
            results.at(1) -= b2_subtract;
            results.at(2) -= b2_subtract;
        }
        /* RSD multipoles: multiply with Legendre polynomials */
        if (rsd) {
            results.at(1) *= 0.5 * (3 * SQUARE(vars.mu_los) - 1);
            results.at(2) *= 0.125 * (
                    35 * POW4(vars.mu_los) - 30 * SQUARE(vars.mu_los) + 3
                    );
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        /* Tell CUBA an error occured */
        return -999;
    }

    for (size_t i = 0; i < n_comp; ++i) {
        ff[i] = results.at(i) * jacobian;
    }

    /* Return success */
    return 0;
}
