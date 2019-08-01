/*
   integrand.c

   Created by Petter Taule on 24.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_spline.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/tables.h"
#include "../include/spt_kernels.h"
#include "../include/diagrams.h"
#include "../include/integrand.h"



static vfloat compute_k1(
        short int m,
        const short int rearrangement[],
        const bool signs[],
        const vfloat bare_scalar_products[][N_COEFFS]
        )
{
    vfloat k1 = bare_scalar_products[N_COEFFS - 1][N_COEFFS - 1];
    for (int i = 2; i <= m; ++i) {
        int index = rearrangement[i-2];
        k1 += bare_scalar_products[index][index];
        // Note that elements in the signs-array correspond to rearranged loop
        // momenta, thus we use 'i-2', not 'index' as index here
        k1 -= 2 * (signs[i-2] ? 1 : -1) * bare_scalar_products[N_COEFFS - 1][index];
    }

    for (int i = 3; i <= m; ++i) {
        for (int j = 2; j < i; ++j) {
            k1 += 2 * (signs[i-2] ? 1 : -1) * (signs[j-2] ? 1 : -1) *
                bare_scalar_products[rearrangement[i-2]][rearrangement[j-2]];
        }
    }

    return sqrt(k1);
}



inline static int heaviside_theta(
        short int m,
        vfloat k1,
        const short int rearrangement[],
        const vfloat Q_magnitudes[]
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
            warning_verbose("Heaviside theta: Q%d < Q%d\n",
                    rearrangement[i-3] + 1,
                    rearrangement[i-2] + 1);
    }
#endif
    return m;
}



static vfloat integrand_term(
        short int diagram_index,
        short int rearrangement_index,
        short int sign_config_index,
        const integration_input_t* input,
        const table_ptrs_t* tables
        )
{
    // Shorthand variables/aliases for convenience
    short int m = input->diagrams[diagram_index].m;
    short int l = input->diagrams[diagram_index].l;
    short int r = input->diagrams[diagram_index].r;

    const short int* const rearrangement =
        input->diagrams[diagram_index].rearrangements[rearrangement_index];
    const bool* const signs =
        input->diagrams[diagram_index].sign_configs[sign_config_index];
    const short int* const arguments_l =
        input->diagrams[diagram_index].argument_configs_l[rearrangement_index][sign_config_index];
    const short int* const arguments_r =
        input->diagrams[diagram_index].argument_configs_r[rearrangement_index][sign_config_index];
    short int kernel_index_l =
        input->diagrams[diagram_index].kernel_indices_l[rearrangement_index][sign_config_index];
    short int kernel_index_r =
        input->diagrams[diagram_index].kernel_indices_r[rearrangement_index][sign_config_index];

#if DEBUG >= 2
    print_integrand_info(m,l,r,arguments_l, arguments_r);
#endif

    vfloat k1 = compute_k1(m, rearrangement, signs,
            tables->bare_scalar_products);
    int h_theta = heaviside_theta(m, k1, rearrangement,
            tables->Q_magnitudes);

    if (h_theta == 0) {
#if DEBUG >= 2
        printf("\t\t=> partial result = 0\n");
#endif
        return 0;
    }

    vfloat result = h_theta * gsl_spline_eval(input->spline,k1,input->acc);

    // If right kernel only has one argument, its value is 1
    if (m == 1 && r == 0) {
        result *= compute_SPT_kernel(arguments_l, kernel_index_l, 2*l + m,
                input->component_a, tables);
    }
    // If there are no "self" loops, and the components to compute are
    // equal, the kernels are equal
    else if (l == 0 && r == 0 && input->component_a == input->component_b) {
        result *= pow(compute_SPT_kernel(arguments_l, kernel_index_l, m,
                    input->component_a, tables) ,2);

    // In DEBUG-mode, check that kernel arguments in fact are equal in this
    // case
#if DEBUG >= 1
        if (kernel_index_l != kernel_index_r) {
            warning("Arguments l & r were wrongly assumed equal.");
        }
        for (int i = 0; i < N_KERNEL_ARGS; ++i) {
            if (arguments_l[i] != arguments_r[i])
                warning("Arguments l & r were wrongly assumed equal.");
        }
#endif
    }
    else {
        result *= compute_SPT_kernel(arguments_l, kernel_index_l, 2*l + m,
                input->component_a, tables)
                * compute_SPT_kernel(arguments_r, kernel_index_r, 2*r + m,
                    input->component_b, tables);
    }

#if DEBUG >= 2
        printf("\t\t=> Partial result = " vfloat_fmt "\n",result);
#endif
    return result;
}



vfloat integrand(
        const integration_input_t* input,
        table_ptrs_t* tables
        )
{
    // Loop over possible diagrams at this loop order
    vfloat result = 0;
    for (int i = 0; i < N_DIAGRAMS; ++i) {
        // Pointer alias for convenience
        const diagram_t* const dg = &(input->diagrams[i]);
        vfloat diagram_result = 0;

        for (int a = 0; a < dg->n_rearrangements; ++a) {
            for (int b = 0; b < dg->n_sign_configs; ++b) {
                diagram_result += integrand_term(i,a,b,input,tables);
            }
        }

        // Multiply and divide by symmetrization & diagram factors
        diagram_result *= dg->diagram_factor;
        diagram_result /= dg->n_rearrangements * dg->n_sign_configs;
        result += diagram_result;
    }

    for (int i = 0; i < LOOPS; ++i) {
        result *= gsl_spline_eval(input->spline, tables->Q_magnitudes[i],
                input->acc);
    }

    return result;
}
