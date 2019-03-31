/*
   integrand.c

   Created by Petter Taule on 24.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_combination.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/kernels.h"
#include "../include/spt_kernels.h"
#include "../include/integrand.h"



inline int diagram_factor(const diagram_t* diagram) {
    short int l = diagram->l;
    short int r = diagram->r;
    short int m = diagram->m;

    int numerator = gsl_sf_fact(2*l + m) * gsl_sf_fact(2*r + m);
    int denominator = pow(2,l+r) * gsl_sf_fact(l) * gsl_sf_fact(r) * gsl_sf_fact(m);

    return numerator/denominator;
}



inline int symmetrization_factor(const diagram_t* diagram) {
    short int l = diagram->l;
    short int r = diagram->r;
    short int m = diagram->m;

    int numerator = gsl_sf_fact(LOOPS) * pow(2,m-1);
    int denominator = gsl_sf_fact(m-1) * gsl_sf_fact(l) * gsl_sf_fact(r);

    return numerator/denominator;
}



// Find (distinct) diagrams for L-loop
// They satisfy: m >= 1; l,r > 0; l + r + m = L + 1
void possible_diagrams(diagram_t diagrams[]) {
    short int m = 0;

    size_t index = 0;

    for (m = 1; m <= LOOPS + 1; ++m) {
        short int l = LOOPS + 1 - m;
        short int r = 0;
        while (l >= r) {
            if (index >= N_DIAGRAMS)
                warning_verbose("Index out of bounds, index = %ld is larger "
                        "than N_DIAGRAMS = %d.", index, N_DIAGRAMS);

            diagrams[index].l = l;
            diagrams[index].r = r;
            diagrams[index].m = m;

            l = LOOPS + 1 - m - (++r);
            index++;
        };
    }
}



void find_kernel_arguments(
        const diagram_t* diagram,        /* in, diagram to compute          */
        const short int rearrangement[], /* in, loop momenta arrangement    */
        const short int signs[],         /* in, signs of "connection" loops */
        short int arguments_l[],         /* out, left kernel arguments      */
        short int arguments_r[]          /* out, right kernel arguments     */
        )
{
    short int l = diagram->l;
    short int r = diagram->r;
    short int m = diagram->m;

    // First argument is on the form k1 = k - k2 - k3 - ... - km
    short int config[N_COEFFS] = {};
    config[N_COEFFS - 1] = 1; // k-coefficient is 1
    for (size_t i = 2; i <= m; ++i) {
        config[rearrangement[i-2]] = - signs[i-2];
    }
    arguments_l[0] = config2label(config,N_COEFFS);
    arguments_r[0] = config2label(config,N_COEFFS);

    // Reset config
    memset(config,0,sizeof(config));

    // Argument indices
    size_t index_l = 1;
    size_t index_r = 1;

    // k2,k3,...,km arguments, the ordering of which is stored by the first
    // (m-1) entries of rearrangement[]
    for (size_t i = 2; i <= m; ++i) {
        config[rearrangement[i-2]] = signs[i-2];
        arguments_l[index_l++] = config2label(config,N_COEFFS);
        arguments_r[index_r++] = config2label(config,N_COEFFS);
        memset(config,0,sizeof(config));
    }

    // l-loop arguments
    for (size_t i = 0; i < l; ++i) {
        short int loop_momentum_index = rearrangement[i + m - 1];
        config[loop_momentum_index] = 1;
        arguments_l[index_l++] = config2label(config,N_COEFFS);
        config[loop_momentum_index] = -1;
        arguments_l[index_l++] = config2label(config,N_COEFFS);

        memset(config,0,sizeof(config));
    }
    // r-loop arguments
    for (size_t i = 0; i < r; ++i) {
        short int loop_momentum_index = rearrangement[i + m - 1 + l];
        config[loop_momentum_index] = 1;
        arguments_r[index_r++] = config2label(config,N_COEFFS);
        config[loop_momentum_index] = -1;
        arguments_r[index_r++] = config2label(config,N_COEFFS);

        memset(config,0,sizeof(config));
    }

    // Fill remaining spots with zero-label
    while (index_l < N_KERNEL_ARGS) {
        arguments_l[index_l++] = ZERO_LABEL;
    }
    while (index_r < N_KERNEL_ARGS) {
        arguments_r[index_r++] = ZERO_LABEL;
    }
}



static vfloat compute_k1(
        short int m,
        const short int rearrangement[],
        const short int signs[],
        const vfloat bare_scalar_products[][N_COEFFS]
        )
{
    vfloat k1 = bare_scalar_products[N_COEFFS - 1][N_COEFFS - 1];
    for (int i = 2; i <= m; ++i) {
        int index = rearrangement[i-2];
        k1 += bare_scalar_products[index][index];
        // Note that elements in the signs-array correspond to rearranged loop
        // momenta, thus we use <i-2>, no <index> as index here
        k1 -= 2 * signs[i-2] * bare_scalar_products[N_COEFFS - 1][index];
    }

    for (int i = 3; i <= m; ++i) {
        for (int j = 2; j < i; ++j) {
            k1 += 2 * signs[i-2] * signs[j-2] *
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



// For debuggin purposes
__attribute__((unused))
static void print_integrand_info(
        const short int arguments_l[N_KERNEL_ARGS],
        const short int arguments_r[N_KERNEL_ARGS],
        const diagram_t* diagram
        )
{
    printf(ANSI_COLOR_MAGENTA "(m,l,r) = (%d,%d,%d)\t",diagram->m,diagram->l,diagram->r);
    printf(ANSI_COLOR_BLUE "F%d(k",2*diagram->l + diagram->m);
    short int config[N_COEFFS];
    label2config(arguments_l[0],config,N_COEFFS);
    for (int i = 0; i < LOOPS; ++i) {
        if (config[i] == 0) continue;
        else if (config[i] == -1) printf("-Q%d",i+1);
        else if (config[i] == 1)  printf("+Q%d",i+1);
    }
    printf(", ");

    for (int i = 1; i < N_KERNEL_ARGS; ++i) {
        if (arguments_l[i] == ZERO_LABEL) break;
        label2config(arguments_l[i],config,N_COEFFS);
        for (int j = 0; j < LOOPS; ++j) {
            if (config[j] == 0) continue;
            else if (config[j] == -1) printf("-Q%d, ",j+1);
            else if (config[j] == 1)  printf("+Q%d, ",j+1);
        }
    }
    printf(") * F%d(k",2*diagram->r + diagram->m);
    label2config(arguments_r[0],config,N_COEFFS);
    for (int i = 0; i < LOOPS; ++i) {
        if (config[i] == 0) continue;
        else if (config[i] == -1) printf("-Q%d",i+1);
        else if (config[i] == 1)  printf("+Q%d",i+1);
    }
    printf(", ");

    for (int i = 1; i < N_KERNEL_ARGS; ++i) {
        if (arguments_l[i] == ZERO_LABEL) break;
        label2config(arguments_r[i],config,N_COEFFS);
        for (int j = 0; j < LOOPS; ++j) {
            if (config[j] == 0) continue;
            else if (config[j] == -1) printf("-Q%d, ",j+1);
            else if (config[j] == 1)  printf("+Q%d, ",j+1);
        }
    }
    printf(")" ANSI_COLOR_RESET);
}



static vfloat integrand_term(
        const short int arguments_l[N_KERNEL_ARGS],
        const short int arguments_r[N_KERNEL_ARGS],
        const diagram_t* diagram,
        const integration_input_t* input,
        const table_pointers_t* data_tables
        )
{
    short int m = diagram->m;
    short int l = diagram->l;
    short int r = diagram->r;

    vfloat result = 1;
    // If right kernel has only one argument, its value is 1
    if (m == 1 && r == 0) {
        result *= compute_SPT_kernel(arguments_l, 2*l + m, input->component_a, data_tables);
    }
    // If there are no "self" loops, and the components to compute are
    // equal, the kernels are equal
    else if (l == 0 && r == 0 && input->component_a == input->component_b) {
        result *= pow(compute_SPT_kernel(arguments_l, m, input->component_a, data_tables) ,2);
    // In DEBUG-mode, check that kernel arguments in fact are equal in this
    // case
#if DEBUG >= 1
        for (int i = 0; i < N_KERNEL_ARGS; ++i) {
            if (arguments_l[i] != arguments_r[i])
                warning("Arguments l & r were wrongly assumed equal.");
        }
#endif
    }
    else {
        result *= compute_SPT_kernel(arguments_l, 2*l + m, input->component_a, data_tables)
                * compute_SPT_kernel(arguments_r, 2*r + m, input->component_b, data_tables);
    }
    return result;
}



vfloat sign_flip_symmetrization(
        const short int rearrangement[],
        const diagram_t* diagram,
        const integration_input_t* input,
        const table_pointers_t* data_tables
        )
{
    short int m = diagram->m;

    // Signs of k1...km momenta (only (m-1) first elements relevant)
    // Start with +1,+1,... signs
    short int signs[LOOPS];
    for (int i = 0; i < LOOPS; ++i) signs[i] = 1;

    // Arguments for "left" and "right" kernels (changed in every
    // symmetrization operation)
    short int arguments_l[N_KERNEL_ARGS];
    short int arguments_r[N_KERNEL_ARGS];

    vfloat result = 0;
    // Loop over possible sign flips, and evaluate integrand
    for (int i = 0; i < pow(2,m-1); ++i) {
        for (int j = 0; j < (m-1); ++j) {
            if ((i/(j+1) % 2) == 1) {
                signs[j] *= -1;
                break;
            }
        }

        find_kernel_arguments(diagram, rearrangement, signs, arguments_l,
                arguments_r);
#if DEBUG >= 2
        print_integrand_info(arguments_l, arguments_r, diagram);
#endif
        vfloat k1 = compute_k1(diagram->m, rearrangement, signs,
                data_tables->bare_scalar_products);
        int h_theta = heaviside_theta(diagram->m, k1, rearrangement,
                data_tables->Q_magnitudes);
        if (h_theta == 0) {
#if DEBUG >= 2
        printf("\t\t=> partial result = 0\n");
#endif
            continue;
        }

        vfloat partial_result = h_theta * gsl_spline_eval(input->spline,k1,input->acc);
        partial_result *= integrand_term(arguments_l, arguments_r, diagram,
                input, data_tables);
#if DEBUG >= 2
        printf("\t\t=> partial result = " vfloat_fmt "\n",partial_result);
#endif
        result += partial_result;
    }
    return result;
}



vfloat loop_momenta_symmetrization(
        const diagram_t* diagram,
        const integration_input_t* input,
        const table_pointers_t* data_tables
        )
{
    short int l = diagram->l;
    short int r = diagram->r;
    short int m = diagram->m;

#if DEBUG >= 1
    if ((m + l + r) != (LOOPS + 1)) warning("m + r + l != LOOPS + 1.");
#endif

    vfloat result = 0;

    // Tag loop momenta by number for symmetrization
    short int loop_momenta[LOOPS];
    for (int i = 0; i < LOOPS; ++i) loop_momenta[i] = i;
    short int rearrangement[LOOPS] = {};

    // Loop momenta combinations:
    // m-1         "connection" loops
    // s = L-(m-1) "self" loops, divided into:
    //      l left "self" loops
    //      r right "self" loops
    gsl_combination *comb_m, *comb_s, *comb_l, *comb_r;

    comb_m = gsl_combination_alloc(LOOPS, (m-1));
    comb_s = gsl_combination_alloc(LOOPS, LOOPS - (m-1));
    gsl_combination_init_first(comb_m);
    gsl_combination_init_last(comb_s);

    // We only need l and r groupings when "self"-loops are present
    if (m < LOOPS + 1) {
        comb_l = gsl_combination_alloc(LOOPS - (m-1), l);
        comb_r = gsl_combination_alloc(LOOPS - (m-1), r);
    }

    // Go through possible (m-1)-groupings
    do {
        for (int i = 0; i < m-1; ++i) {
            rearrangement[i] = loop_momenta[gsl_combination_get(comb_m,i)];
        }
        if (m < LOOPS + 1) {
            gsl_combination_init_first(comb_l);
            gsl_combination_init_last(comb_r);

            // Go through possible l and r groupings
            do {
                for (int i = 0; i < l; ++i) {
                    int j = m - 1 + i;
                    rearrangement[j] = loop_momenta[
                        gsl_combination_get(comb_s,gsl_combination_get(comb_l,i))];
                }
                for (int i = 0; i < r; ++i) {
                    int j = m - 1 + l + i;
                    rearrangement[j] = loop_momenta[
                        gsl_combination_get(comb_s,gsl_combination_get(comb_r,i))];
                }
                result += sign_flip_symmetrization(rearrangement, diagram,
                        input, data_tables);
            } while (
                    gsl_combination_next(comb_l) == GSL_SUCCESS &&
                    gsl_combination_prev(comb_r) == GSL_SUCCESS
                    );
        }
        else {
            result += sign_flip_symmetrization(rearrangement, diagram, input,
                    data_tables);
        }
    } while (
            gsl_combination_next(comb_m) == GSL_SUCCESS &&
            gsl_combination_prev(comb_s) == GSL_SUCCESS
            );

    if (m < LOOPS + 1) {
        gsl_combination_free(comb_l);
        gsl_combination_free(comb_r);
    }
    gsl_combination_free(comb_m);
    gsl_combination_free(comb_s);

    return result;
}



vfloat integrand(
        const integration_input_t* input,
        const integration_variables_t* vars
        )
{
    // Store pointers to the computed tables in struct for convenience
    table_pointers_t data_tables;
    data_tables.Q_magnitudes = vars->magnitudes,
    data_tables.alpha = matrix_alloc(N_CONFIGS,N_CONFIGS);
    data_tables.beta  = matrix_alloc(N_CONFIGS,N_CONFIGS);

    // Allocate space for kernels (calloc also initializes values to 0)
    data_tables.kernels = (kernel_value_t*)
        calloc(COMPONENTS * N_KERNELS, sizeof(kernel_value_t));


    // Initialize sum-, bare_scalar_products-, alpha- and beta-tables
    compute_sum_table(data_tables.sum_table);
    compute_bare_scalar_products(input->k, vars,
            data_tables.bare_scalar_products);
    compute_alpha_beta_tables(data_tables.bare_scalar_products,
            data_tables.alpha, data_tables.beta);

    diagram_t diagrams[N_DIAGRAMS];
    possible_diagrams(diagrams);

    // Loop over possible diagrams at this loop order
    vfloat result = 0;
    for (int i = 0; i < N_DIAGRAMS; ++i) {
        vfloat diagram_result =
            loop_momenta_symmetrization(&(diagrams[i]),input,&data_tables);
        // Multiply and divide by symmetrization & diagram factors
        diagram_result /= symmetrization_factor(&(diagrams[i]));
        diagram_result *= diagram_factor(&(diagrams[i]));
        // If diagram is antisymmetric in l <-> r, multiply by 2 (the algorithm
        // assumes l >= r)
        if (diagrams[i].l != diagrams[i].r) {
            diagram_result *= 2;
        }
        result += diagram_result;
    }

    for (int i = 0; i < LOOPS; ++i) {
        result *= gsl_spline_eval(input->spline, data_tables.Q_magnitudes[i],
                input->acc);
    }

    // Free allocated memory
    free(data_tables.kernels);
    matrix_free(data_tables.alpha);
    matrix_free(data_tables.beta);

    return result;
}
