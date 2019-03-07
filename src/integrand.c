/*
   integrand.c

   Created by Petter Taule on 24.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <math.h>
#include <string.h>

#include <cuba.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_combination.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/kernels.h"
#include "../include/spt_kernels.h"
#include "../include/integrand.h"



int diagram_factor(const diagram_t* diagram) {
    short int l = diagram->l;
    short int r = diagram->r;
    short int m = diagram->m;

    int numerator = gsl_sf_fact(2*l + m) * gsl_sf_fact(2*r + m);
    int denominator = pow(2,l+r) * gsl_sf_fact(l) * gsl_sf_fact(r) * gsl_sf_fact(m);

    return numerator/denominator;
}



int integrand_symmetrization_factor(const diagram_t* diagram) {
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
        config[rearrangement[i-2]] = signs[i-2];
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
    // l-loop arguments
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



vfloat compute_k1(short int m, const vfloat bare_scalar_products[][N_COEFFS]) {
    vfloat k1 = bare_scalar_products[N_COEFFS - 1][N_COEFFS - 1];
    for (int i = 2; i <= m; ++i) {
        k1 += bare_scalar_products[i-2][i-2];
        k1 -= 2 * bare_scalar_products[N_COEFFS - 1][i-2];
    }

    for (int i = 2; i <= m; ++i) {
        for (int j = 2; j < i; ++j) {
            k1 += 2 * bare_scalar_products[i-2][j-2];
        }
    }

    return sqrt(k1);
}



int heaviside_theta(short int m, vfloat k1, const vfloat Q_magnitudes[]) {
    // Heaviside-theta (k1 - k2)
    if (m == 2) {
        if (k1 <= Q_magnitudes[0]) return 0;
        return 2;
    }

    // Heaviside-theta (k2 - k3) etc.
    // Note that (assuming m >= 2), k2 = Q_magnitudes[0] etc.
    for (int i = 3; i <= m; ++i) {
        if (Q_magnitudes[i-3] <= Q_magnitudes[i-2]) return 0;
    }
    return gsl_sf_fact(m);
}



/* vfloat integrand_term( */
/*         const diagram_t* diagram, */
/*         const integration_input_t* data, */
/*         const vfloat Q_magnitudes[], */
/*         const vfloat bare_scalar_products[][N_COEFFS], */
/*         const short int arguments_l[N_KERNEL_ARGS], */
/*         const short int arguments_r[N_KERNEL_ARGS], */
/*         const matrix_vfloat* alpha, */
/*         const matrix_vfloat* beta, */
/*         kernel_value_t* kernels */
/*         ) */
/* { */
/*     short int l = diagram->l; */
/*     short int r = diagram->r; */
/*     short int m = diagram->m; */

/*     vfloat k1 = compute_k1(m,bare_scalar_products); */

/*     gsl_spline* spline    = data->spline; */
/*     gsl_interp_accel* acc = data->acc; */

/*     vfloat result = heaviside_theta(m,k1,Q_magnitudes) * interpolate(k1) */
/*         * compute_SPT_kernel(arguments_l, data->component_a, alpha, beta, kernels) */
/*         * compute_SPT_kernel(arguments_r, data->component_b, alpha, beta, kernels); */

/*     // If diagram is symmetric (l == r), mutliply by 2 */
/*     if (l == r) result *= 2; */

/*     return result; */
/* } */



void symmetrized_integrand(
        const diagram_t* diagram
        /* const integration_input_t* data, */
        /* const vfloat Q_magnitudes[], */
        /* const vfloat bare_scalar_products[][N_COEFFS], */
        /* const matrix_vfloat* alpha, */
        /* const matrix_vfloat* beta, */
        /* kernel_value_t* kernels */
        )
{
    short int l = diagram->l;
    short int r = diagram->r;
    short int m = diagram->m;

#if DEBUG
    if ((m + l + r) != (LOOPS + 1)) error("m + r + l != LOOPS + 1");
#endif

    // Tag loop momenta by number for symmetrization
    short int loop_momenta[LOOPS];
    for (int i = 0; i < LOOPS; ++i) loop_momenta[i] = i;
    short int rearrangement[LOOPS] = {};

    // Signs of k1...km momenta (only (m-1) first elements relevant)
    short int signs[LOOPS];
    for (int i = 0; i < LOOPS; ++i) signs[i] = 1;

    // Arguments for "left" and "right" kernels (changed in every
    // symmetrization operation)
    short int arguments_l[N_KERNEL_ARGS];
    short int arguments_r[N_KERNEL_ARGS];

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

                find_kernel_arguments(diagram,
                        rearrangement, signs, arguments_l, arguments_r);
                printf("L: ");
                for (int i = 0; i < N_KERNEL_ARGS; ++i) {
                    short int config[N_COEFFS] = {};
                    label2config(arguments_l[i], config, N_COEFFS);
                    printf("[");
                    for (int j = 0; j < N_COEFFS; ++j) {
                        printf("%d, ",config[j]);
                    }
                    printf("]");
                }
                printf("\nR: ");
                for (int i = 0; i < N_KERNEL_ARGS; ++i) {
                    short int config[N_COEFFS] = {};
                    label2config(arguments_r[i], config, N_COEFFS);
                    printf("[");
                    for (int j = 0; j < N_COEFFS; ++j) {
                        printf("%d, ",config[j]);
                    }
                    printf("]");
                }
                printf("\n");

            } while (
                    gsl_combination_next(comb_l) == GSL_SUCCESS &&
                    gsl_combination_prev(comb_r) == GSL_SUCCESS
                    );
            printf("\n");
        }
        else {
            find_kernel_arguments(diagram,
                    rearrangement, signs, arguments_l, arguments_r);
            printf("L: ");
            for (int i = 0; i < N_KERNEL_ARGS; ++i) {
                short int config[N_COEFFS] = {};
                label2config(arguments_l[i], config, N_COEFFS);
                printf("[");
                for (int j = 0; j < N_COEFFS; ++j) {
                    printf("%d, ",config[j]);
                }
                printf("]");
            }
            printf("\nR: ");
            for (int i = 0; i < N_KERNEL_ARGS; ++i) {
                short int config[N_COEFFS] = {};
                label2config(arguments_r[i], config, N_COEFFS);
                printf("[");
                for (int j = 0; j < N_COEFFS; ++j) {
                    printf("%d, ",config[j]);
                }
                printf("]");
            }
            printf("\n");
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
}



vfloat integrand(
        const integration_input_t* data,
        const integration_variables_t* vars
        )
{
    // Initialize bare_scalar_products, alpha and beta tables
    vfloat bare_scalar_products[N_COEFFS][N_COEFFS];
    matrix_vfloat* alpha = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);
    matrix_vfloat* beta  = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);
    compute_bare_scalar_products(data->k,vars,bare_scalar_products);
    compute_alpha_beta_tables(bare_scalar_products,alpha,beta);

    // Allocate space for kernels (calloc also initializes values to 0)
    kernel_value_t* kernels =
        (kernel_value_t*)calloc(COMPONENTS * N_KERNELS, sizeof(kernel_value_t));

    diagram_t diagrams[N_DIAGRAMS];
    possible_diagrams(diagrams)

    vfloat result = 0;

    // Free allocated memory
    free(kernels);
    gsl_matrix_free(alpha);
    gsl_matrix_free(beta);

    return result;
}



int cuba_integrand(
        const int *ndim,
        const cubareal xx[],
        const int *ncomp,
        cubareal ff[],
        void *userdata
        )
{
    integration_input_t* data = (integration_input_t*)userdata;

    integration_variables_t vars;

    for (int i = 0; i < LOOPS; ++i) {
        vars.magnitudes[i] = xx[i];
        vars.cos_theta[i]  = xx[i + LOOPS];
    }
    for (int i = 0; i < LOOPS - 1; ++i) {
        vars.phi[i] = xx[i + 2*LOOPS];
    }

    ff[0] = integrand(data,&vars);

    return 0;
}
