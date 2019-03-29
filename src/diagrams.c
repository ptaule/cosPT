/*
   diagrams.c

   Created by Petter Taule on 27.03.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_combination.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/kernels.h"
#include "../include/diagrams.h"



// Find (distinct) diagrams for L-loop
// They satisfy: m >= 1; l,r > 0; l + r + m = L + 1
static void possible_diagrams(diagram_t diagrams[]) {
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



static void integrand_kernel_arguments(
        short int m,
        short int l,
        short int r,
        const short int rearrangement[],
        const bool signs[],
        short int arguments_l[N_KERNEL_ARGS], /* out */
        short int arguments_r[N_KERNEL_ARGS]  /* out */
        )
{
    /* First argument : k1 - (±k2) - (±k3) - ... */
    short int config[N_COEFFS] = {};
    config[N_COEFFS - 1] = 1;
    for (int i = 2; i <= m; ++i) {
        // Note extra minus sign from kernel expression
        config[rearrangement[i-2]] = (signs[i-2] ? - 1 : 1);
    }

    arguments_l[0] = config2label(config,N_COEFFS);
    arguments_r[0] = config2label(config,N_COEFFS);

    // Reset config
    memset(config,0,sizeof(config));

    // Argument indices
    size_t index_l = 1;
    size_t index_r = 1;

    // k2,k3,...,km arguments, the ordering of which is stored by the first
    // (m-1) entries of rearrangement[]. Note that the signs are opposite of
    // corresponding terms in the first argument
    for (size_t i = 2; i <= m; ++i) {
        config[rearrangement[i-2]] = (signs[i-2] ? 1 : -1);
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



static void sign_flips(diagram_t* diagram) {
    short int m = diagram->m;

    // Allocate memory for sign configurations
    diagram->sign_configs = (bool**)
        calloc(diagram->n_sign_configs,  sizeof(bool*));

    for (int j = 0; j < diagram->n_sign_configs; ++j) {
        // Signs of k2,...,km can be flipped, hence allocate (m-1)-length array
        diagram->sign_configs[j] = (bool*)calloc(m - 1, sizeof(bool));
    }

    // Loop over possible sign configurations and store in diagram->sign table
    for (int i = 0; i < diagram->n_sign_configs; ++i) {
        bool* const signs = diagram->sign_configs[i];
        for (int j = 0; j < (m-1); ++j) {
            // Add 1 so that the first sign configuration is +1,+1,...,+1
            signs[j] = (i/(int)pow(2,j) + 1) % 2;
        }
    }
}



static void rearrangements(diagram_t* diagram) {
    short int m = diagram->m;
    short int l = diagram->l;
    short int r = diagram->r;

    // Allocate memory for rearrangements
    diagram->rearrangements = (short int**)
        calloc(diagram->n_rearrangements,sizeof(short int*));
    for (int j = 0; j < diagram->n_rearrangements; ++j) {
        diagram->rearrangements[j] = (short int*)calloc(LOOPS,sizeof(short int));
    }

    // Default loop-momenta ordering (Q1, Q2, Q3, etc.)
    short int loop_momenta[LOOPS];
    for (int i = 0; i < LOOPS; ++i) loop_momenta[i] = i;

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

    // Rearrangement index (should be between 0 and n_rearrangements)
    short int rearrangement_index = 0;
    // Pointer alias for convenience/readability
    short int* rearrangement;

    // Go through possible (m-1)-groupings
    do {

        if (m < LOOPS + 1) {
            gsl_combination_init_first(comb_l);
            gsl_combination_init_last(comb_r);

            // Go through possible l and r groupings
            do {
                rearrangement = diagram->rearrangements[rearrangement_index];
                for (int i = 0; i < m-1; ++i) {
                    rearrangement[i] = loop_momenta[gsl_combination_get(comb_m,i)];
                }

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

                rearrangement_index++;

            } while (
                    gsl_combination_next(comb_l) == GSL_SUCCESS &&
                    gsl_combination_prev(comb_r) == GSL_SUCCESS
                    );
        }
        else {
            rearrangement = diagram->rearrangements[rearrangement_index];
            for (int i = 0; i < m-1; ++i) {
                rearrangement[i] = loop_momenta[gsl_combination_get(comb_m,i)];
            }
            rearrangement_index++;
        }
    } while (
            gsl_combination_next(comb_m) == GSL_SUCCESS &&
            gsl_combination_prev(comb_s) == GSL_SUCCESS
            );

    if (rearrangement_index != diagram->n_rearrangements) {
        warning_verbose("Created %d rearrangements, but n_rearrangements = %d.",
                rearrangement_index, diagram->n_rearrangements);
    }

    if (m < LOOPS + 1) {
        gsl_combination_free(comb_l);
        gsl_combination_free(comb_r);
    }
    gsl_combination_free(comb_m);
    gsl_combination_free(comb_s);
}



void initialize_arguments(diagram_t* diagram) {
    // Allocate memory for argument/kernel_index configurations
    diagram->argument_configs_l = (short int***)
        calloc(diagram->n_rearrangements,sizeof(short int**));
    diagram->argument_configs_r = (short int***)
        calloc(diagram->n_rearrangements,sizeof(short int**));

    for (int a = 0; a < diagram->n_rearrangements; ++a) {
        diagram->argument_configs_l[a] = (short int**)
            calloc(diagram->n_sign_configs, sizeof(short int*));
        diagram->argument_configs_r[a] = (short int**)
            calloc(diagram->n_sign_configs, sizeof(short int*));

        for (int b = 0; b < diagram->n_sign_configs; ++b) {
            diagram->argument_configs_l[a][b] = (short int*)
                calloc(N_KERNEL_ARGS, sizeof(short int));
            diagram->argument_configs_r[a][b] = (short int*)
                calloc(N_KERNEL_ARGS, sizeof(short int));
        }
    }

    // Initialize arguments
    for (int a = 0; a < diagram->n_rearrangements; ++a) {
        for (int b = 0; b < diagram->n_sign_configs; ++b) {
            integrand_kernel_arguments(diagram->m, diagram->l, diagram->r,
                    diagram->rearrangements[a], diagram->sign_configs[b],
                    diagram->argument_configs_l[a][b],
                    diagram->argument_configs_r[a][b]
                    );
        }
    }
}



void initialize_diagrams(diagram_t diagrams[]) {
    possible_diagrams(diagrams);

    for (int i = 0; i < N_DIAGRAMS; ++i) {
        short int m = diagrams[i].m;
        short int l = diagrams[i].l;
        short int r = diagrams[i].r;
#if DEBUG >= 1
        if ((m + l + r) != (LOOPS + 1)) warning("m + r + l != LOOPS + 1.");
#endif

        diagrams[i].diagram_factor = (gsl_sf_fact(2*l + m) * gsl_sf_fact(2*r + m)) /
            (pow(2,l+r) * gsl_sf_fact(l) * gsl_sf_fact(r) * gsl_sf_fact(m));
        // If diagram is antisymmetric in l <-> r, multiply diagram_factor by 2
        // (the possible_diagrams()-function only returns diagrams with l >= r)
        if (l != r) diagrams[i].diagram_factor *= 2;

        diagrams[i].n_rearrangements = gsl_sf_fact(LOOPS) /
            (gsl_sf_fact(m-1) * gsl_sf_fact(l) * gsl_sf_fact(r));
        diagrams[i].n_sign_configs = pow(2,m-1);

        rearrangements(&(diagrams[i]));
        sign_flips(&(diagrams[i]));
        initialize_arguments(&(diagrams[i]));
    }
}



void diagrams_gc(diagram_t diagrams[]) {
    for (int i = 0; i < N_DIAGRAMS; ++i) {
        // Shorthand alias for diagrams[i]
        diagram_t* const dg = &diagrams[i];

        for (int a = 0; a < dg->n_rearrangements; ++a) {
            free(dg->rearrangements[a]);
            for (int b = 0; b < dg->n_sign_configs; ++b) {
                free(dg->argument_configs_l[a][b]);
                free(dg->argument_configs_r[a][b]);
            }
            free(dg->argument_configs_l[a]);
            free(dg->argument_configs_r[a]);
        }
        for (int b = 0; b < dg->n_sign_configs; ++b) {
            free(dg->sign_configs[b]);
        }
        free(dg->rearrangements);
        free(dg->sign_configs);
        free(dg->argument_configs_l);
        free(dg->argument_configs_r);
    }
}



void print_diagram(const diagram_t* diagram) {
    printf(ANSI_COLOR_MAGENTA "(m,l,r) = (%d,%d,%d)\n" ANSI_COLOR_BLUE
            ,diagram->m,diagram->l,diagram->r);

    for (int a = 0; a < diagram->n_rearrangements; ++a) {
        for (int b = 0; b < diagram->n_sign_configs; ++b) {
            printf("\t");
            print_labels(diagram->argument_configs_l[a][b]);
            printf("\t");
            print_labels(diagram->argument_configs_r[a][b]);
            printf("\n");
        }
    }
    printf(ANSI_COLOR_RESET);
}
