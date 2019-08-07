/*
   integrand.h

   Created by Petter Taule on 24.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef INTEGRAND_H
#define INTEGRAND_H

#include <gsl/gsl_spline.h>
#include "constants.h"
#include "tables.h"
#include "diagrams.h"
#include "worker_mem.h"

// Struct for storing overall data to be used in integration
typedef struct {
    vfloat k;
    short int component_a;
    short int component_b;
    gsl_interp_accel* ps_acc;
    gsl_spline* ps_spline;
    const diagram_t* const diagrams;
    const evolution_params_t* params;
    worker_mem_t* worker_mem;
} integration_input_t;

int diagram_factor(const diagram_t* diagram);
int symmetrization_factor(const diagram_t* diagram);

void possible_diagrams(diagram_t diagrams[]);

void find_kernel_arguments(
        const diagram_t* diagram,
        const short int rearrangement[],
        const short int signs[],
        short int arguments_l[],
        short int arguments_r[]
        );

vfloat sign_flip_symmetrization(
        const short int rearrangement[],
        const diagram_t* diagram,
        const integration_input_t* input,
        const table_ptrs_t* tables
        );

vfloat loop_momenta_symmetrization(
        const diagram_t* diagram,
        const integration_input_t* input,
        const table_ptrs_t* tables
        );

vfloat integrand(const integration_input_t* data,
        table_ptrs_t* tables);

#endif /* ifndef INTEGRAND_H */
