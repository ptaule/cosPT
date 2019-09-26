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

// Struct for storing overall data to be used in integration
typedef struct {
    vfloat k;
    short int component_a;
    short int component_b;
    gsl_interp_accel* ps_acc;
    gsl_spline* ps_spline;
    const diagram_t* const diagrams;
    const evolution_params_t* params;
    tables_t* worker_mem;
} integration_input_t;

vfloat integrand(const integration_input_t* data,
        tables_t* tables);

#endif /* ifndef INTEGRAND_H */
