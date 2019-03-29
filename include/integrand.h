/*
   integrand.h

   Created by Petter Taule on 24.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef INTEGRAND_H
#define INTEGRAND_H

#include <gsl/gsl_spline.h>
#include "constants.h"
#include "kernels.h"
#include "diagrams.h"

// Struct for storing overall data to be used in integration
typedef struct {
    vfloat k;
    short int component_a;
    short int component_b;
    const diagram_t* const diagrams;
    gsl_interp_accel* const acc;
    gsl_spline* const spline;
} integration_input_t;

vfloat integrand(const integration_input_t* data, const integration_variables_t* vars);

#endif /* ifndef INTEGRAND_H */
