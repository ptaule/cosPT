/*
   evolve_kernels.h

   Created by Petter Taule on 04.04.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef EVOLVE_KERNELS_H
#define EVOLVE_KERNELS_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "constants.h"
#include "kernels.h"

typedef struct {
    vfloat omega_m0;
    vfloat f2;
    vfloat f_nu;
    vfloat eta_i;
    vfloat eta_f;
    matrix_t* omega;
} parameters_t;


short int compute_time_dependent_kernels(
        const short int arguments[],
        short int n,
        const parameters_t* params,
        const table_pointers_t* data_tables
        );


#endif /* ifndef EVOLVE_KERNELS_H */
