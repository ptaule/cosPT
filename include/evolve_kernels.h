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

short int kernel_evolution(
        const short int arguments[],
        short int index,
        short int n,
        gsl_matrix* omega,
        const parameters_t* params,
        const table_pointers_t* data_tables
        );


#endif /* ifndef EVOLVE_KERNELS_H */
