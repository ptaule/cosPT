/*
   evolve_kernels.h

   Created by Petter Taule on 04.04.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef EVOLVE_KERNELS_H
#define EVOLVE_KERNELS_H

#include "constants.h"
#include "tables.h"

short int kernel_evolution(
        const short int arguments[],
        short int index,
        short int n,
        const evolution_params_t* params,
        tables_t* tables
        );


#endif /* ifndef EVOLVE_KERNELS_H */
