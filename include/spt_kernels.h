/*
   spt_kernels.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef SPT_KERNELS_H
#define SPT_KERNELS_H

#include "constants.h"
#include "tables.h"

vfloat compute_SPT_kernel(
        const short int arguments[],
        short int kernel_index,
        short int n,
        short int component,
        tables_t* tables
        );


#endif /* ifndef SPT_KERNELS_H */
