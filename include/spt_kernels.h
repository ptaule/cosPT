/*
   spt_kernels.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef SPT_KERNELS_H
#define SPT_KERNELS_H

#include "constants.h"
#include "kernels.h"

vfloat partial_SPT_sum(
        const short int arguments[],
        const short int n,
        const short int m,
        const short int a,
        const short int b,
        const table_pointers_t* data_tables
        );

vfloat compute_SPT_kernel(
        const short int arguments[],
        short int n,
        short int component,
        const table_pointers_t* data_tables
        );


#endif /* ifndef SPT_KERNELS_H */
