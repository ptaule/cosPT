/*
   spt_kernels.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef SPT_KERNELS_H
#define SPT_KERNELS_H

#include "constants.h"
#include "tables.h"

short int compute_SPT_kernels(
        const short int arguments[],
        short int n,
        const table_pointers_t* data_tables
        );


#endif /* ifndef SPT_KERNELS_H */
