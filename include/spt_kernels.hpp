/*
   spt_kernels.hpp

   Created by Petter Taule on 03.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef SPT_KERNELS_HPP
#define SPT_KERNELS_HPP

#include "utilities.hpp"
#include "tables.hpp"

short int compute_SPT_kernels(
        const short int arguments[],
        short int kernel_index,
        short int n,
        IntegrandTables& tables
        );

#endif /* ifndef SPT_KERNELS_HPP */
