/*
   spt_kernels.hpp

   Created by Petter Taule on 03.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef SPT_KERNELS_HPP
#define SPT_KERNELS_HPP

class IntegrandTables;

int compute_SPT_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );

#endif /* ifndef SPT_KERNELS_HPP */
