/*
   rsd.hpp

   Created by Petter Taule on 23.09.2022
   Copyright (c) 2022 Petter Taule. All rights reserved.
*/

#ifndef RSD_HPP
#define RSD_HPP

class IntegrandTables;

void compute_rsd_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );

#endif /* !RSD_HPP */
