/*
   rsd.cpp

   Created by Petter Taule on 18.10.2021
   Copyright (c) 2021 Petter Taule. All rights reserved.
*/

#include "../include/combinatorics.hpp"
#include "../include/kernel_evolution.hpp"
#include "../include/rsd.hpp"
#include "../include/spt_kernels.hpp"


double rsd_kernel_unsymmetrized(
        const int arguments[],
        int n,
        IntegrandTables& tables
        )
{
    double result = 0;

    for (int m = 1; m <= n; ++m) {
        double partial_result = 0;

        /* (n-1) choose (m - 1) possible terms */
        Combinations comb(n-1, m-1);
        do {
            Vec1D<int> current = comb.get_current_combination();

            for (auto& el : current) {

            }
        } while(comb.next());

        result += partial_result;
    }

    return result;
}
