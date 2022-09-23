/*
   rsd.cpp

   Created by Petter Taule on 18.10.2021
   Copyright (c) 2021 Petter Taule. All rights reserved.
*/

#include <algorithm>

#include "../include/combinatorics.hpp"
#include "../include/kernel_evolution.hpp"
#include "../include/rsd.hpp"
#include "../include/spt_kernels.hpp"



double rsd_kernel(
        const int arguments[],
        int n,
        IntegrandTables& tables
        )
{
    double result = 0;

    Vec1D<int> current_ordering(n);

    /* m is number of groups that the wavenumbers is divided into. The first
     * group corresponds to (delta + f mu^2 theta), so m >= 1. The expansion of
    * exp(.. theta) gives m = 2,...,n-1 */
    for (int m = 1; m <= n; ++m) {
        /* Go through combinations of "placing (m-1) bars", i.e. divide wavenumbers
        * amoung m groups */
        Combinations comb(n - 1, m - 1);

        /* Current combination and group sizes */
        Vec1D<int> current(m-1);
        Vec1D<int> group_sizes(m);
        do {
            comb.write_current_combination(current);
            std::fill(group_sizes.begin(), group_sizes.end(), 0);

            int elements_left = n;
            if (!current.empty()) {
                /* There is always an element before first bar */
                group_sizes.at(0) = current.at(0) + 1;
                elements_left -= group_sizes.at(0);
            }

            for (size_t i = 1; i < current.size(); ++i) {
                group_sizes.at(i) = current.at(i) - current.at(i - 1);
                elements_left -= group_sizes.at(i);
            }
            group_sizes.back() = elements_left;

            Orderings orderings(n, group_sizes);

            do {
                orderings.write_current(current_ordering );

                size_t i = 0;
                for (auto &el : group_sizes) {
                    for (int j = 0; j < el; ++j) {
                        std::cout << a.at(current_ordering.at(i)) << " ";
                        ++i;
                    }
                    std::cout << "| ";
                }
                std::cout << std::endl;
            } while (orderings.next());

        } while (comb.next());
        std::cout << std::endl;
    }


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
