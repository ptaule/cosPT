/*
   main.cpp

   Created by Petter Taule on 28.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

#include <getopt.h>
#include <gsl/gsl_sf.h>
#include <cuba.h>

#include "include/utilities.hpp"
#include "include/tables.hpp"
#include "include/diagrams.hpp"
#include "include/interpolation.hpp"
#include "include/bispectrum_tree_level.hpp"
#include "include/kernel_evolution.hpp"
#include "include/integrand.hpp"
#include "include/io.hpp"

#include "include/combinatorics.hpp"


void print_help();
void numbers_and_bars();


int main(int argc, char *argv[]) {
    numbers_and_bars();
}



void numbers_and_bars() {
    /* Vec1D<std::string> a = {"k1", "k2", "k3", "k4", "k5"}; */
    // Vec1D<std::string> a = {"k1", "k2", "k3", "k4"};
    Vec1D<std::string> a = {"k1", "k2", "k3"};
    int n = a.size();

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
}



void print_help() {
    const char* help = R"(CosPT

CosPT computes loop corrections to the power spectrum or bispectrum in
cosmological perturbation theory. The program takes a configuration file as
argument, see README.md for how to enter settings in this file. Some settings
may also be set as a command line options, see below.

Usage:
  cosPT config_file
  cosPT config_file [--k_a_idx=IDX] [--n_evals=NUM] [--n_cores=NUM]
  cosPT config_file [--k_a_idx=IDX] [--k_b_idx=IDX] [--k_c_idx=IDX]
                    [--n_evals=NUM] [--n_cores=NUM]
  cosPT -h | --help

Options:
  -h --help     Show this screen.
  -a --k_a_idx  External wavenumber to compute, index in k_a_grid
                file. Overrides index given in configuration file.
  -b --k_b_idx  Similarly for second external wavenumber (bispectrum).
  -c --k_c_idx  Similarly for third external wavenumber (bispectrum).
  -n --n_evals  Number of evaluations in Monte Carlo integration.
  -p --n_cores  Number of cores/threads for CUBA to spawn, in addition to master.
            )";
    std::cout << help;
    return;
}
