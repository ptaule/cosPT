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

using std::pow;

void print_help();

int main(int argc, char* argv[]) {
    /* Vec1D<std::string> a = {"k1", "k2", "k3", "k4", "k5"}; */
    Vec1D<std::string> a = {"k1", "k2", "k3", "k4"};
    /* Vec1D<std::string> a = {"k1", "k2", "k3"}; */
    int n = a.size();

    for (int m = 1; m <= n; ++m) {
        Combinations comb(n-1,m-1);

        do {
            Vec1D<int> current = comb.get_current_combination();

            int i = 0;
            for (auto& el : current) {
                while (i < el + 1) {
                    std::cout << a.at(i) << "\t";
                    ++i;
                }
                std::cout << "|\t";
            }
            for (; i < n; ++i) {
                std::cout << a.at(i) << "\t";
            }
            std::cout << std::endl;
        } while (comb.next());
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
