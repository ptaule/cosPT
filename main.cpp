/*
   main.cpp

   Created by Petter Taule on 28.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <vector>
#include <cmath>

#include <cuba.h>

#include "include/utilities.hpp"
#include "include/tables.hpp"
#include "include/spt_kernels.hpp"
#include "include/kernel_evolution.hpp"

using std::size_t;

int main() {
    /* Settings */
    double k  = 0.1;
    double Q1 = 0.02;
    double Q2 = 0.07;
    double cos_kQ1 = 0.2; /* cosine of angle between k and Q1 */
    double cos_kQ2 = 0.3;
    double phi_2 = 0;     /* Azimuthal angle of Q2 */
    /* (The program uses the freedom to choose coordinate system so that k is
     * aligned with the z-axis and Q1 has no azimuthal angle) */

    size_t time_steps = 100;
    /* Initial/finial time, eta = ln(D) */
    double eta_ini = -3;
    double eta_fin = 0;

    /* zeta = omega_m / f^2 */
    std::string zeta_file =
        "/space/ge52sir/class_public/output/quijote/fiducial/zeta_of_etaD.dat";

    /* Precision settings for ode solver */
    double ode_atol   = 1e-6;
    double ode_rtol   = 1e-4;
    double ode_hstart = 1e-3;

    try {
        /* Setup various table structures */
        int n_loops = 2;

        LoopParameters loop_params(n_loops, POWERSPECTRUM, EVOLVE_IC_EDS);
        int zero_label = loop_params.zero_label();
        SumTable sum_table(loop_params);

        EvolutionParameters ev_params =
            EvolutionParameters(zeta_file, ode_atol, ode_rtol, ode_hstart);
        EtaGrid eta_grid = EtaGrid(time_steps, eta_ini, eta_fin);

        IntegrandTables tables(k, loop_params, sum_table, ev_params, eta_grid);

        tables.vars.magnitudes.at(0) = Q1;
        tables.vars.magnitudes.at(1) = Q2;
        tables.vars.cos_theta.at(0)  = cos_kQ1;
        tables.vars.cos_theta.at(1)  = cos_kQ2;
        tables.vars.phi.at(1)        = phi_2;

        tables.compute_tables();

        /* Compute kernels */

        /* Working at 2-loop, we have three momenta Q1, Q2 (loop) and k
         * (external). Each argument of the kernels can be represented two
         * ways:
         * - as a vector with three elements that are either 0,-1,+1
         *   corresponding to the sign of Q1,Q2, and k.
         * - as an integer label between 0 and 26 (at 2-loop). This
         *   representation is what is used in the code.
         * One can convert the
         * vector representation to the integer label with the function
         * int config2label(const Vec1D<int>& config); */

        /* Example: F4(k-Q1,Q1,Q2,-Q2). */
        int n = 4;
        Vec1D<int> F4_arg_1 = {-1,0,1}; /* -Q1+k */
        Vec1D<int> F4_arg_2 = {1,0,0};  /* +Q1   */
        Vec1D<int> F4_arg_3 = {0,1,0};  /* +Q2   */
        Vec1D<int> F4_arg_4 = {0,-1,0}; /* -Q2   */

        /* Arguments in integer label representation. This vector of arguments
         * must be five elements long (at 2-loop), therefore we use
         * "zero_label" to indicate no argument */
        Vec1D<int> args = {
            config2label(F4_arg_1),
            config2label(F4_arg_2),
            config2label(F4_arg_3),
            config2label(F4_arg_4),
            zero_label
        };

        /* Compute EdS-SPT kernels, needed for the IC. The result(s) are stored
         * in the tables-structure, and the index which is returned can be used
         * to look up the result in this table later. */
        int index = compute_SPT_kernels(args.data(), -1, n, tables);
        /* Compute dynamically evolved kernels. Again the results are stored in
         * the tables-structure, and the index returned should be the same as
         * above (we may give it as an argument here so that the function need
         * not compute it again). */
        kernel_evolution(args.data(), index, n, tables);

        /* Print results to stdout at some given time step t. The print_labels
         * function prints the corresponding k,Q1,Q2 arguments from the labels.
         * Print EdS-SPT results first, then evolved kernels in second column.
         * */
        size_t t = time_steps - 1; /* The last time step */
        std::cout << "F" << n;
        print_labels(args, loop_params.n_coeffs(), loop_params.spectrum(),
                     std::cout);
        std::cout << "\t"
                  << tables.spt_kernels.at(index).values[0] << "\t"
                  << tables.kernels.at(index).values.at(t).at(0) << "\t"
                  << std::endl;
        std::cout << "G" << n;
        print_labels(args, loop_params.n_coeffs(), loop_params.spectrum(),
                     std::cout);
        std::cout << "\t"
                  << tables.spt_kernels.at(index).values[1] << "\t"
                  << tables.kernels.at(index).values.at(t).at(1) << "\t"
                  << std::endl;

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\nExiting." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
