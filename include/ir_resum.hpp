/*
   ir_resum.hpp

   Created by Petter Taule on 15.02.2023
   Copyright (c) 2023 Petter Taule. All rights reserved.
*/

#ifndef IR_RESUM_HPP
#define IR_RESUM_HPP

#include <cstddef>

#include <gsl/gsl_integration.h>

#include "utilities.hpp"

class Interpolation1D;


void DST_II(Vec1D<double>& data);
void DST_III(Vec1D<double>& data);


struct IRresumSettings {
    int n_loops = 1;  /* Loop correction */
    int pt_order = 1; /* At which order are we working? To avoid overcounting
                         of soft loops */

    /* Wiggly/non-wiggly split */
    double k_max = 10;
    int N_power = 16;
    int N_left = 120;
    int N_right = 240;

    /* Damping factor */
    double k_min = 1e-4;
    double k_s = 0.2;
    double k_osc = 1.0/110.0;
    std::size_t integrate_sub_regions = 10000;
    double integrate_atol = 0;
    double integrate_rtol = 1e-6;
    int integrate_key = GSL_INTEG_GAUSS61;
};


/* Remove BAO wiggles. Default parameters taken from CLASS-PT, i.e. 2004.10607 */
void remove_BAO_wiggles(
        const Interpolation1D& ps,
        Interpolation1D& ps_nw,
        const IRresumSettings& settings
        );


void compute_ir_damping(
        const Interpolation1D& ps_nw,
        double& Sigma2,       /* out */
        double& delta_Sigma2, /* out */
        const IRresumSettings& settings
        );


#endif /* ifndef IR_RESUM_HPP */
