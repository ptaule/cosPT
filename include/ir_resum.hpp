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


void compute_ir_damping(
        const Interpolation1D& ps_nw,
        double k_min,
        double k_s,
        double k_osc,
        double& Sigma2,
        double& delta_Sigma2,
        size_t sub_regions = 10000,
        double atol        = 0,
        double rtol        = 1e-6,
        int key            = GSL_INTEG_GAUSS61
        );


#endif /* ifndef IR_RESUM_HPP */
