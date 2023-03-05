/*
   rsd.hpp

   Created by Petter Taule on 23.09.2022
   Copyright (c) 2022 Petter Taule. All rights reserved.
*/

#ifndef RSD_HPP
#define RSD_HPP

#include <cstddef>

#include <gsl/gsl_integration.h>

#include "utilities.hpp"

class IntegrandTables;
class InputPowerSpectrum;


Vec1D<double> rsd_tree_level(
    double k,
    const InputPowerSpectrum& ps,
    std::size_t integrate_sub_regions = 10000,
    double integrate_atol = 0,
    double integrate_rtol = 1e-6,
    int integrate_key = GSL_INTEG_GAUSS61
);


double compute_rsd_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );


#endif /* !RSD_HPP */
