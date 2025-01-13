#ifndef BIASED_TRACERS_HPP
#define BIASED_TRACERS_HPP

#include "utilities.hpp"

class IntegrandTables;

double compute_rsd_biased_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables,
        const Vec1D<double>& biases
        );

#endif /* BIASED_TRACERS_HPP */
