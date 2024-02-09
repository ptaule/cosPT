#ifndef RSD_HPP
#define RSD_HPP

#include <cstddef>

#include <gsl/gsl_integration.h>

#include "utilities.hpp"

class IntegrandTables;


double compute_rsd_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );


#endif /* !RSD_HPP */
