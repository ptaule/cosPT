#ifndef BIASED_TRACERS_HPP
#define BIASED_TRACERS_HPP

class IntegrandTables;

double compute_rsd_biased_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );

#endif /* BIASED_TRACERS_HPP */
