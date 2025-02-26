#ifndef RSD_HPP
#define RSD_HPP


class IntegrandTables;


int compute_rsd_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );


#endif /* !RSD_HPP */
